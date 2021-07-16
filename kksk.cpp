#include "kksk.h"
void kksk::toyb(double y10, double y20, double y30, double y40, double b0)
{
	vector<double>YB(5);
	YB[0] = y10;
	YB[1] = y20;
	YB[2] = y30;
	YB[3] = y40;
	YB[4] = b0;
	yb = YB;
}
//int kksk::round(double a)
//{
//	int b;
//	if (a > 0) {
//		b = (a * 2 + 1) / 2;
//	}
//	else {
//		b = (a * 2 - 1) / 2;
//	}
//	return b;
//}
void kksk::ini(string file)
{
	ifstream infile;
	infile.open(file);   //将文件流对象与文件连接起来 
	assert(infile.is_open());   //若失败,则输出错误消息,并终止程序运行 

	string s;
	int count = -1;
	int flag = false;
	int aa = 0;
	while (getline(infile, s))
	{
		if (s.substr(0, 5) == "    0")
		{
			count = 0;
			int aa = 1;
		}
		if (count == 0 || aa == 1)
		{
			anydata.push_back(s);

			//count = -1;
		}
	}
	infile.close();             //关闭文件输入流 
	pushevery();
}
vector<double> kksk::dealLine(string str)
{
	//以下处理空格切分
	vector<string> res;
	string result;
	stringstream input(str);
	while (input >> result)
		res.push_back(result);
	//字符转化成double数据
	vector<double> data;
	for (size_t i = 0; i < res.size(); i++)
	{
		data.push_back(stod(res[i]));
	}
	return data;
}
void kksk::pushevery()
{
	int k = anydata.size();
	vector<double> data;
	for (int i = 0; i < k; i++)
	{
		data = dealLine(anydata[i]);
		everydata.push_back(data);
	}
}
void kksk::toA()
{
	MatrixXd D(4, 5);
	D << 1, 0, 0, 0, 0,
		0, 1, 0, 0, 0,
		0, 0, 1, 0, 0,
		0, 0, 0, 1, 0;
	A = D;
}
kksk::BC kksk::toBC(double t, double y1, double y2)
{
	MatrixXd C(4, 5);
	C << 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, sin(t),
		0, 0, 0, 0, cos(t);
	double y12;
	y12 = sqrt(pow(y1, 2) + pow(y2, 2));
	MatrixXd B(4, 4);
	B << 0, 0, 1, 0,
		0, 0, 0, 1,
		-1 / pow(y12, 3) + 3 * pow(y1, 2) / pow(y12, 5), 3 * y1 * y2 / pow(y12, 5), 0, 0,
		3 * y1 * y2 / pow(y12, 5), -1 / pow(y12, 3) + 3 * pow(y2, 2) / pow(y12, 5), 0, 0;
	BC bc;
	bc.B = B;
	bc.C = C;
	return bc;
}
void kksk::fun1(double t, double* y, double* g)
{
	g[0] = y[2];
	g[1] = y[3];
	double y12 = y[0] * y[0] + y[1] * y[1];
	y12 = pow(y12, 1.5);
	g[2] = yb[4] * sin(t) - y[0] / y12;
	g[3] = yb[4] * cos(t) - y[1] / y12;
}
void kksk::funy()
{
	double k1[4], k2[4], k3[4], k4[4], y1[4], y1_0[4], y1_1[4], y1_2[4], y1_3[4];
	double h1 = 0.001, e = 0.1;
	double g[4];
	int t;
	vector<double> y(2);
	t = 1200;
	y1[0] = yb[0];
	y1[1] = yb[1];
	y1[2] = yb[2];
	y1[3] = yb[3];
	y[0] = y1[0];
	y[1] = y1[1];
	Y.push_back(y);
	vector<double> ally1;
	int c = round(t / h1);
	int d = 0;
	for (int i = 0; i < c; i++)  // RK4
	{
		double tt;
		tt = i * h1;
		fun1(tt, y1, g);
		for (int j = 0; j < 4; j++)
			k1[j] = g[j];
		for (int j = 0; j < 4; j++)
			y1_1[j] = y1[j] + 0.5 * h1 * g[j];
		tt = i * h1 + 0.5 * h1;
		fun1(tt, y1_1, g);
		for (int j = 0; j < 4; j++)
			y1_2[j] = y1[j] + 0.5 * h1 * g[j];
		for (int j = 0; j < 4; j++)
			k2[j] = g[j];
		fun1(tt, y1_2, g);
		for (int j = 0; j < 4; j++)
			y1_3[j] = y1[j] + h1 * g[j];
		for (int j = 0; j < 4; j++)
			k3[j] = g[j];
		tt = (i + 1.0) * h1;
		fun1(tt, y1_3, g);
		for (int j = 0; j < 4; j++)
			k4[j] = g[j];

		for (int j = 0; j < 4; j++)
		{
			double phi = k1[j] / 6.0 + k2[j] / 3.0 + k3[j] / 3.0 + k4[j] / 6.0;
			y1[j] += h1 * phi;
		}
		y[0] = y1[0];
		y[1] = y1[1];
		Y.push_back(y);
		if ((i + 1) % 20 == 0)
		{
			for (int j = 0; j < 4; j++)
			{
				ally1.push_back(y1[j]);
			}
			ally.push_back(ally1);
			ally1.clear();
		}
	}

}
void kksk::fun2()
{
	int t = 1200;
	double h = 0.001;
	double h1 = 0.001;
	int c = round(t / h);
	for (int i = 0; i < c; i++)  // RK4
	{
		MatrixXd K1, K2, K3, K4, A1, A2, A3;
		double tt;
		vector<double> y(2);
		tt = i * h;
		int j = round(tt / h1);
		y = Y[j];
		BC bc = toBC(tt, y[0], y[1]);
		G = bc.B * A + bc.C;
		/*		for (int m = 0; m < 4; m++)
			{
				for (int n = 0; n < 5; n++)
				{
					sy1.push_back(A(m, n));
				}
				sy.push_back(sy1);
			sy1.clear();
			}*/
		K1 = G;
		A1 = A + 0.5 * h * G;
		tt = i * h + 0.5 * h;
		y[0] = (Y[i][0]+Y[i+1][0])/2;
		y[1] = (Y[i][1] + Y[i + 1][1])/2;
		bc = toBC(tt, y[0], y[1]);
		G = bc.B * A1 + bc.C;
		A2 = A + 0.5 * h * G;
		K2 = G;
		G = bc.B * A2 + bc.C;
		A3 = A + h * G;
		K3 = G;
		tt = (i + 1.0) * h;
		j = round(tt / h1);
		y = Y[j];
		bc = toBC(tt, y[0], y[1]);
		G = bc.B * A3 + bc.C;
		K4 = G;
		MatrixXd phi = K1 / 6.0 + K2 / 3.0 + K3 / 3.0 + K4 / 6.0;
		A = A + h * phi;
		vector<vector<double>> allA1;
		vector<double> allA2;
		if ((i + 1) % 20 == 0)
		{
			for (int m = 0; m < 4; m++)
			{
				for (int n = 0; n < 5; n++)
				{
					allA2.push_back(A(m, n));
				}
				allA1.push_back(allA2);
				allA2.clear();
			}
			allA.push_back(allA1);
			allA1.clear();
		}
	}
}

