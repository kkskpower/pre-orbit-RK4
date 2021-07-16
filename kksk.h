#include <iostream>
#include <Eigen\Dense>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include<iomanip>
#include <vector>
#include <iostream>
#include<cmath>
using namespace std;
using namespace Eigen;
class kksk
{
public:
	vector<vector<vector<double>>> allA;
	vector<vector<double>> ally;
	vector<vector<double>> Y;
	void ini(string file);
	vector<vector<double> > everydata;//读入的数据
	MatrixXd A;
	void toA();//给A赋初值
	void toyb(double y10, double y20, double y30,double y40,double b0);//给y、b赋初值
	void fun2();
	MatrixXd G;
	/*int round(double a);*/
	void funy();//得到300s内的y1 y2
private:
	vector<double>yb;
	vector<string> anydata;
	vector<double> dealLine(string str);
	void pushevery(void);
	struct BC
	{
		MatrixXd C;
		MatrixXd B;
	};
	BC toBC(double t, double y1, double y2);//得到B C
	void fun1(double t, double* y, double* g);
};
