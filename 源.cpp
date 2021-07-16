#include "kksk.h"
int main()
{
    double e = 0.1;
    double y10 = 1 - e;
    double y20 = 0;
    double y30 = 0;
    double y40= sqrt((1.0 + e) / (1.0 - e));
    double b0 = 1.005e-6;
    double rb = 2e-6;
    while (rb>1e-6)
    {
        kksk wow1;
        wow1.toyb(y10, y20, y30, y40, b0);//�����״̬��������ѧģ�Ͳ�������ֵ
        string file = "E:\\360MoveData\\Users\\DELL\\Desktop\\������\\Projecteig\\Obs.dat";
        wow1.ini(file);
        wow1.toA();
        wow1.funy();
        wow1.fun2();
        vector<double> y(2);
        MatrixXd B(60000, 5);
        MatrixXd BTB;
        MatrixXd BTl;
        MatrixXd l(60000, 1);
        for (int i = 0; i < 60000; i++)//���濪ʼ������С���˼���
        {
            y[0] = wow1.ally[i][0];
            y[1] = wow1.ally[i][1];
            double r0 = sqrt(pow(y[0], 2) + pow(y[1], 2));
            l(i, 0) = wow1.everydata[i][1] - r0;
            MatrixXd Hi(1, 4);
            Hi << y[0] / r0, y[1] / r0, 0, 0;
            MatrixXd Ai(4, 5);
            for (int m = 0; m < 4; m++)
            {
                for (int n = 0; n < 5; n++)
                {
                    Ai(m, n) = wow1.allA[i][m][n];
                }
            }
            MatrixXd Bi;
            Bi = Hi * Ai;
            for (int n = 0; n < 5; n++)
            {
                B(i, n) = Bi(0, n);
            }
        }
        /*  vector<vector<double>>showB;
          for (int m = 0; m < 4; m++)
          {
              for (int n = 0; n < 5; n++)
              {

              }
          }*/
        BTB = B.transpose() * B;
        BTl = B.transpose() * l;
        vector<double>showl;
        for (int j = 0; j < 60000; j++)
        {
            showl.push_back(l(j, 0));
        }
        MatrixXd dx;
        dx = BTB.inverse() * BTl;
        double rj0 = sqrt(pow(y10,2) + pow(y20,2));
        y10 = y10 + dx(0, 0);
        y20 = y20 + dx(1, 0);
        y30 = y30 + dx(2, 0);
        y40 = y40 + dx(3, 0);
        b0 = b0 + dx(4, 0);
        rb = fabs(sqrt(pow(y10,2) + pow(y20,2)) - rj0);
        cout << setprecision(10);
        cout << "y1  " << y10 << endl;
        cout << "y2  " << y20 << endl;
        cout << "y3  " << y30 << endl;
        cout << "y4  " << y40 << endl;
        cout<< "b  " << b0 << endl;
        cout << "������ҵҪ�󣬼�ʹ�����޲�Ҫ��ҲҪ���ٵ���һ��" << endl;
    }
    return 1;
}