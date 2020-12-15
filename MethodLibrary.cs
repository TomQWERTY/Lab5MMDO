using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab5MMDO
{
    class MethodLibrary
    {
        public static double[] ReiskiMehod(int n, Function F, double[] x0, double e)
        {
            double[] x = new double[2];
            double[] g = GradF(x0, F);
            if (DoblMod(g) < e)
            {
                do
                {
                    Array.Copy(x0, x, x0.Length);
                    double h = FindH(e, x0, F);
                    for (int i = 0; i < n; i++)
                    {
                        x0[i] = x[i] - h * g[i]; 
                    }
                    g = GradF(x0, F);
                }
                while (DoblMod(x) - DoblMod(x0) >= e || DoblMod(g) >= e);
                
            }
            return x0;
        }

        private static double[] GradF(double[] x, Function F)
        {
            return new double[2] {x[0] / F(x), x[1] / F(x)};
        }

        private static double DoblMod(double[] x)
        {
            return Math.Sqrt(x[0] * x[0] + x[1] * x[1]);
        }

        private static double FindH(double e, double[] start, Function F)
        {
            double h = e * 10;
            double[] x0 = new double[2] {start[0] + h, start[1] + h};
            double[] x1 = new double[2];
            double[] x2 = new double[2];
            double f1 = F(x0);
            double f2;
            do
            {
                h /= 2;
                x2 = new double[2] { start[0] + h, start[1] + h };
                f2 = F(x2);
                if (f1 <= f2)
                {
                    h = -h;
                    x2 = new double[2] { start[0] + h, start[1] + h };
                    f2 = F(x2);
                }                                                                                                 
            }
            while (f1 <= f2 && Math.Abs(h) >= e);

            return h;
        }
    }
}
