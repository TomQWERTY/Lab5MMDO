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
            if (Norm(g) < e)
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
                while (Norm(x) - Norm(x0) >= e || Norm(g) >= e);

            }
            return x0;
        }

        private static double[] GradF(double[] x, Function F)
        {
            int n = x.Length;
            double[] grad = new double[n];
            for (int i = 0; i < n; i++)
            {
                grad[i] = FDiff(i, F)(x);
            }
            return grad;
        }

        private static double Norm(double[] mas)
        {
            double kvSum = 0;
            for (int i = 0; i < mas.Length; i++)
            {
                kvSum += mas[i] * mas[i];
            }
            return Math.Sqrt(kvSum);
        }

        private static double FindH(double e, double[] start, Function F)
        {
            double h = e * 10;
            double[] x0 = new double[2] { start[0] + h, start[1] + h };
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

        public static double[] Newton(int n, Function F, double[] x0, double e)
        {
            double[] g = GradF(x0, F);
            if (Norm(g) > e)
            {
                do
                {
                    double[] x = new double[n];
                    Array.Copy(x0, x, n);
                    double[,] g2 = GessF(x0, F);
                    double[,] g2_1 = Inverse(g2);
                    double[] w = new double[n];
                    for (int i = 0; i < n; i++)
                    {
                        w[i] = 0;
                        for (int j = 0; j < n; j++)
                        {
                            w[i] += g2_1[i, j] * g[j];
                        }
                        x0[i] = x[i] - w[i];
                    }
                    g = GradF(x0, F);
                }
                while (Norm(g) >= e);
            }
            return x0;
        }

        private static double[,] GessF(double[] x, Function F)
        {
            int n = x.Length;
            double[,] gess = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    gess[i, j] = FDiff(j, FDiff(i, F))(x);
                }
            }
            return gess;
        }

        private static Function FDiff(int varInd, Function F_)
        {          
            Function FD = delegate(double[] x)
            {
                int n = x.Length;
                double[] inp = new double[n];
                Array.Copy(x, inp, n);
                double f = F_(inp);
                double dX = 0.0001;
                inp[varInd] += dX;
                double f2 = F_(inp);
                return (f2 - f) / dX;
            };
            return FD;
        }

        private static double[,] Inverse(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            double[,] a1 = new double[n, n];
            Array.Copy(matrix, a1, n * n);
            double[,] b1 = new double[n, n], a2 = new double[n, n], b2 = new double[n, n];
            for (int k = 0; k < n; k++)
            {
                a2[k, k] = 1;
            }
            for (int k = 0; k < n; k++)
            {
                int r = k;
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        b1[i, j] = (a1[r, r] * a1[i, j] - a1[r, j] * a1[i, r]) / a1[r, r];
                        b2[i, j] = (a1[r, r] * a2[i, j] - a2[r, j] * a1[i, r]) / a1[r, r];
                    }
                }
                for (int j = 0; j < n; j++)
                {
                    b1[r, j] = a1[r, j] / a1[r, r];
                    b2[r, j] = a2[r, j] / a1[r, r];
                }
                double[,] temp = a1;
                a1 = b1;
                b1 = temp;
                temp = a2;
                a2 = b2;
                b2 = temp;
            }
            return a2;
        }
    }
}
