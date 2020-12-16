using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab5MMDO
{
    class MethodLibrary
    {
        public static int nI;

        public static double[] TheFastestDownhill(int n, Function F, double[] x0, double h0, double e)
        {
            nI = 0;
            double[] x = new double[n];
            double[] g = GradF(x0, F);
            double[] x_x0 = new double[n];
            if (Norm(g) > e)
            {
                do
                {
                    nI++;
                    Array.Copy(x0, x, n);
                    double h = FindH(n, F, x0, g, h0, e);
                    for (int i = 0; i < n; i++)
                    {
                        x0[i] = x[i] - h * g[i];
                    }
                    g = GradF(x0, F);
                    //Console.WriteLine("\t" + x0[0] + " " + x0[1]);
                    for (int i = 0; i < n; i++)
                    {
                        x_x0[i] = x[i] - x0[i];
                    }
                }
                while (!(Norm(x_x0) < e || Norm(g) < e));

            }
            return x0;
        }

        public static double[] Pokoord(int n, Function F, double[] x0, double h0, double l, double e)
        {
            nI = 0;
            double[] h = new double[n];
            for (int i = 0; i < n; i++)
            {
                h[i] = h0;
            }
            double[] x_int = new double[n];
            double[] x_ext = new double[n];
            Array.Copy(x0, x_int, n);
            double[] i_e = new double[n];
            do
            {
                nI++;
                Array.Copy(x_int, x_ext, n);
                for (int i = 0; i < n; i++)
                {
                    double[] x = new double[n];
                    Array.Copy(x_int, x, n);
                    double fx = F(x);
                    double[] y1 = new double[n];
                    Array.Copy(x, y1, n);
                    y1[i] += 3 * e;
                    double[] y2 = new double[n];
                    Array.Copy(x, y2, n);
                    y2[i] -= 3 * e;
                    double f1 = F(y1);
                    double f2 = F(y2);
                    double z = Math.Sign(f2 - f1);
                    double fx1 = 0;
                    do
                    {
                        x_int[i] = x[i] + h[i] * z;
                        fx1 = F(x_int);
                        if (fx1 >= fx)
                            h[i] *= l;
                    }
                    while (fx1 >= fx && h[i] >= e / 2);
                }
                for (int i = 0; i < n; i++)
                {
                    i_e[i] = x_int[i] - x_ext[i];
                }
            }
            while (Norm(i_e) >= e);
            return x_int;
        }

        public static double[] NewtonClassic(int n, Function F, double[] x0, double e)
        {
            nI = 0;
            double[] g = GradF(x0, F);
            if (Norm(g) > e)
            {
                do
                {
                    nI++;
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

        public static double[] NewtonGeneral(int n, Function F, double[] x0, double h0, double e)
        {
            nI = 0;
            double[] g = GradF(x0, F);
            double h = 0;
            if (Norm(g) > e)
            {
                do
                {
                    nI++;
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
                    }
                    h = FindH(n, F, x0, w, h0, e);
                    for (int i = 0; i < n; i++)
                    {
                        x0[i] = x[i] - h * w[i];
                    }
                    g = GradF(x0, F);
                }
                while (Norm(g) >= e && h >= e);
            }
            return x0;
        }

        private static double FindH(int n, Function F, double[] x0_, double[] g_, double h0, double e)
        {
            double[] x0 = new double[n];
            Array.Copy(x0_, x0, n);
            double[] g = new double[n];
            Array.Copy(g_, g, n);
            double h = 0;
            double f1 = F(x0);
            double f2 = 0;
            double[] x1 = new double[n];
            double[] x2 = new double[n];
            do
            {
                h0 = h0 / 2;
                for (int i = 0; i < n; i++)
                {
                    x2[i] = x0[i] - h0 * g[i];
                }
                f2 = F(x2);
            }
            while (!(f1 > f2 || h0 < e));
            if (h0 > e)
            {
                do
                {
                    Array.Copy(x2, x1, n);
                    f1 = f2;
                    h = h + h0;
                    for (int i = 0; i < n; i++)
                    {
                        x2[i] = x1[i] - h * g[i];
                    }
                    f2 = F(x2);
                }
                while (!(f1 < f2));
                double ha = h - 2 * h0;
                double hb = h;
                double q = e / 3;
                do
                {
                    double h1 = (ha + hb - q) / 2;
                    double h2 = (ha + hb + q) / 2;
                    for (int i = 0; i < n; i++)
                    {
                        x1[i] = x0[i] - h1 * g[i];
                        x2[i] = x0[i] - h2 * g[i];
                    }
                    f1 = F(x1);
                    f2 = F(x2);
                    if (f1 <= f2)
                    {
                        hb = h2;
                    }
                    else
                    {
                        ha = h1;
                    }
                }
                while (!(hb - ha < e));
                h = (ha + hb) / 2;
            }
            return h;
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

        private static Function FDiff(int varInd, Function F_)
        {
            Function FD = delegate(double[] x)
            {
                int n = x.Length;
                double[] inp = new double[n];
                Array.Copy(x, inp, n);
                double f = F_(inp);
                double dX = 0.0000001;
                inp[varInd] += dX;
                double f2 = F_(inp);
                return (f2 - f) / dX;
            };
            return FD;
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

        private static double Norm(double[] mas)
        {
            double kvSum = 0;
            for (int i = 0; i < mas.Length; i++)
            {
                kvSum += mas[i] * mas[i];
            }
            return Math.Sqrt(kvSum);
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
