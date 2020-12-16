using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab5MMDO
{
    delegate double Function(double[] x);
    class Program
    {
        static void Main(string[] args)
        {
            int n = 2;
            double[] x0 = new double[2] { -2, 1 };
            Function f1 = (x) => 5 * x[0] * x[0] + x[0] * x[1] + 25 * x[1] * x[1] - 4 * x[0] + 6 * x[1];
            double e = Math.Pow(10, -4);
            double[] ans = MethodLibrary.ReiskiMehod(n, f1, new double[2] { -2, 1 }, e * 100, Math.Pow(10, -4));
            Console.WriteLine(ans[0] + " " + ans[1]);
            ans = MethodLibrary.PokoordReiski(n, new double[2] { -2, 1 }, e * 100, e * 1000, e, f1);
            Console.WriteLine(ans[0] + " " + ans[1]);
            ans = MethodLibrary.Newton(n, f1, new double[2] { -2, 1 }, Math.Pow(10, -4));
            Console.WriteLine(ans[0] + " " + ans[1]);
            ans = MethodLibrary.Newton2(n, f1, new double[2] { -2, 1 }, Math.Pow(10, -4), e * 100);
            Console.WriteLine(ans[0] + " " + ans[1]);
            Console.ReadKey();
        }
    }
}
