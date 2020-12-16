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
        static int nF;

        static void Main(string[] args)
        {
            Function f1 = delegate(double[] x)
            {
                nF++;
                return
                    //5 * x[0] * x[0] + x[0] * x[1] + 25 * x[1] * x[1] - 4 * x[0] + 6 * x[1];
                    x[0] * x[0] + Math.Exp(2 * x[0] * x[0] + x[1] * x[1]) + x[0] - x[1];
            };
            int n = 2;
            //double x1 = -2, x2 = 1;
            double x1 = 1, x2 = -1;
            double e = Math.Pow(10, -8);

            nF = 0;
            Console.WriteLine("Метод найшвидшого спуску:");
            DateTime time = DateTime.Now;
            double[] ans = MethodLibrary.TheFastestDownhill(n, f1, new double[2] { x1, x2 }, 10, e);
            long t = (DateTime.Now - time).Ticks;
            Console.WriteLine("Точка мiнiмуму: " + ans[0] + ", " + ans[1] + "; Мiнiмум: " + f1(ans) + "; Кiлькiсть iтерацiй: " + MethodLibrary.nI + "; " + "Кiлькiсть обчислень функцiї: " + (nF - 1) + "; Системний час: " + t + ";");

            nF = 0;
            Console.WriteLine("\nМетод покоординатного спуску:");
            time = DateTime.Now;
            ans = MethodLibrary.Pokoord(n, f1, new double[2] { x1, x2 }, 10, 0.1, e);
            t = (DateTime.Now - time).Ticks;
            Console.WriteLine("Точка мiнiмуму: " + ans[0] + ", " + ans[1] + "; Мiнiмум: " + f1(ans) + "; Кiлькiсть iтерацiй: " + MethodLibrary.nI + "; " + "Кiлькiсть обчислень функцiї: " + (nF - 1) + "; Системний час: " + t + ";");

            nF = 0;
            Console.WriteLine("\nКласичний метод Ньютона:");
            time = DateTime.Now;
            ans = MethodLibrary.NewtonClassic(n, f1, new double[2] { x1, x2 }, e);
            t = (DateTime.Now - time).Ticks;
            Console.WriteLine("Точка мiнiмуму: " + ans[0] + ", " + ans[1] + "; Мiнiмум: " + f1(ans) + "; Кiлькiсть iтерацiй: " + MethodLibrary.nI + "; " + "Кiлькiсть обчислень функцiї: " + (nF - 1) + "; Системний час: " + t + ";");

            nF = 0;
            Console.WriteLine("\nУзагальнений метод Ньютона:");
            time = DateTime.Now;
            ans = MethodLibrary.NewtonGeneral(n, f1, new double[2] { x1, x2 }, 10, e);
            t = (DateTime.Now - time).Ticks;
            Console.WriteLine("Точка мiнiмуму: " + ans[0] + ", " + ans[1] + "; Мiнiмум: " + f1(ans) + "; Кiлькiсть iтерацiй: " + MethodLibrary.nI + "; " + "Кiлькiсть обчислень функцiї: " + (nF - 1) + "; Системний час: " + t + ".");
            Console.ReadKey();
        }
    }
}
