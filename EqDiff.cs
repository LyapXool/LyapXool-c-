using System;
namespace LyapXool
{
    public class EqDiff
    {
        static public void eqdff(bool normal, double[] x, double[] f)
        {
            f[0] = -1.0 * x[0] * ((x[0] * x[0]) + (x[1] * x[1]) - (1.0 / 4.0)) * ((x[0] * x[0]) + (x[1] * x[1]) - 1.0) - x[1];
            f[1] = -1.0 * x[1] * ((x[0] * x[0]) + (x[1] * x[1]) - (1.0 / 4.0)) * ((x[0] * x[0]) + (x[1] * x[1]) - 1.0) + x[0];
            int end = x.Length;
            if (normal)
            {
                double norma = Math.Sqrt(ArrayOperations.Dot(f, f));
                for (int i = 0; i < end; ++i)
                    f[i] /= norma;

            }
            ++Instructions.functionodecalls;
        }
    }
}
