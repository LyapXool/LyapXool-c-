using System.IO;

namespace LyapXool
{
    public class Instructions
    {
        public const int ode_dimension = 2;
        public const int c = 1;
        public const int maxmax = 450;
        public const int minmin = -maxmax;
        public const double alpha = 0.09;
        public const bool normal = true;
        public const double cart_grid_density = 0.007;
        public const double radius = 0.49;
        public const double critval = -0.5;
        static public double[] min_geometric_limits = { -1.6, -1.6 };
        static public double[] max_geometric_limits = { 1.6, 1.6 };
        public const int points_directional = 10;
        public const bool printing = true;
        public const int totaliterations = 1;

        /* NO MODIFICAR LAS LÍNEAS SIGUIENTES */
        public static string outputf = "salida.lpx";
        static public StreamWriter woutput = new StreamWriter(outputf);
        static public ulong functionodecalls = 0;
    }
}
