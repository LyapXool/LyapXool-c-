using System;

namespace LyapXool
{
    class Program
    {

        static void Main(string[] args)
        {
            DateTime start = DateTime.Now;

            RBF rbf = new RBF(Instructions.ode_dimension, Instructions.c, Instructions.maxmax, Instructions.minmin, Instructions.alpha, Instructions.points_directional, Instructions.radius, Instructions.cart_grid_density, Instructions.min_geometric_limits, Instructions.max_geometric_limits, Instructions.normal, Instructions.printing);
            Generalities gnl = new Generalities();
            Lyapunov lpv = new Lyapunov(Instructions.ode_dimension, Instructions.c, Instructions.points_directional, Instructions.critval, Instructions.normal, Instructions.printing);
            int totaliterations = Instructions.totaliterations;

            gnl.ExecutionDate(true);

            /* PRIMERO HAZ LA MALLA DE COLOCACIÓN */
            rbf.MakeColGrid();

            /* LUEGO HAZ LA MALLA DIRECCIONAL */
            rbf.DirecGrid();

            /* DESPUÉS LA MALLA CARTESIANA */
            rbf.MakeEvalGrid(ref rbf.cartesianevalgrid);

            /* ENTONCES, COMIENZA CON LA FUNCIÓN DE LIAPUNOV */
            rbf.GetAlphaColMatrix();
            lpv.GetInverseForLyapunov(ref rbf.Amat);

            for (int iteration = 1; iteration <= totaliterations; ++iteration)
            {
                DateTime startite = DateTime.Now;
                Console.WriteLine("\t \t \t La iteración actual es: " + iteration);
                Instructions.woutput.WriteLine("\t \t La iteración actual es: " + iteration);
                lpv.LyapEquation(iteration, rbf);

                lpv.LyapunovFunctions(iteration, true, ref rbf.directgrid, rbf);
                lpv.ChainRecurrentSet(iteration, true, true, ref rbf.directgrid);
                lpv.GetNewAlpha(iteration, rbf);

                lpv.LyapunovFunctions(iteration, false, ref rbf.cartesianevalgrid, rbf);
                lpv.ChainRecurrentSet(iteration, false, true, ref rbf.cartesianevalgrid);
                lpv.FirstDerivative(iteration, false, ref rbf.cartesianevalgrid, rbf);
                DateTime endite = DateTime.Now;
                TimeSpan tsite = (endite - startite);

                Console.WriteLine("=====El tiempo total de la iteración " + iteration + " fue {0:00} días, {1:00} horas, {2:00} minutos, {3:00}.{4} segundos=====", tsite.Days, tsite.Hours, tsite.Minutes, tsite.Seconds, tsite.Milliseconds);
                Instructions.woutput.WriteLine("=====El tiempo total de la iteración " + iteration + " fue {0:00} días, {1:00} horas, {2:00} minutos, {3:00}.{4} segundos=====", tsite.Days, tsite.Hours, tsite.Minutes, tsite.Seconds, tsite.Milliseconds);
                Console.WriteLine("\t \t ============================================");
                Instructions.woutput.WriteLine("\t \t ============================================");
            }
            gnl.ExecutionDate(false);
            DateTime end = DateTime.Now;
            TimeSpan ts = (end - start);
            Console.WriteLine("=====El tiempo total de ejecución de este programa fue {0:00} días, {1:00} horas, {2:00} minutos, {3:00}.{4} segundos=====", ts.Days, ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====El tiempo total de ejecución de este programa fue {0:00} días, {1:00} horas, {2:00} minutos, {3:00}.{4} segundos=====", ts.Days, ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.Close();
        }
    }
}


