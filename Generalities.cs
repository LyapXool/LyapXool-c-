using System;
using System.IO;

namespace LyapXool
{
    class Generalities
    {
        public Generalities()
        {

        }
        public void PrintMatrix(string nameVariable, int interationcount, double[,] matrixToPrint)
        {
            Console.WriteLine(matrixToPrint.Length);

            int minRows = matrixToPrint.GetLength(0);
            int minCols = matrixToPrint.GetLength(1);
            Console.Write(nameVariable + "ite" + interationcount + " = [ ");
            for (int i = 0; i < minRows; ++i)
            {
                for (int j = 0; j < minCols; ++j)
                {
                    Console.Write(matrixToPrint[i, j] + " ");
                }
                Console.Write("\n");
            }
            Console.WriteLine("];");
        }
        public void PrintVector(string nameVariable, int interationcount, double[] vectorToPrint)
        {
            int forLimit = vectorToPrint.Length;
            Console.Write(nameVariable + "ite" + interationcount + " = [ ");
            for (int i = 0; i < forLimit; ++i)
            {
                Console.Write(vectorToPrint[i] + " ");
            }
            Console.WriteLine("];");
        }
        public void PrintColumnsToFile(string nameVariable, int interationcount, double[,] matrixToPrint)
        {
            int minRows = matrixToPrint.GetLength(0);
            int minCols = matrixToPrint.GetLength(1);

            string fileName1 = "s" + nameVariable + "ite" + interationcount + "x.m";
            FileStream fs1 = new FileStream(fileName1, FileMode.OpenOrCreate, FileAccess.Write);
            fs1.Close();
            StreamWriter sw1 = new StreamWriter(fileName1);

            sw1.Write(nameVariable + "ite" + interationcount + "x = [ ");
            for (int i = 0; i < minRows; ++i)
            {
                sw1.Write(matrixToPrint[i, 0] + " ");
            }
            sw1.Write("];");
            sw1.Close();

            string fileName2 = "s" + nameVariable + "ite" + interationcount + "y.m";
            FileStream fs2 = new FileStream(fileName2, FileMode.OpenOrCreate, FileAccess.Write);
            fs2.Close();
            StreamWriter sw2 = new StreamWriter(fileName2);

            sw2.Write(nameVariable + "ite" + interationcount + "y = [ ");
            for (int i = 0; i < minRows; ++i)
            {
                sw2.Write(matrixToPrint[i, 1] + " ");
            }
            sw2.Write("];");
            sw2.Close();
        }
        public void PrintColumnsToFile(string nameVariable, double[,] matrixToPrint)
        {
            int minRows = matrixToPrint.GetLength(0);
            int minCols = matrixToPrint.GetLength(1);

            string fileName1 = "s" + nameVariable + "ite" + "x.m";
            FileStream fs1 = new FileStream(fileName1, FileMode.OpenOrCreate, FileAccess.Write);
            fs1.Close();
            StreamWriter sw1 = new StreamWriter(fileName1);

            sw1.Write(nameVariable + "x = [ ");
            for (int i = 0; i < minRows; ++i)
            {
                sw1.Write(matrixToPrint[i, 0] + " ");
            }
            sw1.Write("];");
            sw1.Close();

            string fileName2 = "s" + nameVariable + "ite" + "y.m";
            FileStream fs2 = new FileStream(fileName2, FileMode.OpenOrCreate, FileAccess.Write);
            fs2.Close();
            StreamWriter sw2 = new StreamWriter(fileName2);

            sw2.Write(nameVariable + "y = [ ");
            for (int i = 0; i < minRows; ++i)
            {
                sw2.Write(matrixToPrint[i, 1] + " ");
            }
            sw2.Write("];");
            sw2.Close();
        }
        public void PrintVectorToFile(string nameVariable, int interationcount, ref double[] vectorToPrint)
        {
            int minRows = vectorToPrint.Length;
            string fileName1 = "s" + nameVariable + "ite" + interationcount + ".m";
            FileStream fs1 = new FileStream(fileName1, FileMode.OpenOrCreate, FileAccess.Write);
            fs1.Close();
            StreamWriter sw1 = new StreamWriter(fileName1);

            sw1.Write(nameVariable + "ite" + interationcount + " = [ ");
            for (int i = 0; i < minRows; ++i)
            {
                sw1.Write(vectorToPrint[i] + " ");
            }
            sw1.Write("];");
            sw1.Close();
        }
        public void ExecutionDate(bool startIt)
        {
            string[] meses = { "enero", "febrero", "marzo", "abril", "mayo", "junio", "julio", "agosto", "septiembre", "octubre", "noviembre", "diciembre" };
            string[] dias = { "lunes", "martes", "miércoles", "jueves", "viernes", "sábado", "domingo" };
            DateTime now = DateTime.Now;
            if (startIt)
            {
                PrintGnlInfo();
                Console.WriteLine("El cálculo de la función de Liapunov comienzó el día {0} {1} de {2} de {3} a la hora {4} (UTC: {5} )", dias[(int)now.DayOfWeek], now.Day, meses[now.Month - 1], now.Year, now.TimeOfDay, DateTime.UtcNow);
                Instructions.woutput.WriteLine("El cálculo de la función de Liapunov comienzó el día {0} {1} de {2} de {3} a la hora {4} (UTC: {5} )", dias[(int)now.DayOfWeek], now.Day, meses[now.Month - 1], now.Year, now.TimeOfDay, DateTime.UtcNow);
            }
            else
            {
                Console.WriteLine("LyapXool# has been properly executed and called {0} times the differential function", Instructions.functionodecalls);
                Console.WriteLine("LyapXool# terminó correctamente y llamó a la función diferencial un total de {0} veces", Instructions.functionodecalls);
                Console.WriteLine("El cálculo de la función de Liapunov terminó el día {0} {1} de {2} de {3} a la hora {4} (UTC: {5} )", dias[(int)now.DayOfWeek], now.Day, meses[now.Month - 1], now.Year, now.TimeOfDay, DateTime.UtcNow);
                Instructions.woutput.WriteLine("LyapXool# has been properly executed and called {0} times the differential function", Instructions.functionodecalls);
                Instructions.woutput.WriteLine("LyapXool# terminó correctamente y llamó a la función diferencial un total de {0} veces", Instructions.functionodecalls);
                Instructions.woutput.WriteLine("El cálculo de la función de Liapunov terminó el día {0} {1} de {2} de {3} a la hora {4} (UTC: {5} )", dias[(int)now.DayOfWeek], now.Day, meses[now.Month - 1], now.Year, now.TimeOfDay, DateTime.UtcNow);
            }
        }
        public void PrintGnlInfo()
        {
            string fileName1 = "data.lpx";
            FileStream fs1 = new FileStream(fileName1, FileMode.OpenOrCreate, FileAccess.Write);
            fs1.Close();
            StreamWriter sw1 = new StreamWriter(fileName1);

            sw1.Write("The parameters used to produced these computations are the following: \n");
            sw1.Write("\n");
            sw1.Write("ode_dimension: " + Instructions.ode_dimension + "\n");
            sw1.Write("c: " + Instructions.c + "\n");
            sw1.Write("maxmax: " + Instructions.maxmax + "\n");
            sw1.Write("minmin: " + Instructions.minmin + "\n");
            sw1.Write("alpha: " + Instructions.alpha + "\n");
            sw1.Write("normal: " + Instructions.normal + "\n");
            sw1.Write("cart_grid_density: " + Instructions.cart_grid_density + "\n");
            sw1.Write("radius: " + Instructions.radius + "\n");
            sw1.Write("critval: " + Instructions.critval + "\n");
            sw1.Write("min_geometric_limits: { ");
            for (int i = 0; i < Instructions.ode_dimension; ++i)
            {
                sw1.Write(Instructions.min_geometric_limits[i] + " ");
            }
            sw1.Write("}" + "\n");
            sw1.Write("max_geometric_limits: { ");
            for (int i = 0; i < Instructions.ode_dimension; ++i)
            {
                sw1.Write(Instructions.max_geometric_limits[i] + " ");
            }
            sw1.Write("}" + "\n");
            sw1.Write("points_directional: " + Instructions.points_directional + "\n");
            sw1.Write("printing: " + Instructions.printing + "\n");
            sw1.Write("\n");
            sw1.Write("\n");
            sw1.Write("\n");
            sw1.Close();
            PrintDisc();
        }
        public void PrintDisc()
        {
            string fileName1 = "gnldisc.lpx";
            FileStream fs1 = new FileStream(fileName1, FileMode.OpenOrCreate, FileAccess.Write);
            fs1.Close();
            StreamWriter sw1 = new StreamWriter(fileName1);
            sw1.Write("LyapXool# has been written by Carlos Argáez from the Marine and Freshwater Research Institute, Iceland.\n" +
                "The code is primarily based on C. Argáez's own code: LyapXool, Argáez et al. DOI: 10.1016/j.softx.2020.100616\n" +
                "However, it does use a great deal of public sources: \n" +
                "\t 1) McCaffrey, James D. Inverting a Matrix using C#.  2015.\n" +
                "\t https://jamesmccaffrey.wordpress.com/2015/03/06/inverting-a-matrix-using-c/ \n" +
                "\t 2) How to resize multidimensional (2D) array in C#? \n" +
                "\t Particular answer provided by Manuel.\n" +
                "\t https://stackoverflow.com/questions/6539571/how-to-resize-multidimensional-2d-array-in-c \n" +
                "\t 3) Cholesky Decomposition : Matrix Decomposition, GeeksforGeeks.\n" +
                "\t https://www.geeksforgeeks.org/cholesky-decomposition-matrix-decomposition/\n" +
                "\t 4) is there in C# a method for List<T> like resize in c++ for vector<T>. \n" +
                "\t Particular answer provided by Jon Hanna.\n" +
                "\t https://stackoverflow.com/questions/12231569/is-there-in-c-sharp-a-method-for-listt-like-resize-in-c-for-vectort \n" +
                "\t 5) ZetCode, https://zetcode.com/csharp/datetime/");
            sw1.Write("\n");
            sw1.Write("\n");
            sw1.Write("\n");
            sw1.Write("Everybody can download, use and modify this code.");
            sw1.Close();
        }
    }
}
