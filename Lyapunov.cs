using System;
using System.Collections.Generic;

namespace LyapXool
{
    class Lyapunov
    {
        Generalities gnl = new Generalities();
        public Lyapunov(int Iode_dimension, int Ic, int Ipoints_directional, double Icritval, bool Inormal, bool Iprinting)
        {
            ode_dimension = Iode_dimension;
            c = Ic;
            points_directional = Ipoints_directional;
            critval = Icritval;
            normal = Inormal;
            printing = Iprinting;
        }
        public void GetInverseForLyapunov(ref double[,] matrixToInvert)
        {
            DateTime start = DateTime.Now;
            int maxbet = (int)matrixToInvert.GetLength(0);
            Array.Resize(ref betaod, maxbet);


            double[][] Amattemp = new double[maxbet][];
            double[][] sds = new double[maxbet][];
            double[,] sss = new double[maxbet, maxbet];
            for (int i = 0; i < maxbet; ++i)
            {
                Amattemp[i] = new double[maxbet];
                sds[i] = new double[maxbet];

            }
            for (int i = 0; i < maxbet; ++i)
            {
                for (int j = 0; j < maxbet; ++j)
                {
                    Amattemp[i][j] = matrixToInvert[i, j];
                }
            }
            Amatperm = ArrayOperations.MatrixInverse(Amattemp);
            DateTime end = DateTime.Now;
            TimeSpan ts = (end - start);
            Console.WriteLine("=====Invertir la matriz tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====Invertir la matriz tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
        }
        public void LyapEquation(int iteration, RBF rbf)
        {
            DateTime start = DateTime.Now;
            int maxbet = (int)rbf.collocationpoints.GetLength(0);
            Array.Resize(ref betaod, maxbet);
            betaod = ArrayOperations.MatrixVectorProduct(ArrayOperations.MatrixInverse(Amatperm), rbf.alphavector);
            gnl.PrintVectorToFile("vecbet", iteration, ref betaod);
            DateTime end = DateTime.Now;
            TimeSpan ts = (end - start);
            Console.WriteLine("=====Resolver la ecuación de Liapunov tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====Resolver la ecuación de Liapunov tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
        }

        public void LyapunovFunctions(int iteration, bool type_of_grid, ref double[,] evalcoordinates, RBF rbf)
        {
            DateTime start = DateTime.Now;
            int i = 0, k = 0;
            int maxite = (int)evalcoordinates.GetLength(0);
            int maxbet = (int)rbf.collocationpoints.GetLength(0);
            int pointdim = (int)rbf.collocationpoints.GetLength(1);
            {
                Array.Resize(ref lyapfunc, maxite);
                Array.Resize(ref orbder, maxite);
                Array.Fill(lyapfunc, 0.0);
                Array.Fill(orbder, 0.0);

                double[] diffpoints = new double[pointdim];
                double[] diffpointski = new double[pointdim];
                double[] diffpointskineg = new double[pointdim];
                double[] resulti = new double[pointdim];
                double[] resultk = new double[pointdim];
                double[] saving = new double[pointdim];
                double[] savingdomain = new double[pointdim];

                Array.Fill(diffpoints, 0.0);
                Array.Fill(diffpointski, 0.0);
                Array.Fill(diffpointskineg, 0.0);
                Array.Fill(resulti, 0.0);
                Array.Fill(resultk, 0.0);
                Array.Fill(saving, 0.0);
                Array.Fill(savingdomain, 0.0);

                double proctk = 0.0;
                double producting = 0.0;
                double twopointsdistance = 0.0;
                double wdlfvalue1 = 0.0;
                double wdlfvalue2 = 0.0;
                double checking = 0.0;

                for (i = 0; i < maxite; ++i)
                {
                    for (int all = 0; all < pointdim; ++all)
                        savingdomain[all] = evalcoordinates[i, all];

                    EqDiff.eqdff(normal, savingdomain, resulti);
                    for (k = 0; k < maxbet; ++k)
                    {
                        for (int all = 0; all < pointdim; ++all)
                            saving[all] = rbf.collocationpoints[k, all];

                        diffpoints = ArrayOperations.VDiff(savingdomain, saving);
                        twopointsdistance = Math.Sqrt(ArrayOperations.Dot(diffpoints, diffpoints));
                        checking = 1.0 - c * twopointsdistance;
                        if (checking > 0.0)
                        {
                            EqDiff.eqdff(normal, saving, resultk);
                            wdlfvalue1 = Wendland.WndlndFnctnFirst(twopointsdistance, c);
                            wdlfvalue2 = Wendland.WndlndFnctnSecond(twopointsdistance, c);
                            diffpointski = ArrayOperations.VDiff(saving, savingdomain);
                            proctk = ArrayOperations.Dot(diffpointski, resultk);
                            producting = betaod[k] * proctk;
                            lyapfunc[i] += producting * wdlfvalue1;
                            orbder[i] += -wdlfvalue2 * producting * ArrayOperations.Dot(diffpointski, resulti) - betaod[k] * wdlfvalue1 * ArrayOperations.Dot(resulti, resultk);
                        }
                    }
                }
                if (printing)
                {
                    if (type_of_grid)
                    {
                        gnl.PrintVectorToFile("lyapfuncd", iteration, ref lyapfunc);
                        gnl.PrintVectorToFile("orbderd", iteration, ref orbder);
                    }
                    else
                    {
                        gnl.PrintVectorToFile("lyapfuncc", iteration, ref lyapfunc);
                        gnl.PrintVectorToFile("orbderc", iteration, ref orbder);
                    }
                }
            }
            DateTime end = DateTime.Now;
            TimeSpan ts = (end - start);
            Console.WriteLine("=====Obtener la función de Liapunov tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====Obtener la función de Liapunov tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
        }
        public void ChainRecurrentSet(int currentiteration, bool type_of_grid, bool with_orbder, ref double[,] evalcoordinates)
        {
            DateTime start = DateTime.Now;
            List<int> counterzero = new List<int>();
            int maxlength = (int)evalcoordinates.GetLength(0);
            int maxwidth = (int)evalcoordinates.GetLength(1);
            for (int j = 0; j < maxlength; ++j)
            {
                if (with_orbder)
                {
                    if (orbder[j] > critval)
                    {
                        counterzero.Add(j);
                    }
                }
                else
                {
                    if (-normed[j] > -critval)
                    {
                        counterzero.Add(j);
                    }
                }
            }
            int faillength = (int)counterzero.Count;

            double[] crslyapun = new double[faillength];
            double[] crsorbder = new double[faillength];
            double[,] failinggrid = new double[faillength, maxwidth];
            double[] failinglyapunov = new double[faillength];
            double[] failingorbder = new double[faillength];
            int end = counterzero.Count;
            int m = 0;
            {
                for (int ii = 0; ii < end; ++ii)
                {
                    crslyapun[m] = lyapfunc[counterzero[ii]];
                    crsorbder[m] = orbder[counterzero[ii]];
                    failinglyapunov[m] = lyapfunc[counterzero[ii]];
                    failingorbder[m] = orbder[counterzero[ii]];

                    for (int All = 0; All < maxwidth; ++All)
                        failinggrid[m, All] = evalcoordinates[counterzero[ii], All];

                    m++;
                }
            }
            counterzero.Clear();
            if (printing)
            {
                if (type_of_grid)
                {
                    gnl.PrintColumnsToFile("fdirecgrid", currentiteration, failinggrid);
                    gnl.PrintVectorToFile("flfdirecgrid", currentiteration, ref failinglyapunov);
                    gnl.PrintVectorToFile("flfpdirecgrid", currentiteration, ref failingorbder);
                }
                else
                {
                    gnl.PrintColumnsToFile("fcartesian", currentiteration, failinggrid);
                    gnl.PrintVectorToFile("flfcartesian", currentiteration, ref failinglyapunov);
                    gnl.PrintVectorToFile("flfpcartesian", currentiteration, ref failingorbder);
                }
            }
            DateTime end1 = DateTime.Now;
            TimeSpan ts = (end1 - start);
            Console.WriteLine("=====Obtener la recurrencia tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====Obtener la recurrencia tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
        }
        public void FirstDerivative(int currentiteration, bool type_of_grid, ref double[,] evalcoordinates, RBF rbf)
        {
            DateTime start = DateTime.Now;
            int i = 0, j = 0, k = 0;
            int evaldim = (int)evalcoordinates.GetLength(0);
            ArrayOperations.ResizeArray<double>(ref fdvector, evaldim, ode_dimension);
            ArrayOperations.FillMatrixConst(ref fdvector, 0.0);

            double checking = 0.0;
            double twopointsdistance = 0.0;
            double wdlfvalue1 = 0.0;
            double wdlfvalue2 = 0.0;
            {
                double[] saving = new double[ode_dimension];
                double[] savingdomain = new double[ode_dimension];
                double[] diffpoints = new double[ode_dimension];
                double[] resultk = new double[ode_dimension];

                int maxite = (int)betaod.Length;

                for (j = 0; j < evaldim; ++j)
                {
                    for (int All = 0; All < ode_dimension; ++All)
                        saving[All] = evalcoordinates[j, All];

                    for (i = 0; i < ode_dimension; ++i)
                    {
                        for (k = 0; k < maxite; ++k)
                        {
                            for (int All = 0; All < ode_dimension; ++All)
                                savingdomain[All] = rbf.collocationpoints[k, All];

                            EqDiff.eqdff(normal, savingdomain, resultk);
                            diffpoints = ArrayOperations.VDiff(saving, savingdomain);
                            twopointsdistance = Math.Sqrt(ArrayOperations.Dot(diffpoints, diffpoints));
                            checking = 1.0 - c * twopointsdistance;
                            if (checking > 0.0)
                            {
                                wdlfvalue1 = Wendland.WndlndFnctnFirst(twopointsdistance, c);
                                wdlfvalue2 = Wendland.WndlndFnctnSecond(twopointsdistance, c);
                                fdvector[j, i] += betaod[k] * (-resultk[i] * wdlfvalue1
                                                          - diffpoints[i]
                                                          * ArrayOperations.Dot(diffpoints, resultk)
                                                          * wdlfvalue2);
                            }
                        }
                    }
                }
            }
            double numbernormsquare = 0.0;
            Array.Resize(ref normed, evaldim);
            double[] temp1 = new double[ode_dimension];
            for (int p = 0; p < evaldim; ++p)
            {
                for (int All = 0; All < ode_dimension; ++All)
                {
                    temp1[All] = fdvector[p, All];
                }
                numbernormsquare = ArrayOperations.Dot(temp1, temp1);
                normed[p] = Math.Sqrt(numbernormsquare);
            }
            if (printing)
            {
                if (type_of_grid)
                {
                    gnl.PrintColumnsToFile("lyapprimexdir", currentiteration, fdvector);
                    gnl.PrintVectorToFile("normeddire", currentiteration, ref normed);
                }
                else
                {
                    gnl.PrintColumnsToFile("lyapprimexcar", currentiteration, fdvector);
                    gnl.PrintVectorToFile("normedcar", currentiteration, ref normed);
                }
            }
            DateTime end = DateTime.Now;
            TimeSpan ts = (end - start);
            Console.WriteLine("=====Obtener el gradiente y su norma tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====Obtener el gradiente y su norma tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
        }
        public void GetNewAlpha(int currentiteration, RBF rbf)
        {
            DateTime start = DateTime.Now;
            double summing = 0.0;
            double normalizationfactor = 0.0;
            int evaldim = rbf.alphavector.Length;
            Array.Resize(ref rbf.alphavector, evaldim);
            for (int iii = 0; iii < evaldim; ++iii)
            {
                summing = 0.0;
                for (int j = 0; j < 2 * points_directional; ++j)
                {
                    summing += orbder[(2 * points_directional) * iii + j];
                }
                if (summing > 0.0)
                {
                    summing = 0.0;
                }
                rbf.alphavector[iii] = summing / ((double)(2 * points_directional));
                normalizationfactor += rbf.alphavector[iii];
            }
            for (int All = 0; All < evaldim; ++All)
                rbf.alphavector[All] = Math.Abs(evaldim / normalizationfactor) * rbf.alphavector[All];

            if (printing)
            {
                gnl.PrintVectorToFile("alphavector", currentiteration, ref rbf.alphavector);
            }
            DateTime end = DateTime.Now;
            TimeSpan ts = (end - start);
            Console.WriteLine("=====Obtener el nuevo vector alfa tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====Obtener el nuevo vector alfa tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
        }
        public int ode_dimension;
        public int c;
        public int points_directional;
        public double critval;
        public bool normal;
        public bool printing;
        public double[] betaod = new double[] { };
        public double[] lyapfunc = new double[] { };
        public double[] orbder = new double[] { };
        double[] normed = new double[] { };
        double[,] fdvector = new double[,] { };
        double[][] Amatperm = new double[][] { };
    }
}
