using System;
using System.Collections.Generic;
using System.Linq;

namespace LyapXool
{
    class RBF
    {
        Generalities gnl = new Generalities();
        public RBF(int Iode_dimension, int Ic, int Imaxpos, int Imaxneg, double Ialpha, int Ipoints_directional, double Iradius, double Icart_grid_density, double[] Imin_geometric_limits, double[] Imax_geometric_limits, bool Inormal, bool Iprinting)
        {
            ode_dimension = Iode_dimension;
            c = Ic;
            maxpos = Imaxpos;
            maxneg = Imaxneg;
            alpha = Ialpha;
            points_directional = Ipoints_directional;
            radius = Iradius;
            cart_grid_density = Icart_grid_density;
            min_geometric_limits = Imin_geometric_limits;
            max_geometric_limits = Imax_geometric_limits;
            normal = Inormal;
            printing = Iprinting;
        }

        public void WBase()
        {
            ArrayOperations.ResizeArray<double>(ref rbfbasis, (int)ode_dimension, (int)ode_dimension);
            ArrayOperations.FillMatrixConst(ref rbfbasis, 0.0);

            double[] ek = new double[ode_dimension];

            Array.Fill(ek, 0.0);

            for (int i = 1; i <= ode_dimension; ++i)
            {
                for (int j = 1; j <= i; ++j)
                {
                    if (j == i)
                    {
                        rbfbasis[i - 1, j - 1] = (i + 1) * Math.Sqrt(1.0 / (2.0 * j * (j + 1.0)));
                    }
                    else
                    {
                        rbfbasis[i - 1, j - 1] = Math.Sqrt(1.0 / (2.0 * j * (j + 1.0)));
                    }
                }
            }

        }
        public void Grid()
        {
            int i = 0, j = 0, k = 0;

            int elements;
            int mn = 1;
            ArrayOperations.FillMatrixConst(ref coord, 0.0);
            elements = Math.Abs(maxneg) + Math.Abs(maxpos) + 1;

            int[] m = new int[ode_dimension];
            List<List<int>> v = new List<List<int>>();
            v.Resize(ode_dimension);
            {
                for (i = 0; i < ode_dimension; ++i)
                {
                    v[i].Resize(elements);
                    for (j = 0; j < elements; ++j)
                    {
                        v[i][j] = maxneg + j;
                    }
                }
                for (i = 0; i < ode_dimension; ++i)
                {
                    m[i] = (int)v[i].Count;
                    mn *= m[i];
                }
                ArrayOperations.ResizeArray<double>(ref coord, (int)mn, (int)ode_dimension);

                for (i = 0; i < mn; ++i)
                {
                    k = i;
                    for (j = ode_dimension - 1; j >= 0; --j)
                    {
                        coord[i, j] = v[j][k % m[j]];
                        k /= m[j];
                    }
                }
            }
            v.Clear();
        }
        public void EvaluatingGrid()
        {
            int elements;
            int mn = 1;
            double maxmax = double.NegativeInfinity;
            double minmin = double.PositiveInfinity;

            for (int jc = 0; jc < ode_dimension; ++jc)
            {
                if (min_geometric_limits[jc] <= minmin)
                {
                    minmin = min_geometric_limits[jc];
                }
                if (max_geometric_limits[jc] >= maxmax)
                {
                    maxmax = max_geometric_limits[jc];
                }
            }
            int dimension = 1 + (int)((Math.Abs(maxmax) + Math.Abs(minmin)) / cart_grid_density);
            double[] evaluatingpoints = new double[dimension];
            Array.Resize(ref evaluatingpoints, dimension);
            Array.Fill(evaluatingpoints, 0.0);
            for (int i = 0; i < dimension; ++i)
            {
                evaluatingpoints[i] = minmin + i * cart_grid_density;
            }
            elements = (int)evaluatingpoints.Length;

            int[] m = new int[ode_dimension];

            List<List<double>> v = new List<List<double>>();
            v.Resize(ode_dimension);
            for (int i = 0; i != v.Count; ++i)
            {
                v[i].Clear();
                v[i].Resize(elements);
                int kj = 0;
                for (int j = 0; j != v[i].Count; ++j)
                {
                    v[i][j] = evaluatingpoints[kj];
                    ++kj;
                }
            }
            for (int i = 0; i < ode_dimension; ++i)
            {
                m[i] = (int)v[i].Count;
                mn *= m[i];
            }

            ArrayOperations.ResizeArray<double>(ref cuadricula, (int)mn, (int)ode_dimension);
            for (int i = 0; i < mn; ++i)
            {
                int k = i;
                for (int j = ode_dimension - 1; j >= 0; --j)
                {
                    cuadricula[i, j] = v[j][k % m[j]];
                    k /= m[j];
                }
            }
            v.Clear();
        }
        public void EffectiveGrid(ref double[,] gridtobeclean, ref double[,] cleanedgrid)
        {
            for (int jc = 0; jc < ode_dimension; ++jc)
            {
                if (max_geometric_limits[jc] <= min_geometric_limits[jc])
                {
                    Console.WriteLine("ERROR: Maximum should be larger than Minimum");
                    Console.WriteLine("Entry: " + jc + " value: " + max_geometric_limits[jc]);
                    System.Environment.Exit(9);
                }
            }
            int dim1 = (int)gridtobeclean.GetLength(0);//longitud
            int dim2 = (int)gridtobeclean.GetLength(1);//anchura
            List<int> counter = new List<int>();
            for (int i = 0; i < dim1; ++i)
            {
                int inside = 0;
                for (int jc = 0; jc < ode_dimension; ++jc)
                {
                    if ((gridtobeclean[i, jc] <= max_geometric_limits[jc]
                        && gridtobeclean[i, jc] >= min_geometric_limits[jc]))
                    {
                        ++inside;
                    }
                    else
                    {
                        break;
                    }
                }
                if (inside == ode_dimension)
                {
                    counter.Add(i);
                }
            }

            int fin = (int)counter.Count;

            ArrayOperations.ResizeArray<double>(ref cleanedgrid, fin, dim2);
            int n = 0;
            for (var i = 0; i < fin; ++i)
            {
                for (int j = 0; j < dim2; ++j)
                {
                    cleanedgrid[n, j] = gridtobeclean[counter[i], j];
                }
                n++;
            }
            counter.Clear();
        }
        public void RBFGrid()
        {
            int i = 0, j = 0;
            int length = (int)Math.Pow((Math.Abs(maxneg) + Math.Abs(maxpos) + 1), ode_dimension);

            ArrayOperations.ResizeArray<double>(ref gridrbf, length, ode_dimension);
            ArrayOperations.FillMatrixConst(ref gridrbf, 0.0);

            int dimgpoint1 = (int)coord.GetLength(0);
            int dimgpoint2 = (int)coord.GetLength(1);
            for (i = 0; i < dimgpoint1; ++i)
            {
                for (j = 0; j < dimgpoint2; ++j)
                {
                    for (int All = 0; All < ode_dimension; ++All)
                    {
                        gridrbf[i, All] += alpha * coord[i, j] * rbfbasis[j, All];
                    }
                }
            }


            {
                for (i = 0; i < dimgpoint1; ++i)
                {
                    for (j = 0; j < dimgpoint2; ++j)
                    {
                        for (int All = 0; All < ode_dimension; ++All)
                        {
                            gridrbf[i, All] += alpha * coord[i, j] * rbfbasis[j, All];
                        }
                    }
                    if (dimgpoint2 == 2)
                    {
                        for (int All = 0; All < ode_dimension; ++All)
                        {
                            gridrbf[i, All] += 0.5 * alpha * (rbfbasis[0, All] + rbfbasis[1, All]);
                        }
                    }
                    if (dimgpoint2 == 3)
                    {
                        for (int All = 0; All < ode_dimension; ++All)
                        {
                            gridrbf[i, All] += 0.5 * alpha * (rbfbasis[0, All] + rbfbasis[1, All] + rbfbasis[2, All]);
                        }
                    }
                }

            }
        }
        public void AlphaFunction()
        {
            int N = collocationpoints.GetLength(0);
            Array.Resize(ref alphavector, N);
            Array.Fill(alphavector, -1.0);
        }
        public void InterpolationMatrixA()
        {
            DateTime start = DateTime.Now;

            int j = 0, k = 0;


            double atzero;

            int dimA = (int)collocationpoints.GetLength(0);
            int dimAc = (int)collocationpoints.GetLength(1);

            ArrayOperations.ResizeArray<double>(ref Amat, dimA, dimA);
            ArrayOperations.FillMatrixConst(ref Amat, 0.0);


            atzero = Wendland.WndlndFnctnFirst(0.0, c);


            {
                double[] diffsave = new double[dimAc];
                double[] savingcallj = new double[dimAc];
                double[] savingcallk = new double[dimAc];
                double[] resultj = new double[dimAc];
                double[] resultk = new double[dimAc];

                Array.Fill(diffsave, 0.0);
                Array.Fill(savingcallj, 0.0);
                Array.Fill(savingcallk, 0.0);
                Array.Fill(resultj, 0.0);
                Array.Fill(resultk, 0.0);

                double twopointsdistance = 0.0;

                double wdlfvalue1 = 0.0;
                double wdlfvalue2 = 0.0;
                double checking = 0.0;


                for (j = 0; j < dimA; ++j)
                {
                    for (int all = 0; all < dimAc; ++all)
                        savingcallj[all] = collocationpoints[j, all];


                    EqDiff.eqdff(normal, savingcallj, resultj);
                    for (k = 0; k < dimA; ++k)
                    {
                        for (int all = 0; all < dimAc; ++all)
                            savingcallk[all] = collocationpoints[k, all];

                        diffsave = ArrayOperations.VDiff(savingcallj, savingcallk);
                        if (k == j)
                        {
                            Amat[j, k] = -atzero * ArrayOperations.Dot(resultj, resultj);
                        }
                        else
                        {
                            twopointsdistance = Math.Sqrt(ArrayOperations.Dot(diffsave, diffsave));
                            checking = 1.0 - c * twopointsdistance;
                            if (checking > 0.0)
                            {
                                EqDiff.eqdff(normal, savingcallk, resultk);
                                wdlfvalue1 = Wendland.WndlndFnctnFirst(twopointsdistance, c);
                                wdlfvalue2 = Wendland.WndlndFnctnSecond(twopointsdistance, c);
                                Amat[j, k] = -wdlfvalue2 * ArrayOperations.Dot(diffsave, resultj) * ArrayOperations.Dot(diffsave, resultk)
                                             - wdlfvalue1 * ArrayOperations.Dot(resultj, resultk);
                            }
                        }
                    }
                }
            }
            DateTime end = DateTime.Now;
            TimeSpan ts = (end - start);
            Console.WriteLine("=====Calcular la matriz de colocación tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Days, ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====Calcular la matriz de colocación tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Days, ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
        }
        public void DirecGrid()
        {
            int lrows = (int)collocationpoints.GetLength(0);
            int lcols = (int)collocationpoints.GetLength(1);


            int newlenght = (int)(points_directional * 2 * lrows);

            int j, jd;
            stride = points_directional * 2 + 1;

            double norm;
            ArrayOperations.ResizeArray<double>(ref coldirectgrid, lrows * stride, lcols);
            ArrayOperations.ResizeArray<double>(ref directgrid, newlenght, lcols);


            double[,] domain = new double[newlenght, lcols];

            double[] savingdomain = new double[lcols];
            double[] evaldfunction = new double[lcols];
            {
                for (int i = 0; i < lrows; ++i)
                {
                    j = stride * i;
                    jd = (stride - 1) * i;
                    for (int All = 0; All < lcols; ++All)
                    {
                        coldirectgrid[j, All] = collocationpoints[i, All];
                        savingdomain[All] = collocationpoints[i, All];
                    }

                    EqDiff.eqdff(normal, savingdomain, evaldfunction);
                    norm = Math.Sqrt(ArrayOperations.Dot(evaldfunction, evaldfunction));
                    int kp = 0;
                    for (int kd = 0; kd < points_directional; kd += 1)
                    {
                        for (int All = 0; All < lcols; ++All)
                        {
                            directgrid[jd + kp, All] = collocationpoints[i, All] + (radius / points_directional) * (kd + 1) * alpha * (evaldfunction[All] / norm);
                            directgrid[jd + kp + 1, All] = collocationpoints[i, All] - (radius / points_directional) * (kd + 1) * alpha * (evaldfunction[All] / norm);
                            coldirectgrid[j + kp + 1, All] = directgrid[jd + kp, All];
                            coldirectgrid[j + kp + 2, All] = directgrid[jd + kp + 1, All];
                        }
                        kp += 2;
                    }
                }
            }


            List<int> counter = new List<int>();
            List<int> counterf = new List<int>();
            {
                int cdrows = (int)coldirectgrid.GetLength(0);
                int drows = (int)directgrid.GetLength(0);
                for (int i = 0; i < cdrows; ++i)
                {
                    for (int jc = 0; jc < ode_dimension; ++jc)
                    {
                        if ((coldirectgrid[i, jc] <= max_geometric_limits[jc]) && (coldirectgrid[i, jc] >= min_geometric_limits[jc]))
                        {
                            counter.Add(i);
                        }

                    }
                }
                for (int ii = 0; ii < drows; ++ii)
                {
                    for (int jc = 0; jc < ode_dimension; ++jc)
                    {
                        if ((directgrid[ii, jc] <= max_geometric_limits[jc]) && (directgrid[ii, jc] >= min_geometric_limits[jc]))
                        {
                            counterf.Add(ii);
                        }
                    }
                }
            }

            int ana = (int)counter.Count;
            int flo = (int)counterf.Count;

            double[,] cleanbigag = new double[,] { };
            double[,] cleanbigfg = new double[,] { };
            var boolcoldirectgrid = Enumerable.Repeat<bool>(false, stride * lrows).ToArray(); ;
            var booldirectgrid = Enumerable.Repeat<bool>(false, lrows * (stride - 1)).ToArray();


            ArrayOperations.ResizeArray<double>(ref cleanbigag, ana, lcols);
            ArrayOperations.ResizeArray<double>(ref cleanbigfg, flo, lcols);

            int n = 0;
            int m = 0;
            int end = counter.Count;
            for (int i = 0; i < end; ++i)
            {
                boolcoldirectgrid[counter[i]] = true;
                for (int All = 0; All < lcols; ++All)
                    cleanbigag[n, All] = coldirectgrid[counter[i], All];
                n++;
            }
            end = counterf.Count;
            for (int i = 0; i < end; ++i)
            {
                booldirectgrid[counterf[i]] = true;
                for (int All = 0; All < lcols; ++All)
                    cleanbigfg[m, All] = directgrid[counterf[i], All];
                m++;
            }

            counter.Clear();
            counterf.Clear();

            if (printing)
            {
                gnl.PrintColumnsToFile("direcgrid", 0, directgrid);
            }
            Console.WriteLine("La cantidad de puntos en la malla direccional de evalucación es: " + directgrid.GetLength(0));
            Instructions.woutput.WriteLine("La cantidad de puntos en la malla direccional de evalucación es: " + directgrid.GetLength(0));
        }
        public void GetAlphaColMatrix()
        {
            AlphaFunction();
            InterpolationMatrixA();
        }
    public void MakeColGrid()
        {
            DateTime start = DateTime.Now;
            WBase();
            Grid();
            RBFGrid();
            EffectiveGrid(ref gridrbf, ref collocationpoints);
            Console.WriteLine("La cantidad de puntos en la malla de colocación es: " + collocationpoints.GetLength(0));
            Instructions.woutput.WriteLine("La cantidad de puntos en la malla de colocación es: " + collocationpoints.GetLength(0));
            gnl.PrintColumnsToFile("col", collocationpoints);

            DateTime end = DateTime.Now;
            TimeSpan ts = (end - start);
            Console.WriteLine("=====Construir la malla de colocación tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Days, ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====Construir la malla de colocación tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Days, ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
        }

        public void MakeEvalGrid(ref double[,] evalcoordinates)
        {
            DateTime start = DateTime.Now;
            EvaluatingGrid();
            EffectiveGrid(ref cuadricula, ref evalcoordinates);
            gnl.PrintColumnsToFile("cg", cartesianevalgrid);
            Console.WriteLine("La cantidad de puntos en la malla cartesiana de evalucación es: " + cartesianevalgrid.GetLength(0));
            Instructions.woutput.WriteLine("La cantidad de puntos en la malla cartesiana de evalucación es: " + cartesianevalgrid.GetLength(0));

            DateTime end = DateTime.Now;
            TimeSpan ts = (end - start);
            Console.WriteLine("=====Construir la malla cartesiana tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Days, ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
            Instructions.woutput.WriteLine("=====Construir la malla cartesiana tomó: {0:00} horas, {1:00} minutos, {2:00}.{3} segundos=====", ts.Days, ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds);
        }
        static public int ode_dimension;
        public int c;
        public int maxpos;
        public int maxneg;
        public double alpha;
        public int points_directional;
        public double radius;
        double cart_grid_density;
        public double[] min_geometric_limits;
        public double[] max_geometric_limits;
        public bool normal;
        public bool printing;
        public int stride;
        public double[] alphavector = new double[] { };
        public double[,] coldirectgrid = new double[,] { };
        public double[,] gridrbf = new double[,] { };
        public double[,] directgrid = new double[,] { };
        public double[,] rbfbasis = new double[,] { };
        public double[,] coord = new double[,] { };
        public double[,] cuadricula = new double[,] { };
        public double[,] collocationpoints = new double[,] { };
        public double[,] Amat = new double[,] { };
        public double[,] cartesianevalgrid = new double[,] { };
    }
}
