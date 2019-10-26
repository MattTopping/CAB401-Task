﻿using System;
using System.Numerics;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Collections.Generic;
using System.Linq;

namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        public float[][] timeFreqData;
        public int wSamp;
        public Complex[] twiddles;

        public timefreq(float[] x, int windowSamp)
        {
            int ii;
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];
            for (ii = 0; ii < wSamp; ii++)
            {
                double a = 2 * pi * ii / (double)wSamp;
                twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            }

            timeFreqData = new float[wSamp/2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            Complex[] compX = new Complex[nearest];
            for (int kk = 0; kk < nearest; kk++)
            {
                if (kk < x.Length)
                {
                    compX[kk] = x[kk];
                }
                else
                {
                    compX[kk] = Complex.Zero;
                }
            }


            int cols = 2 * nearest /wSamp;

            // Overhead too large to make parallel viable
            //Parallel.For(0, wSamp / 2, jj =>
            //{
            //    timeFreqData[jj] = new float[cols];
            //});

            for (int jj = 0; jj < wSamp / 2; jj++)
            {
                timeFreqData[jj] = new float[cols];
            }

            timeFreqData = stft(compX, wSamp);
	
        }

        float[][] stft(Complex[] x, int wSamp)
        {
            int N = x.Length;
            float fftMax = 0;
            
            float[][] Y = new float[wSamp / 2][];

            //Overhead too great, original execution more effective
            //Parallel.For(0, wSamp / 2, ll =>
            //{
            //   Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            //});

            for (int ll = 0; ll < wSamp / 2; ll++)
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            }

            Complex[] temp = new Complex[wSamp];
            Complex[] tempFFT = new Complex[wSamp];

            //List<float> threadMaxes = new List<float>();
            //float fftMax_thread = 0;

            //Parallel.For(0, numThreads, iterator =>
            //{
            //    int guard = 2 * (int)Math.Floor(N / (double)wSamp) - 1;
            //    int chunk_size = (guard + (numThreads - 1)) / numThreads;
            //    int start = chunk_size * iterator;
            //    int end = Math.Min(start + chunk_size, guard);

            //    Debug.WriteLine("\n \nStart chunk: " + start + "\n  End chunk: " + end);

            //    for (int ii = start; ii < end; ii++)
            //    {
            //        for (int jj = 0; jj < wSamp; jj++)
            //        {
            //            temp[jj] = x[ii * (wSamp / 2) + jj];
            //        }

            //        tempFFT = fft(temp);

            //        for (int kk = 0; kk < wSamp / 2; kk++)
            //        {
            //            Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

            //            if (Y[kk][ii] > fftMax_thread)
            //            {
            //                fftMax_thread = Y[kk][ii];
            //            }
            //        }
            //    }
            //    Debug.WriteLine(fftMax_thread + "\n");
            //    threadMaxes.Add(fftMax_thread);
            //});

            //Debug.WriteLine("Max value: " + threadMaxes.Max());
            //fftMax = threadMaxes.Max();

            //////////////////////////////////////////////////////////////////////////////////////////////

            //Parallel.For(0, (int)(2 * Math.Floor(N / (double)wSamp) - 1), ii =>
            //{
            //    for (int jj = 0; jj < wSamp; jj++)
            //    {
            //        temp[jj] = x[ii * (wSamp / 2) + jj];
            //    }

            //    tempFFT = fft(temp);

            //    for (int kk = 0; kk < wSamp / 2; kk++)
            //    {
            //        Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

            //        if (Y[kk][ii] > fftMax)
            //        {
            //            fftMax = Y[kk][ii];
            //        }
            //    }
            //}
            //);

            //Debug.WriteLine("\n\n fftMax: " + fftMax);

            /////////////////////////////////////////////////////////////////////////////////////////////////

            for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (int jj = 0; jj < wSamp; jj++)
                {
                    temp[jj] = x[ii * (wSamp / 2) + jj];
                }

                tempFFT = fft(temp);

                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

                    if (Y[kk][ii] > fftMax)
                    {
                        fftMax = Y[kk][ii];
                    }
                }
            }



            for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            }


            return Y;
        }

        Complex[] fft(Complex[] x)
        {
            int ii = 0;
            int kk = 0;
            int N = x.Length;

            Complex[] Y = new Complex[N];

            // NEED TO MEMSET TO ZERO?

            if (N == 1)
            {
                Y[0] = x[0];
            }
            else{

                Complex[] E = new Complex[N/2];
                Complex[] O = new Complex[N/2];
                Complex[] even = new Complex[N/2];
                Complex[] odd = new Complex[N/2];

                for (ii = 0; ii < N; ii++)
                {

                    if (ii % 2 == 0)
                    {
                        even[ii / 2] = x[ii];
                    }
                    if (ii % 2 == 1)
                    {
                        odd[(ii - 1) / 2] = x[ii];
                    }
                }

                E = fft(even);
                O = fft(odd);

                for (kk = 0; kk < N; kk++)
                {
                    Y[kk] = E[(kk % (N / 2))] + O[(kk % (N / 2))] * twiddles[kk * wSamp / N];
                }
            }

           return Y;
        }
        
    }
}
