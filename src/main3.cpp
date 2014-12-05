#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>

#define MAX_POWER 4
// Filter Length
#define Nf 256
// Length of Signal
#define N 496125
// Sampling frequency
const double Fs = 11025;

const std::complex<double> WN[] = {1.0, 1.0/2.0+sqrt(3.0)/2.0i, -1.0/2.0+sqrt(3.0)/2.0i,
                                  -1.0, -1.0/2.0-sqrt(3.0)/2.0i, 1.0/2.0-sqrt(3.0)/2.0i};

void fft6(int, int, std::complex<double>*);
void twiddle(std::complex<double>*, int, double);
void bit_reorder(std::complex<double>*, int);

int main()
{
  // Input stream for filter
  std::ifstream filterIn("../data/LowPassFilter.dat");
  // Input stream for signal x[n]
  std::ifstream xIn("../data/flute.dat");

  // filter of length Nf = 256
  // Nf*4 for zero padding
  std::complex<double> h[4*Nf];

  // input vairable
  double in;

  // Read in the filter data
  for(int i = 0; i < 4*Nf; ++i)
  {
    h[i] = 0;
  }
  for(int n = 0; n < Nf; ++n)
  {
    filterIn >> in;
    h[n] = in;
  }

  // Output streams for the input x int argc, char** argvsignal
  // and the output y signalt
  std::ofstream x_dat("../data/x3.dat");
  std::ofstream y_dat("../data/y3.dat");
  std::ofstream H_dat("../data/H3.dat");

  
  // input x signal of length N = 25600
  std::complex<double>* x;
  x = (std::complex<double>*)malloc(sizeof(std::complex<double>)*N);

  // output y signal of Length N = 25600
  std::complex<double>* y;
  y = (std::complex<double>*)malloc(sizeof(std::complex<double>)*(N+Nf-1));

  // Generate input signal x[n]
  for(int n = 0; n < N; ++n)
  {
    xIn >> in;
    x[n].real(in);
    x[n].imag(0);
    x_dat << x[n].real() << std::endl;
  }

  //const int nfft = 1024;
  
  int M = 256;
  int overlap = M-1;
  int nfft = 1296; 
  int stepsize = nfft - overlap;

  std::complex<double> H[nfft];
  memcpy(H, h, sizeof(h)); 
  // generate fft
  fft6(0, nfft, H);
  bit_reorder(H, nfft);
  for(int i = 0; i < nfft; ++i)
  {
    H_dat << H[i].real() << "\t" << H[i].imag() <<std::endl;
  }

  std::complex<double> yt[nfft];
  std::complex<double> xt[nfft];

  int position = 0;
  while(position + nfft <= N)
  {
    for(int j = 0; j < nfft; ++j)
    {
      xt[j] = x[j + position];
    }

    fft6(0, nfft, xt);
    bit_reorder(xt, nfft);

    for(int k = 0; k < nfft; ++k)
    {
      yt[k] = xt[k] * H[k];
    }
    fft6(1, nfft, yt);
    bit_reorder(yt, nfft);
    for(int j = M-1; j < nfft; ++j)
    {
      y[j-M+position] = yt[j];
    }
    position += stepsize;
  }
  for(int n = 0; n < N; ++n)
  {
    y_dat << y[n].real() << std::endl;
  }
  free(x);
  //free(y);
  return 0;
}

void bit_reorder(std::complex<double>* x, int n)
{
  int power, N1, N2, N3;
  int N4 = 1;
  std::complex<double> temp[n];
  for(int i = 0; i < n; ++i)
  {
    temp[i] = x[i];
  }

  for(int i = 0; i <= MAX_POWER; ++i)
  {
    if(pow(6,i) == n)
    {
      power = i;
      N1 = pow(6,i)/6.0;
      if(N1 > 1)
      {
        N2 = N1/6.0;
      }
      else
      {
        N2 = 1;
      }
      if(N2 > 1)
      {
        N3 = N2/6.0;
      }
      else
      {
        N3 = 1;
      }
    }
  }

  int index = 0;
  for(int i = 0; i < 6; i++)
  {
    for(int j = 0; j < N1/N2; j++)
    {
      for(int k = 0; k < N2/N3; k++)
      {
        for(int l = 0; l < N3/N4; l++)
        {
          if(index > N)
            break;

          if(N1 == 1)
          {
            x[i] = temp[index++];
          }
          else if(N2 == 1)
          {
            x[i+j*N1] = temp[index++];
          }
          else if(N3 == 1)
          {
            x[i*N3 + j*N2 + k*N1] = temp[index++];
          }
          else
          {
            x[i*N4 + j*N3 + k*N2 + l*N1] = temp[index++];
          }
        }
      }
    }
  }
}

void twiddle(std::complex<double>* W, int n, double k)
{
  W->real(cos(k*2*M_PI/(double)n));
  W->imag(-sin(k*2*M_PI/(double)n));
}

void fft6(int in, int M, std::complex<double>* x)
{
  std::complex<double> W, butterfly[6];

  int N1 = 6;
  int N2 = M/6;

  if(in == 1)
  {
    for(int i = 0; i < M; i++)
    {
      x[i] = std::conj(x[i]);
    }
  }
  for(int n = 0; n < N2; n++)
  {
    butterfly[0] = (WN[0]*x[n] + WN[0]*x[N2+n] + WN[0]*x[2*N2+n] + WN[0]*x[3*N2+n] + WN[0]*x[4*N2+n] + WN[0]*x[5*N2+n]);
    butterfly[1] = (WN[0]*x[n] + WN[1]*x[N2+n] + WN[2]*x[2*N2+n] + WN[3]*x[3*N2+n] + WN[4]*x[4*N2+n] + WN[5]*x[5*N2+n]);
    butterfly[2] = (WN[0]*x[n] + WN[2]*x[N2+n] + WN[4]*x[2*N2+n] + WN[0]*x[3*N2+n] + WN[2]*x[4*N2+n] + WN[4]*x[5*N2+n]);
    butterfly[3] = (WN[0]*x[n] + WN[3]*x[N2+n] + WN[0]*x[2*N2+n] + WN[3]*x[3*N2+n] + WN[0]*x[4*N2+n] + WN[3]*x[5*N2+n]);
    butterfly[4] = (WN[0]*x[n] + WN[4]*x[N2+n] + WN[2]*x[2*N2+n] + WN[0]*x[3*N2+n] + WN[4]*x[4*N2+n] + WN[2]*x[5*N2+n]);
    butterfly[5] = (WN[0]*x[n] + WN[5]*x[N2+n] + WN[4]*x[2*N2+n] + WN[3]*x[3*N2+n] + WN[2]*x[4*N2+n] + WN[1]*x[5*N2+n]);
    for(int k = 0; k < N1; ++k)
    {
      twiddle(&W, M, (double)k*(double)n);
      x[n + N2*k] = butterfly[k]*W;
    }
  }
  if(N2 != 1)
  {
    for(int k = 0; k < N1; k++)
    {
      fft6(2, N2, &x[N2*k]);
    }
  }
  if(in == 1)
  {
    for(int i = 0; i < M; i++)
    {
      x[i] /= M;
    }
  }
  if(in == 1)
  {
    for(int i = 0; i < M; i++)
    {
      x[i] = std::conj(x[i]);
    }
  }
}
