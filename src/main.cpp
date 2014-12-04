// Tyler Travis A01519795
#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#define MAX_POWER 4

void fft6(int, int, std::complex<double>*);
void twiddle(std::complex<double>*, int, double);
void bit_reorder(std::complex<double>*, int);

const std::complex<double> WN[] = {1.0, 1.0/2.0+sqrt(3.0)/2.0i, -1.0/2.0+sqrt(3.0)/2.0i,
                                  -1.0, -1.0/2.0-sqrt(3.0)/2.0i, 1.0/2.0-sqrt(3.0)/2.0i};

int main(int argc, char** argv)
{
  const int N = atoi(argv[1]);
  std::ofstream x_dat("../data/x.dat");
  std::ofstream y_dat("../data/y.dat");
  std::complex<double> x[N];
  double freq1 = 17.01/11025;
  double freq2 = 297.74/11025;
  double freq3 = 425.35/11025;
  double freq4 = 2637/11025;
  for(int n = 0; n < N; ++n)
  {
    x[n] = cos(2*M_PI*freq1*n) + cos(2*M_PI*freq2*n) + cos(2*M_PI*freq3*n) + cos(2*M_PI*freq4*n);
  }
  for(int i = 0; i < N; ++i)
  {
    //std::cout << x[i] << std::endl;
    x_dat << x[i].real() << '\t' << x[i].imag() << std::endl;
  }
  std::cout << "FFT" << std::endl;
  fft6(0, N, x);
  if(N != 6)
    bit_reorder(x, N);
  for(int i = 0; i < N; ++i)
  {
    //std::cout << x[i] << std::endl;
    y_dat << x[i].real() << '\t' << x[i].imag() << std::endl;
  }
  return 0;
}

void bit_reorder(std::complex<double>* x, int N)
{
  int power, N1, N2, N3;
  std::cout << "N: " << N << std::endl;
  std::complex<double> temp[N];
  for(int i = 0; i < N; ++i)
  {
    temp[i] = x[i];
  }

  for(int i = 0; i <= MAX_POWER; ++i)
  {
    if(pow(6,i) == N)
    {
      power = i;
      N1 = pow(6,i)/6.0;
      N2 = N1/6.0;
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

  std::cout << "power: " << power << std::endl;
  std::cout << "N1: " << N1 << std::endl;
  std::cout << "N2: " << N2 << std::endl;
  int index = 0;
  for(int i = 0; i < 6; i++)
  {
    for(int j = 0; j < N1/N2*N3; j++)
    {
      for(int k = 0; k < N2/N3; k++)
      {
        if(index > N)
          break;

        if(i + j*N1/N2 + k*N1 < N)
        {
          std::cout << i << ',' << j << ',' << k << '\t';
          if(N2 == 1)
          {
            std::cout << index << "->" << i+j*6+k*N1 << std::endl;
            x[i+j*6+k*N1] = temp[index++];
          }
          else
          {
            std::cout << index << "->" << i+j*N2+k*N1 << std::endl;
            x[i+j*N2+k*N1] = temp[index++];
          }
        }
      }
    }
  }
}

void twiddle(std::complex<double>* W, int N, double k)
{
  W->real(cos(k*2*M_PI/(double)N));
  W->imag(-sin(k*2*M_PI/(double)N));
}

void fft6(int in, int N, std::complex<double>* x)
{
  std::complex<double> W, butterfly[6];

  int N1 = 6;
  int N2 = N/6;

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
      twiddle(&W, N, (double)k*(double)n);
      x[n + N2*k] = butterfly[k]*W;
    }
  }
  if(N2 != 1)
  {
    for(int k = 0; k < N1; k++)
    {
      fft6(in, N2, &x[N2*k]);
    }
  }
}
