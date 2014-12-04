// Tyler Travis A01519795
#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>


void fft6(int, int, std::complex<double>*);
void twiddle(std::complex<double>*, int, double);
void bit_reorder(std::complex<double>*, int);

const std::complex<double> WN[] = {1.0, 1.0/2.0+sqrt(3.0)/2.0i, -1.0/2.0+sqrt(3.0)/2.0i,
                                  -1.0, -1.0/2.0-sqrt(3.0)/2.0i, 1.0/2.0-sqrt(3.0)/2.0i};

int main(int argc, char** argv)
{
  const int N = 216;
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
  bit_reorder(x, N);
  for(int i = 0; i < N; ++i)
  {
    //std::cout << x[i] << std::endl;
    x_dat << x[i].real() << '\t' << x[i].imag() << std::endl;
  }
  std::cout << "FFT" << std::endl;
  fft6(0, N, x);
  for(int i = 0; i < N; ++i)
  {
    //std::cout << x[i] << std::endl;
    y_dat << x[i].real() << '\t' << x[i].imag() << std::endl;
  }
  return 0;
}

void twiddle(std::complex<double>* W, int N, double k)
{
  W->real(cos(k*2*M_PI/(double)N));
  W->imag(-sin(k*2*M_PI/(double)N));
}

void bit_reorder(std::complex<double>* x, int N)
{
  int power, N1;
  std::complex<double> temp[N];
  for(int i = 0; i < N; ++i)
  {
    temp[i] = x[i];
  }

  for(int i = 0; i < 10; ++i)
  {
    if(pow(6,i) == N)
    {
      power = i;
      N1 = pow(6,i)/6.0;
    }
  }
  std::cout << power << std::endl;
  for(int i = 0; i < N; i+=N1)
  {
    std::cout << "\ti: " << i << std::endl;
    for(int j = 0; j < N1; j++)
    {
      if(i+N1*j < N)
      {
        std::cout << "temp at " << i+N1*j << std::endl;
        x[i+j] = temp[i+N1*j];
      }
    }
  }
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
