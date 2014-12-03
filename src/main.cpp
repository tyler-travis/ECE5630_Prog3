// Tyler Travis A01519795
#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>

const int N = 36;

void fft6(int, int, std::complex<double>*);
void radix6(std::complex<double>*);

const std::complex<double> WN[] = {1.0, 1.0/2.0+sqrt(3.0)/2.0i, -1.0/2.0+sqrt(3.0)/2.0i,
                                  -1.0, -1.0/2.0-sqrt(3.0)/2.0i, 1.0/2.0-sqrt(3.0)/2.0i};

int main(int argc, char** argv)
{
  std::ofstream x_dat("../data/x.dat");
  std::complex<double> x[6];
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
    std::cout << x[i] << std::endl;
    x_dat << x[i] << std::endl;
  }
  std::cout << "FFT" << std::endl;
  fft6(0, N, x);
  for(int i = 0; i < N; ++i)
  {
    std::cout << x[i] << std::endl;
  }
  return 0;
}

void radix6(std::complex<double>* x)
{
  std::complex<double> temp[6];
  for(int i = 0; i < 6; ++i)
  {
    temp[i] = x[i];
  }

  x[0] = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5];
  x[1] = temp[0] + WN[1]*temp[1] + WN[2]*temp[2] + WN[3]*temp[3] + WN[4]*temp[4] + WN[5]*temp[5];
  x[2] = temp[0] + WN[2]*temp[1] + WN[4]*temp[2] + WN[0]*temp[3] + WN[2]*temp[4] + WN[4]*temp[5];
  x[3] = temp[0] + WN[3]*temp[1] + WN[0]*temp[2] + WN[3]*temp[3] + WN[0]*temp[4] + WN[3]*temp[5];
  x[4] = temp[0] + WN[4]*temp[1] + WN[2]*temp[2] + WN[0]*temp[3] + WN[4]*temp[4] + WN[2]*temp[5];
  x[5] = temp[0] + WN[5]*temp[1] + WN[4]*temp[2] + WN[3]*temp[3] + WN[2]*temp[4] + WN[1]*temp[5];
}

void fft6(int in, int n, std::complex<double>* x)
{
  int power = 0;
  for(int i = 0; i < 100; ++i)
  {
    if(pow(6,i) == n)
    {
      power = i;
    }
  }
  std::complex<double> xi[n/6][6];
  for(int i = 0; i < n/6; ++i)
  {
    for(int j = 0; j < 6; ++j)
    {
      xi[i][j] = x[i*6+j];
    }
  }
  for(int i = 0; i < power; ++i)
  {
    
  }
}
