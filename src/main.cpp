// Tyler Travis A01519795
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <complex>
#include <iomanip>


int main(int argc, char** argv)
{
  //using namespace std::literals;
  std::cout << std::fixed << std::setprecision(1);

  std::complex<double> z1 = 1i * 1i;
  std::cout << "i * i = " << z1 << std::endl;


  //std::complex<double> z2 = std::pow(1i, 2);
  //std::cout << "pow(i, 2) = " << z2 << std::endl;

  /*double PI = std::acos(-1);
  std::complex<double> z3 = std::exp(1i * PI);
  std::cout << "exp(i, pi) = " z3 << std::endl;*/

  std::complex<double> z4 = 1. + 2i, z5 = 1. -2i;
  std::cout << "(1+2i)*(1-2i) = " << z4*z5 << std::endl;
}
