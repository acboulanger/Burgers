#include<stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI atan(1)

double* GenerateGaussianWhiteNoise(double sigma, int size)
{
  double* a = (double*)malloc(size*sizeof(double));
  double* b = (double*)malloc(size*sizeof(double));
  double* res = (double*)malloc(size*sizeof(double));
  int i;
  
  //Generate first uniform law
  for(i=0;i<size;i++)
  {
    a[i] = 1.0*rand()/RAND_MAX;
  }
  //Generate second uniform law
  for(i=0;i<size;i++)
  {
    b[i] = 1.0*rand()/RAND_MAX;
  }
  //Generate Gaussian noise, Box Muller method
  for(i=0;i<size;i++)
  {
    res[i] = sigma*sqrt(-2*log(a[i]))*cos(2*PI*b[i]);
  }
  return res;
}

int main()
{
  srand(time(NULL));
  int size = 10;
  double sigma = 0.1;
  int i;
  double* bruit = GenerateGaussianWhiteNoise(sigma,size);
  for(i = 0;i< size;i++)
  {
    printf("test %f \n", bruit[i]);
  }
  
  return EXIT_SUCCESS;
}