/********************
  Modified NL means filter ODCT

  Michael Eager
  Zhaolin Chen

  Monash Biomedical Imaging
  Monash University, 2015
*******************/

/***************************************************************************                  
/* Jose V. Manjon - jmanjon@fis.upv.es                                     */
/* Universidad Politecnica de Valencia, Spain                              */
/* Pierrick Coupe - pierrick.coupe@gmail.com                               */
/* Brain Imaging Center, Montreal Neurological Institute.                  */
/* Mc Gill University                                                      */
/*                                                                         */
/* Copyright (C) 2010 Jose V. Manjon and Pierrick Coupe                    */

/***************************************************************************
 * New Methods for MRI Denoising based on Sparseness and Self-Similarity    *
 * Jos� V. Manj�n, Pierrick Coup�, Antonio Buades,D. Louis Collins 
 * and Montserrat Robles
 ***************************************************************************/

#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include "matrix.h"
#ifndef PI
#define PI 3.14159265358979323846
#endif
#define B 4
#define BB 16

/* Multithreading stuff*/
#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

typedef struct{
  int rows;
  int cols;
  int slices;
  float * in_image;
  float * out_image;
  float * out_image2;
  float * acu_image;
  float * B1_image;
  int ini;
  int fin;    
  float sigma;
  double * mibias;
}myargument;

double W1,W2; 
double *tabla1,*tabla2; /* tablas de cosenos*/
int rician;



/*
  https://dst.lbl.gov/ACSSoftware/colt/
*/
	
/****************************************
 *    COEFFICIENTS FOR METHODS i1, i1e  *
 ****************************************/

static double chbevl( double x, double coef[], int N ) {
  double b0, b1, b2;

  int p = 0;
  int i;

  b0 = coef[p++];
  b1 = 0.0;
  i = N - 1;

  do {
    b2 = b1;
    b1 = b0;
    b0 = x * b1  -  b2  + coef[p++];
  } while( --i > 0);

  return( 0.5*(b0-b2) );
}
/**
 * Returns the modified Bessel function of order 0 of the
 * argument.
 * <p>
 * The function is defined as <tt>i0(x) = j0( ix )</tt>.
 * <p>
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 * @param x the value to compute the bessel function of.
 */
static  double i0(double x) {
  double y;

  static  double A_i0[] = {
    -4.41534164647933937950E-18,
    3.33079451882223809783E-17,
    -2.43127984654795469359E-16,
    1.71539128555513303061E-15,
    -1.16853328779934516808E-14,
    7.67618549860493561688E-14,
    -4.85644678311192946090E-13,
    2.95505266312963983461E-12,
    -1.72682629144155570723E-11,
    9.67580903537323691224E-11,
    -5.18979560163526290666E-10,
    2.65982372468238665035E-9,
    -1.30002500998624804212E-8,
    6.04699502254191894932E-8,
    -2.67079385394061173391E-7,
    1.11738753912010371815E-6,
    -4.41673835845875056359E-6,
    1.64484480707288970893E-5,
    -5.75419501008210370398E-5,
    1.88502885095841655729E-4,
    -5.76375574538582365885E-4,
    1.63947561694133579842E-3,
    -4.32430999505057594430E-3,
    1.05464603945949983183E-2,
    -2.37374148058994688156E-2,
    4.93052842396707084878E-2,
    -9.49010970480476444210E-2,
    1.71620901522208775349E-1,
    -3.04682672343198398683E-1,
    6.76795274409476084995E-1
  };



  /**
   * Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
   * in the inverted interval [8,infinity].
   *
   * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
   */
  static  double B_i0[] = {
    -7.23318048787475395456E-18,
    -4.83050448594418207126E-18,
    4.46562142029675999901E-17,
    3.46122286769746109310E-17,
    -2.82762398051658348494E-16,
    -3.42548561967721913462E-16,
    1.77256013305652638360E-15,
    3.81168066935262242075E-15,
    -9.55484669882830764870E-15,
    -4.15056934728722208663E-14,
    1.54008621752140982691E-14,
    3.85277838274214270114E-13,
    7.18012445138366623367E-13,
    -1.79417853150680611778E-12,
    -1.32158118404477131188E-11,
    -3.14991652796324136454E-11,
    1.18891471078464383424E-11,
    4.94060238822496958910E-10,
    3.39623202570838634515E-9,
    2.26666899049817806459E-8,
    2.04891858946906374183E-7,
    2.89137052083475648297E-6,
    6.88975834691682398426E-5,
    3.36911647825569408990E-3,
    8.04490411014108831608E-1
  };


  if( x < 0 ) x = -x;
  if( x <= 8.0 ) {
    y = (x/2.0) - 2.0;
    return( exp(x) * chbevl( y, A_i0, 30 ) );
  }

  return(  exp(x) * chbevl( 32.0/x - 2.0, B_i0, 25 ) / sqrt(x) );
}
/**
 * Returns the exponentially scaled modified Bessel function
 * of order 0 of the argument.
 * <p>
 * The function is defined as <tt>i0e(x) = exp(-|x|) j0( ix )</tt>.
 *
 *
 * @param x the value to compute the bessel function of.
 *
 static  double i0e(double x){
 double y;
	
 if( x < 0 ) x = -x;
 if( x <= 8.0 ) {
 y = (x/2.0) - 2.0;
 return( chbevl( y, A_i0, 30 ) );
 }

 return( chbevl( 32.0/x - 2.0, B_i0, 25 ) / sqrt(x) );
 }*/
/**
 * Returns the modified Bessel function of order 1 of the
 * argument.
 * <p>
 * The function is defined as <tt>i1(x) = -i j1( ix )</tt>.
 * <p>
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 * @param x the value to compute the bessel function of.
 */
static  double i1(double x)  {
  double y, z;
  /**
   * Chebyshev coefficients for exp(-x) I1(x) / x
   * in the interval [0,8].
   *
   * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
   */
  static  double A_i1[] = {
    2.77791411276104639959E-18,
    -2.11142121435816608115E-17,
    1.55363195773620046921E-16,
    -1.10559694773538630805E-15,
    7.60068429473540693410E-15,
    -5.04218550472791168711E-14,
    3.22379336594557470981E-13,
    -1.98397439776494371520E-12,
    1.17361862988909016308E-11,
    -6.66348972350202774223E-11,
    3.62559028155211703701E-10,
    -1.88724975172282928790E-9,
    9.38153738649577178388E-9,
    -4.44505912879632808065E-8,
    2.00329475355213526229E-7,
    -8.56872026469545474066E-7,
    3.47025130813767847674E-6,
    -1.32731636560394358279E-5,
    4.78156510755005422638E-5,
    -1.61760815825896745588E-4,
    5.12285956168575772895E-4,
    -1.51357245063125314899E-3,
    4.15642294431288815669E-3,
    -1.05640848946261981558E-2,
    2.47264490306265168283E-2,
    -5.29459812080949914269E-2,
    1.02643658689847095384E-1,
    -1.76416518357834055153E-1,
    2.52587186443633654823E-1
  };

  /*
   * Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
   * in the inverted interval [8,infinity].
   *
   * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
   */
  static  double B_i1[] = {
    7.51729631084210481353E-18,
    4.41434832307170791151E-18,
    -4.65030536848935832153E-17,
    -3.20952592199342395980E-17,
    2.96262899764595013876E-16,
    3.30820231092092828324E-16,
    -1.88035477551078244854E-15,
    -3.81440307243700780478E-15,
    1.04202769841288027642E-14,
    4.27244001671195135429E-14,
    -2.10154184277266431302E-14,
    -4.08355111109219731823E-13,
    -7.19855177624590851209E-13,
    2.03562854414708950722E-12,
    1.41258074366137813316E-11,
    3.25260358301548823856E-11,
    -1.89749581235054123450E-11,
    -5.58974346219658380687E-10,
    -3.83538038596423702205E-9,
    -2.63146884688951950684E-8,
    -2.51223623787020892529E-7,
    -3.88256480887769039346E-6,
    -1.10588938762623716291E-4,
    -9.76109749136146840777E-3,
    7.78576235018280120474E-1
  };

  z = fabs(x);
  if( z <= 8.0 )
    {
      y = (z/2.0) - 2.0;
      z = chbevl( y, A_i1, 29 ) * z * exp(z);
    }
  else
    {
      z = exp(z) * chbevl( 32.0/z - 2.0, B_i1, 25 ) / sqrt(z);
    }
  if( x < 0.0 )
    z = -z;
  return( z );
}
/**
 * Returns the exponentially scaled modified Bessel function
 * of order 1 of the argument.
 * <p>
 * The function is defined as <tt>i1(x) = -i exp(-|x|) j1( ix )</tt>.
 * 
 * @param x the value to compute the bessel function of.
 *
 static  double i1e(double x) { 
 double y, z;

 z = fabs(x);
 if( z <= 8.0 )
 {
 y = (z/2.0) - 2.0;
 z = chbevl( y, A_i1, 29 ) * z;
 }
 else
 {
 z = chbevl( 32.0/z - 2.0, B_i1, 25 ) / sqrt(z);
 }
 if( x < 0.0 )
 z = -z;
 return( z );
 }*/



/*Returns the modified Bessel function I0(x) for any real x.*/
double Bessel0(double x)
{
  /*  double ax,ans,a;
      double y; 
      if ((ax=fabs(x)) < 3.75) 
      { 
      y=x/3.75;
      y*=y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
      } 
      else 
      {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax));
      a=y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))));
      ans=ans*(0.39894228 + y*(0.1328592e-1 +y*(0.225319e-2+y*(-0.157565e-2+a))));    
      }
      return ans;
  */
  return i0(x);
}

/*Returns the modified Bessel function I1(x) for any real x.*/
double Bessel1(double x)
{
  /*
    double ax,ans;
    double y; 
    if ((ax=fabs(x)) < 3.75)
    { 
    y=x/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
    } 
    else 
    {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
  */
  return i1(x);
}

void RicianBias(float level,int max,double * mibias, int N)
{
  int cz,j;
  double c;
  float z,i;          
    
  for(i=0;i<max*N;i=i+1) 
    {
      mibias[(int)i]=0;         
    }
    
  for(i=0;i<max;i=i+0.1)
    {      
      c=double((i*i))/(double((level*level)));
      z=(float) ((PI/2.0f)*exp(-(c/2.0f))*((1+c/2.0f)*Bessel0(c/4.0f)+(c/2.0f)*Bessel1(c/4.0f))*((1+c/2.0f)*Bessel0(c/4.0f)+(c/2.0f)*Bessel1(c/4.0f)));
      z=(sqrt(z*(level*level)));
      cz=(int)z;
      if(cz<=0) cz=(int)i;
      if(cz>=0 && cz<max)
	{ 
          if(mibias[cz]==0) 
	    {  
              for(j=0;j<N;j++) 
		{
                  mibias[j*max+cz]+=i;
		}              
	    }
	}
    }
}



/*http://epubs.siam.org/doi/pdf/10.1137/070693370*/
void FastDCT3D(float *in, float *out, int dim)
{
  float opa[BLK][MAX2][MAX2][MAX2]; float opb[BLK][MAX2][MAX2][MAX2];
  float h1,h2,h3,h4,h5,h6,h7,h8; float h1a,h2a,h3a,h4a,h5a,h6a,h7a,h8a;
  float h1b,h2b,h3b,h4b,h5b,h6b,h7b,h8b;
  register int i,j,k,l; register int dim2 = dim/2; register int index=dim/4-1;
  if (dim==2)
    {
      h1 = *(in); h2 = *(in+hc1); h3 = *(in+hc2); h4 = *(in+hc3);
      h5 = *(in+hc4); h6 = *(in+hc5); h7 = *(in+hc6); h8 = *(in+hc7);
      h1a=h1+h2; h2a=h1-h2; h3a=h3+h4; h4a=h3-h4;
      h5a=h5+h6; h6a=h5-h6; h7a=h7+h8; h8a=h7-h8;
      h1b=h1a+h3a; h2b=h2a+h4a; h3b=h1a-h3a; h4b=h2a-h4a;
      h5b=h5a+h7a; h6b=h6a+h8a; h7b=h5a-h7a; h8b=h6a-h8a;
      *(out+ad0) = (h1b+h5b); *(out+ad1) = (h2b+h6b)*C14a;
      *(out+ad2) = (h3b+h7b)*C14a; *(out+ad3) = (h4b+h8b)*C14b;
      *(out+ad4) = (h1b-h5b)*C14a; *(out+ad5) = (h2b-h6b)*C14b;
      *(out+ad6) = (h3b-h7b)*C14b; *(out+ad7) = (h4b-h8b)*C14c;
      return ;
    }
  for (i=0;i<dim2;i++)
    for (j=0;j<dim2;j++)
      for (k=0;k<dim2;k++)
	for (l=0;l<BLK;l++)
	  opa[l][i][j][k] = *(in+adrev[index][l][i][j][k]);
  for (i=0;i<dim2;i++)
    for (j=0;j<dim2;j++)
      for (k=0;k<dim2;k++)
	{
	  h1a=opa[0][i][j][k]+opa[1][i][j][k]; h2a=opa[0][i][j][k]-opa[1][i][j][k];
	  h3a=opa[2][i][j][k]+opa[3][i][j][k]; h4a=opa[2][i][j][k]-opa[3][i][j][k];
	  h5a=opa[4][i][j][k]+opa[5][i][j][k]; h6a=opa[4][i][j][k]-opa[5][i][j][k];
	  h7a=opa[6][i][j][k]+opa[7][i][j][k]; h8a=opa[6][i][j][k]-opa[7][i][j][k];
	  h1b=h1a+h3a; h2b=h2a+h4a;h3b=h1a-h3a;
	  h4b=h2a-h4a; h5b=h5a+h7a;h6b=h6a+h8a; h7b=h5a-h7a; h8b=h6a-h8a;
	  opb[0][i][j][k] = h1b+h5b; opb[1][i][j][k] = h2b+h6b;
	  opb[2][i][j][k] = h3b+h7b; opb[3][i][j][k] = h4b+h8b;
	  opb[4][i][j][k] = h1b-h5b; opb[5][i][j][k] = h2b-h6b;
	  opb[6][i][j][k] = h3b-h7b; opb[7][i][j][k] = h4b-h8b;
	}
  for (i=0;i<dim2;i++)
    for (j=0;j<dim2;j++)
      for (k=0;k<dim2;k++)
	{
	  opa[0][i][j][k]=opb[0][i][j][k];
	  for (l=1;l<BLK;l++)
	    opa[l][i][j][k]=opb[l][i][j][k]*cosmul[index][l][i][j][k];
	}
  for (l=0;l<BLK;l++)
    DCT3(&opa[l][0][0][0],&opb[l][0][0][0],dim2);
  for (i=0;i<dim2;i++)
    for (j=0;j<dim2;j++)
      for (k=0;k<dim2;k++)
	{
	  opa[0][i][j][k]=opb[0][i][j][k];
	  if (k<dim2-1)
	    opa[1][i][j][k]=opb[1][i][j][k]+opb[1][i][j][k+1];
	  else
	    opa[1][i][j][k]=opb[1][i][j][k];
	  if (j<dim2-1)
	    opa[2][i][j][k]=opb[2][i][j][k]+opb[2][i][j+1][k];
	  else
	    opa[2][i][j][k]=opb[2][i][j][k];
	  if (j<dim2-1 && k<dim2-1)
	    opa[3][i][j][k]=opb[3][i][j][k]+opb[3][i][j][k+1]+
	      opb[3][i][j+1][k]+opb[3][i][j+1][k+1];}
  else if (j<dim2-1 && k==dim2-1)
    opa[3][i][j][k]=opb[3][i][j][k]+opb[3][i][j+1][k];
  else if (j==dim2-1 && k<dim2-1)
    opa[3][i][j][k]=opb[3][i][j][k]+opb[3][i][j][k+1];
  else
    opa[3][i][j][k]=opb[3][i][j][k];
  if (i<dim2-1)
    opa[4][i][j][k]=opb[4][i][j][k]+opb[4][i+1][j][k];
  else
    opa[4][i][j][k]=opb[4][i][j][k];
  if (i<dim2-1 && k<dim2-1)
    opa[5][i][j][k]=opb[5][i][j][k]+opb[5][i][j][k+1]+
      opb[5][i+1][j][k]+opb[5][i+1][j][k+1];
  else if (i<dim2-1 && k==dim2-1)
    opa[5][i][j][k]=opb[5][i][j][k]+opb[5][i+1][j][k];
  else if (i==dim2-1 && k<dim2-1)
    opa[5][i][j][k]=opb[5][i][j][k]+opb[5][i][j][k+1];
  else
    opa[5][i][j][k]=opb[5][i][j][k];
  if (i<dim2-1 && j<dim2-1)
    opa[6][i][j][k]=opb[6][i][j][k]+opb[6][i][j+1][k]+
      opb[6][i+1][j][k]+opb[6][i+1][j+1][k];
  else if (i<dim2-1 && j==dim2-1)
    opa[6][i][j][k]=opb[6][i][j][k]+opb[6][i+1][j][k];
  else if (i==dim2-1 && j<dim2-1)
    opa[6][i][j][k]=opb[6][i][j][k]+opb[6][i][j+1][k];
  else
    opa[6][i][j][k]=opb[6][i][j][k];
  if (i<dim2-1 && j<dim2-1 && k<dim2-1)
    opa[7][i][j][k]=opb[7][i][j][k]+opb[7][i][j][k+1]+opb[7][i][j+1][k]+opb[7][i][j+
										  1][k+1]+opb[7][i+1][j][k]+opb[7][i+1][j][k+1]+opb[7][i+1][j+1][k]+opb[7][i+1][j+1][k+1];
  else if (i<dim2-1 && j<dim2-1 && k==dim2-1)
    opa[7][i][j][k]=opb[7][i][j][k]+opb[7][i][j+1][k]+opb[7][i+1][j][k]+opb[7][i+1][j+1][k];
  else if (i<dim2-1 && j==dim2-1 && k<dim2-1)
    opa[7][i][j][k]=opb[7][i][j][k]+opb[7][i][j][k+1]+opb[7][i+1][j][k]+opb[7][i+1][j][k+1];
  else if (i<dim2-1&&j==dim2-1 && k==dim2-1)
    opa[7][i][j][k]=opb[7][i][j][k]+opb[7][i+1][j][k];
  else if (i==dim2-1 && j<dim2-1 && k<dim2-1)
    opa[7][i][j][k]=opb[7][i][j][k]+opb[7][i][j][k+1]+opb[7][i][j+1][k]+opb[7][i][j+1][k+1];
  else if (i==dim2-1 && j<dim2-1 && k==dim2-1)
    opa[7][i][j][k]=opb[7][i][j][k]+opb[7][i][j+1][k];
  else if (i==dim2-1 && j==dim2-1 && k<dim2-1)
    opa[7][i][j][k]=opb[7][i][j][k]+opb[7][i][j][k+1];
  else
    opa[7][i][j][k]=opb[7][i][j][k];
}
for (i=0;i<dim2;i++)
  for (j=0;j<dim2;j++)
    for (k=0;k<dim2;k++)
      for (l=0;l<BLK;l++)
	*(out+admul[index][l][i][j][k])=opa[l][i][j][k];
}



static void dct(float *data,float * vol)
{
  int r,k;
  double spec;
  extern double *tabla1;

  vol[0] = (float) (W1*((double)data[0]*tabla1[0] + (double)data[1]*tabla1[1] + (double)data[2]*tabla1[2] + (double)data[3]*tabla1[3]));   
  for(r=1;r<B;r++) 
    {
      spec= (double)data[0]*tabla1[r*B] + (double)data[1]*tabla1[r*B+1] + (double)data[2]*tabla1[r*B+2] + data[3]*tabla1[r*B+3];
      vol[r] = (float)(W2*spec);
    }
}

static void dct3(float *data,float * vol)
{
  int r1,r2,r3;
  float v2[4];
  float v3[4];    
  float s2[4];
  float s3[4];
  float tmp[64];    
      
  /* do the stuff     */
  for(r3=0;r3<B;r3++) 
    for(r2=0;r2<B;r2++) 
      {      
	dct(&data[r3*BB+r2*B],&tmp[r3*BB+r2*B]);           
      }
    
  for(r1=0;r1<B;r1++) 
    for(r3=0;r3<B;r3++) 
      {      
	v2[0]=tmp[r3*BB+r1]; 
	v2[1]=tmp[r3*BB+B+r1]; 
	v2[2]=tmp[r3*BB+2*B+r1]; 
	v2[3]=tmp[r3*BB+3*B+r1]; 
	dct(v2,s2);
	tmp[r3*BB+r1]=s2[0]; 
	tmp[r3*BB+B+r1]=s2[1]; 
	tmp[r3*BB+2*B+r1]=s2[2]; 
	tmp[r3*BB+3*B+r1]=s2[3];      
      }
    
  for(r1=0;r1<B;r1++) 
    for(r2=0;r2<B;r2++) 
      {      
	v3[0]=tmp[r2*B+r1];                      
	v3[1]=tmp[BB+r2*B+r1];                     
	v3[2]=tmp[2*BB+r2*B+r1];                     
	v3[3]=tmp[3*BB+r2*B+r1];                     
	dct(v3,s3);
	vol[r2*B+r1]=s3[0];                      
	vol[BB+r2*B+r1]=s3[1];                     
	vol[2*BB+r2*B+r1]=s3[2];                     
	vol[3*BB+r2*B+r1]=s3[3];                                
      }
}

static void idct(float *data,float * vol)
{
  int r,k;
  double spec;
  extern double *tabla2;

  /* do the stuff         */
  for(r=0;r<B;r++) 
    {
      vol[r] = (float)( W1*(double)data[0]*tabla2[r*B] + W2*(double)data[1]*tabla2[r*B+1] + W2*(double)data[2]*tabla2[r*B+2] + W2*(double)data[3]*tabla2[r*B+3] );                  
    }
}

static void idct3(float *data,float * vol)
{
  int r1,r2,r3;    
  float v2[4];
  float v3[4];    
  float s2[4];
  float s3[4];
  float tmp[64];    
      
  /* do the stuff     */
  for(r3=0;r3<B;r3++) 
    for(r2=0;r2<B;r2++) 
      {      
	idct(&data[r3*BB+r2*B],&tmp[r3*BB+r2*B]);     
      }
    
  for(r1=0;r1<B;r1++) 
    for(r3=0;r3<B;r3++) 
      {      
	v2[0]=tmp[r3*BB+r1]; 
	v2[1]=tmp[r3*BB+B+r1]; 
	v2[2]=tmp[r3*BB+2*B+r1]; 
	v2[3]=tmp[r3*BB+3*B+r1]; 
	idct(v2,s2);
	tmp[r3*BB+r1]=s2[0]; 
	tmp[r3*BB+B+r1]=s2[1]; 
	tmp[r3*BB+2*B+r1]=s2[2]; 
	tmp[r3*BB+3*B+r1]=s2[3];      
      }
    
  for(r1=0;r1<B;r1++) 
    for(r2=0;r2<B;r2++) 
      {      
	v3[0]=tmp[r2*B+r1];                      
	v3[1]=tmp[BB+r2*B+r1];                     
	v3[2]=tmp[2*BB+r2*B+r1];                     
	v3[3]=tmp[3*BB+r2*B+r1];     
	idct(v3,s3);
	vol[r2*B+r1]=s3[0];                      
	vol[BB+r2*B+r1]=s3[1];                     
	vol[2*BB+r2*B+r1]=s3[2];                     
	vol[3*BB+r2*B+r1]=s3[3];             
      }      
}

void* ThreadFunc( void* pArguments )
{
  float *ima,*fima,*fimaB1,*acu,sigma,T,z,zB1,fac,max;
  int ind,indB1,ii,jj,kk,i,j,k,ini,fin,rows,cols,slices,p,p1;   
  float b[64]; 
  float cb[64]; 
  float b3[64]; 
  float cb3[64];
  myargument arg;
  arg=*(myargument *) pArguments;

  rows=arg.rows;    
  cols=arg.cols;
  slices=arg.slices;
  ini=arg.ini;    
  fin=arg.fin;
  ima=arg.in_image;
  fima=arg.out_image;
  acu=arg.acu_image;    
  sigma=arg.sigma; 
            
  T=2.7*sigma;
  /*ZC for B1 estimate*/
  int thresholdB1 = 1; 
  float TB1 = sigma*0.15;
  fimaB1 = arg.B1_image;

  for(k=ini;k<=fin;k=k+1)
    {  
      for(j=0;j<=rows-B;j=j+1)
	{
	  for(i=0;i<=cols-B;i=i+1)
	    {	            
	      max=0;
	      for(kk=0;kk<B;kk++)
                for(jj=0;jj<B;jj++)	
		  for(ii=0;ii<B;ii++)
                    {	
                      p=kk*B*B+jj*B+ii;             
                      b[p] = ima[(k+kk)*rows*cols+(j+jj)*cols+(i+ii)];
		      b3[p] = fimaB1[(k+kk)*rows*cols+(j+jj)*cols+(i+ii)];
                      if(b[p]>max) max=b[p];
                    }         
	      if(max==0) continue;
                                              
	      dct3(b,cb);
	      dct3(b3,cb3); 
	      ind=0; indB1=0;
	      for(kk=0;kk<B;kk++)
                for(jj=0;jj<B;jj++)	
		  for(ii=0;ii<B;ii++)
                    {
		      p=kk*B*B+jj*B+ii;
		      
		      if(abs(cb[p])<T) cb[p]=0; /* hardthreshold(cb,T);       */
		      else ind++;     
		      if(abs(cb3[p])<TB1) cb3[p]=0; 
		      else indB1++;  
		      /*
			if (kk>thresholdB1 || jj>thresholdB1 || ii>thresholdB1 ) {
			cb3[p]=0;
			}else{ indB1++;}*/ 
                    }
	      z=1.0/(ind+1); zB1=1.0/(indB1+1);

	      idct3(cb,b);
	      idct3(cb3,b3); 

	      for(kk=0;kk<B;kk++)
                for(jj=0;jj<B;jj++)	
		  for(ii=0;ii<B;ii++)
                    {	
                      p1=(k+kk)*rows*cols+(j+jj)*cols+(i+ii);                

                      /* agregation                  */
		      /*fimaB1[p1] += (1.0f/(thresholdB1*thresholdB1*thresholdB1))**/
		      fimaB1[p1] += zB1 * cb3[kk*B*B+jj*B+ii];
                      fima[p1] += zB1 * b3[kk*B*B+jj*B+ii];  
                      acu[p1]+=z;     

                    }  
            
	    }
	} 
    }   
  
#ifdef _WIN32
  _endthreadex(0);    
#else
  pthread_exit(0);    
#endif

  return 0;
} 

void* ThreadFunc2( void* pArguments )
{
  double *mibias;
  float *ima,*fima,*fima2,*fimaB1,*acu,sigma,T,z,fac,val,v1,max;
  int ind,indB1,ii,jj,kk,i,j,k,ini,fin,rows,cols,slices,p,p1,iii;
  float b[64]; 
  float cb[64]; 
  float b2[64]; 
  float cb2[64]; 
  float b3[64]; 
  float cb3[64]; 
  extern int rician;    
 
  myargument arg;
  arg=*(myargument *) pArguments;

  /*ZC for B1 estimate*/
  int thresholdB1 = 1;float TB1=0.5*sigma;
  fimaB1 = arg.B1_image;

  rows=arg.rows;    
  cols=arg.cols;
  slices=arg.slices;
  ini=arg.ini;    
  fin=arg.fin;
  ima=arg.in_image;
  fima=arg.out_image;
  fima2=arg.out_image2;
  acu=arg.acu_image;    
  sigma=arg.sigma; 
  mibias=arg.mibias;
    
  /*filter*/
    
  T=sigma;
  

  for(k=ini;k<=fin;k=k+1)
    {  
      for(j=0;j<=rows-B;j=j+1)
	{
	  for(i=0;i<=cols-B;i=i+1)
	    {	                                   
	      max=0; 
	      for(kk=0;kk<B;kk++)
	        for(jj=0;jj<B;jj++)
		  for(ii=0;ii<B;ii++)
		    {	
		      p=kk*B*B+jj*B+ii;
		      p1=(k+kk)*rows*cols+(j+jj)*cols+(i+ii);
		      b[p] = ima[p1];                     
		      b2[p] = fima[p1];   
		      b3[p] = fimaB1[p1];
		      if(b[p]>max) max=b[p];
		    } 
	      if(max==0) continue;
            
	      dct3(b,cb);
	      dct3(b2,cb2);
	      dct3(b3,cb3); /*ZC*/

	      /*/////////////////////////////////////////////////*/
                       
	      ind=0;           
	      for(kk=0;kk<B;kk++)
	        for(jj=0;jj<B;jj++)	
		  for(ii=0;ii<B;ii++)
		    {
		      p=kk*B*B+jj*B+ii;
		      if(abs(cb2[p])<T) cb[p]=0; /* hardthreshold(cb,T); */
		      else ind++;  
		      if(abs(cb3[p])<TB1) cb3[p]=0; 
		      else indB1++; 
		    }
	      z=1.0/(ind+1); 

	      /*ZC DCT threshold at rank 2 for B1*/
	      /*for(kk=thresholdB1+1;kk<B;kk++)
	        for(jj=thresholdB1+1;jj<B;jj++)	
		for(ii=thresholdB1+1;ii<B;ii++)
		{
		p=kk*B*B+jj*B+ii;              
		cb3[p]=0;                
		}
	      */
	      idct3(cb,b);
	      idct3(cb3,b3); /*ZC*/
                                  
	      for(kk=0;kk<B;kk++)
	        for(jj=0;jj<B;jj++)	
		  for(ii=0;ii<B;ii++)
		    {	
		      p1=(k+kk)*rows*cols+(j+jj)*cols+(i+ii);  
		      p=kk*B*B+jj*B+ii;
		      /*		      fimaB1[p1]+=(1/(indB1+1))*b3[p]; /*ZC  thresholdB1*thresholdB1*thresholdB1*/
		      if(rician>0)
			{                 
			  iii=(int)(b[p]);                 
			  val=mibias[iii];
			  if(val>2*b[p]) val=b[p];                
			  fima2[p1]+=z*val;
			}
              
		      else 
			{
			  fima2[p1]+=z*b[p];  /* agregation     */
			}
		      acu[p1]+=z;                                               
		    }
	    }
	} 
    }
 
#ifdef _WIN32
  _endthreadex(0);    
#else
  pthread_exit(0);    
#endif

  return 0;
} 


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*Declarations*/
  mxArray *xData,*pv;
  double *mibias;
  float *ima,*fima,*fima2,*acu,*fimaB1;
  float sigma,val,v1,max;
  int i,ndim,Nthreads,ini,fin,N,r,k,ind,imax;
  const int *dims;
  mwSize *bdims;

  myargument *ThreadArgs; 
 
#ifdef _WIN32
  HANDLE *ThreadList; /* Handles to the worker threads*/
#else
  pthread_t * ThreadList;
  int *ThreadIndex;
#endif

  mexPrintf("cM_ODCT3D: \n");
  if (nrhs != 3) {
    mexErrMsgIdAndTxt("MyODCT3d:nrhs",
                      "Two inputs required.");
    return;
  }
  if (nlhs > 5 || nlhs < 1) {
    mexErrMsgIdAndTxt("MyODCT3d:nlhs",
                      "One to Five outputs required.");
    return;
  }
  /* Parse input args */
  /* make sure the first input argument is a matrix */
  if ( mxIsSparse(prhs[0]) || 
       mxIsDouble(prhs[0]) || 
       mxIsComplex(prhs[0]) ||
       mxGetNumberOfElements(prhs[0]) == 1 ||
       mxGetNumberOfDimensions(prhs[0]) != 3) {
    mexErrMsgIdAndTxt("MyODCT3d:not3DImage",
                      "Input image must be a float matrix with 3 dimensions.");
  } 
  
  ima = (float*)mxGetPr(prhs[0]);
  ndim = mxGetNumberOfDimensions(prhs[0]);
  dims= mxGetDimensions(prhs[0]);
  
  mexPrintf("cM_ODCT3D: ndims %d\n",ndim);

  /* Copy input parameters */
  if( mxIsComplex(prhs[1]) ||
      mxGetNumberOfElements(prhs[1]) != 1 ) {
    mexErrMsgIdAndTxt("MyODCT3d:notScalar",
                      "Input sigma must be a scalar.");
  }  
  sigma = (float)(mxGetScalar(prhs[1]));
  mexPrintf("cM_ODCT3D: sigma %f\n",sigma);
  

  if( mxIsComplex(prhs[2]) ||
      mxGetNumberOfElements(prhs[2]) != 1 ) {
    mexErrMsgIdAndTxt("MyODCT3d:notIntScalar",
                      "Input rician must be a scalar integer.");
  } 
  rician = (int)(mxGetScalar(prhs[2]));
  mexPrintf("cM_ODCT3D: rician %d\n",rician);
  

  /* Allocate memory and assign output pointer */
  plhs[0] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
  fima = (float*) mxGetPr(plhs[0]);
  if (nlhs > 1) { 
    plhs[1] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
    fimaB1 = (float*) mxGetPr(plhs[1]); 
  }else{
    xData = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
    fimaB1 = (float*) mxGetPr(xData);
  }
  if (nlhs > 2) {
    plhs[2] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
    acu = (float*) mxGetPr(plhs[2]);
  }else{
    xData = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
    acu = (float*) mxGetPr(xData);
  }
  if (nlhs > 3) {
    plhs[3] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
    fima2 = (float*) mxGetPr( plhs[3]);
  }else{
    xData = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
    fima2 = (float*) mxGetPr(xData);
  }
  max=0;
  for(i=0;i<dims[0]*dims[1]*dims[2];i++)
    {        
      fima[i]=0;    
      acu[i]=0;
      fimaB1[i]=ima[i]; 
      if(ima[i]>max) max=ima[i];
    }					
  imax=(int)ceil(2*max);

  Nthreads = 16;/*floor((float)dims[2]/(float)B);*/

  N=B;
  tabla1 = (double*) mxMalloc(N*N*sizeof(double));
  tabla2 = (double*) mxMalloc(N*N*sizeof(double));

  for(r=0;r<N;r++) 
    for(k=0;k<N;k++) 
      tabla1[r*N + k] = cos((double)(PI*(2*k+1)*r)/(double)(2*N));               
  for(r=0;r<N;r++) 
    for(k=0;k<N;k++) 
      tabla2[r*N + k] = cos((double)(PI*(2*r+1)*k)/(double)(2*N));      
  W1=sqrt(1/(double)N);
  W2=sqrt(2/(double)N);	


  /*mibias=mxMalloc(max*sizeof(float));*/
  bdims = (mwSize *) mxMalloc (2 * sizeof(mwSize));
  bdims[0] = 1;
  bdims[1] = imax*Nthreads;
  if (nlhs > 4){
    plhs[4] = mxCreateNumericArray (2, bdims, mxDOUBLE_CLASS, mxREAL);
    mibias = (double*)mxGetPr(plhs[4]);
  }else{
    pv = mxCreateNumericArray (2, bdims, mxDOUBLE_CLASS, mxREAL);
    mibias = (double*)mxGetPr(pv);
  }
  mexPrintf("cM_ODCT3D: Rician Bias \n");

  if(rician>0) RicianBias(sigma,imax,mibias,Nthreads);


#ifdef _WIN32

  /* Reserve room for handles of threads in ThreadList*/
  ThreadList = (HANDLE*) mxMalloc(Nthreads*sizeof( HANDLE ));
  ThreadArgs = (myargument*) mxMalloc(Nthreads*sizeof(myargument));

  for (i=0; i<Nthreads; i++)
    {       
      /* Make Thread Structure*/
      ini=ceil((i*(dims[2]-B))/Nthreads);
      if(i==Nthreads-1) fin=floor(((i+1)*(dims[2]-B))/Nthreads);
      else fin=floor(((i+1)*(dims[2]-B))/Nthreads)-1;   
       
      ThreadArgs[i].rows=dims[1];
      ThreadArgs[i].cols=dims[0];
      ThreadArgs[i].slices=dims[2];
      ThreadArgs[i].in_image=ima;
      ThreadArgs[i].out_image=fima;
      ThreadArgs[i].acu_image=acu;
      ThreadArgs[i].B1_image=fimaB1;
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;    
      ThreadArgs[i].sigma=sigma;        
	
      ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &ThreadFunc, &ThreadArgs[i] , 0, NULL );    
    }
    
  for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
  for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }

  for(i=0;i<dims[0]*dims[1]*dims[2];i++)
    {
      if(acu[i]==0) acu[i]=1; 
      fima[i]/=acu[i];
    }

  /*//////////////////////////////////////////////*/
  /* second pass*/
  /*//////////////////////////////////////////////*/

  for(i=0;i<dims[0]*dims[1]*dims[2];i++)
    {        
      fima2[i]=0;    
      acu[i]=0;
    }					
  for (i=0; i<Nthreads; i++)
    {       
      /* Make Thread Structure*/
      ini=ceil((i*(dims[2]-B))/Nthreads);
      if(i==Nthreads-1) fin=floor(((i+1)*(dims[2]-B))/Nthreads);
      else fin=floor(((i+1)*(dims[2]-B))/Nthreads)-1;    
       
      ThreadArgs[i].rows=dims[1];
      ThreadArgs[i].cols=dims[0];
      ThreadArgs[i].slices=dims[2];
      ThreadArgs[i].in_image=ima;
      ThreadArgs[i].out_image=fima;
      ThreadArgs[i].out_image2=fima2;
      ThreadArgs[i].acu_image=acu;
      ThreadArgs[i].B1_image=fimaB1;
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;   
      ThreadArgs[i].sigma=sigma;      
      ThreadArgs[i].mibias=&mibias[i*imax];      
	
      ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &ThreadFunc2, &ThreadArgs[i] , 0, NULL );    
    }
    
  for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
  for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }


#else
  mexPrintf("cM_ODCT3D: Pthreading \n");

  /* 
     Reserve room for handles of threads in ThreadList
  */
  ThreadList = (pthread_t *) calloc(Nthreads,sizeof(pthread_t));
  ThreadArgs = (myargument*) calloc( Nthreads,sizeof(myargument));
  ThreadIndex = (int *)malloc(Nthreads* sizeof( int ));
    
  for (i=0; i<Nthreads; i++)
    {       
      /* 
	 Make Thread Structure 
      */
      ini=ceil((i*(dims[2]-B))/Nthreads);
      if(i==Nthreads-1) fin=floor(((i+1)*(dims[2]-B))/Nthreads);
      else fin=floor(((i+1)*(dims[2]-B))/Nthreads)-1;   
       
      ThreadArgs[i].rows=dims[1];
      ThreadArgs[i].cols=dims[0];
      ThreadArgs[i].slices=dims[2];
      ThreadArgs[i].in_image=ima;
      ThreadArgs[i].out_image=fima;
      ThreadArgs[i].acu_image=acu;
      ThreadArgs[i].B1_image=fimaB1;
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;    
      ThreadArgs[i].sigma=sigma;        
    }

  for (i=0; i<Nthreads; i++)
    {
      if(ThreadIndex[i] = pthread_create(&ThreadList[i], NULL, ThreadFunc, (void*) &ThreadArgs[i]))
        /*pthread_create( &ThreadList[i], NULL, ThreadFunc, (void*) &ThreadArgs[i]);*/
	{
	  printf("Threads cannot be created\n");
	  exit(1);
	}        
    }
   
  for (i=0; i<Nthreads; i++)
    {
      pthread_join(ThreadList[i],NULL);
    }

  for(i=0;i<dims[0]*dims[1]*dims[2];i++)
    {
      if(acu[i]==0) acu[i]=1; 
      fima[i]/=acu[i];
    }

  /*//////////////////////////////////////////////*/
  /* second pass*/
  /*//////////////////////////////////////////////*/
  mexPrintf("cM_ODCT3D: Second pass \n");
  for(i=0;i<dims[0]*dims[1]*dims[2];i++)
    {        
      fima2[i]=0;    
      acu[i]=0;
    }	
				
  for (i=0; i<Nthreads; i++)
    {       
      /* Make Thread Structure*/
      ini=ceil((i*(dims[2]-B))/Nthreads);
      if(i==Nthreads-1) fin=floor(((i+1)*(dims[2]-B))/Nthreads);
      else fin=floor(((i+1)*(dims[2]-B))/Nthreads)-1;    
       
      ThreadArgs[i].rows=dims[1];
      ThreadArgs[i].cols=dims[0];
      ThreadArgs[i].slices=dims[2];
      ThreadArgs[i].in_image=ima;
      ThreadArgs[i].out_image=fima;
      ThreadArgs[i].out_image2=fima2;
      ThreadArgs[i].acu_image=acu;
      ThreadArgs[i].B1_image=fimaB1;
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;   
      ThreadArgs[i].sigma=sigma;      
      ThreadArgs[i].mibias=&mibias[i*imax];     
    }

	
  for (i=0; i<Nthreads; i++)
    {
      if(ThreadIndex[i] = pthread_create(&ThreadList[i], NULL, ThreadFunc2,&ThreadArgs[i]))
	{
	  printf("Threads cannot be created\n");
	  exit(1);
	}        
    }
   
  for (i=0; i<Nthreads; i++)
    {
      pthread_join(ThreadList[i],NULL);
    }

  mexPrintf("cM_ODCT3D: Threads Done \n");
#endif

 
  for(i=0;i<dims[0]*dims[1]*dims[2];i++) 
    {
      if(acu[i]==0) acu[i]=1;
      fima[i]=fima2[i]/acu[i];   
      if(rician) if(fima[i]<0) fima[i]=0; /* background*/
    }
  mexPrintf("cM_ODCT3D: arrays \n");
  mxFree(tabla1);
  mxFree(tabla2);   
  mxFree(bdims);
  mexPrintf("cM_ODCT3D: Freeing threads \n");



#ifdef _WIN32
  mxFree(ThreadArgs); 
  mxFree(ThreadList);
#else
  free(ThreadArgs); 
  free(ThreadList);

#endif
  mexPrintf("cM_ODCT3D: Done \n");

  return;

}

