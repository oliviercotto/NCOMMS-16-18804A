/**  $Id: Uniform.h,v 1.9.2.3 2014-05-01 07:01:43 fred Exp $
*
*  @file Uniform.h
*  Nemo2
*
*   Copyright (C) 2006-2011 Frederic Guillaume
*   frederic.guillaume@env.ethz.ch
*
*   This file is part of Nemo
*
*   Nemo is free software; you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation; either version 2 of the License, or
*   (at your option) any later version.
*
*   Nemo is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program; if not, write to the Free Software
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
*  @author fred
*/

#ifndef __UNIFORM_H
#define __UNIFORM_H
#include <cmath>
#include <limits>
#include <assert.h>
#include "output.h"
#include <iostream>

#ifdef HAS_SPRNG
 #define SIMPLE_SPRNG
 #include <sprng.h>
#endif

#ifdef HAS_GSL
 #include <gsl/gsl_rng.h>
 #include <gsl/gsl_randist.h>
 #include <gsl/gsl_blas.h>
 #include <gsl/gsl_math.h>
 #include <gsl/gsl_vector.h>
 #include <gsl/gsl_matrix.h>
 #include <gsl/gsl_permutation.h>
#endif

/**Random number generation class, uses various types of random generators depending on the implementation.*/
class RAND {
private:
  
  RAND();
  
public:
  
#if defined(HAS_GSL) && !defined(HAS_SPRNG)

  static gsl_rng * r;
  
  static void init_gsl(const gsl_rng_type* T, unsigned long seed) 
  {
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,seed);
  }
  
  static void free() 
  {
    gsl_rng_free(r);
  }

#elif !defined(HAS_GSL) && !defined(HAS_SPRNG)

  static long Seed1,Seed2;
  
#endif
  /**Initialize the random generator's seed.*/
  static void init (unsigned long seed) {

#if defined(HAS_SPRNG)
    
    init_sprng( SPRNG_LFG, seed, SPRNG_DEFAULT );
  
#elif defined(HAS_GSL)
    
    init_gsl( gsl_rng_mt19937, seed );
  
#else
    
    Seed1 = seed;
    
#endif
  }
  /**Generates a random number from [0.0, 1.0[ uniformly distributed.
   * If SPRNG or GSL libraries are not used, implement a random generator from:
   * L'Ecuyer, 1988, "Efficient and Portable Combined Random Number Generators",
   *                      Communication of the ACM, 31(6):742-774.
    **/
  static inline double Uniform () {

#ifdef HAS_SPRNG
    return sprng();
#elif defined(HAS_GSL)
    return gsl_rng_uniform(r);
#else
    register long z, w;

    do{
      w = Seed1 / 53668;

      Seed1 = 40014 * (Seed1 - w * 53668) - w * 12211;

      if (Seed1 < 0) Seed1 += 2147483563;

      w = (Seed2 / 52774);

      Seed2 = 40692 * (Seed2 - w * 52774) - w * 3791;

      if (Seed2 < 0) Seed2 += 2147483399;

      z = Seed1 - Seed2;

      if (z < 1) z += 2147483562;

    }while (!((z * 4.656613e-10) < 1.0));

    return (z * 4.656613e-10);
#endif
  }
  /**Returns a uniformly distributed random number from [0.0, max[.*/
  static inline unsigned int Uniform (unsigned int max) 
  {
    return (unsigned int)(Uniform() * max);
  }
  
  /**Returns a random boolean.*/
  static inline bool RandBool() {
    //generate a random int (or long)
#ifdef HAS_SPRNG
    static int intrand = isprng();
#elif defined(HAS_GSL)
    static unsigned long int intrand =  gsl_rng_get(r);
#else
    static int intrand = (int)(Uniform()*(double)std::numeric_limits<int>::max());
#endif
    //read up to the first 16 bits
    static unsigned int num = 0;
    //redraw an number after that limit
    if ( ++num > 16 ) {
      num = 0;

#ifdef HAS_SPRNG
      intrand = isprng();
#elif defined(HAS_GSL)
      intrand =  gsl_rng_get(r);
#else
      intrand = (int)(Uniform()*(double)std::numeric_limits<int>::max());
#endif

    }
    return (intrand & (1 << num) );

  }
  
  /**Return a random unsigned long, from uniform distribution.*/
  static inline unsigned long RandULong() {
#if defined(HAS_GSL) && !defined(HAS_SPRNG)
    return gsl_rng_get(r);
#else
    unsigned long rnd, limit = 0x10000000; //=2^28 this gives ~7% redraws
    do{
      rnd = (unsigned long)(Uniform()*ULONG_MAX);
    }while(rnd < limit);

    return rnd;
#endif
  }
  /**From the Numerical Recieps.*/
  static inline double gammln (double xx) {
    double x,y,tmp,ser=1.000000000190015;
    static double cof[6]={76.18009172947146,-86.50532032941677,
      24.01409824083091,-1.231739572450155,
      0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    for (j = 0; j < 6; ++j) ser += cof[j]/++y;

    return -tmp+log(2.5066282746310005*ser/x);

  }
 /**From the Numerical Recieps.*/
  static inline double Poisson (double mean) {
    
#if defined(HAS_GSL) && !defined(HAS_SPRNG)
    return gsl_ran_poisson(r, mean);
#else
    static double sq,alxm,g,oldm=(-1.0);
    double em,t,y;

    if (mean < 12.0)
    {
      if (mean != oldm){
        oldm=mean;
        g=exp(-mean);
      }
      em = -1;
      t=1.0;
      do {
        ++em;
        t *= Uniform();
      } while (t > g);
    }else
    {
      if (mean != oldm)
      {
        oldm=mean;
        sq=sqrt(2.0*mean);
        alxm=log(mean);
        g=mean*alxm-gammln(mean+1.0);
      }
      do {
        do {
          y=tan(M_PI*Uniform());
          em=sq*y+mean;
        } while (em < 0.0);
        em=floor(em);
        t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (Uniform() > t);
    }
    return em;
#endif
  }
  
  static inline double Gaussian(double sigma)
  {
#if defined(HAS_GSL) && !defined(HAS_SPRNG)
    
    return gsl_ran_gaussian_ratio_method (r, sigma);
    
#else
  /**From the GSL.*/
//    double x, y, r2;
//    
//    do
//    {
//      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
//      
//      x = -1 + 2 * Uniform ( );
//      y = -1 + 2 * Uniform ( );
//      
//      /* see if it is in the unit circle */
//      r2 = x * x + y * y;
//    }
//    while (r2 > 1.0 || r2 == 0);
//    
//    /* Box-Muller transform */
//    return sigma * y * sqrt (-2.0 * log (r2) / r2);
    
  /*trying the ratio method from GSL*/
  /* see code in gsl/randist/gauss.c */  
    double u, v, x, y, Q;
    const double s = 0.449871;    /* Constants from Leva */
    const double t = -0.386595;
    const double a = 0.19600;
    const double b = 0.25472;
    const double r1 = 0.27597;
    const double r2 = 0.27846;
    
    do                            /* This loop is executed 1.369 times on average  */
    {
      /* Generate a point P = (u, v) uniform in a rectangle enclosing
       the K+M region v^2 <= - 4 u^2 log(u). */
      
      /* u in (0, 1] to avoid singularity at u = 0 */
      u = 1 - RAND::Uniform();
      
      /* v is in the asymmetric interval [-0.5, 0.5).  However v = -0.5
       is rejected in the last part of the while clause.  The
       resulting normal deviate is strictly symmetric about 0
       (provided that v is symmetric once v = -0.5 is excluded). */
      v = RAND::Uniform() - 0.5;
      
      /* Constant 1.7156 > sqrt(8/e) (for accuracy); but not by too
       much (for efficiency). */
      v *= 1.7156;
      
      /* Compute Leva's quadratic form Q */
      x = u - s;
      y = fabs (v) - t;
      Q = x * x + y * (a * y - b * x);
      
      /* Accept P if Q < r1 (Leva) */
      /* Reject P if Q > r2 (Leva) */
      /* Accept if v^2 <= -4 u^2 log(u) (K+M) */
      /* This final test is executed 0.012 times on average. */
    }
    while (Q >= r1 && (Q > r2 || v * v > -4 * u * u * log (u)));
    
    return sigma * (v / u);
#endif
  }
  
  static inline void BivariateGaussian(double sigma1,double sigma2, double rho, double *out1,double *out2) {
#if defined(HAS_GSL) && !defined(HAS_SPRNG)
    gsl_ran_bivariate_gaussian(r,sigma1,sigma2,rho,out1,out2);
#else
    //gsl code:
    double u, v, r2, scale;
   //double *x = out, *y = (++out);
    do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
      
      u = -1 + 2 * Uniform ();
      v = -1 + 2 * Uniform ();
      
      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
    while (r2 > 1.0 || r2 == 0);
    
    scale = sqrt (-2.0 * log (r2) / r2);
    
    *out1 = sigma1 * u * scale;
    *out2 = sigma2 * (rho * u + sqrt(1 - rho*rho) * v) * scale;
    
#endif
  }
  
  static inline double LogNormal (double zeta, double sigma) {
#if defined(HAS_GSL) && !defined(HAS_SPRNG)
    return gsl_ran_lognormal(r,zeta,sigma);
#else
    //this is the GSL code:
    double u, v, r2, normal, z;
    
    do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
      
      u = -1 + 2 * Uniform();
      v = -1 + 2 * Uniform();
      
      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
    while (r2 > 1.0 || r2 == 0);
    
    normal = u * sqrt (-2.0 * log (r2) / r2);
    
    z =  exp (sigma * normal + zeta);
    
    return z;
    
#endif
  }
  
  static inline double Gamma (double a, double b) {
#if defined(HAS_GSL) && !defined(HAS_SPRNG)
    return gsl_ran_gamma(r, a, b);
#else
    /*pasting GSL code in here*/
    
    /* assume a > 0 */
    
    if (a < 1)
    {
      double u = Uniform();//might give singularity at 0....
      return Gamma(1.0 + a, b) * pow (u, 1.0 / a);
    }
    
    {
      double x, v, u;
      double d = a - 1.0 / 3.0;
      double c = (1.0 / 3.0) / sqrt (d);
      
      while (1)
      {
        do
        {
          x = Gaussian(1.0); //should use the Ziggurat method?
          v = 1.0 + c * x;
        }
        while (v <= 0);
        
        v = v * v * v;
        u = Uniform();//might give singularity at 0....
        
        if (u < 1 - 0.0331 * x * x * x * x) 
          break;
        
        if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
          break;
      }
      
      return b * d * v;
    }
    
#endif
  }
  
  static inline double Bernoulli (double p)
  {
#if defined(HAS_GSL) && !defined(HAS_SPRNG)
    return gsl_ran_bernoulli(r, p); //supposed to return 1 with proba p
#else    
    double u = RAND::Uniform() ;
    
    if (u < p)
    {
      return 1 ;
    }
    else
    {
      return 0 ;
    }
#endif
  }
  // renaming by Olivier
  static inline double Bernoulli2 (double p) {return Bernoulli(p);}

  static inline double Exponential (double mu) {
    return -mu * log(RAND::Uniform());
  }
    
#ifdef HAS_GSL
  static inline void MultivariateGaussian(gsl_vector *eval, gsl_matrix *evec,
                                          gsl_vector *workspace, gsl_vector *px)
  {
    size_t i;
    
    for (i=0; i<eval->size; i++)
      gsl_vector_set (workspace, i, 
					  Gaussian(gsl_vector_get (eval, i))); 
    
    gsl_blas_dgemv (CblasNoTrans, 1.0,
                    evec, workspace, 0.0, px);   /* px = evec * px */
                    
  }

  static inline void MultivariateGaussianCholesky(gsl_vector *sigma, gsl_matrix *M, gsl_vector *px)
  {//generate n i.i.d. random normal variates
    for (size_t i = 0; i < sigma->size; i++)
      gsl_vector_set (px, i, Gaussian(gsl_vector_get(sigma, i)));

    gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit, M, px);
  }

  static inline long double factoriel(unsigned int x)   {
    long double f = 1;
    for(unsigned int i=1;i<=x;++i)
      f *= i;
    return f;
  }

//  static inline double Binomial (double p,unsigned int k,long double N,
//                                const unsigned int n) {
//    register long double bincoeff = N/(factoriel(k)*factoriel(n - k));
//    return (bincoeff*pow(p,(int)k)*pow(1-p,(int)(n - k)));
//  }
  
  static inline double Binomial(double p, unsigned int n)
  {
    /*implements Knuth method*/
    unsigned int i, a, b, k = 0;
    
    while (n > 10)        /* This parameter is tunable */
    {
      double X;
      a = 1 + (n / 2);
      b = 1 + n - a;
      
      X = RAND::Beta((double) a, (double) b);
      
      if (X >= p)
      {
        n = a - 1;
        p /= X;
      }
      else
      {
        k += a;
        n = b - 1;
        p = (p - X) / (1 - X);
      }
    }
    
    for (i = 0; i < n; i++)
    {
      double u = RAND::Uniform();
      if (u < p)
        k++;
    }
    
    return k;    
    
  }
  
  static inline double Beta (const double a, const double b)
  {
    /*from Knuth*/
    double x1 = RAND::Gamma(a, 1.0);
    double x2 = RAND::Gamma(b, 1.0);
    
    return x1 / (x1 + x2);
  }
#endif
  
  static inline unsigned int Binomial2 (double p, unsigned int n){
    
#if defined(HAS_GSL) && !defined(HAS_SPRNG)
    return gsl_ran_binomial (r,  p,  n);
#else
    
    cout<<"Problem: need GSL for Binomial"<<endl;
    return n;
#endif
  }
  
  static inline void SampleSeq(int from, int to, int by, unsigned int num, int* result, bool replace = false)
  {
    assert(from < to && by < to && by > 0);
    
    unsigned int seq_length = (int)((to - from) / by); //don't include last (to)
    
    assert(num <= seq_length);
    
    int *seq = new int [seq_length];
    
    seq[0] = from;
    
    for(unsigned int i = 1; seq[i-1] + by < to && i < seq_length; ++i) seq[i] = seq[i-1] + by;
     
    if(!replace) { //without replacement
      
      unsigned int size = seq_length, pos, last = seq_length - 1;
      
      for (unsigned int i = 0; i < num; i++) {
        pos = RAND::Uniform(size);
        result[i] = seq[pos];
        seq[pos] = seq[last];
        seq[last] = result[i];
        size--; last--;
      }
      
    } else { //with replacement
      for (unsigned int i = 0; i < num; i++) result[i] = seq[ RAND::Uniform(seq_length) ];
    }

  }
  
  static inline void SampleSeqWithReciprocal(int from, int to, int by, unsigned int num1, int* result1, unsigned int num2, int* result2)
  {
    assert(from < to && by < to && by > 0);
    
    unsigned int seq_length = (int)((to - from) / by); //don't include last (to)
    
    assert(num1 + num2 == seq_length);
    
    int *seq = new int [seq_length];
    
    seq[0] = from;
    
    for(unsigned int i = 1; seq[i-1] + by < to && i < seq_length; ++i) seq[i] = seq[i-1] + by;

    unsigned int size = seq_length, pos, last = seq_length - 1;
    
    for (unsigned int i = 0; i < num1; i++) {
      pos = RAND::Uniform(size);
      result1[i] = seq[pos];
      seq[pos] = seq[last];
      seq[last] = result1[i];
      size--; last--;
    }
    
    for (unsigned int i = 0; i < num2 && i < size; i++)
      result2[i] = seq[i];
  }
  
  };

#endif
