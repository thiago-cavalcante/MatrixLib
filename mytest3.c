#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#include "seahorn/seahorn.h"

#define MAX_SIZE 10

#define MACHEPS 2.22045e-16

#define SQRT_MAGIC_F 0x5f3759df

#define DEBUG

#ifdef DEBUG
#define	m_output(mat) m_foutput(stdout, mat)
#endif

static const char *format = "%14.9g ";

#define	v_chk_idx(x, dim, i) ((i)>=0 && (i)<dim)

#define v_get_val(x, dim, i) (v_chk_idx(x, dim, i) ? (x)[(i)] : (0))

#define	m_chk_idx(A, A_m, A_n, i, j) ((i)>=0 && (i)<A_m && (j)>=0 && (j)<=A_n)

#define	m_get_val(A, A_m, A_n, i, j) (m_chk_idx(A, A_m, A_n, i, j) ? \
    (A)[(i)][(j)] : (0))

#define	m_entry(A, A_m, A_n, i, j) m_get_val(A, A_m, A_n, i, j)

#ifdef DEBUG
#define printfc(c) printf("%f%c%fi\n", c.real, (c.imag>=0.0f)? '+':'\0', c.imag)
#endif

/* type independent min and max operations */
#ifndef max
#define	max(a, b) ((a) > (b) ? (a) : (b))
#endif /* max */
#ifndef min
#define	min(a, b) ((a) > (b) ? (b) : (a))
#endif /* min */

/* complex number definition */
typedef	struct
{
  double real, imag;
}
CMPLX;

/* peak value definition */
typedef	struct
{
  double mp;
  int kp;
}
PKVL;

/* last state vector definition */
typedef	struct
{
  double xk[MAX_SIZE][MAX_SIZE];
  int lastState;
}
LST;

/* Global variables */
// xk stores the last computed state of the system
LST xk;

/******************************math.h Functions******************************/

/* sp_fabs -- absolute value of floating-point number */
double sp_fabs(double n)
{
  if(n >= 0)
    return n; //if positive, return without ant change
  else
    return (-n); //if negative, return a positive version
}

/* sp_ceil -- the smallest integer value greater than or equal to x */
double sp_ceil(double x)
{
  union
  {
    float f;
    int i;
  }float_int;
  float_int.f = x;
  // Extract sign, exponent and mantissa
  // Bias is removed from exponent
  int sign=float_int.i >> 31;
  int exponent=((float_int.i & 0x7fffffff) >> 23) - 127;
  int mantissa=float_int.i & 0x7fffff;
  // Is the exponent less than zero?
  if(exponent < 0)
  {
    // In this case, x is in the open interval (-1, 1)
    if(x <= 0.0f)
      return 0.0f;
    else
      return 1.0f;
  }
  else
  {
    // Construct a bit mask that will mask off the
    // fractional part of the mantissa
    int mask = 0x7fffff >> exponent;
    // Is x already an integer (i.e. are all the
    // fractional bits zero?)
    if((mantissa & mask) == 0)
      return x;
    else
    {
      // If x is positive, we need to add 1 to it
      // before clearing the fractional bits
      if(!sign)
      {
        mantissa += 1 << (23-exponent);

        // Did the mantissa overflow?
        if(mantissa & 0x800000)
        {
          // The mantissa can only overflow if all the
          // integer bits were previously 1 -- so we can
          // just clear out the mantissa and increment
          // the exponent
          mantissa = 0;
          exponent++;
        }
      }
      // Clear the fractional bits
      mantissa &= ~mask;
    }
  }
  // Put sign, exponent and mantissa together again
  float_int.i = (sign << 31) | ((exponent+127) << 23) | mantissa;
  return (double)float_int.f;
}

/* sp_pow -- returns a raised to the power of n i.e. a^n */
double sp_pow(double a, int n)
{
  double r = 1;
  while(n > 0)
  {
    if(n & 1)
      r *= a;
    a *= a;
    n >>= 1;
  }
  return r;
}

/**
 * Calculate ln logarithm using integers with 16 bit precision
 * min: sp_fxp_ln(0.000015259<<16)
 * max: sp_fxp_ln(32767<<16)
 */
int sp_fxp_ln(int x)
{
  int t, y;
  y = 0xa65af;
  if(x < 0x00008000)
    x <<= 16, y -= 0xb1721;
  if(x < 0x00800000)
    x <<= 8, y -= 0x58b91;
  if(x < 0x08000000)
    x <<= 4, y -= 0x2c5c8;
  if(x < 0x20000000)
    x <<= 2, y -= 0x162e4;
  if(x < 0x40000000)
    x <<= 1, y -= 0x0b172;
  t = x + (x >> 1);
  if((t & 0x80000000) == 0)
    x = t, y -= 0x067cd;
  t = x + (x >> 2);
  if((t & 0x80000000) == 0)
    x = t, y -= 0x03920;
  t = x + (x >> 3);
  if((t & 0x80000000) == 0)
    x = t, y -= 0x01e27;
  t = x + (x >> 4);
  if((t & 0x80000000) == 0)
    x = t, y -= 0x00f85;
  t = x + (x >> 5);
  if((t & 0x80000000) == 0)
    x = t, y -= 0x007e1;
  t = x + (x >> 6);
  if((t & 0x80000000) == 0)
    x = t, y -= 0x003f8;
  t = x + (x >> 7);
  if((t & 0x80000000) == 0)
    x = t, y -= 0x001fe;
  x = 0x80000000 - x;
  y -= x >> 15;
  return y;
}

/**
 * Calculate log10 logarithm using 16 bit precision
 * min: sp_log10_2(0.000015259)
 * max: sp_log10_2(32767.0)
 */
double sp_log10_2(double x)
{
  int xint = (int) (x * 65536.0 + 0.5);
  int lnum = sp_fxp_ln(xint);
  int lden = sp_fxp_ln(655360);
  return ((double) lnum / (double) lden);
}

/* sp_floor -- returns the largest integer value less than or equal to num */
double sp_floor(double num)
{
  long long n;
  double d;
  if(num >= 9.2234e+18 || num <= -9.2234e+18 || num != num)
  {
  /* handle large values, infinities and nan */
    return num;
  }
  n = (long long)num;
  d = (double)n;
  if (d == num || num >= 0)
    return d;
  else
    return (d - 1);
}

///* sp_sqrt -- returns the square root of fg */
//double sp_sqrt(double fg)
//{
//  double n, lstX;
//  n = fg / 2;
//  lstX = 0;
//  while(n != lstX)
//  {
//    lstX = n;
//    n = (n + (fg/n))/2;
//  }
//  return n;
//}

/* sp_sqrt -- returns the square root of fg */
double sp_sqrt(double x)
{
	union {
		int i;
		float x;
	} u;
	u.x = x;
	u.i = (1 << 29) + (u.i >> 1) - (1 << 22);

	// Two Babylonian Steps (simplified from:)
	// u.x = 0.5f * (u.x + x/u.x);
	// u.x = 0.5f * (u.x + x/u.x);
	u.x = u.x + x / u.x;
	u.x = 0.25f * u.x + x / u.x;

	return (double)u.x;
}

///* sp_sqrt -- returns the square root of fg */
//double sp_sqrt(double x)
//{
//	float xhalf = 0.5f*x;
//
//	  union // get bits for floating value
//	  {
//	    float x;
//	    int i;
//	  } u;
//	  u.x = x;
//	  u.i = SQRT_MAGIC_F - (u.i >> 1);  // gives initial guess y0
//	  return (double) (x*u.x*(1.5f - xhalf*u.x*u.x));// Newton step, repeating increases accuracy
//}

/******************************Matrix Functions******************************/

/* m_add -- matrix addition -- may be in-situ */
void m_add(double m1[MAX_SIZE][MAX_SIZE], double m2[MAX_SIZE][MAX_SIZE],
           double m3[MAX_SIZE][MAX_SIZE], int row, int col)
{
  unsigned int i, j;
  for(i = 0; i < row; i++)
  {
    for(j = 0; j < col; j++)
	{
	  m3[i][j] =  (m1[i][j] + m2[i][j]);
	}
  }
}

/* m_sub -- matrix subtraction -- may be in-situ */
void m_sub(double m1[MAX_SIZE][MAX_SIZE], double m2[MAX_SIZE][MAX_SIZE],
           double m3[MAX_SIZE][MAX_SIZE], int row, int col)
{
  unsigned int i, j;
  for(i = 0; i < row; i++)
  {
    for(j = 0; j < col; j++)
	{
	  m3[i][j] =  (m1[i][j] - m2[i][j]);
	}
  }
}

/* m_zero -- zero the matrix A */
void m_zero(double A[MAX_SIZE][MAX_SIZE], int m, int n)
{
  int i, j, A_m, A_n;
  for(i = 0; i < m; i++)
    for(j = 0; j < n; j++)
      A[i][j] = 0.0;
}

/* __mltadd__ -- scalar multiply and add c.f. v_mltadd() */
void __mltadd__(double *dp1, const double *dp2, double s, int len)
{
  register int i;
//#ifdef VUNROLL
//  register int len4;
//  len4 = len / 4;
//  len  = len % 4;
//  for(i = 0; i < len4; i++)
//  {
//    dp1[4*i]   += s*dp2[4*i];
//    dp1[4*i+1] += s*dp2[4*i+1];
//    dp1[4*i+2] += s*dp2[4*i+2];
//    dp1[4*i+3] += s*dp2[4*i+3];
//  }
//  dp1 += 4*len4;    dp2 += 4*len4;
//#endif
  for(i = 0; i < len; i++)
    dp1[i] += s*dp2[i];
}

/* m_copy -- copies matrix into new area
  	-- B <- A */
void m_copy(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE],
            int m, int n)
{
  int i, j;
  for(i = 0; i < m; i++)
  {
    for(j = 0; j < n; j++)
	{
      B[i][j] = A[i][j];
	}
  }
}

/* m_mlt -- matrix-matrix multiplication */
void m_mlt(double m1[MAX_SIZE][MAX_SIZE], int m1_m, int m1_n,
           double m2[MAX_SIZE][MAX_SIZE], int m2_m, int m2_n,
           double m3[MAX_SIZE][MAX_SIZE])
{
  unsigned int i, j, k;
  double m4[MAX_SIZE][MAX_SIZE];
  m_zero(m4, m1_m, m2_n);
  m_zero(m3, m1_m, m2_n);
  if(m1_n == m2_m)
  {
    double mult;
    // Checking if the multiplication is possible
    // Initialising Matrix 3
    // Calculating multiplication result
    for(i = 0; i < m1_m; i++)
    {
      for(j = 0; j < m2_n; j++)
      {
        for(k = 0; k < m1_n; k++)
        {
          mult = (m1[i][k] * m2[k][j]);
          m4[i][j] = m4[i][j] + (m1[i][k] * m2[k][j]);
        }
      }
    }
    m_copy(m4, m3, m1_m, m2_n);
  }
  else
  {
    #ifdef DEBUG
	printf("\nError! Operation invalid, please enter with valid matrices.\n");
    #endif
  }
}


/* sm_mlt -- scalar-matrix multiply -- may be in-situ */
void sm_mlt(double scalar, double matrix[MAX_SIZE][MAX_SIZE],
            double out[MAX_SIZE][MAX_SIZE], int m, int n)
{
  unsigned int i, j;
  for(i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      out[i][j] = scalar*matrix[i][j];
}

/* v_zero -- zero the vector x */
void v_zero(double x[MAX_SIZE], int dim)
{
  int i;
  for(i = 0; i < dim; i++)
    x[i] = 0.0;
}

/* set_col -- sets column of matrix to values given in vec (in situ)
	-- that is, mat(i0:lim,col) <- vec(i0:lim) */
void set_col(double mat[MAX_SIZE][MAX_SIZE], unsigned int dim, unsigned int col,
             double vec[MAX_SIZE], double matOut[MAX_SIZE][MAX_SIZE])
{
  unsigned int i;
  m_copy(mat, matOut, dim, dim);
  for(i = 0; i < dim; i++)
    matOut[i][col] = vec[i];
}

/* v_copy -- copies vector into new area
	-- B <- A */
void v_copy(double A[MAX_SIZE], double B[MAX_SIZE], int dim)
{
  int i;
  for(i = 0; i < dim; i++)
    B[i] = A[i];
}

/* __ip__ -- inner product */
double __ip__(const double *dp1, const double *dp2, int len)
{
  register int i;
  register double sum;
  sum = 0.0;
//#ifdef VUNROLL
//  register int len4;
//  register double sum1, sum2, sum3;
//  sum1 = sum2 = sum3 = 0.0;
//  len4 = len / 4;
//  len  = len % 4;
//
//  for(i = 0; i < len4; i++)
//  {
//    sum  += dp1[4*i]*dp2[4*i];
//    sum1 += dp1[4*i+1]*dp2[4*i+1];
//    sum2 += dp1[4*i+2]*dp2[4*i+2];
//    sum3 += dp1[4*i+3]*dp2[4*i+3];
//  }
//  sum  += sum1 + sum2 + sum3;
//  dp1 += 4*len4;	dp2 += 4*len4;
//#endif
  for(i = 0; i < len; i++)
  {
    sum  += dp1[i]*dp2[i];
  }
  return sum;
}

/* m_inverse -- returns inverse of A, provided A is not too rank deficient
  -- uses Gauss - Jordan */
void m_inverse(double A[MAX_SIZE][MAX_SIZE], double out[MAX_SIZE][MAX_SIZE],
               int dim)
{
  int i, j, k, matsize;
  double temp;
  double AUX[MAX_SIZE][MAX_SIZE];
  m_copy(A, AUX, dim, dim);
  matsize = dim;
  // automatically initialize the unit matrix, e.g.
  for(i = 0; i < matsize; i++)
  {
    for(j = 0; j < matsize; j++)
    {
      if(i == j)
      {
        out[i][j] = 1;
      }
      else
        out[i][j] = 0;
    }
  }
/*---------------Logic starts here------------------*/
  /* procedure to make the matrix A to unit matrix
   --by some row operations,and the same row operations of
   --Unit mat. I gives the inverse of matrix A
   --'temp' stores the A[k][k] value so that A[k][k] will not change
   --during the operation A[i][j]/=A[k][k] when i=j=k
  --*/
  for(k = 0; k < matsize; k++)
  {
    // it performs the following row operations to make A to unit matrix
    // R0=R0/A[0][0],similarly for I also R0=R0/A[0][0]
    // R1=R1-R0*A[1][0] similarly for I
    // R2=R2-R0*A[2][0]
    temp = AUX[k][k];
    for(j = 0; j < matsize; j++)
    {
      AUX[k][j] /= temp;
      out[k][j] /= temp;
    }
    for(i = 0; i < matsize; i++)
    {
      // R1=R1/A[1][1]
      // R0=R0-R1*A[0][1]
      // R2=R2-R1*A[2][1]
      temp = AUX[i][k];
      for(j = 0; j < matsize; j++)
      {
        if(i == k)
          break;
        // R2=R2/A[2][2]
        // R0=R0-R2*A[0][2]
        // R1=R1-R2*A[1][2]
        AUX[i][j] -= AUX[k][j]*temp;
        out[i][j] -= out[k][j]*temp;
      }
    }
  }
/*---------------Logic ends here--------------------*/
}

/* m_ident -- set A to being closest to identity matrix as possible
  -- i.e. A[i][j] == 1 if i == j and 0 otherwise */
void m_ident(double A[MAX_SIZE][MAX_SIZE], int dim)
{
  int i, size;
  m_zero(A, dim, dim);
  size = dim;
  for(i = 0; i < size; i++)
    A[i][i] = 1.0;
}

void print_arr(double m[MAX_SIZE][MAX_SIZE], int row, int col)
{
  int i, j;
  for (i = 0; i < row; i++) {
    for (j = 0; j < col; j++) {
	  printf("%f ", m[i][j]);
    }
    printf("\n");
  }
}

/* fast_m_pow -- auxiliary function to compute integer powers of a
 * square matrix M, M^n */
void fast_m_pow(double m[MAX_SIZE][MAX_SIZE], int power,
                double result[MAX_SIZE][MAX_SIZE], int dim)
{
  double out[MAX_SIZE][MAX_SIZE], out2[MAX_SIZE][MAX_SIZE];
  double temp_m[MAX_SIZE][MAX_SIZE], aux_m[MAX_SIZE][MAX_SIZE];
  m_ident(result, dim);
//  print_arr(result, 5, 5);
  m_copy(m, temp_m, dim, dim);
  if(power == 0)
  {
	  // do nothing
  }
  else if(power == 1)
	m_copy(m, result, dim, dim);
  else
  {
    while(power)
    {
      if(power & 1)
      {
        m_mlt(result, dim, dim, temp_m, dim, dim, out);
        m_copy(out, result, dim, dim);
        m_zero(out, dim, dim);
      }
      m_mlt(temp_m, dim, dim, temp_m, dim, dim, out2);
      m_copy(out2, temp_m, dim, dim);
      m_zero(out2, dim, dim);
      power >>= 1;
    }
  }
}

/* smart_m_pow -- auxiliary function to compute integer powers of a
 * square matrix M, M^n */
void smart_m_pow(double m[MAX_SIZE][MAX_SIZE], int power,
                 double result[MAX_SIZE][MAX_SIZE], int dim)
{
  double out[MAX_SIZE][MAX_SIZE], out2[MAX_SIZE][MAX_SIZE];
  double temp_m[MAX_SIZE][MAX_SIZE];
  m_ident(result, dim);
  if(power == 0)
  {
	  // do nothing
  }
  else if(power == 1)
	m_copy(m, result, dim, dim);
  else
  {
    while(power > 0)
    {
      if((power % 2) == 1)
      {
        m_mlt(result, dim, dim, m, dim, dim, out);
        m_copy(out, result, dim, dim);
      }
      m_copy(m, temp_m, dim, dim);
      m_mlt(m, dim, dim, temp_m, dim, dim, out2);
      m_copy(out2, temp_m, dim, dim);
      m_copy(out2, m, dim, dim);
      power = power/2;
    }
  }
}

/* m_pow -- computes integer powers of a square matrix A, A^p */
void m_pow(double A[MAX_SIZE][MAX_SIZE], int p, double out[MAX_SIZE][MAX_SIZE],
           int dim)
{
  double tmp[MAX_SIZE][MAX_SIZE];
  if(p < 0)
  {
//	printf("we're actually here!\n");
    m_inverse(A, tmp, dim);
//    print_arr(tmp, dim, dim);
    fast_m_pow(tmp, -p, out, dim);
  }
  else
  {
    fast_m_pow(A, p, out, dim);
  }
}

/* get_col -- gets a specified column of a matrix and retruns it as a vector */
void get_col(double mat[MAX_SIZE][MAX_SIZE], unsigned int dim,
             unsigned int col, double vec[MAX_SIZE])
{
  unsigned int i;
  for(i = 0; i < dim; i++)
    vec[i] = mat[i][col];
}

/* _in_prod -- inner product of two vectors from i0 downwards
   -- that is, returns a(i0:dim)^T.b(i0:dim) */
double _in_prod(double a[MAX_SIZE], double b[MAX_SIZE], int dim, unsigned int i0)
{
  return __ip__(&(a[i0]), &(b[i0]), (int)(dim-i0));
}

/* hhvec -- calculates Householder vector to eliminate all entries after the
   i0 entry of the vector vec. It is returned as out. May be in-situ */
void hhvec(double vec[MAX_SIZE], unsigned int i0, double beta, double newval,
           double out[MAX_SIZE],  int dim)
{
  double norm, temp;
  v_copy(vec, out, dim);
  temp = (double)_in_prod(out, out, dim, i0);
  norm = sp_sqrt(temp);
  if(norm <= 0.0)
  {
    beta = 0.0;
//    return (out);
  }
  else
  {
    beta = 1.0/(norm * (norm + sp_fabs(out[i0])));
    if(out[i0] > 0.0)
      newval = -norm;
    else
      newval = norm;
    out[i0] -= newval;
  }
//  return (out);
}

/* _hhtrcols -- transform a matrix by a Householder vector by columns
    starting at row i0 from column j0
    -- that is, M(i0:m,j0:n) <- (I-beta.hh(i0:m).hh(i0:m)^T)M(i0:m,j0:n)
    -- in-situ
    -- scratch vector w passed as argument
    -- raises error if w == NULL
*/
void _hhtrcols(double M[MAX_SIZE][MAX_SIZE], int dim, unsigned int i0,
               unsigned int j0, double hh[MAX_SIZE], double beta,
               double w[MAX_SIZE], double out[MAX_SIZE][MAX_SIZE])
{
  int i;
  m_copy(M, out, dim, dim);
  if(beta == 0.0)
  {
//    do nothing -> return (M);
  }
  else
  {
    v_zero(w, dim);
    for(i = i0; i < dim; i++)
      if(hh[i] != 0.0)
        __mltadd__(&(w[j0]), &(out[i][j0]), hh[i],
                   (int)(dim - j0));
    for(i = i0; i < dim; i++)
      if(hh[i] != 0.0)
        __mltadd__(&(out[i][j0]), &(w[j0]), -beta*hh[i],
                   (int)(dim - j0));
//    return (M);
  }
}

/* hhtrrows -- transform a matrix by a Householder vector by rows
    starting at row i0 from column j0 -- in-situ
    -- that is, M(i0:m,j0:n) <- M(i0:m,j0:n)(I-beta.hh(j0:n).hh(j0:n)^T) */
void hhtrrows(double M[MAX_SIZE][MAX_SIZE], int dim, unsigned int i0,
              unsigned int j0, double hh[MAX_SIZE], double beta,
              double out[MAX_SIZE][MAX_SIZE])
{
  double ip, scale;
  int i;
  m_copy(M, out, dim, dim);
  if(beta == 0.0)
  {
//    do nothing -> return (M);
  }
  else
  {
    /* for each row ... */
    for(i = i0; i < dim; i++)
    { /* compute inner product */
      ip = __ip__(&(out[i][j0]), &(hh[j0]), (int)(dim - j0));
      scale = beta*ip;
      if(scale == 0.0)
        continue;
      /* do operation */
      __mltadd__(&(out[i][j0]), &(hh[j0]), -scale,
                 (int)(dim - j0));
    }
//  return (M);
  }
}

/* Hfactor -- compute Hessenberg factorization in compact form.
    -- factorization performed in situ
*/
void Hfactor(double A[MAX_SIZE][MAX_SIZE], int A_dim, double diag[MAX_SIZE],
             double beta[MAX_SIZE], double out[MAX_SIZE][MAX_SIZE])
{
  static double hh[MAX_SIZE], w[MAX_SIZE], hh_tmp[MAX_SIZE];
  int k, limit;
  double b;
  double A_temp[MAX_SIZE][MAX_SIZE], out2[MAX_SIZE][MAX_SIZE];
  m_copy(A, A_temp, A_dim, A_dim);
  limit = A_dim - 1;
  for(k = 0; k < limit; k++)
  {
    /* compute the Householder vector hh */
	get_col(A_temp, (unsigned int)A_dim, (unsigned int)k, hh);
	hhvec(hh, k+1, (beta[k]), (A_temp[k+1][k]), hh_tmp, A_dim);
	v_zero(hh, A_dim);
	v_copy(hh_tmp, hh, A_dim);
	v_zero(hh_tmp, A_dim);
	if(((k+1) >= 0) && ((k+1) < A_dim))
	{
	  diag[k] = hh[k+1];
	}
	else
	{
	  diag[k] = 0;
	}
    /* apply Householder operation symmetrically to A */
	if((k >= 0) && (k < A_dim))
	{
	  b = beta[k];
	}
    _hhtrcols(A_temp, A_dim, k+1, k+1, hh, b, w, out2);
    hhtrrows(out2, A_dim, 0, k+1, hh, b, out);
    m_copy(out, A_temp, A_dim, A_dim);
    m_zero(out2, A_dim, A_dim);
    m_zero(out, A_dim, A_dim);
    v_zero(hh, A_dim);
  }
  m_copy(A_temp, out, A_dim, A_dim);
//  return (A);
}

/* hhtrvec -- apply Householder transformation to vector
    -- that is, out <- (I-beta.hh(i0:n).hh(i0:n)^T).in
    -- may be in-situ */
void hhtrvec(double hh[MAX_SIZE], int dim, double beta, unsigned int i0,
             double in[MAX_SIZE], double out[MAX_SIZE])
{
  double scale, temp;
  temp = (double)_in_prod(hh, in, dim, i0);
  scale = beta*temp;
  v_copy(in, out, dim);
  __mltadd__(&(out[i0]), &(hh[i0]), -scale, (int)(dim - i0));
//  return (out);
}

/* makeHQ -- construct the Hessenberg orthogonalising matrix Q;
    -- i.e. Hess M = Q.M.Q'	*/
void makeHQ(double H[MAX_SIZE][MAX_SIZE], int dim, double diag[MAX_SIZE],
            double beta[MAX_SIZE], double Qout[MAX_SIZE][MAX_SIZE])
{
  int i, j, limit;
  double tmp1[MAX_SIZE], tmp2[MAX_SIZE], tmp3[MAX_SIZE];
  double tmp4[MAX_SIZE][MAX_SIZE];
  limit = dim - 1;
  for(i = 0; i < dim; i++)
  {
    /* tmp1 = i'th basis vector */
    for(j = 0; j < dim; j++)
      tmp1[j] = 0.0;
    tmp1[i] = 1.0;
    /* apply H/h transforms in reverse order */
    for(j = limit-1; j >= 0; j--)
    {
      get_col(H, dim, (unsigned int)j, tmp2);
      tmp2[j+1] = (((j)>=0 && (j)<dim) ? (diag)[(j)] : 0);
      hhtrvec(tmp2, dim, beta[j], j+1, tmp1, tmp3);
      v_copy(tmp3, tmp1, dim);
      v_zero(tmp3, dim);
    }
    /* insert into Qout */
    set_col(Qout, (unsigned int)dim, (unsigned int)i, tmp1, tmp4);
    m_zero(Qout, dim, dim);
    m_copy(tmp4, Qout, dim, dim);
    m_zero(tmp4, dim, dim);
  }
//  return (Qout);
}

/* makeH -- construct actual Hessenberg matrix */
void makeH(double H[MAX_SIZE][MAX_SIZE], int dim, double Hout[MAX_SIZE][MAX_SIZE])
{
  int i, j, limit;
  m_copy(H, Hout, dim, dim);
  limit = dim;
  for(i = 1; i < limit; i++)
    for(j = 0; j < i-1; j++)
    	Hout[i][j] = 0.0;
//  return (Hout);
}

/* rot_cols -- postmultiply mat by givens rotation described by c, s */
void rot_cols(double mat[MAX_SIZE][MAX_SIZE], unsigned int dim, unsigned int i,
              unsigned int k, double c, double s, double out[MAX_SIZE][MAX_SIZE])
{
  unsigned int j;
  double temp, temp1, temp2;
  m_copy(mat, out, dim, dim);
  for(j = 0; j < dim; j++)
  {
    temp = c*m_entry(out, dim, dim, j, i) + s*m_entry(out, dim, dim, j, k);
    temp1 = (((j)<dim && (i)<=dim) ? (out)[(j)][(i)] : (0));
    temp2 = (((j)<dim && (k)<=dim) ? (out)[(j)][(k)] : (0));
    out[j][k] = (-s*temp1 +	c*temp2);
    out[j][i] = temp;
  }
//  return (out);
}

/* rot_rows -- pre-multiply mat by givens rotation described by c, s */
void rot_rows(double mat[MAX_SIZE][MAX_SIZE], unsigned int dim, unsigned int i,
              unsigned int k, double c, double s, double out[MAX_SIZE][MAX_SIZE])
{
  unsigned int j;
  double temp, temp1, temp2, temp3;
  m_copy(mat, out, dim, dim);
  for(j = 0; j < dim; j++)
  {
	temp = c*m_entry(out, dim, dim, i, j) + s*m_entry(out, dim, dim, k, j);
	temp1 = (((i)<dim && (j)<=dim) ? (out)[(i)][(j)] : (0));
	temp2 = (((k)<dim && (j)<=dim) ? (out)[(k)][(j)] : (0));
    out[k][j] = (-s*temp1 + c*temp2);
    out[i][j] = temp;
  }
//  return (out);
}

/* hhldr3 -- computes */
static void hhldr3(double x, double y, double z, double *nu1,
                   double *beta, double *newval)
{
  double alpha;
  if(x >= 0.0)
    alpha = sp_sqrt(x*x+y*y+z*z);
  else
    alpha = -sp_sqrt(x*x+y*y+z*z);
  *nu1 = x + alpha;
  *beta = 1.0/(alpha*(*nu1));
  *newval = alpha;
}

/* hhldr3rows */
static void hhldr3rows(double A[MAX_SIZE][MAX_SIZE], int dim, int k, int i0,
                       double beta, double nu1, double nu2, double nu3,
                       double A_temp[MAX_SIZE][MAX_SIZE])
{
  double ip, prod, temp1, temp2, temp3;
  int i;
  m_copy(A, A_temp, dim, dim);
  i0 = min(i0, dim-1);
  for(i = 0; i <= i0; i++)
  {
	temp1 = nu1*m_entry(A_temp, dim, dim, i, k);
	temp2 = nu2*m_entry(A_temp, dim, dim, i, k+1);
	temp3 = nu3*m_entry(A_temp, dim, dim, i, k+2);
    ip = temp1 + temp2 + temp3;
    prod = ip*beta;
    A_temp[i][k] += (-prod*nu1);
    A_temp[i][k+1] += (-prod*nu2);
    A_temp[i][k+2] += (-prod*nu3);
  }
//  return A;
}

/* givens -- returns c,s parameters for Givens rotation to
       eliminate y in the vector [ x y ]' */
void givens(double x, double y, double *c, double *s)
{
  double norm;
  norm = sp_sqrt(x*x+y*y);
  if(norm == 0.0)
  {
    *c = 1.0;
    *s = 0.0;
  }	/* identity */
  else
  {
    *c = x/norm;
    *s = y/norm;
  }
}

/* schur -- computes the Schur decomposition of the matrix A in situ
    -- optionally, gives Q matrix such that Q^T.A.Q is upper triangular
    -- returns upper triangular Schur matrix */
void schur(double A[MAX_SIZE][MAX_SIZE], int dim, double Q[MAX_SIZE][MAX_SIZE],
           double A_out[MAX_SIZE][MAX_SIZE])
{
  int i, j, iter, k, k_min, k_max, k_tmp, split;
  double beta2, c, discrim, dummy, nu1, s, tmp, x, y, z;
  double A_me[MAX_SIZE][MAX_SIZE], A_temp[MAX_SIZE][MAX_SIZE];
  double Q_temp[MAX_SIZE][MAX_SIZE];
  double sqrt_macheps;
  static double diag[MAX_SIZE], beta[MAX_SIZE];
  double a00, a01, a10, a11;
  double scale, t, numer, denom;
  /* compute Hessenberg form */
  Hfactor(A, dim, diag, beta, A_temp);
  /* save Q if necessary */
  makeHQ(A_temp, dim, diag, beta, Q_temp);
  makeH(A_temp, dim, A_out);
  m_zero(A_temp, dim, dim);
  sqrt_macheps = sp_sqrt(MACHEPS);
  k_min = 0;
//  A_me = A_out;
//  memcpy(A_me, A.me, MAX_SIZE*MAX_SIZE*8);
//  for(i = 0; i < dim; i++)
//  {
//    for(j = 0; j < dim; j++)
//    {
//    	A_me[i][j] = A_out[i][j];
//    }
//  }
  while(k_min < dim)
  {
    /* find k_max to suit:
       submatrix k_min..k_max should be irreducible */
    k_max = dim-1;
    for(k = k_min; k < k_max; k++)
      if(m_entry(A_out, dim, dim, k+1, k) == 0.0)
      {
        k_max = k;
        break;
      }
    if(k_max <= k_min)
    {
      k_min = k_max + 1;
      continue;      /* outer loop */
    }
    /* check to see if we have a 2 x 2 block
       with complex eigenvalues */
    if(k_max == (k_min + 1))
    {
      a00 = m_entry(A_out, dim, dim, k_min, k_min);
      a01 = m_entry(A_out, dim, dim, k_min, k_max);
      a10 = m_entry(A_out, dim, dim, k_max, k_min);
      a11 = m_entry(A_out, dim, dim, k_max, k_max);
      tmp = a00 - a11;
      discrim = tmp*tmp + 4*a01*a10;
      if(discrim < 0.0)
      {
        /* yes -- e-vals are complex
               -- put 2 x 2 block in form [a b; c a];
        then eigenvalues have real part a & imag part sp_sqrt(|bc|) */
        numer = - tmp;
        denom = (a01+a10 >= 0.0) ?
                (a01+a10) + sp_sqrt((a01+a10)*(a01+a10)+tmp*tmp) :
                (a01+a10) - sp_sqrt((a01+a10)*(a01+a10)+tmp*tmp);
        if(denom != 0.0)
        {    /* t = s/c = numer/denom */
          t = numer/denom;
          scale = c = 1.0/sp_sqrt(1+t*t);
          s = c*t;
        }
        else
        {
          c = 1.0;
          s = 0.0;
        }
        rot_cols(A_out, dim, k_min, k_max, c, s, A_temp);
        m_zero(A_out, dim, dim);
        rot_rows(A_temp, dim, k_min, k_max, c, s, A_out);
        m_zero(A_temp, dim, dim);
        rot_cols(Q_temp, dim, k_min, k_max, c, s, A_temp);
        m_zero(Q_temp, dim, dim);
        m_copy(A_temp, Q_temp, dim, dim);
        m_zero(A_temp, dim, dim);
        k_min = k_max + 1;
        continue;
      }
      else
      {
        /* discrim >= 0; i.e. block has two real eigenvalues */
        /* no -- e-vals are not complex;
         split 2 x 2 block and continue */
        /* s/c = numer/denom */
        numer = (tmp >= 0.0) ?
              - tmp - sp_sqrt(discrim) : - tmp + sp_sqrt(discrim);
        denom = 2*a01;
        if(sp_fabs(numer) < sp_fabs(denom))
        {    /* t = s/c = numer/denom */
          t = numer/denom;
          scale = c = 1.0/sp_sqrt(1+t*t);
          s = c*t;
        }
        else if(numer != 0.0)
        {    /* t = c/s = denom/numer */
          t = denom/numer;
          scale = 1.0/sp_sqrt(1+t*t);
          c = sp_fabs(t)*scale;
          s = (t >= 0.0) ? scale : -scale;
        }
        else /* numer == denom == 0 */
        {
          c = 0.0;
          s = 1.0;
        }
        rot_cols(A_out, dim, k_min, k_max, c, s, A_temp);
        m_zero(A_out, dim, dim);
        rot_rows(A_temp, dim, k_min, k_max, c, s, A_out);
        m_zero(A_temp, dim, dim);
        rot_cols(Q_temp, dim, k_min, k_max, c, s, A_temp);
        m_zero(Q_temp, dim, dim);
        m_copy(A_temp, Q_temp, dim, dim);
        m_zero(A_temp, dim, dim);
        k_min = k_max + 1;  /* go to next block */
        continue;
      }
    }
    /* now have r x r block with r >= 2:
     apply Francis QR step until block splits */
    split = 0;
    iter = 0;
    while(!split)
    {
      iter++;
      /* set up Wilkinson/Francis complex shift */
      k_tmp = k_max - 1;
      a00 = m_entry(A_out, dim, dim, k_tmp, k_tmp);
      a01 = m_entry(A_out, dim, dim, k_tmp, k_max);
      a10 = m_entry(A_out, dim, dim, k_max, k_tmp);
      a11 = m_entry(A_out, dim, dim, k_max, k_max);
      /* treat degenerate cases differently
         -- if there are still no splits after five iterations
            and the bottom 2 x 2 looks degenerate, force it to
         split */
      #ifdef DEBUG
        printf("# schur: bottom 2 x 2 = [%lg, %lg; %lg, %lg]\n",
               a00, a01, a10, a11);
      #endif
      if(iter >= 5 &&
         sp_fabs(a00-a11) < sqrt_macheps*(sp_fabs(a00)+sp_fabs(a11)) &&
         (sp_fabs(a01) < sqrt_macheps*(sp_fabs(a00)+sp_fabs(a11)) ||
          sp_fabs(a10) < sqrt_macheps*(sp_fabs(a00)+sp_fabs(a11))) )
      {
        if(sp_fabs(a01) < sqrt_macheps*(sp_fabs(a00)+sp_fabs(a11)))
        {
          A_out[k_tmp][k_max] = 0.0;
        }
        if(sp_fabs(a10) < sqrt_macheps*(sp_fabs(a00)+sp_fabs(a11)))
        {
          A_out[k_max][k_tmp] = 0.0;
          split = 1;
          continue;
        }
      }
      s = a00 + a11;
      t = a00*a11 - a01*a10;
      /* break loop if a 2 x 2 complex block */
      if(k_max == k_min + 1 && s*s < 4.0*t)
      {
        split = 1;
        continue;
      }
      /* perturb shift if convergence is slow */
      if((iter % 10) == 0)
      {
        s += iter*0.02;
        t += iter*0.02;
      }
      /* set up Householder transformations */
      k_tmp = k_min + 1;
      a00 = m_entry(A_out, dim, dim, k_min, k_min);
      a01 = m_entry(A_out, dim, dim, k_min, k_tmp);
      a10 = m_entry(A_out, dim, dim, k_tmp, k_min);
      a11 = m_entry(A_out, dim, dim, k_tmp, k_tmp);
      x = a00*a00 + a01*a10 - s*a00 + t;
      y = a10*(a00+a11-s);
      if(k_min + 2 <= k_max)
        z = a10*A_out[k_min+2][k_tmp];
      else
        z = 0.0;
      for(k = k_min; k <= k_max-1; k++)
      {
        if(k < k_max - 1)
        {
          hhldr3(x, y, z, &nu1, &beta2, &dummy);
          hhldr3rows(Q_temp, dim, k, dim-1, beta2, nu1, y, z, A_temp);
          m_zero(Q_temp, dim, dim);
          m_copy(A_temp, Q_temp, dim, dim);
          m_zero(A_temp, dim, dim);
        }
        else
        {
          givens(x, y, &c, &s);
          rot_cols(A_out, dim, k, k+1, c, s, A_temp);
          m_zero(A_out, dim, dim);
          rot_rows(A_temp, dim, k, k+1, c, s, A_out);
          m_zero(A_temp, dim, dim);
          rot_cols(Q_temp, dim, k, k+1, c, s, A_temp);
          m_zero(Q_temp, dim, dim);
          m_copy(A_temp, Q_temp, dim, dim);
          m_zero(A_temp, dim, dim);
        }
        x = m_entry(A_out, dim, dim, k+1, k);
        if(k <= k_max - 2)
          y = m_entry(A_out, dim, dim, k+2, k);
        else
          y = 0.0;
        if(k <= k_max - 3)
          z = m_entry(A_out, dim, dim, k+3, k);
        else
          z = 0.0;
      }
	  for(k = k_min; k <= k_max-2; k++)
	  {
        /* zero appropriate sub-diagonals */
		A_out[k+2][k] = 0.0;
        if(k < k_max-2)
        {
          A_out[k+3][k] = 0.0;
        }
      }

      /* test to see if matrix should split */
      for(k = k_min; k < k_max; k++)
        if(sp_fabs(A_out[k+1][k]) < MACHEPS*
          (sp_fabs(A_out[k][k])+sp_fabs(A_out[k+1][k+1])))
        {
          A_out[k+1][k] = 0.0;
          split = 1;
        }
	}
  }
  /* polish up A by zeroing strictly lower triangular elements
     and small sub-diagonal elements */
  for(i = 0; i < dim; i++)
    for(j = 0; j < i-1; j++)
      A_out[i][j] = 0.0;
    for(i = 0; i < dim - 1; i++)
      if(sp_fabs(A_out[i+1][i]) < MACHEPS*
         (sp_fabs(A_out[i][i])+sp_fabs(A_out[i+1][i+1])))
        A_out[i+1][i] = 0.0;
//  return A;
}

/* schur_vals -- compute real & imaginary parts of eigenvalues
	-- assumes T contains a block upper triangular matrix
		as produced by schur()
	-- real parts stored in real_pt, imaginary parts in imag_pt */
void schur_evals(double T[MAX_SIZE][MAX_SIZE], int dim, double real_pt[MAX_SIZE],
                 double imag_pt[MAX_SIZE])
{
  int i, j, n;
  double discrim, T_me[MAX_SIZE][MAX_SIZE];
  double diff, sum, tmp;
  n = dim;
  m_copy(T, T_me, dim, dim);
  i = 0;
  while(i < n)
  {
    if(i < n-1 && T_me[i+1][i] != 0.0)
    {   /* should be a complex eigenvalue */
      sum  = 0.5*(T_me[i][i]+T_me[i+1][i+1]);
      diff = 0.5*(T_me[i][i]-T_me[i+1][i+1]);
      discrim = diff*diff + T_me[i][i+1]*T_me[i+1][i];
      if(discrim < 0.0)
      { /* yes -- complex e-vals */
        real_pt[i] = real_pt[i+1] = sum;
        imag_pt[i] = sp_sqrt(-discrim);
        imag_pt[i+1] = -imag_pt[i];
      }
      else
      { /* no -- actually both real */
        tmp = sp_sqrt(discrim);
        real_pt[i]   = sum + tmp;
        real_pt[i+1] = sum - tmp;
        imag_pt[i]   = imag_pt[i+1] = 0.0;
      }
      i += 2;
    }
    else
    {   /* real eigenvalue */
      real_pt[i] = T_me[i][i];
      imag_pt[i] = 0.0;
      i++;
    }
  }
}

/******************************Verification Functions******************************/

/* m_get_eigenvalues -- get the eigenvalues of a matrix A
	-- */
CMPLX *m_get_eigenvalues(double A[MAX_SIZE][MAX_SIZE], int dim)
{
  double T[MAX_SIZE][MAX_SIZE], Q[MAX_SIZE][MAX_SIZE];
  double evals_re[MAX_SIZE], evals_im[MAX_SIZE];
  static CMPLX z[MAX_SIZE];
//  m_copy(A, T, dim, dim);
  /* compute Schur form: A = Q.T.Q^T */
  schur(A, dim, Q, T);
  /* extract eigenvalues */
  schur_evals(T, dim, evals_re, evals_im);
  for(int i = 0; i < dim; i++)
  {
    z[i].real = evals_re[i];
    z[i].imag = evals_im[i];
  }
  return z;
}

/* cmplx_mag -- get the magnitude of a complex number taking its real
 * and imaginary parts */
double cmplx_mag(double real, double imag)
{
  return sp_sqrt(real * real + imag * imag);
}

/* is_same_sign -- check if a has the same sign as b */
int is_same_sign(double a, double b)
{
  if(((a >= 0) && (b >= 0)) || ((a <= 0) && (b <= 0)))
    return 1;
  else
    return 0;
}

/* x_k -- computes the state signal in the k-th sample */
void x_k(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE], int dim, double u, int k)
{
  double AUX[MAX_SIZE][MAX_SIZE], AUX2[MAX_SIZE][MAX_SIZE], AUX3[MAX_SIZE][MAX_SIZE];
  if(xk.lastState == (k - 1))
  {
    m_mlt(A, dim, dim, xk.xk, dim, 1, AUX);
    sm_mlt(u, B, AUX2, dim, 1);
    m_add(AUX, AUX2, AUX3, dim, 1);
    xk.lastState = k;
//    m_zero(xk.xk, dim, 1);
    m_copy(AUX3, xk.xk, dim, 1);
//    xk.xk = AUX3;
  }
//  for(int i=0;i<dim;i++){
//  	printf("x(%d)=%f\n", i, xk.xk[i][0]);
//  }
}

/* y_k2 -- computes the output signal in the k-th sample */
double y_k2(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE],
            double C[MAX_SIZE][MAX_SIZE], double D[MAX_SIZE][MAX_SIZE],
            double u, int k, int dim)
{
  double Ak[MAX_SIZE][MAX_SIZE], AUX[MAX_SIZE][MAX_SIZE], AUX2[MAX_SIZE][MAX_SIZE];
  double y, temp;
  // y[k]=Cx[k]+Du[k]
  x_k(A, B, dim, u, k);
//  print_arr(C, 1, dim);
//  print_arr(xk.xk, dim, 1);
  m_mlt(C, 1, 5, xk.xk, 5, 1, AUX);
//  printf("AUX[0][0]=%f\n", AUX[0][0]);
//  print_arr(AUX, 1, 1);
  temp = D[0][0] * u;
  y = AUX[0][0] + temp;
//  printf("D=%f\n", D[0][0]);
//  printf("u=%f\n", u);
//  printf("temp=%f\n", temp);
  return y;
}

/* x_k2 -- computes the state signal in the k-th sample */
void x_k2(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE],
         double C[MAX_SIZE][MAX_SIZE], double D[MAX_SIZE][MAX_SIZE],
         double u, int k, double X0[MAX_SIZE][MAX_SIZE], int dim,
         double out[MAX_SIZE][MAX_SIZE])
{
  double x[MAX_SIZE][MAX_SIZE], Ak[MAX_SIZE][MAX_SIZE], AUX[MAX_SIZE][MAX_SIZE];
  double AUX3[MAX_SIZE][MAX_SIZE], AUX4[MAX_SIZE][MAX_SIZE], x_tmp[MAX_SIZE][MAX_SIZE];
  double AUX2[MAX_SIZE][MAX_SIZE], AUX5[MAX_SIZE][MAX_SIZE];
  int m;
  // y = C * A.pow(k) * X0;
  m_pow(A, k, Ak, dim);
  m_mlt(Ak, dim, dim, X0, dim, 1, x);
  for(m = 0; m <= (k - 1); m++)
  {
    // y += (C * A.pow(k - m - 1) * B * u) + D * u;
    m_pow(A, (k-m-1), Ak, dim);
    m_mlt(Ak, dim, dim, B, dim, 1, AUX);
    sm_mlt(u, AUX, AUX2, dim, 1);
    m_add(x, AUX2, x_tmp, dim, 1);
    m_zero(x, dim, 1);
    m_copy(x_tmp, x, dim, 1);
    m_zero(x_tmp, dim, 1);
  }
  xk.lastState = k;
  m_copy(x, xk.xk, dim, 1);
  m_copy(x, out, dim, 1);
}

/* y_k3 -- computes the output signal in the k-th sample */
double y_k3(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE],
            double C[MAX_SIZE][MAX_SIZE], double D[MAX_SIZE][MAX_SIZE],
            double u, int k, double X0[MAX_SIZE][MAX_SIZE], int dim)
{
  double Ak[MAX_SIZE][MAX_SIZE], Xk[MAX_SIZE][MAX_SIZE];
  double AUX[MAX_SIZE][MAX_SIZE], AUX2[MAX_SIZE][MAX_SIZE];
  double y, temp;
  x_k2(A, B, C, D, u, k, X0, dim, Xk);
  // y[k]=Cx[k]+Du[k]
//  x_k(A, B, u, k);
  m_mlt(C, 1, dim, Xk, dim, 1, AUX);
  sm_mlt(u, D, AUX2, 1, 1);
  y = AUX[0][0] + AUX2[0][0];
  return y;
}

/* y_k -- computes the output signal in the k-th sample */
double y_k(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE],
           double C[MAX_SIZE][MAX_SIZE], double D[MAX_SIZE][MAX_SIZE],
           double u, int k, double X0[MAX_SIZE][MAX_SIZE], int dim)
{
  double y[MAX_SIZE][MAX_SIZE], Ak[MAX_SIZE][MAX_SIZE], AUX[MAX_SIZE][MAX_SIZE];
  double AUX3[MAX_SIZE][MAX_SIZE], AUX4[MAX_SIZE][MAX_SIZE], y_tmp[MAX_SIZE][MAX_SIZE];
  double AUX2[MAX_SIZE][MAX_SIZE], AUX5[MAX_SIZE][MAX_SIZE];
  int m;
  // y = C * A.pow(k) * X0;
  m_pow(A, k, Ak, dim);
  m_mlt(C, 1, dim, Ak, dim, dim, AUX);
  m_mlt(AUX, 1, dim, X0, dim, 1, y);
  for(m = 0; m <= (k - 1); m++)
  {
    // y += (C * A.pow(k - m - 1) * B * u) + D * u;
    m_pow(A, (k-m-1), Ak, dim);
    m_mlt(C, 1, dim, Ak, dim, dim, AUX);
    m_mlt(AUX, 1, dim, B, dim, 1, AUX2);
    m_add(AUX2, D, AUX5, 1, 1);
    m_add(y, AUX5, y_tmp, 1, 1);
    m_copy(y_tmp, y, 1, 1);
    m_zero(y_tmp, 1, 1);
  }
  return y[0][0]*u;
}

/* y_ss -- computes steady-state output value of a given system */
double y_ss(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE],
            double C[MAX_SIZE][MAX_SIZE], double D[MAX_SIZE][MAX_SIZE],
            double u, int dim)
{
  double yss;
  double AUX[MAX_SIZE][MAX_SIZE], AUX2[MAX_SIZE][MAX_SIZE], AUX3[MAX_SIZE][MAX_SIZE];
  double AUX4[MAX_SIZE][MAX_SIZE], AUX5[MAX_SIZE][MAX_SIZE], Id[MAX_SIZE][MAX_SIZE];
  // get the expression y_ss=(C(I-A)^(-1)B+D)u
  m_ident(Id, dim);
  // Id - A
  m_sub(Id, A, AUX, dim, dim);
  m_inverse(AUX, AUX2, dim);
  m_mlt(C, 1, dim, AUX2, dim, dim, AUX3);
  m_mlt(AUX3, 1, dim, B, dim, 1, AUX4);
  m_add(AUX4, D, AUX5, 1, 1);
  yss = AUX5[0][0] * u;
  return yss;
}

/* peak_output -- computes the biggest peak value of a signal (Mp) */
PKVL peak_output(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE],
                 double C[MAX_SIZE][MAX_SIZE], double D[MAX_SIZE][MAX_SIZE],
                 double X0[MAX_SIZE][MAX_SIZE], double yss, double u, int dim)
{
  PKVL out;
  double greater, cmp, o;
  int i = 0;
  greater = sp_fabs(y_k2(A, B, C, D, u, i, dim));
  o = y_k2(A, B, C, D, u, i+1, dim);
  cmp = sp_fabs(o);
//  printf("greater=%f\n", y_k2(A, B, C, D, u, 1, dim));
  while((cmp >= greater))
  {
//    printf("(%f >= %f)\n", cmp, greater);
    if(greater < cmp)
    {
//      printf("(%f < %f)\n", greater, cmp);
      greater = cmp;
      out.mp = o;
      out.kp = i+2;
//      printf("greater1=%f\n", y_k2(A, B, C, D, u, i+1, dim));
//      j++;
    }
    else
    {
      out.mp = o;
      out.kp = i+2;
//      printf("greater2=%f\n", y_k2(A, B, C, D, u, i+1, dim));
//      j++;
    }
    if(!is_same_sign(yss, out.mp))
    {
      greater = 0;
    }
    i++;
    o = y_k2(A, B, C, D, u, i+1, dim);
    cmp = sp_fabs(o);
  }
//  printf("j=%d\n", j);
  return out;
}

/* c_bar -- computes an auxiliary variable to calculate k_bar */
double c_bar(double mp, double yss, double lambmax, int kp)
{
  double cbar;
  cbar = (mp-yss)/(sp_pow(lambmax, kp));
  return cbar;
}

/* log_b -- computes the log of x in the base 'base' */
double log_b(double base, double x)
{
  return (double) (sp_log10_2(x) / sp_log10_2(base));
}

/* k_bar -- computes instant in which the system enters in the settling
 * -time region */
int k_bar(double lambdaMax, double p, double cbar, double yss, int order)
{
  double k_ss, x;
  x = (p * yss) / (100 * cbar);
  k_ss = log_b(lambdaMax, x);
  return sp_ceil(k_ss)+order;
}

/* max_mag_eigenvalue -- computes biggest magnitude among the eigenvalues */
double max_mag_eigenvalue(double A[MAX_SIZE][MAX_SIZE], int dim)
{
  double maximum = 0, aux;
  CMPLX *z;
  int i;
  z = m_get_eigenvalues(A, dim);
  for(i = 0; i < dim; i++)
  {
    aux = cmplx_mag(z[i].real, z[i].imag);
    if(aux > maximum)
    {
      maximum = aux;
    }
  }
  return maximum;
}

/* check_settling_time -- check if a given settling time satisfies to
 * a given system */
int check_settling_time(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE],
                        double C[MAX_SIZE][MAX_SIZE], double D[MAX_SIZE][MAX_SIZE],
                        double X0[MAX_SIZE][MAX_SIZE], double u, double tsr,
                        double p, double ts, int dim)
{
  double yss, mp, lambMax, cbar, output, inf, sup;
  PKVL out;
  int kbar, kp, i;
//  xk.xk = m_get(A.m, 1);
  xk.lastState = 0;
  yss = y_ss(A, B, C, D, u, dim);
  out = peak_output(A, B, C, D, X0, yss, u, dim);
  mp = out.mp;
  kp = out.kp;
  lambMax = max_mag_eigenvalue(A, dim);
  cbar = c_bar(mp, yss, lambMax, kp);
  kbar = k_bar(lambMax, p, cbar, yss, dim);
  #ifdef DEBUG
  printf("Mp=%f\n", mp);
  printf("yss=%f\n", yss);
  printf("lambMax=%f\n", lambMax);
  printf("kp=%d\n", kp);
  printf("cbar=%f\n", cbar);
  #endif
  if(kbar * ts < tsr)
  {
    #ifdef DEBUG
    printf("kbar=%d\n", kbar);
    #endif
    return 1;
  }
  i = (int)sp_ceil(tsr / ts)-1;
  output = y_k3(A, B, C, D, u, i, X0, dim);
  while(i <= kbar)
  {
    if(yss > 0)
    {
      inf = (yss - (yss * (p/100)));
      sup = (yss * (p/100) + yss);
    }
    else
    {
      sup = (yss - (yss * (p/100)));
      inf = (yss * (p/100) + yss);
    }
    if(!(output < sup && (output > inf)))
    {
      #ifdef DEBUG
      printf("kbar=%d\n", kbar);
      #endif
      return 0;
    }
    i++;
    output = y_k2(A, B, C, D, u, i, dim);
  }
  return 1;
}

int main() {
//	MAT A, B, C, D, X0;
	CMPLX *z;
//	CMPLX my_z[MAX_SIZE];
//	double u, tsr, p, ts;
//	int res;
	double u = 1.0;
//	ts = 0.5;
//	p = 5;
//
//	A = m_get(5, 5);
//	A.me[0][0]=-0.5000;A.me[0][1]=1.000;A.me[0][2]=0;A.me[0][3]=0;A.me[0][4]=0;
//	A.me[1][0]=0.0000;A.me[1][1]=-0.5000;A.me[1][2]=1.0;A.me[1][3]=0;A.me[1][4]=0;
//	A.me[2][0]=0;A.me[2][1]=0;A.me[2][2]=-0.5000;A.me[2][3]=1.0000;A.me[2][4]=0;
//	A.me[3][0]=0;A.me[3][1]=0;A.me[3][2]=0.0000;A.me[3][3]=-0.5000;A.me[3][4]=1.0;
//	A.me[4][0]=0;A.me[4][1]=0;A.me[4][2]=0;A.me[4][3]=0;A.me[4][4]=-0.5;
//
//
//	B = m_get(5,1);
//	B.me[0][0]=0;
//	B.me[1][0]=0;
//	B.me[2][0]=2.5;
//	B.me[3][0]=1;
//	B.me[4][0]=0;
//
//	C = m_get(1,5);
//	C.me[0][0]=0;C.me[0][1]=2.6;C.me[0][2]=0.5;C.me[0][3]=1.2;C.me[0][4]=0;
//
//	D=m_get(1,1);
//
//	X0=m_get(5,1);
//	int tvar = 5;
//	xk.xk = m_get(tvar, 1);

	double myA[MAX_SIZE][MAX_SIZE];
	myA[0][0]=-0.5000;myA[0][1]=1.000;myA[0][2]=0;myA[0][3]=0;myA[0][4]=0;
	myA[1][0]=0.0000;myA[1][1]=-0.5000;myA[1][2]=1.0;myA[1][3]=0;myA[1][4]=0;
	myA[2][0]=0;myA[2][1]=0;myA[2][2]=-0.5000;myA[2][3]=1.0000;myA[2][4]=0;
	myA[3][0]=0;myA[3][1]=0;myA[3][2]=0.0000;myA[3][3]=-0.5000;myA[3][4]=1.0;
	myA[4][0]=0;myA[4][1]=0;myA[4][2]=0;myA[4][3]=0;myA[4][4]=-0.5000;

	double myB[MAX_SIZE][MAX_SIZE];
	myB[0][0]=-0.4;
	myB[1][0]=-0.6;
	myB[2][0]=5.5;
	myB[3][0]=1;
	myB[4][0]=-0.3;

	double myC[MAX_SIZE][MAX_SIZE];
	myC[0][0]=0;myC[0][1]=2.6;myC[0][2]=0.5;myC[0][3]=1.2;myC[0][4]=0;

	double myD[MAX_SIZE][MAX_SIZE];
	m_zero(myD, 1, 1);

	double myX0[MAX_SIZE][MAX_SIZE];

//	double myT[MAX_SIZE][MAX_SIZE];
//	m_add(myA,myA,myT,5,5);
//	double myS[MAX_SIZE][MAX_SIZE];
//	m_sub(myA,myA,myS,5,5);
//	m_zero(myS, 5, 5);
//	double w[MAX_SIZE];
//	__mltadd__(&(w[0]), &(myA[0][0]), w[0],
//	                   (int)(5-0));
//	m_mlt(myA, 5, 5, myA, 5, 5, myS);
//	m_copy(myA, 5, 5, myS, 5, 5);
//	sm_mlt(2.0, myA, myS, 5, 5);
//	m_inverse(myA, myS, 5);
//	m_ident(myS, 5);
//	m_pow(myA,-4, myS, 5);
//	m_mlt(myA, 5, 5, myA, 5, 5, myS);
//	print_arr(myS,5,5);
//	m_mlt(myS, 5, 5, myS, 5, 5, myT);
//	smart_m_pow(myA, 4, myS, 5);
//	MAT t;
//	t = m_add(A, A);
//	print_arr(myA,5,5);
//	double vec[MAX_SIZE], out[MAX_SIZE];
//	get_col(myA, 5, 2, vec);
//	for(int i = 0;i<5;i++){
//		printf("vec[%d]=%f | ", i, vec[i]);
//	}
//	printf("\n");
//
//	printf("_in_prod=%f\n", _in_prod(vec, 5, vec, 5, 2));
//	hhvec(vec, 2, 1.5, 2, out, 5);
//	for(int i = 0;i<5;i++){
//	  printf("out[%d]=%f | ", i, out[i]);
//	}
//	printf("\n");

//    z = m_get_eigenvalues(myA, 5);
//	  double T[MAX_SIZE][MAX_SIZE], Q[MAX_SIZE][MAX_SIZE];
	  /* compute Schur form: A = Q.T.Q^T */
//	  schur(myA, 5, Q, T);
//	  double diag[MAX_SIZE], beta[MAX_SIZE];
//	  Hfactor(myA, 5, diag, beta, myT);
//	double evals_re[MAX_SIZE], evals_im[MAX_SIZE];
//	  /* extract eigenvalues */
//	  schur_evals(myA, 5, evals_re, evals_im);
//    int size = 5;
//    for(int i = 0;i < size;i++)
//    {
////      printf("%f+%f i", z[i].real, z[i].imag);
////      printfc(z[i]);
////      assert(z[i].real==-0.5);
//    	assert(-0.5==-0.5);
//    }
//	print_arr(myS,5,5);
//    double yss = y_ss(myA, myB, myC, myD, u, 5);
//    printf("yss = %f\n", yss);
//    assert(yss!=yss);
//    PKVL p = peak_output(myA, myB, myC, myD, myX0, yss, u, 5);
//    printf("Mp=%f e kp=%d\n", p.mp, p.kp);
    double tsr = 7*0.5;
    int res = check_settling_time(myA, myB, myC, myD, myX0, u, tsr, 5, 0.5, 5);
    assert(res != 1);
//    for(int i=0;i<5;i++){
//    	printf("x(%d)=%f\n", i, xk.xk[i][0]);
//    }
	return 0;
}
