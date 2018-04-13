#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#include "seahorn/seahorn.h"

#define MAX_SIZE 10

#define MACHEPS 2.22045e-16

#define DEBUG

#ifdef DEBUG
#define	m_output(mat) m_foutput(stdout, mat)
#endif

static const char *format = "%14.9g ";

#ifndef ANSI_C
#define ANSI_C 1
#endif

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

/* matrix definition */
typedef	struct
{
  unsigned int m, n;
  unsigned int max_m, max_n, max_size;
  double me[MAX_SIZE][MAX_SIZE];
}
MAT;

/* vector definition */
typedef	struct
{
  unsigned int dim, max_dim;
  double ve[MAX_SIZE];
}
VEC;

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
  MAT xk;
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
    return d - 1;
}

/* sp_sqrt -- returns the square root of fg */
double sp_sqrt(double fg)
{
  double n = fg / 2.0;
  double lstX = 0.0;
  while(n != lstX)
  {
    lstX = n;
    n = (n + fg/n) / 2.0;
  }
  return n;
}

/******************************Matrix Functions******************************/

/* m_add -- matrix addition -- may be in-situ */
#ifndef ANSI_C
void m_add(m1, m2, m3, row, col)
double m1[MAX_SIZE][MAX_SIZE], m2[MAX_SIZE][MAX_SIZE], m3[MAX_SIZE][MAX_SIZE];
int row, col;
#else
void m_add(double m1[MAX_SIZE][MAX_SIZE], double m2[MAX_SIZE][MAX_SIZE], double m3[MAX_SIZE][MAX_SIZE], int row, int col)
#endif
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
#ifndef ANSI_C
void m_sub(m1, m2, m3, row, col)
double m1[MAX_SIZE][MAX_SIZE], m2[MAX_SIZE][MAX_SIZE], m3[MAX_SIZE][MAX_SIZE];
int row, col;
#else
void m_sub(double m1[MAX_SIZE][MAX_SIZE], double m2[MAX_SIZE][MAX_SIZE], double m3[MAX_SIZE][MAX_SIZE], int row, int col)
#endif
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
#ifndef ANSI_C
MAT m_zero(A, m, n)
double A[MAX_SIZE][MAX_SIZE];
int m, n;
#else
void m_zero(double A[MAX_SIZE][MAX_SIZE], int m, int n)
#endif
{
  int i, j, A_m, A_n;
  for(i = 0; i < m; i++)
    for(j = 0; j < n; j++)
      A[i][j] = 0.0;
}

/* __mltadd__ -- scalar multiply and add c.f. v_mltadd() */
#ifndef ANSI_C
void __mltadd__(dp1, dp2, s, len)
register double	*dp1, *dp2;
register double s;
register int len;
#else
void __mltadd__(double *dp1, const double *dp2, double s, int len)
#endif
{
  register int i;
#ifdef VUNROLL
  register int len4;
  len4 = len / 4;
  len  = len % 4;
  for(i = 0; i < len4; i++)
  {
    dp1[4*i]   += s*dp2[4*i];
    dp1[4*i+1] += s*dp2[4*i+1];
    dp1[4*i+2] += s*dp2[4*i+2];
    dp1[4*i+3] += s*dp2[4*i+3];
  }
  dp1 += 4*len4;    dp2 += 4*len4;
#endif
  for(i = 0; i < len; i++)
    dp1[i] += s*dp2[i];
}

/* m_mlt -- matrix-matrix multiplication */
#ifndef ANSI_C
void m_mlt(m1, m1_m, m1_n, m2, m2_m, m2_n, m3)
const double m1[MAX_SIZE][MAX_SIZE], m2[MAX_SIZE][MAX_SIZE], m3[MAX_SIZE][MAX_SIZE];
int m1_m, m1_n, m2_m, m2_n;
#else
void m_mlt(const double m1[MAX_SIZE][MAX_SIZE], int m1_m, int m1_n, const double m2[MAX_SIZE][MAX_SIZE], int m2_m, int m2_n, double m3[MAX_SIZE][MAX_SIZE])
#endif
{
  unsigned int i1, j1, i2, j2;
  unsigned int i, j, k;
  i1 = m1_m;
  j1 = m1_n;
  i2 = m2_m;
  j2 = m2_n;
  if(j1 == i2)
  {
    double mult;
    // Checking if the multiplication is possible
    // Initialising Matrix 3
//    m3 = m_get(i1, j2);
    // Calculating multiplication result
    for(i = 0; i < i1; i++)
    {
      for(j = 0; j < j2; j++)
      {
        for(k = 0; k < j1; k++)
        {
          mult = (m1[i][k] * m2[k][j]);
          m3[i][j] = m3[i][j] + (m1[i][k] * m2[k][j]);
        }
      }
    }
  }
  else
  {
    #ifdef DEBUG
	printf("\nError! Operation invalid, please enter with valid matrices.\n");
    #endif
  }
}

/* m_copy -- copies matrix into new area
  	-- B <- A */
void m_copy(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE], int m, int n)
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

/* sm_mlt -- scalar-matrix multiply -- may be in-situ */
#ifndef ANSI_C
void sm_mlt(scalar, matrix, out, m, n)
double scalar;
double matrix[MAX_SIZE][MAX_SIZE], out[MAX_SIZE][MAX_SIZE];
int m, n;
#else
void sm_mlt(double scalar, double matrix[MAX_SIZE][MAX_SIZE], double out[MAX_SIZE][MAX_SIZE], int m, int n)
#endif
{
  unsigned int i, j;
  for(i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      out[i][j] = scalar*matrix[i][j];
}

/* v_zero -- zero the vector x */
#ifndef ANSI_C
void v_zero(x, dim)
double x[MAX_SIZE];
int dim;
#else
void v_zero(double x[MAX_SIZE], int dim)
#endif
{
  int i;
  for(i = 0; i < dim; i++)
    x[i] = 0.0;
}

/* set_col -- sets column of matrix to values given in vec (in situ)
	-- that is, mat(i0:lim,col) <- vec(i0:lim) */
#ifndef ANSI_C
void set_col(mat, mat_m, col, vec, dim)
double mat[MAX_SIZE][MAX_SIZE];
double vec[MAX_SIZE];
unsigned int col, mat_m, dim;
#else
void set_col(double mat[MAX_SIZE][MAX_SIZE], unsigned int mat_m, unsigned int col, double vec[MAX_SIZE], unsigned int dim)
#endif
{
  unsigned int i, lim, i0;
  lim = min(mat_m, dim);
  for(i = i0; i < lim; i++)
    mat[i][col] = vec[i];
}

/* v_copy -- copies vector into new area
	-- B <- A */
#ifndef ANSI_C
void v_copy(A, B)
const double A[MAX_SIZE], B[MAX_SIZE];
#else
void v_copy(double A[MAX_SIZE], double B[MAX_SIZE], int dim)
#endif
{
  int i;
  for(i = 0; i < dim; i++)
    B[i] = A[i];
}

/* __ip__ -- inner product */
#ifndef ANSI_C
double __ip__(dp1, dp2, len)
register double	*dp1, *dp2;
int len;
#else
  double __ip__(const double *dp1, const double *dp2, int len)
#endif
{
#ifdef VUNROLL
  register int len4;
  register double sum1, sum2, sum3;
#endif
  register int i;
  register double sum;
  sum = 0.0;
#ifdef VUNROLL
  sum1 = sum2 = sum3 = 0.0;
  len4 = len / 4;
  len  = len % 4;

  for(i = 0; i < len4; i++)
  {
    sum  += dp1[4*i]*dp2[4*i];
    sum1 += dp1[4*i+1]*dp2[4*i+1];
    sum2 += dp1[4*i+2]*dp2[4*i+2];
    sum3 += dp1[4*i+3]*dp2[4*i+3];
  }
  sum  += sum1 + sum2 + sum3;
  dp1 += 4*len4;	dp2 += 4*len4;
#endif
  for(i = 0; i < len; i++)
  {
    sum  += dp1[i]*dp2[i];
  }
  return sum;
}

/* m_inverse -- returns inverse of A, provided A is not too rank deficient
  -- uses Gauss - Jordan */
#ifndef ANSI_C
void m_inverse(A, out, dim)
double A[MAX_SIZE][MAX_SIZE];
double out[MAX_SIZE][MAX_SIZE];
int dim;
#else
void m_inverse(double A[MAX_SIZE][MAX_SIZE], double out[MAX_SIZE][MAX_SIZE], int dim)
#endif
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
#ifndef ANSI_C
void m_ident(A, dim)
double A[MAX_SIZE];
int dim;
#else
void m_ident(double A[MAX_SIZE][MAX_SIZE], int dim)
#endif
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

/* fast_m_pow -- auxiliary function to compute integer powers of a square matrix M, M^n */
void fast_m_pow(double m[MAX_SIZE][MAX_SIZE], int power, double result[MAX_SIZE][MAX_SIZE], int dim)
{
  double out[MAX_SIZE][MAX_SIZE], out2[MAX_SIZE][MAX_SIZE], temp_m[MAX_SIZE][MAX_SIZE], aux_m[MAX_SIZE][MAX_SIZE];
  m_ident(result, dim);
  print_arr(result, 5, 5);
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
    	printf("aqui!\n");
    	print_arr(result, 5, 5);
        m_mlt(result, dim, dim, temp_m, dim, dim, out);
        print_arr(out, 5, 5);
        m_copy(out, result, dim, dim);
        m_zero(out, dim, dim);
        print_arr(result, 5, 5);
      }
      printf("here2!\n");
//      print_arr(m, 5, 5);
//      m_copy(m, temp_m, dim, dim);
//      print_arr(temp_m, 5, 5);
      m_mlt(temp_m, dim, dim, temp_m, dim, dim, out2);
      print_arr(out2, 5, 5);
      m_copy(out2, temp_m, dim, dim);
      print_arr(temp_m, 5, 5);
      m_zero(out2, dim, dim);
      power >>= 1;
    }
  }
}

/* smart_m_pow -- auxiliary function to compute integer powers of a square matrix M, M^n */
void smart_m_pow(double m[MAX_SIZE][MAX_SIZE], int power, double result[MAX_SIZE][MAX_SIZE], int dim)
{
  double out[MAX_SIZE][MAX_SIZE], out2[MAX_SIZE][MAX_SIZE], temp_m[MAX_SIZE][MAX_SIZE];
  m_ident(result, dim);
  if(power == 0)
  {
	  // do nothing
  }
  else if(power == 1)
	m_copy(m, result, dim, dim);
  else
  {
    while(power>0)
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
#ifndef ANSI_C
void m_pow(A, p, out, dim)
double A[MAX_SIZE][MAX_SIZE], out[MAX_SIZE][MAX_SIZE];
int	p, dim;
#else
MAT m_pow(double A[MAX_SIZE][MAX_SIZE], int p, double out[MAX_SIZE][MAX_SIZE], int dim)
#endif
{
  double tmp[MAX_SIZE][MAX_SIZE];
  if(p < 0)
  {
	printf("we're actually here!\n");
    m_inverse(A, tmp, dim);
    print_arr(tmp, dim, dim);
    fast_m_pow(tmp, -p, out, dim);
  }
  else
  {
    fast_m_pow(A, p, out, dim);
  }
}

/* get_col -- gets a specified column of a matrix and retruns it as a vector */
#ifndef ANSI_C
void get_col(mat, m_dim, col, vec)
unsigned int col;
double mat[MAX_SIZE][MAX_SIZE];
double vec[MAX_SIZE];
#else
void get_col(double mat[MAX_SIZE][MAX_SIZE], unsigned int m_dim, unsigned int col, double vec[MAX_SIZE])
#endif
{
  unsigned int i;
  for(i = 0; i < m_dim; i++)
    vec[i] = mat[i][col];
}

/* _in_prod -- inner product of two vectors from i0 downwards
   -- that is, returns a(i0:dim)^T.b(i0:dim) */
#ifndef ANSI_C
double _in_prod(a, a_dim, b, b_dim, i0)
double a[MAX_SIZE], b[MAX_SIZE];
unsigned int i0;
int a_dim, b_dim;
#else
double _in_prod(double a[MAX_SIZE], int a_dim, double b[MAX_SIZE], int b_dim, unsigned int i0)
#endif
{
  unsigned int limit;
  limit = min(a_dim, b_dim);
  return __ip__(&(a[i0]), &(b[i0]), (int)(limit-i0));
}

/* hhvec -- calculates Householder vector to eliminate all entries after the
   i0 entry of the vector vec. It is returned as out. May be in-situ */
#ifndef ANSI_C
void hhvec(vec, i0, beta, newval, out, dim)
double vec[MAX_SIZE], out[MAX_SIZE];
unsigned int i0;
double beta, newval;
int dim;
#else
void hhvec(double vec[MAX_SIZE], unsigned int i0, double beta, double newval, double out[MAX_SIZE],  int dim)
#endif
{
  double norm, temp;
  v_copy(vec, out, dim);
  temp = (double)_in_prod(out, dim, out, dim, i0);
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
#ifndef ANSI_C
void _hhtrcols(M, i0, j0, hh, beta, w, out)
double M[MAX_SIZE][MAX_SIZE], out[MAX_SIZE][MAX_SIZE];
unsigned int i0, j0;
double hh[MAX_SIZE];
double	beta;
double w[MAX_SIZE];
#else
void _hhtrcols(double M[MAX_SIZE][MAX_SIZE], int dim, unsigned int i0, unsigned int j0,
               double hh[MAX_SIZE], double beta, double w[MAX_SIZE], double out[MAX_SIZE][MAX_SIZE])
#endif
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
#ifndef ANSI_C
void hhtrrows(M, dim, i0, j0, hh, beta)
double M[MAX_SIZE][MAX_SIZE];
unsigned int i0, j0;
double hh[MAX_SIZE];
double beta;
int dim;
#else
void hhtrrows(double M[MAX_SIZE][MAX_SIZE], int dim, unsigned int i0, unsigned int j0,
              double hh[MAX_SIZE], double beta, double out[MAX_SIZE][MAX_SIZE])
#endif
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
#ifndef ANSI_C
void Hfactor(A, A_dim, diag, beta, out)
double A[MAX_SIZE][MAX_SIZE], out[MAX_SIZE][MAX_SIZE];
double diag[MAX_SIZE], beta[MAX_SIZE];
int A_dim;
#else
void Hfactor(double A[MAX_SIZE][MAX_SIZE], int A_dim, double diag[MAX_SIZE], double beta[MAX_SIZE], double out[MAX_SIZE][MAX_SIZE])
#endif
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
	hhvec(hh, k+1, (beta[k]), (A[k+1][k]), hh_tmp, A_dim);
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
  }
//  return (A);
}

/* hhtrvec -- apply Householder transformation to vector
    -- that is, out <- (I-beta.hh(i0:n).hh(i0:n)^T).in
    -- may be in-situ */
#ifndef ANSI_C
void hhtrvec(hh, dim, beta, i0, in, out)
double hh[MAX_SIZE], in[MAX_SIZE], out[MAX_SIZE];	/* hh = Householder vector */
unsigned int i0;
double beta;
int dim;
#else
void hhtrvec(double hh[MAX_SIZE], int dim, double beta, unsigned int i0, double in[MAX_SIZE], double out[MAX_SIZE])
#endif
{
  double scale, temp;
  temp = (double)_in_prod(hh, dim, in, dim, i0);
  scale = beta*temp;
  v_copy(in, out, dim);
  __mltadd__(&(out[i0]), &(hh[i0]), -scale, (int)(dim - i0));
//  return (out);
}

/* makeHQ -- construct the Hessenberg orthogonalising matrix Q;
    -- i.e. Hess M = Q.M.Q'	*/
#ifndef ANSI_C
void makeHQ(H, dim, diag, beta, Qout)
double H[MAX_SIZE][MAX_SIZE], Qout[MAX_SIZE][MAX_SIZE];
double diag[MAX_SIZE], beta[MAX_SIZE];
int dim;
#else
void makeHQ(double H[MAX_SIZE][MAX_SIZE], int dim, double diag[MAX_SIZE], double beta[MAX_SIZE], double Qout[MAX_SIZE][MAX_SIZE])
#endif
{
  int i, j, limit;
  static double tmp1[MAX_SIZE], tmp2[MAX_SIZE], tmp3[MAX_SIZE];
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
    set_col(Qout, (unsigned int)dim, (unsigned int)i, tmp1, (unsigned int)dim);
  }
//  return (Qout);
}

/* makeH -- construct actual Hessenberg matrix */
#ifndef ANSI_C
void makeH(H, dim, Hout)
double H[MAX_SIZE][MAX_SIZE], Hout[MAX_SIZE][MAX_SIZE];
int dim;
#else
void makeH(double H[MAX_SIZE][MAX_SIZE], int dim, double Hout[MAX_SIZE][MAX_SIZE])
#endif
{
  int i, j, limit;
  m_copy(H, Hout, dim, dim);
  limit = dim;
  for(i = 1; i < limit; i++)
    for(j = 0; j < i-1; j++)
    	Hout[i][j] = 0.0;
//  return (Hout);
}

// TODO
/* rot_cols -- postmultiply mat by givens rotation described by c, s */
#ifndef ANSI_C
void rot_cols(mat, dim, i, k, c, s, out)
double mat[MAX_SIZE][MAX_SIZE], out[MAX_SIZE][MAX_SIZE];
unsigned int i, k, dim;
double c, s;
#else
void rot_cols(double mat[MAX_SIZE][MAX_SIZE], unsigned int dim, unsigned int i, unsigned int k, double c, double s, double out[MAX_SIZE][MAX_SIZE])
#endif
{
  unsigned int j;
  double temp;
  m_copy(mat, out, dim, dim);
  for(j = 0; j < dim; j++)
  {
    temp = c*m_entry(out, dim, dim, j, i) + s*m_entry(out, dim, dim, j, k);
    out[j][k] = (-s*(((j)<dim && (i)<=dim) ? \
            (out)[(j)][(i)] : (0)) +
        	c*(((j)<dim && (k)<=dim) ? \
            (out)[(j)][(k)] : (0)));
    out[j][i] = temp;
  }
//  return (out);
}

/* rot_rows -- pre-multiply mat by givens rotation described by c, s */
#ifndef ANSI_C
void rot_rows(mat, i, k, c, s)
double mat[MAX_SIZE][MAX_SIZE], out[MAX_SIZE][MAX_SIZE];
unsigned int i, k, dim;
double c, s;
#else
void rot_rows(double mat[MAX_SIZE][MAX_SIZE], unsigned int dim, unsigned int i, unsigned int k, double c, double s, double out[MAX_SIZE][MAX_SIZE])
#endif
{
  unsigned int j;
  double temp;
  m_copy(mat, out, dim, dim);
  for(j = 0; j < dim; j++)
  {
    temp = c*m_entry(out, dim, dim, i, j) + s*m_entry(out, dim, dim, k, j);
    out[k][j] = (-s*(((i)<dim && (j)<=dim) ? \
                       (out)[(i)][(j)] : (0)) +
        		       c*(((k)<dim && (j)<=dim) ? \
                       (out)[(k)][(j)] : (0)));
    out[i][j] = temp;
  }
//  return (out);
}

/* hhldr3 -- computes */
#ifndef ANSI_C
static void hhldr3(x, y, z, nu1, beta, newval)
double x, y, z;
double *nu1, *beta, *newval;
#else
static void hhldr3(double x, double y, double z, double *nu1,
                   double *beta, double *newval)
#endif
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
#ifndef ANSI_C
static void hhldr3rows(A, dim, k, i0, beta, nu1, nu2, nu3, A_temp)
double A[MAX_SIZE][MAX_SIZE], A_temp[MAX_SIZE][MAX_SIZE];
int	dim, k, i0;
double beta, nu1, nu2, nu3;
#else
static void hhldr3rows(double A[MAX_SIZE][MAX_SIZE], int dim, int k, int i0, double beta,
                       double nu1, double nu2, double nu3, double A_temp[MAX_SIZE][MAX_SIZE])
#endif
{
  double ip, prod;
  int i;
  m_copy(A, A_temp, dim, dim);
  i0 = min(i0, dim-1);
  for(i = 0; i <= i0; i++)
  {
    ip = nu1*m_entry(A_temp, dim, dim, i, k) + nu2*m_entry(A_temp, dim, dim, i, k+1)+nu3 * m_entry(A_temp, dim, dim, i, k+2);
    prod = ip*beta;
    A_temp[i][k] += (-prod*nu1);
    A_temp[i][k+1] += (-prod*nu2);
    A_temp[i][k+2] += (-prod*nu3);
  }
//  return A;
}

/* givens -- returns c,s parameters for Givens rotation to
       eliminate y in the vector [ x y ]' */
#ifndef ANSI_C
void givens(x, y, c, s)
double x, y;
double *c, *s;
#else
void givens(double x, double y, double *c, double *s)
#endif
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
#ifndef ANSI_C
void schur(A, Q)
double A[MAX_SIZE][MAX_SIZE], Q[MAX_SIZE][MAX_SIZE];
#else
void schur(double A[MAX_SIZE][MAX_SIZE], int dim, double Q[MAX_SIZE][MAX_SIZE], double A_out[MAX_SIZE][MAX_SIZE])
#endif
{
  int i, j, iter, k, k_min, k_max, k_tmp, n, split;
  double beta2, c, discrim, dummy, nu1, s, tmp, x, y, z;
  double A_me[MAX_SIZE][MAX_SIZE], A_temp[MAX_SIZE][MAX_SIZE], Q_temp[MAX_SIZE][MAX_SIZE];
  double sqrt_macheps;
  static double diag[MAX_SIZE], beta[MAX_SIZE];
  n = dim;

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
  double a00, a01, a10, a11;
  double scale, t, numer, denom;
  while(k_min < n)
  {
    /* find k_max to suit:
       submatrix k_min..k_max should be irreducible */
    k_max = n-1;
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
          hhldr3rows(Q_temp, dim, k, n-1, beta2, nu1, y, z, A_temp);
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
#ifndef ANSI_C
void schur_evals(T, real_pt, imag_pt)
MAT *T;
VEC *real_pt, *imag_pt;
#else
void schur_evals(double T[MAX_SIZE][MAX_SIZE], int dim, double real_pt[MAX_SIZE], double imag_pt[MAX_SIZE])
#endif
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

int main() {
//	MAT A, B, C, D, X0;
	CMPLX *z;
//	CMPLX my_z[MAX_SIZE];
//	double u, tsr, p, ts;
//	int res;
//	u = 1.0;
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
	myA[4][0]=0;myA[4][1]=0;myA[4][2]=0;myA[4][3]=0;myA[4][4]=-0.5;

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

    z = m_get_eigenvalues(myA, 5);
//    int size = 5;
//    for(int i = 0;i < size;i++)
//    {
//    //		printf("%f+%f i", z[i].real, z[i].imag);
//      printfc(z[i]);
//    }
//	print_arr(myS,5,5);

//	MAT T, Q;
//	static VEC diag, beta;
//	VEC evals_re, evals_im;
//	Q = m_get(A.m, A.n);
//	T = m_copy(A);
////	T = schur(T, Q);
//	diag = v_get(A.n);
//	beta = v_get(A.n);
//	A = Hfactor(A, diag, beta);
//	evals_re = v_get(A.m);
//	evals_im = v_get(A.m);

//	static VEC hh, w;
//	VEC diag = v_get(5);
//	VEC beta = v_get(A.m);
//	int k=0, limit;
//	double b, temp;
//	limit = A.m - 1;
//	hh = v_get(A.m);
//	w  = v_get(A.n);
//
//	hh = get_col(A, (unsigned int)k);
////	hh = hhvec(hh, k+1, (beta.ve[k]), (A.me[k+1][k]));
//	temp = (double)_in_prod(hh, hh, k+1);
////	schur_evals(&T, &evals_re, &evals_im);
//
//	assert(sp_fabs(-5) == 5);
//	printf("ceil(%f)=%f\n", 5.89, sp_ceil(5.89));
//	assert(sp_ceil(5.89) == 6.0);
//	assert(sp_floor(4.3) == 4);
//	assert(sp_pow(2, 2) == 4);
//	assert(sp_sqrt(4) == 2);
	return 0;
}
