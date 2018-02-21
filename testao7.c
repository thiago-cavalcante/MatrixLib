#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define MAX_SIZE 10

/* miscellaneous constants */
#define	VNULL ((VEC *)NULL)
#define	MNULL ((MAT *)NULL)
#define	MNULL2 ((MAT2 *)NULL)

/* macros that also check types and sets pointers to NULL */
#define	M_FREE(mat) (m_free(mat), (mat)=(MAT *)NULL)
#define V_FREE(vec)	(v_free(vec), (vec)=(VEC *)NULL)
#define	M_FREE2(mat) (m_free2(mat), (mat)=(MAT2 *)NULL)
//#define V_FREE2(vec)	(v_free2(vec), (&vec)=(VEC2 *)NULL)

#define MACHEPS 2.22045e-16

#define	m_output(mat) m_foutput(stdout, mat)
#define	m_output2(mat) m_foutput2(stdout, mat)

static const char *format = "%14.9g ";

#ifndef ANSI_C
#define ANSI_C 1
#endif

#define SEGMENTED

#ifndef THREADSAFE /* for use as a shared library */
#define THREADSAFE 1
#endif

#define TYPE_MAT mem_bytes(0, 0, sizeof(MAT))
#define TYPE_VEC mem_bytes(0, 0, sizeof(VEC))

#define	v_chk_idx(x, i) ((i)>=0 && (i)<(x)->dim)

#define	v_chk_idx2(x, i) ((i)>=0 && (i)<(x).dim)

#define v_get_val(x, i) (v_chk_idx(x, i) ? (x)->ve[(i)] : \
    (printf("Error!\n")))

#define v_get_val2(x, i) (v_chk_idx2(x, i) ? (x).ve[(i)] : \
    (printf("Error!\n")))

#define	v_entry(x, i) v_get_val(x, i)

#define	v_entry2(x, i) v_get_val2(x, i)

#define	v_set_val(x, i, val) (x->ve[i] = val)

#define	v_set_val2(x, i, val) (x.ve[i] = val)

#define	m_set_val(A, i, j, val) ((A)->me[(i)][(j)] = (val))

#define	m_set_val2(A, i, j, val) ((A).me[(i)][(j)] = (val))

#define	m_add_val(A, i, j, val)	((A)->me[(i)][(j)] += (val))

#define	m_add_val2(A, i, j, val)	((A).me[(i)][(j)] += (val))

#define	m_chk_idx(A, i, j) ((i)>=0 && (i)<(A)->m && (j)>=0 && (j)<=(A)->n)

#define	m_chk_idx2(A, i, j) ((i)>=0 && (i)<(A).m && (j)>=0 && (j)<=(A).n)

#define	m_get_val(A, i, j) (m_chk_idx(A, i, j) ? \
    (A)->me[(i)][(j)] : (printf("Error!")))

#define	m_get_val2(A, i, j) (m_chk_idx2(A, i, j) ? \
    (A).me[(i)][(j)] : (printf("Error!")))

#define	m_entry(A, i, j) m_get_val(A, i, j)

#define	m_entry2(A, i, j) m_get_val2(A, i, j)

#define printfc(c) printf("%f%c%fi\n", c.real, (c.imag>=0.0f)? '+':'\0', c.imag)

/* standard copy & zero functions */
#define	MEM_COPY(from, to, size) memmove((to), (from), (size))
#define	MEM_ZERO(where, size) memset((where), '\0', (size))

/* allocate one object of given type */
#define	NEW(type) ((type *)calloc((size_t)1, (size_t)sizeof(type)))

/* allocate num objects of given type */
#define	NEW_A(num, type)  ((type *)calloc((size_t)(num), (size_t)sizeof(type)))

#define	MEMCOPY(from, to, n_items, type) \
  MEM_COPY((char *)(from), (char *)(to), (unsigned)(n_items)*sizeof(type))

/* type independent min and max operations */
#ifndef max
#define	max(a, b) ((a) > (b) ? (a) : (b))
#endif /* max */
#ifndef min
#define	min(a, b) ((a) > (b) ? (b) : (a))
#endif /* min */

#ifndef THREADSAFE
#define MEM_STAT_REG(var, type) mem_stat_reg_list((void **)&(var),
                     type, 0, __FILE__, __LINE__)
#else
#define MEM_STAT_REG(var, type)
#endif

/* matrix definition */
typedef	struct
{
  unsigned int m, n;
  unsigned int max_m, max_n, max_size;
  double **me, *base;   /* base is base of alloc'd mem */
}
MAT;

/* matrix definition */
typedef	struct
{
  unsigned int m, n;
  unsigned int max_m, max_n, max_size;
  double me[MAX_SIZE][MAX_SIZE];   /* base is base of alloc'd mem */
}
MAT2;

/* vector definition */
typedef	struct
{
  unsigned int dim, max_dim;
  double *ve;
}
VEC;

/* vector definition */
typedef	struct
{
  unsigned int dim, max_dim;
  double ve[MAX_SIZE];
}
VEC2;

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

/******************************Matrix Functions******************************/


MAT2 m_get2(int m, int n)
{
  MAT2 A;
  A.m = m;A.n = n;
  A.max_m = m;A.max_n = n;
  A.max_size = A.max_m * A.max_n;
  for(int i = 0; i < m; i++)
  {
    for(int j = 0; j < n; j++){
      A.me[i][j] = 0;
    }
  }
  return A;
}

/* m_add2 -- matrix addition -- may be in-situ */
#ifndef ANSI_C
MAT2 m_add2(mat1, mat2)
MAT2 mat1, mat2;
#else
MAT2 m_add2(MAT2 mat1, MAT2 mat2)
#endif
{
  MAT2 out;
  unsigned int i, j, mat1_m, mat1_n, mat2_m, mat2_n;
  mat1_m = mat1.m; mat1_n = mat1.n;
  mat2_m = mat2.m; mat2_n = mat2.n;

//  mat2 = m_resize(mat2, mat1_m, mat1_n);
  out = m_get2(mat1_m, mat1_n);
  for(i=0; i<mat1_m; i++)
  {
    for(j = 0; j < mat1_n; j++)
	  out.me[i][j] = mat1.me[i][j]+mat2.me[i][j];
  }
  return (out);
}

/* m_sub2 -- matrix subtraction -- may be in-situ */
#ifndef ANSI_C
MAT2 m_sub2(mat1, mat2)
MAT2 mat1, mat2;
#else
MAT2 m_sub2(MAT2 mat1, MAT2 mat2)
#endif
{
  MAT2 out;
  unsigned int i, j, mat1_m, mat1_n, mat2_m, mat2_n;
  mat1_m = mat1.m; mat1_n = mat1.n;
  mat2_m = mat2.m; mat2_n = mat2.n;

//  mat2 = m_resize(mat2, mat1_m, mat1_n);
  out = m_get2(mat1_m, mat1_n);
  for(i=0; i<mat1_m; i++)
  {
    for(j = 0; j < mat1_n; j++)
	  out.me[i][j] = mat1.me[i][j]-mat2.me[i][j];
  }
  return (out);
}

/* m_zero2 -- zero the matrix A */
#ifndef ANSI_C
MAT2 m_zero2(m, n)
int m, n;
#else
MAT2 m_zero2(int m, int n)
#endif
{
  MAT2 A;
  int i, j, A_m, A_n;
//  double **A_me;
  A = m_get2(m, n);
  A_m = A.m;	A_n = A.n;/*	A_me = A.me;*/
  for(i = 0; i < A_m; i++)
    for(j = 0; j < A_n; j++)
      A.me[i][j] = 0.0;
  return A;
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

/* m_mlt2 -- matrix-matrix multiplication */
#ifndef ANSI_C
MAT2 m_mlt2(A, B)
MAT2 A, B;
#else
MAT2 m_mlt2(const MAT2 A, const MAT2 B)
#endif
{
  MAT2 OUT;
  unsigned int i, /* j, */ k, m, n, p;
  double A_v[MAX_SIZE][MAX_SIZE], B_v[MAX_SIZE][MAX_SIZE], sum, tmp;

  m = A.m;	n = A.n; p = B.n;
  memcpy(A_v, A.me, sizeof (A_v));
  memcpy(B_v, B.me, sizeof (B_v));
//  memcpy(A_v, A.me, min(sizeof(A_v), sizeof(A.me)));
//  memcpy(B_v, B.me, min(sizeof(B_v), sizeof(B.me)));
//  A_v = A.me; B_v = B.me;

  if(OUT.m != A.m || OUT.n != B.n)
    OUT = m_get2(A.m, B.n);
//  m_zero(OUT);
  for(i=0; i<m; i++)
    for( k=0; k<n; k++)
    {
      if(A_v[i][k] != 0.0)
        __mltadd__(OUT.me[i], B_v[k], A_v[i][k], (int)p);
    }

  return OUT;
}

/* m_foutput2 -- prints a representation of the matrix a onto file/stream fp */
#ifndef ANSI_C
void m_foutput2(fp, a)
FILE *fp;
MAT2 a;
#else
void m_foutput2(FILE *fp, const MAT2 a)
#endif
{
  unsigned int i, j, tmp;

  fprintf(fp, "Matrix: %d by %d\n", a.m, a.n);
  for(i = 0; i < a.m; i++)   /* for each row... */
  {
    fprintf(fp, "row %u: ", i);
    for(j = 0, tmp = 2; j < a.n; j++, tmp++)
    {             /* for each col in row... */
      fprintf(fp, format, a.me[i][j]);
        if(!(tmp % 5))
          putc('\n', fp);
    }
    if(tmp % 5 != 1)
      putc('\n', fp);
  }
}

/* v_zero2 -- zero the vector x */
#ifndef ANSI_C
VEC2 v_zero2(x)
VEC2 x;
#else
VEC2 v_zero2(VEC2 x)
#endif
{
  for(int i = 0; i < x.dim; i++)
    x.ve[i] = 0.0;
  return x;
}

/* set_col2 -- sets column of matrix to values given in vec (in situ)
	-- that is, mat(i0:lim,col) <- vec(i0:lim) */
#ifndef ANSI_C
MAT2 set_col2(mat, col, vec)
MAT2 mat;
VEC2 vec;
unsigned int col;
#else
MAT2 set_col2(MAT2 mat, unsigned int col, const VEC2 vec/*, unsigned int i0*/)
#endif
{
  unsigned int i, lim, i0;

  lim = min(mat.m, vec.dim);
  for(i=i0; i<lim; i++)
    mat.me[i][col] = vec.ve[i];

  return (mat);
}

VEC2 v_get2(int dim)
{
  VEC2 V;
  V.dim = dim;
  V.max_dim = dim;
  for(int i = 0; i < dim; i++)
  {
    V.ve[i] = 0;
  }
  return V;
}

/* v_copy2 -- copies vector into new area
	-- out(i0:dim) <- in(i0:dim) */
#ifndef ANSI_C
VEC2 v_copy2(in)
VEC2 in;
#else
VEC2 v_copy2(const VEC2 in)
#endif
{
  VEC2 out;
  unsigned int i0 = 0;

  if(out.dim < in.dim)
    out = v_get2(in.dim);

  MEM_COPY(&(in.ve[i0]), &(out.ve[i0]), (in.dim - i0)*sizeof(double));

  return (out);
}

VEC2 v_copy22(VEC2 A)
{
  VEC2 B;
  B.dim = A.dim;
  B.max_dim = A.max_dim;
  memcpy(B.ve, A.ve, sizeof (B.ve));
  return B;
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

MAT2 m_copy2(MAT2 A)
{
  MAT2 B;
  B.m = A.m;
  B.max_m = A.max_m;
  B.max_n = A.max_n;
  B.max_size = A.max_size;
  memcpy(B.me, A.me, sizeof (B.me));
  B.n = A.n;
  return B;
}

/* m_inverse2 -- returns inverse of A, provided A is not too rank deficient
-- uses Gauss - Jordan */
#ifndef ANSI_C
MAT2 m_inverse2(A)
MAT2 A;
#else
MAT2 m_inverse2(const MAT2 A)
#endif
{
  MAT2 out = m_get2(A.m, A.n);
  int i, j, k, matsize;
  double temp;
  MAT2 AUX = m_copy2(A);
  matsize = AUX.m;
  // automatically initialize the unit matrix, e.g.
  for(i = 0; i < matsize; i++)
  {
    for(j = 0; j < matsize; j++)
    {
      if(i == j)
      {
        out.me[i][j]=1;
      }
      else
        out.me[i][j]=0;
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
    temp = AUX.me[k][k];
    for(j = 0; j < matsize; j++)
    {
      AUX.me[k][j] /= temp;
      out.me[k][j] /= temp;
    }
    for(i = 0; i < matsize; i++)
    {
      // R1=R1/A[1][1]
      // R0=R0-R1*A[0][1]
      // R2=R2-R1*A[2][1]
      temp = AUX.me[i][k];
      for(j = 0; j < matsize; j++)
      {
        if(i == k)
          break;
        // R2=R2/A[2][2]
        // R0=R0-R2*A[0][2]
        // R1=R1-R2*A[1][2]
        AUX.me[i][j] -= AUX.me[k][j]*temp;
        out.me[i][j] -= out.me[k][j]*temp;
      }
    }
  }
/*---------------Logic ends here--------------------*/
  return out;
}

/* mat_ident2 -- set A to being closest to identity matrix as possible
  -- i.e. A[i][j] == 1 if i == j and 0 otherwise */
#ifndef ANSI_C
MAT2 m_ident2(dim)
int dim;
#else
MAT2 m_ident2(int dim)
#endif
{
  MAT2 A = m_get2(dim, dim);
  int i, size;
  A = m_zero2(dim, dim);
  size = min(A.m, A.n);
  for(i = 0; i < size; i++)
    A.me[i][i] = 1.0;
  return A;
}

MAT2 matpow(MAT2 m, int p)
{
  MAT2 tmp;
  if(p == 0)
    return m_ident2(m.m);
  else if(p % 2 == 0)
  {
    tmp = matpow(m, p / 2);
    return m_mlt2(tmp, tmp);
  }
  else
    return m_mlt2(m, matpow(m, p - 1));
}

/* m_pow3 -- computes integer powers of a square matrix A, A^p */
#ifndef ANSI_C
MAT2 m_pow3(A, p)
MAT2 A;
int	p;
#else
MAT2 m_pow3(const MAT2 A, int p)
#endif
{
  MAT2 out;
  static MAT2 wkspace, tmp;

  wkspace = m_get2(A.m, A.n);
  out = m_get2(A.m, A.n);
  if(p < 0)
  {
    tmp = m_get2(A.m, A.n);
    tmp = m_inverse2(A);
    out = matpow(tmp, -p);
  }
  else
  {
    out = matpow(A, p);
  }
  return out;
}

/* get_col2 -- gets a specified column of a matrix and retruns it as a vector */
#ifndef ANSI_C
VEC2 get_col2(mat, col, vec)
unsigned int col;
MAT2 mat;
VEC2 vec;
#else
VEC2 get_col2(const MAT2 mat, unsigned int col)
#endif
{
  VEC2 vec;
  unsigned int i;

  if(vec.dim < mat.m)
    vec = v_get2(mat.m);

  for(i = 0; i < mat.m; i++)
    vec.ve[i] = mat.me[i][col];

  return (vec);
}

/* _v_copy2 -- copies vector into new area
	-- out(i0:dim) <- in(i0:dim) */
#ifndef ANSI_C
VEC2 _v_copy2(in, out, i0)
VEC2 in, out;
unsigned int i0;
#else
VEC2 _v_copy2(const VEC2 in, unsigned int i0)
#endif
{
  VEC2 out;
  if(out.dim < in.dim)
    out = v_get2(in.dim);
  MEM_COPY(&(in.ve[i0]), &(out.ve[i0]), (in.dim - i0)*sizeof(double));
  return (out);
}

/* _in_prod2 -- inner product of two vectors from i0 downwards
   -- that is, returns a(i0:dim)^T.b(i0:dim) */
#ifndef ANSI_C
double	_in_prod2(a, b, i0)
VEC2 a, b;
unsigned int i0;
#else
double _in_prod2(const VEC2 a, const VEC2 b, unsigned int i0)
#endif
{
  unsigned int limit;
  limit = min(a.dim,b.dim);
  return __ip__(&(a.ve[i0]), &(b.ve[i0]), (int)(limit-i0));
}

#define SQRT_MAGIC_F 0x5f3759df
double sqrt2(const float x)
{
  const float xhalf = 0.5f*x;

  union // get bits for floating value
  {
    float x;
    int i;
  } u;
  u.x = x;
  u.i = SQRT_MAGIC_F - (u.i >> 1);  // gives initial guess y0
  return (double)(x*u.x*(1.5f - xhalf*u.x*u.x));// Newton step, repeating increases accuracy
}

double sqrt5(const float m)
{
   float i=0;
   float x1,x2;
   while( (i*i) <= m )
          i+=0.1f;
   x1=i;
   for(int j=0;j<10;j++)
   {
       x2=m;
      x2/=x1;
      x2+=x1;
      x2/=2;
      x1=x2;
   }
   return x2;
}

double sqrt11(const double number)
{
const double ACCURACY=0.001;
double lower, upper, guess;

 if (number < 1)
 {
  lower = number;
  upper = 1;
 }
 else
 {
  lower = 1;
  upper = number;
 }

 while ((upper-lower) > ACCURACY)
 {
  guess = (lower + upper)/2;
  if(guess*guess > number)
   upper =guess;
  else
   lower = guess;
 }
 return (lower + upper)/2;

}

double sqrt9(const double fg)
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

double fabs2(double n)
{
  if(n >= 0)
    return n; //if positive, return without ant change
  else
    return (-n); //if negative, return a positive version
}

/* hhvec2 -- calulates Householder vector to eliminate all entries after the
   i0 entry of the vector vec. It is returned as out. May be in-situ */
#ifndef ANSI_C
VEC2 hhvec2(vec, i0, beta, out, newval)
VEC2 vec, out;
unsigned int i0;
double *beta, *newval;
#else
VEC2 hhvec2(const VEC2 vec, unsigned int i0, double *beta, double *newval)
#endif
{
  VEC2 out;
  double norm, temp;
//  out = _v_copy2(vec, i0);
  out = v_copy22(vec);
  temp = (double)_in_prod2(out, out, i0);
  norm = sqrt9(temp);
  if(norm <= 0.0)
  {
    *beta = 0.0;
    return (out);
  }
  *beta = 1.0/(norm * (norm+fabs2(out.ve[i0])));
  if(out.ve[i0] > 0.0)
    *newval = -norm;
  else
    *newval = norm;
  out.ve[i0] -= *newval;

  return (out);
}

/* _hhtrcols2 -- transform a matrix by a Householder vector by columns
    starting at row i0 from column j0
    -- that is, M(i0:m,j0:n) <- (I-beta.hh(i0:m).hh(i0:m)^T)M(i0:m,j0:n)
    -- in-situ
    -- scratch vector w passed as argument
    -- raises error if w == NULL
*/
#ifndef ANSI_C
MAT2 _hhtrcols2(M, i0, j0, hh, beta, w)
MAT2 M;
unsigned int i0, j0;
VEC2 hh;
double	beta;
VEC2 w;
#else
MAT2 _hhtrcols2(MAT2 M, unsigned int i0, unsigned int j0,
               const VEC2 hh, double beta, VEC2 w)
#endif
{
  int i;

  if(beta == 0.0)
    return (M);

  w = v_zero2(w);

  for(i = i0; i < M.m; i++)
    if(hh.ve[i] != 0.0)
      __mltadd__(&(w.ve[j0]), &(M.me[i][j0]), hh.ve[i],
                 (int)(M.n-j0));
  for(i = i0; i < M.m; i++)
    if(hh.ve[i] != 0.0)
      __mltadd__(&(M.me[i][j0]), &(w.ve[j0]), -beta*hh.ve[i],
                 (int)(M.n-j0));
  return (M);
}

/* hhtrrows2 -- transform a matrix by a Householder vector by rows
    starting at row i0 from column j0 -- in-situ
    -- that is, M(i0:m,j0:n) <- M(i0:m,j0:n)(I-beta.hh(j0:n).hh(j0:n)^T) */
#ifndef ANSI_C
MAT2 hhtrrows2(M, i0, j0, hh, beta)
MAT2 M;
unsigned int i0, j0;
VEC2 hh;
double beta;
#else
MAT2 hhtrrows2(MAT2 M, unsigned int i0, unsigned int j0,
              const VEC2 hh, double beta)
#endif
{
  double ip, scale;
  int i;

  if(beta == 0.0)
    return (M);

  /* for each row ... */
  for(i = i0; i < M.m; i++)
  { /* compute inner product */
    ip = __ip__(&(M.me[i][j0]), &(hh.ve[j0]), (int)(M.n-j0));
    scale = beta*ip;
    if(scale == 0.0)
      continue;

    /* do operation */
    __mltadd__(&(M.me[i][j0]), &(hh.ve[j0]), -scale,
               (int)(M.n-j0));
  }
  return (M);
}

/* Hfactor2 -- compute Hessenberg factorization in compact form.
    -- factorization performed in situ
*/
#ifndef ANSI_C
MAT2 Hfactor2(A, diag, beta)
MAT2 A;
VEC2 diag, beta;
#else
MAT2 Hfactor2(MAT2 A, VEC2 diag, VEC2 beta)
#endif
{
  static VEC2 hh, w;
  int k, limit;
  double b;

  limit = A.m - 1;
  hh = v_get2(A.m);
  w  = v_get2(A.n);

  for(k = 0; k < limit; k++)
  {
    /* compute the Householder vector hh */
	hh = get_col2(A, (unsigned int)k);
	hh = hhvec2(hh, k+1, &beta.ve[k], &A.me[k+1][k]);
    diag.ve[k] = (((k+1)>=0 && (k+1)<(hh).dim) ? (hh).ve[(k+1)] : \
        0);
//    v_set_val2(diag, k, v_entry2(hh, k+1));

    /* apply Householder operation symmetrically to A */
    b = v_entry2(beta, k);
    A = _hhtrcols2(A, k+1, k+1, hh, b, w);
    A = hhtrrows2(A, 0, k+1, hh, b);
  }
  return (A);
}

/* hhtrvec2 -- apply Householder transformation to vector
    -- that is, out <- (I-beta.hh(i0:n).hh(i0:n)^T).in
    -- may be in-situ */
#ifndef ANSI_C
VEC2 hhtrvec2(hh, beta, i0, in, out)
VEC2 hh, in, out;	/* hh = Householder vector */
unsigned int i0;
double beta;
#else
VEC2 hhtrvec2(const VEC2 hh, double beta, unsigned int i0,
             const VEC2 in)
#endif
{
  VEC2 out;
  double scale, temp;
  temp = (double)_in_prod2(hh, in, i0);
  scale = beta*temp;
  out = v_copy2(in);
  __mltadd__(&(out.ve[i0]), &(hh.ve[i0]), -scale, (int)(in.dim-i0));

  return (out);
}

/* makeHQ2 -- construct the Hessenberg orthogonalising matrix Q;
    -- i.e. Hess M = Q.M.Q'	*/
#ifndef ANSI_C
MAT2 makeHQ2(H, diag, beta)
MAT2 H;
VEC2 diag, beta;
#else
MAT2 makeHQ2(MAT2 H, VEC2 diag, VEC2 beta)
#endif
{
  MAT2 Qout;
  int i, j, limit;
  static VEC2 tmp1, tmp2;
  Qout = m_get2(H.m, H.m);
  tmp1 = v_get2(H.m);
  tmp2 = v_get2(H.m);;
  for(i = 0; i < H.m; i++)
  {
    /* tmp1 = i'th basis vector */
    for(j = 0; j < H.m; j++)
      tmp1.ve[j] = 0.0;
    tmp1.ve[i] = 1.0;

    /* apply H/h transforms in reverse order */
    for(j = limit-1; j >= 0; j--)
    {
      tmp2 = get_col2(H, (unsigned int)j);
      tmp2.ve[j+1] = (((j)>=0 && (j)<(diag).dim) ? (diag).ve[(j)] : \
          0);
//      v_set_val2(tmp2, j+1, v_entry2(diag, j));
      tmp1 = hhtrvec2(tmp2, beta.ve[j], j+1, tmp1);
    }

    /* insert into Qout */
    Qout = set_col2(Qout, (unsigned int)i, tmp1);
  }
  return (Qout);
}

/* makeH2 -- construct actual Hessenberg matrix */
#ifndef ANSI_C
MAT2 makeH2(H)
MAT2 H;
#else
MAT2 makeH2(const MAT2 H)
#endif
{
  MAT2 Hout;
  int i, j, limit;

  Hout = m_get2(H.m, H.m);
  Hout = m_copy2(H);

  limit = H.m;
  for(i = 1; i < limit; i++)
    for(j = 0; j < i-1; j++)
    	Hout.me[i][j] = 0.0;

  return (Hout);
}

/* rot_cols2 -- postmultiply mat by givens rotation described by c, s */
#ifndef ANSI_C
MAT2 rot_cols2(mat, i, k, c, s)
MAT2 mat;
unsigned int i, k;
double c, s;
#else
MAT2 rot_cols2(const MAT2 mat, unsigned int i, unsigned int k,
              double c, double s)
#endif
{
  MAT2 out;
  unsigned int j;
  double temp;

  out = m_copy2(mat);

  for(j=0; j<mat.m; j++)
  {
    temp = c*m_entry2(out, j, i) + s*m_entry2(out, j, k);
//    m_set_val2(out, j, k, -s*m_entry2(out, j, i) + c*m_entry2(out, j, k));
    out.me[j][k] = (-s*(((j)>=0 && (j)<(out).m && (i)>=0 && (i)<=(out).n) ? \
        (out).me[(j)][(i)] : (printf("Error!"))) + c*(((j)>=0 && (j)<(out).m && (k)>=0 && (k)<=(out).n) ? \
        (out).me[(j)][(k)] : (printf("Error!"))));
//    m_set_val2(out, j, i, temp);
    out.me[j][i] = temp;
  }

  return (out);
}

/* rot_rows2 -- premultiply mat by givens rotation described by c, s */
#ifndef ANSI_C
MAT2 rot_rows2(mat, i, k, c, s)
MAT2 mat;
unsigned int i, k;
double c, s;
#else
MAT2 rot_rows2(const MAT2 mat, unsigned int i, unsigned int k,
              double c, double s)
#endif
{
  MAT2 out;
  unsigned int j;
  double temp;

  out = m_copy2(mat);

  for(j=0; j<mat.n; j++)
  {
    temp = c*m_entry2(out, i, j) + s*m_entry2(out, k, j);
//    m_set_val(out, k, j, -s*m_entry2(out, i, j) + c*m_entry2(out, k, j));
    out.me[k][j] = (-s*(((i)>=0 && (i)<(out).m && (j)>=0 && (j)<=(out).n) ? \
        (out).me[(i)][(j)] : (printf("Error!"))) + c*(((k)>=0 && (k)<(out).m && (j)>=0 && (j)<=(out).n) ? \
        (out).me[(k)][(j)] : (printf("Error!"))));
//    m_set_val(out, i, j, temp);
    out.me[i][j] = temp;
  }

  return (out);
}

/* hhldr3 -- computes */
#ifndef ANSI_C
static void hhldr3(x, y, z, nu1, beta, newval)
double x, y, z;
double *nu1, *beta, *newval;
#else
static void hhldr3(double x, double y, double z,
                   double *nu1, double *beta, double *newval)
#endif
{
  double alpha;

  if(x >= 0.0)
    alpha = sqrt9(x*x+y*y+z*z);
  else
    alpha = -sqrt9(x*x+y*y+z*z);
  *nu1 = x + alpha;
  *beta = 1.0/(alpha*(*nu1));
  *newval = alpha;
}

/*hhldr3rows2 */
#ifndef ANSI_C
static MAT2 hhldr3rows2(A, k, i0, beta, nu1, nu2, nu3)
MAT2 A;
int	k, i0;
double beta, nu1, nu2, nu3;
#else
static MAT2 hhldr3rows2(MAT2 A, int k, int i0, double beta,
                       double nu1, double nu2, double nu3)
#endif
{
  double ip, prod;
  int i, m;
  m = A.m;
  i0 = min(i0, m-1);

  for(i = 0; i <= i0; i++)
  {
    ip = nu1*m_entry2(A, i, k) + nu2*m_entry2(A, i, k+1)+nu3 * m_entry2(A, i, k+2);
    prod = ip*beta;
//    m_add_val(A, i, k, -prod*nu1);
    A.me[i][k] += (-prod*nu1);
//    m_add_val(A, i, k+1, -prod*nu2);
    A.me[i][k+1] += (-prod*nu2);
//    m_add_val(A, i, k+2, -prod*nu3);
    A.me[i][k+2] += (-prod*nu3);
  }
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

  norm = sqrt9(x*x+y*y);
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

/* schur2 -- computes the Schur decomposition of the matrix A in situ
    -- optionally, gives Q matrix such that Q^T.A.Q is upper triangular
    -- returns upper triangular Schur matrix */
#ifndef ANSI_C
MAT2 schur2(A, Q)
MAT2 A, Q;
#else
MAT2 schur2(MAT2 A, MAT2 Q)
#endif
{
  int i, j, iter, k, k_min, k_max, k_tmp, n, split;
  double beta2, c, discrim, dummy, nu1, s, t, tmp, x, y, z;
  double A_me[MAX_SIZE][MAX_SIZE];
  double sqrt_macheps;
  static VEC2 diag, beta;
  n = A.n;
  diag = v_get2(A.n);
  beta = v_get2(A.n);
  /* compute Hessenberg form */
  A = Hfactor2(A, diag, beta);
  /* save Q if necessary */
  if(&Q)
    Q = makeHQ2(A, diag, beta);
  A = makeH2(A);
  sqrt_macheps = sqrt9(MACHEPS);

  k_min = 0;
  memcpy(A_me, A.me, sizeof(A.me));
//  A_me = A.me;
  while(k_min < n)
  {
    double a00, a01, a10, a11;
    double scale, t, numer, denom;

    /* find k_max to suit:
       submatrix k_min..k_max should be irreducible */
    k_max = n-1;
    for(k = k_min; k < k_max; k++)
      if(m_entry2(A, k+1, k) == 0.0)
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
    if(k_max == k_min + 1)
    {
      a00 = m_entry2(A, k_min, k_min);
      a01 = m_entry2(A, k_min, k_max);
      a10 = m_entry2(A, k_max, k_min);
      a11 = m_entry2(A, k_max, k_max);
      tmp = a00 - a11;
      discrim = tmp*tmp + 4*a01*a10;
      if(discrim < 0.0)
      {
        /* yes -- e-vals are complex
               -- put 2 x 2 block in form [a b; c a];
        then eigenvalues have real part a & imag part sqrt(|bc|) */
        numer = - tmp;
        denom = (a01+a10 >= 0.0) ?
                (a01+a10) + sqrt9((a01+a10)*(a01+a10)+tmp*tmp) :
                (a01+a10) - sqrt9((a01+a10)*(a01+a10)+tmp*tmp);
        if(denom != 0.0)
        {    /* t = s/c = numer/denom */
          t = numer/denom;
          scale = c = 1.0/sqrt9(1+t*t);
          s = c*t;
        }
        else
        {
          c = 1.0;
          s = 0.0;
        }
        A = rot_cols2(A, k_min, k_max, c, s);
        A = rot_rows2(A, k_min, k_max, c, s);
//        if(&Q != MNULL2)
          Q = rot_cols2(Q, k_min, k_max, c, s);
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
              - tmp - sqrt9(discrim) : - tmp + sqrt9(discrim);
        denom = 2*a01;
        if(fabs2(numer) < fabs2(denom))
        {    /* t = s/c = numer/denom */
          t = numer/denom;
          scale = c = 1.0/sqrt9(1+t*t);
          s = c*t;
        }
        else if(numer != 0.0)
        {    /* t = c/s = denom/numer */
          t = denom/numer;
          scale = 1.0/sqrt9(1+t*t);
          c = fabs2(t)*scale;
          s = (t >= 0.0) ? scale : -scale;
        }
        else /* numer == denom == 0 */
        {
          c = 0.0;
          s = 1.0;
        }
        A = rot_cols2(A, k_min, k_max, c, s);
        A = rot_rows2(A, k_min, k_max, c, s);
//        if(&Q != MNULL2)
          Q = rot_cols2(Q, k_min, k_max, c, s);
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

      a00 = m_entry2(A, k_tmp, k_tmp);
      a01 = m_entry2(A, k_tmp, k_max);
      a10 = m_entry2(A, k_max, k_tmp);
      a11 = m_entry2(A, k_max, k_max);

      /* treat degenerate cases differently
         -- if there are still no splits after five iterations
            and the bottom 2 x 2 looks degenerate, force it to
         split */
      #ifdef DEBUG
        printf("# schur: bottom 2 x 2 = [%lg, %lg; %lg, %lg]\n",
               a00, a01, a10, a11);
      #endif
      if(iter >= 5 &&
         fabs2(a00-a11) < sqrt_macheps*(fabs2(a00)+fabs2(a11)) &&
         (fabs2(a01) < sqrt_macheps*(fabs2(a00)+fabs2(a11)) ||
          fabs2(a10) < sqrt_macheps*(fabs2(a00)+fabs2(a11))) )
      {
        if(fabs2(a01) < sqrt_macheps*(fabs2(a00)+fabs2(a11)))
        {
//          m_set_val2(A, k_tmp, k_max, 0.0);
          A.me[k_tmp][k_max] = 0.0;
        }
        if(fabs2(a10) < sqrt_macheps*(fabs2(a00)+fabs2(a11)))
        {
//          m_set_val2(A, k_max, k_tmp, 0.0);
          A.me[k_max][k_tmp] = 0.0;
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

      a00 = m_entry2(A, k_min, k_min);
      a01 = m_entry2(A, k_min, k_tmp);
      a10 = m_entry2(A, k_tmp, k_min);
      a11 = m_entry2(A, k_tmp, k_tmp);

      x = a00*a00 + a01*a10 - s*a00 + t;
      y = a10*(a00+a11-s);
      if(k_min + 2 <= k_max)
        z = a10*A.me[k_min+2][k_tmp];
      else
        z = 0.0;

      for(k = k_min; k <= k_max-1; k++)
      {
        if(k < k_max - 1)
        {
          hhldr3(x, y, z, &nu1, &beta2, &dummy);
//          if(&Q != MNULL2)
            Q = hhldr3rows2(Q, k, n-1, beta2, nu1, y, z);
        }
        else
        {
          givens(x, y, &c, &s);
          A = rot_cols2(A, k, k+1, c, s);
          A = rot_rows2(A, k, k+1, c, s);
          if(&Q)
            Q = rot_cols2(Q, k, k+1, c, s);
        }
        x = m_entry2(A, k+1, k);
        if(k <= k_max - 2)
          y = m_entry2(A, k+2, k);
        else
          y = 0.0;
        if(k <= k_max - 3)
          z = m_entry2(A, k+3, k);
        else
          z = 0.0;
      }
	  for(k = k_min; k <= k_max-2; k++)
	  {
        /* zero appropriate sub-diagonals */
//        m_set_val(A, k+2, k, 0.0);
		A.me[k+2][k] = 0.0;
        if(k < k_max-2)
        {
//	      m_set_val(A, k+3, k, 0.0);
          A.me[k+3][k] = 0.0;
        }
      }

      /* test to see if matrix should split */
      for(k = k_min; k < k_max; k++)
        if(fabs2(A_me[k+1][k]) < MACHEPS*
          (fabs2(A_me[k][k])+fabs2(A_me[k+1][k+1])))
        {
          A_me[k+1][k] = 0.0;
          split = 1;
        }
	}
  }

  /* polish up A by zeroing strictly lower triangular elements
     and small sub-diagonal elements */
  for(i = 0; i < A.m; i++)
    for(j = 0; j < i-1; j++)
      A_me[i][j] = 0.0;
    for(i = 0; i < A.m - 1; i++)
      if(fabs2(A_me[i+1][i]) < MACHEPS*
         (fabs2(A_me[i][i])+fabs2(A_me[i+1][i+1])))
        A_me[i+1][i] = 0.0;

//#ifdef THREADSAFE
//  v_free2(&diag);  v_free2(&beta);
//#endif

  return A;
}

/* schur_vals2 -- compute real & imaginary parts of eigenvalues
	-- assumes T contains a block upper triangular matrix
		as produced by schur()
	-- real parts stored in real_pt, imaginary parts in imag_pt */
#ifndef ANSI_C
void schur_evals2(T, real_pt, imag_pt)
MAT2 *T;
VEC2 *real_pt, *imag_pt;
#else
void schur_evals2(MAT2 *T, VEC2 *real_pt, VEC2 *imag_pt)
#endif
{
  int i, n;
  double discrim, T_me[MAX_SIZE][MAX_SIZE];
  double diff, sum, tmp;

  n = T->n;
  memcpy(T_me, T->me, sizeof(T->me));
//  T_me = T->me;
//  real_pt = v_resize(real_pt, (unsigned int)n);
//  imag_pt = v_resize(imag_pt, (unsigned int)n);

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
        real_pt->ve[i] = real_pt->ve[i+1] = sum;
        imag_pt->ve[i] = sqrt9(-discrim);
        imag_pt->ve[i+1] = - imag_pt->ve[i];
      }
      else
      { /* no -- actually both real */
        tmp = sqrt9(discrim);
        real_pt->ve[i]   = sum + tmp;
        real_pt->ve[i+1] = sum - tmp;
        imag_pt->ve[i]   = imag_pt->ve[i+1] = 0.0;
      }
      i += 2;
    }
    else
    {   /* real eigenvalue */
      real_pt->ve[i] = T_me[i][i];
      imag_pt->ve[i] = 0.0;
      i++;
    }
  }
}

/* m_get_eigenvalues2 -- get the eigenvalues of a matrix A
	-- */
CMPLX *m_get_eigenvalues2(MAT2 A)
{
  MAT2 T, Q;
  VEC2 evals_re, evals_im;

  static CMPLX z[MAX_SIZE];
  Q = m_get2(A.m, A.n);
  T = m_copy2(A);
  /* compute Schur form: A = Q.T.Q^T */
  T = schur2(T, Q);
  /* extract eigenvalues */
  evals_re = v_get2(A.m);
  evals_im = v_get2(A.m);
  schur_evals2(&T, &evals_re, &evals_im);
//  z = malloc(evals_re->dim*sizeof(CMPLX));
  for(int i = 0; i < evals_re.dim; i++)
  {
//    z[i] = evals_re->ve[i] + I*evals_im->ve[i];
    z[i].real = evals_re.ve[i];
    z[i].imag = evals_im.ve[i];
  }
  return z;
}

double cmplx_mag(double real, double imag)
{
  return sqrt9(real * real + imag * imag);
}

/* max_mag_eigenvalue -- extracts the magnitude of the maximum eigenvalue
	-- */
double max_mag_eigenvalue(CMPLX *z, int size)
{
  double maximum = 0, aux;
  for(int c = 1; c < size; c++)
  {
    aux = cmplx_mag(z[c].real, z[c].imag);
    if(aux > maximum)
    {
      maximum  = aux;
    }
  }
  return (double)maximum;
}

/* is_same_sign -- check if a has the same sign as b
	-- */
int is_same_sign(double a, double b)
{
  if(((a >= 0) && (b >= 0)) || ((a <= 0) && (b <= 0)))
    return 1;
  else
    return 0;
}

/* y_k2 -- computes the output signal in the k-th sample
	-- */
double y_k2(MAT2 A, MAT2 B, MAT2 C, MAT2 D, double u, int k, MAT2 x0)
{
  MAT2 y, Ak, AUX, AUX2;
  MAT2 AUX3, AUX4, AUX5;
//  U = m_get(A->m, A->n);
//  U = m_same_elements(U,u);
  // y = C * A.pow(k) * x0;
  Ak = m_get2(A.m, A.n);
  Ak = m_pow3(A, k);
  AUX = m_get2(A.m, A.n);
  AUX = m_mlt2(C, Ak);
  y = m_get2(A.m, A.n);
  y = m_mlt2(AUX, x0);

  AUX2 = m_get2(A.m, A.n);
  for(int m = 0; m <= (k - 1); m++)
  {
    // y += (C * A.pow(k - m - 1) * B * u) + D * u;
    Ak = m_pow3(A, (k-m-1));
    AUX = m_mlt2(C, Ak);
    AUX2 = m_mlt2(AUX, B);
//    AUX3 = m_mlt(AUX2, U, MNULL);
//    AUX4 = m_mlt(D, U, MNULL);
    AUX5 = m_add2(AUX2, D);
    y = m_add2(y, AUX5);
  }
  return y.me[0][0]*u;
}

/* peak_output2 -- computes the biggest peak value of a signal (Mp)
	-- */
PKVL peak_output2(MAT2 A, MAT2 B, MAT2 C, MAT2 D, MAT2 x0, double yss, double u)
{
  PKVL out;
  double greater;
  int i = 1;
  greater = fabs2(y_k2(A, B, C, D, u, i, x0));
  while((fabs2(y_k2(A, B, C, D, u, i+1, x0)) >= fabs2(yss)))
  {
    if(greater < fabs2(y_k2(A, B, C, D, u, i+1, x0)))
    {
      greater = fabs2(y_k2(A, B, C, D, u, i+1, x0));
      out.mp = y_k2(A, B, C, D, u, i+1, x0);
      out.kp = i+2;
    }
    if(!is_same_sign(yss, out.mp))
    {
      greater = 0;
    }
    i++;
  }
//  printf("Mp=%f e kp=%d\n", out.mp, out.kp);
  return out;
}

double y_ss2(MAT2 A, MAT2 B, MAT2 C, MAT2 D, double u)
{
  double yss;
  MAT2 AUX, AUX2, AUX3, AUX4, AUX5;
  MAT2 Id;
  // get the expression y_ss=(C(I-A)^(-1)B+D)u
  Id = m_get2(A.m, A.n);
  Id = m_ident2(Id.m);/*printf("I\n");m_output(Id);*/
  AUX = m_get2(A.m, A.n);
  // Id - A
  AUX = m_sub2(Id, A);/*printf("I-A\n");m_output(AUX);*/
  AUX2 = m_get2(A.m, A.n);
  AUX2 = m_inverse2(AUX);/*printf("(I-A)^(-1)\n");m_output(AUX2);*/
  AUX3 = m_get2(A.m, A.n);
  AUX3 = m_mlt2(C, AUX2);/*printf("C(I-A)^(-1))\n");m_output(AUX3);*/
  AUX4 = m_get2(A.m, A.n);
  AUX4 = m_mlt2(AUX3, B);/*printf("C(I-A)^(-1)B\n");m_output(AUX4);*/
  AUX5 = m_get2(A.m, A.n);
  AUX5 = m_add2(AUX4, D);/*printf("(C(I-A)^(-1)B+D)\n");m_output(AUX5);*/
  yss = AUX5.me[0][0] * u;/*printf("yss=\n");*/

  return yss;
}

double pow2(double a, int n)
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

double c_bar(double mp, double yss, double lambmax, int kp)
{
  double cbar;
  cbar = (mp-yss)/(pow2(lambmax, kp));
  return cbar;
}

/**
 * Calculate ln logarithm using integers with 16 bit precision
 * min: fxp_ln(0.000015259<<16)
 * max: fxp_ln(32767<<16)
 */
int fxp_ln(int x)
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
 * min: fxp_log10(0.000015259)
 * max: fxp_log10(32767.0)
 */
double fxp_log10_low(double x)
{
  int xint = (int) (x * 65536.0 + 0.5);
  int lnum = fxp_ln(xint);
  int lden = fxp_ln(655360);
  return ((double) lnum / (double) lden);
}

/**
 * Calculate log10 logarithm using 16 bit precision
 * min: fxp_log10(0.000015259)
 * max: fxp_log10(2147483647.0)
 */
double fxp_log10(double x)
{
  if(x > 32767.0)
  {
    if(x > 1073676289.0)
    {
      x = x / 1073676289.0;
      return fxp_log10_low(x) + 9.030873362;
    }
    x = x / 32767.0;
    return fxp_log10_low(x) + 4.515436681;
  }
  return fxp_log10_low(x);
}

double log_b(double base, double x)
{
  return (double) (fxp_log10(x) / fxp_log10(base));
}

double ceil2(double x)
{
  union
  {
    float f;
    int i;
  } float_int;

//    float_int val;
  float_int.f=x;

  // Extract sign, exponent and mantissa
  // Bias is removed from exponent
  int sign=float_int.i >> 31;
  int exponent=((float_int.i & 0x7fffffff) >> 23) - 127;
  int mantissa=float_int.i & 0x7fffff;

  // Is the exponent less than zero?
  if(exponent<0)
  {
    // In this case, x is in the open interval (-1, 1)
    if(x<=0.0f)
      return 0.0f;
    else
      return 1.0f;
  }
  else
  {
    // Construct a bit mask that will mask off the
    // fractional part of the mantissa
    int mask=0x7fffff >> exponent;

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
        mantissa+=1 << (23-exponent);

        // Did the mantissa overflow?
        if(mantissa & 0x800000)
        {
          // The mantissa can only overflow if all the
          // integer bits were previously 1 -- so we can
          // just clear out the mantissa and increment
          // the exponent
          mantissa=0;
          exponent++;
        }
      }

      // Clear the fractional bits
      mantissa&=~mask;
    }
  }

  // Put sign, exponent and mantissa together again
  float_int.i=(sign << 31) | ((exponent+127) << 23) | mantissa;

  return (double)float_int.f;
}

int k_bar(double lambdaMax, double p, double cbar, double yss, int order)
{
  double k_ss, x;
  x = (p * yss) / (100 * cbar);
  k_ss = log_b(lambdaMax, x);
  return ceil2(k_ss)+order;
}

double max_mag_eigenvalue22(MAT2 A)
{
  double maximum = 0, aux;
  CMPLX *z;
  z = m_get_eigenvalues2(A);
  for(int i = 0; i < A.m; i++)
  {
    aux = cmplx_mag(z[i].real, z[i].imag);
    if(aux > maximum)
    {
      maximum = aux;
    }
  }
  return maximum;
}

int check_settling_time2(MAT2 A, MAT2 B, MAT2 C, MAT2 D, MAT2 x0,
                        double u, double tsr, double p, double ts)
{
  double yss, mp, lambMax, cbar, output;
  PKVL out;
  int kbar, kp, i;
  yss = y_ss2(A, B, C, D, u);
  out = peak_output2(A, B, C, D, x0, yss, u);
  mp = out.mp;
  kp = out.kp;
  lambMax = max_mag_eigenvalue22(A);
  printf("Mp=%f", mp);
  printf("yss=%f", yss);
  printf("lambMax=%f", lambMax);
  printf("kp=%d", kp);

  cbar = c_bar(mp, yss, lambMax, kp);

  kbar = k_bar(lambMax, p, cbar, yss, A.m);
  printf("cbar=%f", cbar);
  if(kbar * ts < tsr)
  {
    printf("kbar=%i", kbar);
    return 1;
  }

  i = ceil2(tsr / ts);
  while(i <= kbar)
  {
    output = y_k2(A, B, C, D, u, i, x0);
    if(!(output > (yss - (yss * (p/100))) && (output < (yss * (p/100) + yss))))
    {
      printf("kbar=%i", kbar);
      return 0;
    }
    i++;
  }
  printf("kbar=%i", kbar);
  return 1;
}

int main(){
	CMPLX *z, *mz;
	CMPLX my_z[MAX_SIZE];

    MAT2 T = m_get2(5,5);
    T.me[0][0]=-0.5000;T.me[0][1]=0.6000;T.me[0][2]=0;T.me[0][3]=0;T.me[0][4]=0;
    T.me[1][0]=-0.6000;T.me[1][1]=-0.5000;T.me[1][2]=0;T.me[1][3]=0;T.me[1][4]=0;
    T.me[2][0]=0;T.me[2][1]=0;T.me[2][2]=0.2000;T.me[2][3]=0.8000;T.me[2][4]=0;
    T.me[3][0]=0;T.me[3][1]=0;T.me[3][2]=-0.8000;T.me[3][3]=0.2000;T.me[3][4]=0;
    T.me[4][0]=0;T.me[4][1]=0;T.me[4][2]=0;T.me[4][3]=0;T.me[4][4]=0.6;printf("A ");
    m_output2(T);
    printf("A+A ");m_output2(m_add2(T, T));
    printf("A-A ");m_output2(m_sub2(T, T));
    printf("A*A=\n");m_output2(m_mlt2(T, T));
    printf("inv(A)=\n");m_output2(m_inverse2(T));
    printf("pow(A, 4)=\n");m_output2(m_pow2(T, 4));
    printf("pow2(A, 4)=\n");m_output2(m_power_opt(T, 4));
    printf("pow3(A, 4)=\n");m_output2(matpow(T, 4));
    printf("m_pow3(A, 4)=\n");m_output2(m_pow3(T, 4));

   //setting up A matrix
//    A=m_get(4,4);
//    A->me[0][0]=-0.5000;A->me[0][1]=0.6000;A->me[0][2]=0;A->me[0][3]=0;
//    A->me[1][0]=-0.6000;A->me[1][1]=-0.5000;A->me[1][2]=0;A->me[1][3]=0;
//    A->me[2][0]=0;A->me[2][1]=0;A->me[2][2]=0.2000;A->me[2][3]=0.8000;
//    A->me[3][0]=0;A->me[3][1]=0;A->me[3][2]=-0.8000;A->me[3][3]=0.2000;
//    printf("A ");m_output(A);
//    A=m_get(5,5);
//    A->me[0][0]=-0.5000;A->me[0][1]=0.6000;A->me[0][2]=0;A->me[0][3]=0;A->me[0][4]=0;
//    A->me[1][0]=-0.6000;A->me[1][1]=-0.5000;A->me[1][2]=0;A->me[1][3]=0;A->me[1][4]=0;
//    A->me[2][0]=0;A->me[2][1]=0;A->me[2][2]=0.2000;A->me[2][3]=0.8000;A->me[2][4]=0;
//    A->me[3][0]=0;A->me[3][1]=0;A->me[3][2]=-0.8000;A->me[3][3]=0.2000;A->me[3][4]=0;
//    A->me[4][0]=0;A->me[4][1]=0;A->me[4][2]=0;A->me[4][3]=0;A->me[4][4]=0.6;printf("A ");
//    m_output(A);

    //setting up B matrix
//    B=m_get(4,1);
//    B->me[0][0]=0;
//    B->me[1][0]=0;
//    B->me[2][0]=2.5;
//    B->me[3][0]=1;printf("B ");m_output(B);
    B=m_get(5,1);
    B->me[0][0]=0;
    B->me[1][0]=0;
    B->me[2][0]=2.5;
    B->me[3][0]=1;
    B->me[4][0]=0;printf("B ");m_output(B);
    //setting up C matrix
//    C=m_get(1,4);
//    C->me[0][0]=0;C->me[0][1]=2.6;C->me[0][2]=0.5;C->me[0][3]=1.2;
//    printf("C ");m_output(C);
    C=m_get(1,5);
    C->me[0][0]=0;C->me[0][1]=2.6;C->me[0][2]=0.5;C->me[0][3]=1.2;C->me[0][4]=0;
    printf("C ");m_output(C);
    //setting up D matrix
    D=m_get(1,1);
    m_zero(D);printf("D ");m_output(D);
//    D->me[0][0] = 0;printf("D ");m_output(D);
    x0=m_get(5,1);
    m_zero(x0);printf("x0 ");m_output(x0);
    printf("-----------------------------------------------------------\n");
    A2=m_get(5,5);
    A2 = m_add(A, A, A2);printf("A+A=\n");
    m_output(A2);
    A3=m_get(5,5);
    A3 = m_sub(A, A, A3);printf("A-A=\n");
    m_output(A3);
    A4=m_get(5,5);
    A4 = m_mlt(A, A, A4);printf("A*A=\n");
    m_output(A4);
    A5=m_get(5,5);
    A5 = m_inverse(A,A5);printf("inv(A)=\n");
    m_output(A5);
    A6=m_get(5,5);
    A6 = m_pow(A,4,A6);printf("pow(A, 4)=\n");
    m_output(A6);

    z = m_get_eigenvalues(A);
    printf("test\n");
    mz = m_get_eigenvalues2(T);
    printf("test2\n");
    int size = A->m, s = T.m;
    for(int i=0;i<size;i++){
    //		printf("%f+%f i", z[i].real, z[i].imag);
      printfc(z[i]);
    }
    double lambmax = max_mag_eigenvalue(z, size);
    printf("Maximum:%f\n", lambmax);
    double lambmax2 = max_mag_eigenvalue(z, s);
    printf("Maximum2:%f\n", lambmax2);

//    PKVL out;
//    double yss = y_ss(A, B, C, D, 1.0);
//    out = peak_output(A, B, C, D, x0, yss, 1.0);
//    double mp = out.mp;
//    int kp = out.kp;
//    printf("Mp=%f\n", mp);
//    printf("kp=%d\n", out.kp);
//    double cbar = c_bar(mp, yss, lambmax, kp);
//    printf("c_bar=%f\n", cbar);
//    int kbar = k_bar(lambmax, 5, cbar, yss, A->m);
//    printf("k_bar=%d\n", kbar);
//    printf("y_ss=%f\n", yss);
//    printf("result log2: %f\n", sqrt2(2));
////    printf("result log: %f\n", sqrt(2));
//    printf("result log5: %f\n", sqrt5(2));
//    printf("result log11: %f\n", sqrt11(2));
//    printf("result log9: %f\n", sqrt9(2));
//    printf("result ceil2: %f\n", ceil2(-2.34));
//    // Testing
//    assert(lambmax == 0.824621);
//    assert(mp == -1.330000);
//    assert(kp == 4);
//    assert(cbar == -2.808715397923877);
//    assert(kbar == 44);
//    assert(yss == -0.031250000000000);
    return 0;
}
