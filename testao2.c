#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <strings.h>
void bcopy(const void *src, void *dest, size_t n);


#define	MEM_COPY(from,to,size)	bcopy((char *)(from),(char *)(to),(int)(size))

/* miscellaneous constants */
#define	VNULL	((VEC *)NULL)
#define	MNULL	((MAT *)NULL)
#define	PNULL	((PERM *)NULL)

/* available standard types */
#define TYPE_NULL              (-1)
#define TYPE_MAT    	        0
#define TYPE_BAND               1
#define TYPE_PERM		2
#define TYPE_VEC		3
#define TYPE_IVEC		4

/* max number of lists of types */
#define MEM_CONNECT_MAX_LISTS    3

#define	v_chk_idx(x,i)		((i)>=0 && (i)<(x)->dim)

#define	v_get_val(x,i)	( v_chk_idx(x,i) ? (x)->ve[(i)] : \
	(printf("Error!")))

#define	v_entry(x,i)		v_get_val(x,i)

#define CHAR0ISDBL0 1
#define ANSI_C 1

/* standard copy & zero functions */
#define	MEM_COPY(from,to,size)	memmove((to),(from),(size))
#define	MEM_ZERO(where,size)	memset((where),'\0',(size))

/* allocate num objects of given type */
#define	NEW_A(num,type)	((type *)calloc((size_t)(num),(size_t)sizeof(type)))

 /* re-allocate arry to have num objects of the given type */
#define	RENEW(var,num,type) \
    ((var)=(type *)((var) ? \
		    realloc((char *)(var),(size_t)((num)*sizeof(type))) : \
		    calloc((size_t)(num),(size_t)sizeof(type))))

/* type independent min and max operations */
#ifndef max
#define	max(a,b)	((a) > (b) ? (a) : (b))
#endif /* max */
#ifndef min
#define	min(a,b)	((a) > (b) ? (b) : (a))
#endif /* min */

#define mem_bytes(type,old_size,new_size)  \
  mem_bytes_list(type,old_size,new_size,0)

#define	v_set_val(x,i,val)	(x->ve[i] = val)

/* macros */

#define mem_info()   mem_info_file(stdout,0)

#ifndef THREADSAFE
#define mem_stat_reg(var,type)  mem_stat_reg_list((void **)var,type,0,__FILE__,__LINE__)
#define MEM_STAT_REG(var,type)  mem_stat_reg_list((void **)&(var),type,0,__FILE__,__LINE__)
#define mem_stat_free(mark)   mem_stat_free_list(mark,0)
#else
#define mem_stat_reg(var,type)
#define MEM_STAT_REG(var,type)
#define mem_stat_free(mark)
#endif

#define mem_bytes(type,old_size,new_size)  \
  mem_bytes_list(type,old_size,new_size,0)

#define mem_numvar(type,num) mem_numvar_list(type,num,0)

/* MACROS */
#define	m_copy(in,out)	_m_copy(in,out,0,0)
#define	v_copy(in,out)	_v_copy(in,out,0)

#define	m_set_val(A,i,j,val)	((A)->me[(i)][(j)] = (val))
#define	m_get_val(A,i,j)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] : (printf("Error!")))

#define	m_entry(A,i,j)		m_get_val(A,i,j)
#define	v_entry(x,i)		v_get_val(x,i)

#undef TRUE
#define	TRUE	1
#undef FALSE
#define	FALSE	0

#define	m_add_val(A,i,j,val)	((A)->me[(i)][(j)] += (val))

#define MACHEPS 2.22045e-16

/* macros that also check types and sets pointers to NULL */
//#define	M_FREE(mat)	( m_free(mat),	(mat)=(MAT *)NULL )
//#define V_FREE(vec)	( v_free(vec),	(vec)=(VEC *)NULL )
//#define	PX_FREE(px)	( px_free(px),	(px)=(PERM *)NULL )

#define	m_output(mat)	m_foutput(stdout,mat)

#define NOT_SEGMENTED 1

#ifndef NOT_SEGMENTED
#define	SEGMENTED
#endif

#ifndef THREADSAFE	/* for use as a shared library */
#define	THREADSAFE 1
#endif

/* vector definition */
typedef	struct	{
		unsigned int	dim, max_dim;
		double	*ve;
} VEC;

/* matrix definition */
typedef	struct	{
		unsigned int	m, n;
		unsigned int	max_m, max_n, max_size;
		double	**me,*base;	/* base is base of alloc'd mem */
} MAT;

/* band matrix definition */
typedef struct {
               MAT   *mat;       /* matrix */
               int   lb,ub;    /* lower and upper bandwidth */
} BAND;


/* permutation definition */
typedef	struct	{
		unsigned int	size, max_size, *pe;
} PERM;

/* integer vector definition */
typedef struct	{
		unsigned int	dim, max_dim;
		int	*ive;
} IVEC;

/* structure for memory information */
typedef struct {
   long bytes;       /* # of allocated bytes for each type (summary) */
   int  numvar;      /* # of allocated variables for each type */
} MEM_ARRAY;

typedef struct {
   char **type_names;        /* array of names of types (strings) */
   int  (**free_funcs)();    /* array of functions for releasing types */
   unsigned ntypes;          /* max number of types */
   MEM_ARRAY *info_sum;      /* local array for keeping track of memory */
} MEM_CONNECT;

/* __zero__ -- zeros an array of floating point numbers */
#ifndef ANSI_C
void	__zero__(dp,len)
register double	*dp;
register int	len;
#else
void	__zero__(double *dp, int len)
#endif
{
#ifdef CHAR0ISDBL0
    /* if a floating point zero is equivalent to a string of nulls */
    MEM_ZERO((char *)dp,len*sizeof(double));
#else
    /* else, need to zero the array entry by entry */
    int	i;
    for ( i = 0; i < len; i++ )
	  dp[i] = 0.0;
#endif
}

/* names of types */
static char *mem_type_names[] = {
   "MAT",
   "PERM",
   "VEC"};

#define MEM_NUM_STD_TYPES  (sizeof(mem_type_names)/sizeof(mem_type_names[0]))

/* for freeing various types */
/*static int (*mem_free_funcs[MEM_NUM_STD_TYPES])() = {
   m_free,
   v_free};*/

/* local array for keeping track of memory */
static MEM_ARRAY   mem_info_sum[MEM_NUM_STD_TYPES];  


/* it is a global variable for passing 
   pointers to local arrays defined here */
MEM_CONNECT mem_connect[MEM_CONNECT_MAX_LISTS];
mem_connect[0] = mem_type_names;
mem_connect[1] = MEM_NUM_STD_TYPES;
mem_connect[2] = mem_info_sum;

/* mem_bytes_list
   
   Arguments:
   type - the number of type;
   old_size - old size of allocated memory (in bytes);
   new_size - new size of allocated memory (in bytes);
   list - list of types
   */
void mem_bytes_list(int type, int old_size, int new_size, int list)
{
   MEM_CONNECT *mlist;
   
   if ( list < 0 || list >= MEM_CONNECT_MAX_LISTS )
     return;
   
   mlist = &mem_connect[list];
   if (  type < 0 || type >= mlist->ntypes
       || mlist->free_funcs[type] == NULL )
     return;

   mlist->info_sum[type].bytes += new_size - old_size;
}


/* m_get -- gets an mxn matrix (in MAT form) by dynamic memory allocation
	-- normally ALL matrices should be obtained this way
	-- if either m or n is negative this will raise an error
	-- note that 0 x n and m x 0 matrices can be created */
/*MAT	*m_get(int m, int n)
{
   MAT	*matrix;
   int	i;

   mem_bytes(TYPE_MAT,0,sizeof(MAT));
   mem_numvar(TYPE_MAT,1);
   
   matrix->m = m;		matrix->n = matrix->max_n = n;
   matrix->max_m = m;	matrix->max_size = m*n;
   
   return (matrix);
}*/

/* m_get -- gets an mxn matrix (in MAT form) by dynamic memory allocation
	-- normally ALL matrices should be obtained this way
	-- if either m or n is negative this will raise an error
	-- note that 0 x n and m x 0 matrices can be created */
#ifndef ANSI_C
MAT	*m_get(m,n)
int	m,n;
#else
MAT	*m_get(int m, int n)
#endif
{
   MAT	*matrix;
   int	i;
      mem_bytes(TYPE_MAT,0,sizeof(MAT));
      mem_numvar(TYPE_MAT,1);
   
   matrix->m = m;		matrix->n = matrix->max_n = n;
   matrix->max_m = m;	matrix->max_size = m*n;
#ifndef SEGMENTED
   if ((matrix->base = NEW_A(m*n,double)) == (double *)NULL )
   {
      free(matrix);
   }
   else {
      mem_bytes(TYPE_MAT,0,m*n*sizeof(double));
   }
#else
   matrix->base = (double *)NULL;
#endif
    mem_bytes(TYPE_MAT,0,m*sizeof(double *));
   
#ifndef SEGMENTED
   /* set up pointers */
   for ( i=0; i<m; i++ )
     matrix->me[i] = &(matrix->base[i*n]);
#else
   for ( i = 0; i < m; i++ )
	mem_bytes(TYPE_MAT,0,n*sizeof(double));
#endif
   
   return (matrix);
}

/* m_resize -- returns the matrix A of size new_m x new_n; A is zeroed
   -- if A == NULL on entry then the effect is equivalent to m_get() */
/*MAT	*m_resize(MAT *A,int new_m, int new_n)
{
   int	i;
   int	new_max_m, new_max_n, old_m, old_n;
   
   if (new_m < 0 || new_n < 0)
     printf("The size must be positive!");

   if ( ! A )
     return m_get(new_m,new_n);

   // nothing was changed 
   if (new_m == A->m && new_n == A->n)
     return A;

   old_m = A->m;	old_n = A->n;
   if ( new_m > A->max_m )
   {	// re-allocate A->me
	 mem_bytes(TYPE_MAT,A->max_m*sizeof(double *),
		      new_m*sizeof(double *));
      A->me = RENEW(A->me,new_m,double *);
   }
   new_max_m = max(new_m,A->max_m);
   new_max_n = max(new_n,A->max_n);
      
   A->max_m = new_max_m;
   A->max_n = new_max_n;
   A->max_size = A->max_m*A->max_n;
   A->m = new_m;	A->n = new_n;
   
   return A;
}*/

/* m_resize -- returns the matrix A of size new_m x new_n; A is zeroed
   -- if A == NULL on entry then the effect is equivalent to m_get() */
#ifndef ANSI_C
MAT	*m_resize(A,new_m,new_n)
MAT	*A;
int	new_m, new_n;
#else
MAT	*m_resize(MAT *A,int new_m, int new_n)
#endif
{
   int	i;
   int	new_max_m, new_max_n, new_size, old_m, old_n;

   if ( ! A )
     return m_get(new_m,new_n);

   /* nothing was changed */
   if (new_m == A->m && new_n == A->n)
     return A;

   old_m = A->m;	old_n = A->n;
   if ( new_m > A->max_m )
   {	/* re-allocate A->me */
	 mem_bytes(TYPE_MAT,A->max_m*sizeof(double *),
		      new_m*sizeof(double_t *));

      A->me = RENEW(A->me,new_m,double *);
   }
   new_max_m = max(new_m,A->max_m);
   new_max_n = max(new_n,A->max_n);
   
#ifndef SEGMENTED
   new_size = new_max_m*new_max_n;
   if ( new_size > A->max_size )
   {	/* re-allocate A->base */
	 mem_bytes(TYPE_MAT,A->max_m*A->max_n*sizeof(double),
		      new_size*sizeof(double));

      A->base = RENEW(A->base,new_size,double);
      A->max_size = new_size;
   }
   
   /* now set up A->me[i] */
   for ( i = 0; i < new_m; i++ )
     A->me[i] = &(A->base[i*new_n]);
   
   /* now shift data in matrix */
   if ( old_n > new_n )
   {
      for ( i = 1; i < min(old_m,new_m); i++ )
	MEM_COPY((char *)&(A->base[i*old_n]),
		 (char *)&(A->base[i*new_n]),
		 sizeof(double)*new_n);
   }
   else if ( old_n < new_n )
   {
      for ( i = (int)(min(old_m,new_m))-1; i > 0; i-- )
      {   /* copy & then zero extra space */
	 MEM_COPY((char *)&(A->base[i*old_n]),
		  (char *)&(A->base[i*new_n]),
		  sizeof(double)*old_n);
	 __zero__(&(A->base[i*new_n+old_n]),(new_n-old_n));
      }
      __zero__(&(A->base[old_n]),(new_n-old_n));
      A->max_n = new_n;
   }
   /* zero out the new rows.. */
   for ( i = old_m; i < new_m; i++ )
     __zero__(&(A->base[i*new_n]),new_n);
#else
   if ( A->max_n < new_n )
   {
      double	*tmp;
      
      for ( i = 0; i < A->max_m; i++ )
      {
	    mem_bytes(TYPE_MAT,A->max_n*sizeof(double),
			 new_max_n*sizeof(double));

	    A->me[i] = tmp;
      }
      for ( i = A->max_m; i < new_max_m; i++ )
      {
	    A->me[i] = tmp;
	    mem_bytes(TYPE_MAT,0,new_max_n*sizeof(double));
      }
   }
   else if ( A->max_m < new_m )
   {
      for ( i = A->max_m; i < new_m; i++ )
	   mem_bytes(TYPE_MAT,0,new_max_n*sizeof(double));
      
   }
   
   if ( old_n < new_n )
   {
      for ( i = 0; i < old_m; i++ )
	__zero__(&(A->me[i][old_n]),new_n-old_n);
   }
   
   /* zero out the new rows.. */
   for ( i = old_m; i < new_m; i++ )
     __zero__(A->me[i],new_n);
#endif
   
   A->max_m = new_max_m;
   A->max_n = new_max_n;
   A->max_size = A->max_m*A->max_n;
   A->m = new_m;	A->n = new_n;
   
   return A;
}

/* m_copy -- copies matrix into new area
	-- out(i0:m,j0:n) <- in(i0:m,j0:n) */
/*MAT	*m_copy(const MAT *in, MAT *out, unsigned int i0, unsigned int j0)
{
	unsigned int	i, j;

	if ( in==MNULL )
		printf("Error in is MNULL");
	if ( in==out )
		return (out);
	if ( out==MNULL || out->m < in->m || out->n < in->n )
		out = m_resize(out,in->m,in->n);

	for ( i=i0; i < in->m; i++ )
		//MEM_COPY(&(in->me[i][j0]),&(out->me[i][j0]),
		//		(in->n - j0)*sizeof(double));
		for ( j=j0; j < in->n; j++ )
			out->me[i][j] = in->me[i][j];

	return (out);
}*/

/* _m_copy -- copies matrix into new area
	-- out(i0:m,j0:n) <- in(i0:m,j0:n) */
#ifndef ANSI_C
MAT	*_m_copy(in,out,i0,j0)
MAT	*in,*out;
unsigned int	i0,j0;
#else
MAT	*_m_copy(const MAT *in, MAT *out, unsigned int i0, unsigned int j0)
#endif
{
	unsigned int	i /* ,j */;

	if ( in==out )
		return (out);
	if ( out==MNULL || out->m < in->m || out->n < in->n )
		out = m_resize(out,in->m,in->n);

	for ( i=i0; i < in->m; i++ )
		MEM_COPY(&(in->me[i][j0]),&(out->me[i][j0]),
				(in->n - j0)*sizeof(double));
		/* for ( j=j0; j < in->n; j++ )
			out->me[i][j] = in->me[i][j]; */

	return (out);
}

/* v_resize -- returns the vector x with dim new_dim
   -- x is set to the zero vector */
/*#ifndef ANSI_C
VEC	*v_resize(x,new_dim)
VEC	*x;
int	new_dim;
#else
VEC	*v_resize(VEC *x, int new_dim)
#endif
{
   if ( ! x )
     return v_get(new_dim);

   // nothing is changed
   if (new_dim == x->dim)
     return x;

   if ( x->max_dim == 0 )	// assume that it's from sub_vec
     return v_get(new_dim);
   
   if ( new_dim > x->max_dim )
   {
	 mem_bytes(TYPE_VEC,x->max_dim*sizeof(double),
			 new_dim*sizeof(double));

      x->ve = RENEW(x->ve,new_dim,double);
      x->max_dim = new_dim;
   }
   
   if ( new_dim > x->dim )
     __zero__(&(x->ve[x->dim]),new_dim - x->dim);
   x->dim = new_dim;
   
   return x;
}*/

/* v_get -- gets a VEC of dimension 'size'
   -- Note: initialized to zero */
#ifndef ANSI_C
VEC	*v_get(size)
int	size;
#else
VEC	*v_get(int size)
#endif
{
   VEC	*vector;
      mem_bytes(TYPE_VEC,0,sizeof(VEC));
      mem_numvar(TYPE_VEC,1);
   
   vector->dim = vector->max_dim = size;
   if ((vector->ve=NEW_A(size,double)) == (double *)NULL )
   {
      free(vector);
   }
   else {
      mem_bytes(TYPE_VEC,0,size*sizeof(double));
   }
   
   return (vector);
}

/* v_resize -- returns the vector x with dim new_dim
   -- x is set to the zero vector */
#ifndef ANSI_C
VEC	*v_resize(x,new_dim)
VEC	*x;
int	new_dim;
#else
VEC	*v_resize(VEC *x, int new_dim)
#endif
{
   if ( ! x )
     return v_get(new_dim);

   /* nothing is changed */
   if (new_dim == x->dim)
     return x;

   if ( x->max_dim == 0 )	/* assume that it's from sub_vec */
     return v_get(new_dim);
   
   if ( new_dim > x->max_dim )
   {
	 mem_bytes(TYPE_VEC,x->max_dim*sizeof(double),
			 new_dim*sizeof(double));

      x->ve = RENEW(x->ve,new_dim,double);
      x->max_dim = new_dim;
   }
   
   if ( new_dim > x->dim )
     __zero__(&(x->ve[x->dim]),new_dim - x->dim);
   x->dim = new_dim;
   
   return x;
}

/* get_col -- gets a specified column of a matrix and retruns it as a vector */
#ifndef ANSI_C
VEC	*get_col(mat,col,vec)
unsigned int	col;
MAT	*mat;
VEC	*vec;
#else
VEC	*get_col(const MAT *mat, unsigned int col, VEC *vec)
#endif
{
   unsigned int	i;
   
   if ( vec==(VEC *)NULL || vec->dim<mat->m )
     vec = v_resize(vec,mat->m);
   
   for ( i=0; i<mat->m; i++ )
     vec->ve[i] = mat->me[i][col];
   
   return (vec);
}

/* get_row -- gets a specified row of a matrix and retruns it as a vector */
#ifndef ANSI_C
VEC	*get_row(mat,row,vec)
unsigned int	row;
MAT	*mat;
VEC	*vec;
#else
VEC	*get_row(const MAT *mat, unsigned int row, VEC *vec)
#endif
{
   unsigned int	i;
   
   if ( vec==(VEC *)NULL || vec->dim<mat->n )
     vec = v_resize(vec,mat->n);
   
   for ( i=0; i<mat->n; i++ )
     vec->ve[i] = mat->me[row][i];
   
   return (vec);
}

/* _set_col -- sets column of matrix to values given in vec (in situ)
	-- that is, mat(i0:lim,col) <- vec(i0:lim) */
#ifndef ANSI_C
MAT	*_set_col(mat,col,vec,i0)
MAT	*mat;
VEC	*vec;
unsigned int	col,i0;
#else
MAT	*_set_col(MAT *mat, unsigned int col, const VEC *vec, unsigned int i0)
#endif
{
   unsigned int	i,lim;
   
   lim = min(mat->m,vec->dim);
   for ( i=i0; i<lim; i++ )
     mat->me[i][col] = vec->ve[i];
   
   return (mat);
}

/* _set_row -- sets row of matrix to values given in vec (in situ) */
#ifndef ANSI_C
MAT	*_set_row(mat,row,vec,j0)
MAT	*mat;
VEC	*vec;
unsigned int	row,j0;
#else
MAT	*_set_row(MAT *mat, unsigned int row, const VEC *vec, unsigned int j0)
#endif
{
   unsigned int	j,lim;

   lim = min(mat->n,vec->dim);
   for ( j=j0; j<lim; j++ )
     mat->me[row][j] = vec->ve[j];
   
   return (mat);
}

/* _v_copy -- copies vector into new area
	-- out(i0:dim) <- in(i0:dim) */
#ifndef ANSI_C
VEC	*_v_copy(in,out,i0)
VEC	*in,*out;
unsigned int	i0;
#else
VEC	*_v_copy(const VEC *in, VEC *out, unsigned int i0)
#endif
{
	/* unsigned int	i,j; */

	if ( in==out )
		return (out);
	if ( out==VNULL || out->dim < in->dim )
		out = v_resize(out,in->dim);

	MEM_COPY(&(in->ve[i0]),&(out->ve[i0]),(in->dim - i0)*sizeof(double));
	/* for ( i=i0; i < in->dim; i++ )
		out->ve[i] = in->ve[i]; */

	return (out);
}

/* __ip__ -- inner product */
#ifndef ANSI_C
double	__ip__(dp1,dp2,len)
register double	*dp1, *dp2;
int	len;
#else
double	__ip__(const double *dp1, const double *dp2, int len)
#endif
{
#ifdef VUNROLL
    register int	len4;
    register double	sum1, sum2, sum3;
#endif
    register int	i;
    register double     sum;

    sum = 0.0;
#ifdef VUNROLL
    sum1 = sum2 = sum3 = 0.0;
    
    len4 = len / 4;
    len  = len % 4;
    
    for ( i = 0; i < len4; i++ )
    {
	sum  += dp1[4*i]*dp2[4*i];
	sum1 += dp1[4*i+1]*dp2[4*i+1];
	sum2 += dp1[4*i+2]*dp2[4*i+2];
	sum3 += dp1[4*i+3]*dp2[4*i+3];
    }
    sum  += sum1 + sum2 + sum3;
    dp1 += 4*len4;	dp2 += 4*len4;
#endif
    
    for ( i = 0; i < len; i++ )
	sum  += dp1[i]*dp2[i];
    
    return sum;
}

/* _in_prod -- inner product of two vectors from i0 downwards
	-- that is, returns a(i0:dim)^T.b(i0:dim) */
#ifndef ANSI_C
double	_in_prod(a,b,i0)
VEC	*a,*b;
unsigned int	i0;
#else
double	_in_prod(const VEC *a, const VEC *b, unsigned int i0)
#endif
{
	unsigned int	limit;
	/* double	*a_v, *b_v; */
	/* register double	sum; */

	return __ip__(&(a->ve[i0]),&(b->ve[i0]),(int)(limit-i0));
	/*****************************************
	a_v = &(a->ve[i0]);		b_v = &(b->ve[i0]);
	for ( i=i0; i<limit; i++ )
		sum += a_v[i]*b_v[i];
		sum += (*a_v++)*(*b_v++);

	return (double)sum;
	******************************************/
}

/* hhvec -- calulates Householder vector to eliminate all entries after the
	i0 entry of the vector vec. It is returned as out. May be in-situ */
#ifndef ANSI_C
VEC	*hhvec(vec,i0,beta,out,newval)
VEC	*vec,*out;
unsigned int	i0;
double	*beta,*newval;
#else
VEC	*hhvec(const VEC *vec, unsigned int i0, double *beta,
	       VEC *out, double *newval)
#endif
{
	double	norm,temp;

	out = _v_copy(vec,out,i0);
	temp = (double)_in_prod(out,out,i0);
	norm = sqrt(temp);
	if ( norm <= 0.0 )
	{
		*beta = 0.0;
		return (out);
	}
	*beta = 1.0/(norm * (norm+fabs(out->ve[i0])));
	if ( out->ve[i0] > 0.0 )
		*newval = -norm;
	else
		*newval = norm;
	out->ve[i0] -= *newval;

	return (out);
}

/* v_zero -- zero the vector x */
#ifndef ANSI_C
VEC	*v_zero(x)
VEC	*x;
#else
VEC	*v_zero(VEC *x)
#endif
{

	__zero__(x->ve,x->dim);
	/* for ( i = 0; i < x->dim; i++ )
		x->ve[i] = 0.0; */

	return x;
}

/* __mltadd__ -- scalar multiply and add c.f. v_mltadd() */
#ifndef ANSI_C
void	__mltadd__(dp1,dp2,s,len)
register double	*dp1, *dp2;
register double s;
register int	len;
#else
void	__mltadd__(double *dp1, const double *dp2, double s, int len)
#endif
{
    register int	i;
#ifdef VUNROLL
    register int        len4;
    
    len4 = len / 4;
    len  = len % 4;
    for ( i = 0; i < len4; i++ )
    {
	dp1[4*i]   += s*dp2[4*i];
	dp1[4*i+1] += s*dp2[4*i+1];
	dp1[4*i+2] += s*dp2[4*i+2];
	dp1[4*i+3] += s*dp2[4*i+3];
    }
    dp1 += 4*len4;	dp2 += 4*len4;
#endif
    
    for ( i = 0; i < len; i++ )
	dp1[i] += s*dp2[i];
}

/* _hhtrcols -- transform a matrix by a Householder vector by columns
	starting at row i0 from column j0 
	-- that is, M(i0:m,j0:n) <- (I-beta.hh(i0:m).hh(i0:m)^T)M(i0:m,j0:n)
	-- in-situ
	-- scratch vector w passed as argument
	-- raises error if w == NULL
*/
#ifndef ANSI_C
MAT	*_hhtrcols(M,i0,j0,hh,beta,w)
MAT	*M;
unsigned int	i0, j0;
VEC	*hh;
double	beta;
VEC	*w;
#else
MAT	*_hhtrcols(MAT *M, unsigned int i0, unsigned int j0,
		   const VEC *hh, double beta, VEC *w)
#endif
{
	/* double	ip, scale; */
	int	i /*, k */;
	/*  static	VEC	*w = VNULL; */

	if ( beta == 0.0 )	return (M);

	if ( w->dim < M->n )
	  w = v_resize(w,M->n);
	/*  MEM_STAT_REG(w,TYPE_VEC); */
	v_zero(w);

	for ( i = i0; i < M->m; i++ )
	    if ( hh->ve[i] != 0.0 )
		__mltadd__(&(w->ve[j0]),&(M->me[i][j0]),hh->ve[i],
							(int)(M->n-j0));
	for ( i = i0; i < M->m; i++ )
	    if ( hh->ve[i] != 0.0 )
		__mltadd__(&(M->me[i][j0]),&(w->ve[j0]),-beta*hh->ve[i],
							(int)(M->n-j0));
	return (M);
}

/* __ip__ -- inner product */
/*#ifndef ANSI_C
double	__ip__(dp1,dp2,len)
register double	*dp1, *dp2;
int	len;
#else
double	__ip__(const double *dp1, const double *dp2, int len)
#endif
{
#ifdef VUNROLL
    register int	len4;
    register double	sum1, sum2, sum3;
#endif
    register int	i;
    register double     sum;

    sum = 0.0;
#ifdef VUNROLL
    sum1 = sum2 = sum3 = 0.0;
    
    len4 = len / 4;
    len  = len % 4;
    
    for ( i = 0; i < len4; i++ )
    {
	sum  += dp1[4*i]*dp2[4*i];
	sum1 += dp1[4*i+1]*dp2[4*i+1];
	sum2 += dp1[4*i+2]*dp2[4*i+2];
	sum3 += dp1[4*i+3]*dp2[4*i+3];
    }
    sum  += sum1 + sum2 + sum3;
    dp1 += 4*len4;	dp2 += 4*len4;
#endif
    
    for ( i = 0; i < len; i++ )
	sum  += dp1[i]*dp2[i];
    
    return sum;
}*/

/* hhtrrows -- transform a matrix by a Householder vector by rows
	starting at row i0 from column j0 -- in-situ
	-- that is, M(i0:m,j0:n) <- M(i0:m,j0:n)(I-beta.hh(j0:n).hh(j0:n)^T) */
#ifndef ANSI_C
MAT	*hhtrrows(M,i0,j0,hh,beta)
MAT	*M;
unsigned int	i0, j0;
VEC	*hh;
double	beta;
#else
MAT	*hhtrrows(MAT *M, unsigned int i0, unsigned int j0,
		  const VEC *hh, double beta)
#endif
{
	double	ip, scale;
	int	i /*, j */;

	if ( beta == 0.0 )	return (M);

	/* for each row ... */
	for ( i = i0; i < M->m; i++ )
	{	/* compute inner product */
		ip = __ip__(&(M->me[i][j0]),&(hh->ve[j0]),(int)(M->n-j0));
		/**************************************************
		ip = 0.0;
		for ( j = j0; j < M->n; j++ )
			ip += M->me[i][j]*hh->ve[j];
		**************************************************/
		scale = beta*ip;
		if ( scale == 0.0 )
		    continue;

		/* do operation */
		__mltadd__(&(M->me[i][j0]),&(hh->ve[j0]),-scale,
							(int)(M->n-j0));
		/**************************************************
		for ( j = j0; j < M->n; j++ )
			M->me[i][j] -= scale*hh->ve[j];
		**************************************************/
	}

	return (M);
}

/* Hfactor -- compute Hessenberg factorisation in compact form.
	-- factorisation performed in situ
	-- for details of the compact form see QRfactor.c and matrix2.doc */
#ifndef ANSI_C
MAT	*Hfactor(A, diag, beta)
MAT	*A;
VEC	*diag, *beta;
#else
MAT	*Hfactor(MAT *A, VEC *diag, VEC *beta)
#endif
{
	static	VEC	*hh = VNULL, *w = VNULL;
	int	k, limit;
	double b;

	limit = A->m - 1;

	hh = v_resize(hh,A->m);
	w  = v_resize(w,A->n);
	MEM_STAT_REG(hh,TYPE_VEC);
	MEM_STAT_REG(w, TYPE_VEC);

	for ( k = 0; k < limit; k++ )
	  {
	    /* compute the Householder vector hh */
	    get_col(A,(unsigned int)k,hh);
	    /* printf("the %d'th column = ");	v_output(hh); */
	    hhvec(hh,k+1,&beta->ve[k],hh,&A->me[k+1][k]);
	    /* diag->ve[k] = hh->ve[k+1]; */
	    v_set_val(diag,k,v_entry(hh,k+1));
	    /* printf("H/h vector = ");	v_output(hh); */
	    /* printf("from the %d'th entry\n",k+1); */
	    /* printf("beta = %g\n",beta->ve[k]); */

	    /* apply Householder operation symmetrically to A */
		b = v_entry(beta,k);
	    _hhtrcols(A,k+1,k+1,hh,b,w);
	    hhtrrows(A,0  ,k+1,hh,b);
	    /* printf("A = ");		m_output(A); */
	  }

#ifdef THREADSAFE
	V_FREE(hh);	V_FREE(w);
#endif

	return (A);
}

/* hhtrvec -- apply Householder transformation to vector 
	-- that is, out <- (I-beta.hh(i0:n).hh(i0:n)^T).in
	-- may be in-situ */
#ifndef ANSI_C
VEC	*hhtrvec(hh,beta,i0,in,out)
VEC	*hh,*in,*out;	/* hh = Householder vector */
unsigned int	i0;
double	beta;
#else
VEC	*hhtrvec(const VEC *hh, double beta, unsigned int i0,
		 const VEC *in, VEC *out)
#endif
{
	double	scale,temp;
	/* unsigned int	i; */

	temp = (double)_in_prod(hh,in,i0);
	scale = beta*temp;
	out = v_copy(in,out);
	__mltadd__(&(out->ve[i0]),&(hh->ve[i0]),-scale,(int)(in->dim-i0));
	/************************************************************
	for ( i=i0; i<in->dim; i++ )
		out->ve[i] = in->ve[i] - scale*hh->ve[i];
	************************************************************/

	return (out);
}

/* makeHQ -- construct the Hessenberg orthogonalising matrix Q;
	-- i.e. Hess M = Q.M.Q'	*/
#ifndef ANSI_C
MAT	*makeHQ(H, diag, beta, Qout)
MAT	*H, *Qout;
VEC	*diag, *beta;
#else
MAT	*makeHQ(MAT *H, VEC *diag, VEC *beta, MAT *Qout)
#endif
{
	int	i, j, limit;
	static	VEC	*tmp1 = VNULL, *tmp2 = VNULL;

	Qout = m_resize(Qout,H->m,H->m);

	tmp1 = v_resize(tmp1,H->m);
	tmp2 = v_resize(tmp2,H->m);
	MEM_STAT_REG(tmp1,TYPE_VEC);
	MEM_STAT_REG(tmp2,TYPE_VEC);

	for ( i = 0; i < H->m; i++ )
	{
		/* tmp1 = i'th basis vector */
		for ( j = 0; j < H->m; j++ )
			/* tmp1->ve[j] = 0.0; */
		    v_set_val(tmp1,j,0.0);
		/* tmp1->ve[i] = 1.0; */
		v_set_val(tmp1,i,1.0);

		/* apply H/h transforms in reverse order */
		for ( j = limit-1; j >= 0; j-- )
		{
			get_col(H,(unsigned int)j,tmp2);
			/* tmp2->ve[j+1] = diag->ve[j]; */
			v_set_val(tmp2,j+1,v_entry(diag,j));
			hhtrvec(tmp2,beta->ve[j],j+1,tmp1,tmp1);
		}

		/* insert into Qout */
		set_col(Qout,(unsigned int)i,tmp1);
	}

#ifdef THREADSAFE
	V_FREE(tmp1);	V_FREE(tmp2);
#endif

	return (Qout);
}

/* makeH -- construct actual Hessenberg matrix */
#ifndef ANSI_C
MAT	*makeH(H,Hout)
MAT	*H, *Hout;
#else
MAT	*makeH(const MAT *H, MAT *Hout)
#endif
{
	int	i, j, limit;

	Hout = m_resize(Hout,H->m,H->m);
	Hout = m_copy(H,Hout);

	limit = H->m;
	for ( i = 1; i < limit; i++ )
		for ( j = 0; j < i-1; j++ )
			/* Hout->me[i][j] = 0.0;*/
		    m_set_val(Hout,i,j,0.0);

	return (Hout);
}

/* rot_rows -- premultiply mat by givens rotation described by c,s */
#ifndef ANSI_C
MAT	*rot_rows(mat,i,k,c,s,out)
MAT	*mat,*out;
unsigned int	i,k;
double	c,s;
#else
MAT	*rot_rows(const MAT *mat, unsigned int i, unsigned int k,
		  double c, double s, MAT *out)
#endif
{
	unsigned int	j;
	double	temp;

	if ( mat != out )
		out = m_copy(mat,m_resize(out,mat->m,mat->n));

	for ( j=0; j<mat->n; j++ )
	{
		/* temp = c*out->me[i][j] + s*out->me[k][j]; */
		temp = c*m_entry(out,i,j) + s*m_entry(out,k,j);
		/* out->me[k][j] = -s*out->me[i][j] + c*out->me[k][j]; */
		m_set_val(out,k,j, -s*m_entry(out,i,j) + c*m_entry(out,k,j));
		/* out->me[i][j] = temp; */
		m_set_val(out,i,j, temp);
	}

	return (out);
}

/* rot_cols -- postmultiply mat by givens rotation described by c,s */
#ifndef ANSI_C
MAT	*rot_cols(mat,i,k,c,s,out)
MAT	*mat,*out;
unsigned int	i,k;
double	c,s;
#else
MAT	*rot_cols(const MAT *mat,unsigned int i,unsigned int k,
		  double c, double s, MAT *out)
#endif
{
	unsigned int	j;
	double	temp;

	if ( mat != out )
		out = m_copy(mat,m_resize(out,mat->m,mat->n));

	for ( j=0; j<mat->m; j++ )
	{
		/* temp = c*out->me[j][i] + s*out->me[j][k]; */
		temp = c*m_entry(out,j,i) + s*m_entry(out,j,k);
		/* out->me[j][k] = -s*out->me[j][i] + c*out->me[j][k]; */
		m_set_val(out,j,k, -s*m_entry(out,j,i) + c*m_entry(out,j,k));
		/* out->me[j][i] = temp; */
		m_set_val(out,j,i,temp);
	}

	return (out);
}

/* hhldr3 -- computes */
#ifndef ANSI_C
static	void	hhldr3(x,y,z,nu1,beta,newval)
double	x, y, z;
double	*nu1, *beta, *newval;
#else
static	void	hhldr3(double x, double y, double z,
		       double *nu1, double *beta, double *newval)
#endif
{
	double	alpha;

	if ( x >= 0.0 )
		alpha = sqrt(x*x+y*y+z*z);
	else
		alpha = -sqrt(x*x+y*y+z*z);
	*nu1 = x + alpha;
	*beta = 1.0/(alpha*(*nu1));
	*newval = alpha;
}


/*hhldr3rows */
#ifndef ANSI_C
static	void	hhldr3rows(A,k,i0,beta,nu1,nu2,nu3)
MAT	*A;
int	k, i0;
double	beta, nu1, nu2, nu3;
#else
static	void	hhldr3rows(MAT *A, int k, int i0, double beta,
			   double nu1, double nu2, double nu3)
#endif
{
	double	**A_me, ip, prod;
	int	i, m;

	/* printf("hhldr3rows:(l.%d) A at 0x%lx\n", __LINE__, (long)A); */
	/* printf("hhldr3rows: k = %d\n", k); */
	A_me = A->me;		m = A->m;
	i0 = min(i0,m-1);

	for ( i = 0; i <= i0; i++ )
	{
	    /****
	    ip = nu1*A_me[i][k] + nu2*A_me[i][k+1] + nu3*A_me[i][k+2];
	    prod = ip*beta;
	    A_me[i][k]   -= prod*nu1;
	    A_me[i][k+1] -= prod*nu2;
	    A_me[i][k+2] -= prod*nu3;
	    ****/

	    ip = nu1*m_entry(A,i,k)+nu2*m_entry(A,i,k+1)+nu3*m_entry(A,i,k+2);
	    prod = ip*beta;
	    m_add_val(A,i,k  , - prod*nu1);
	    m_add_val(A,i,k+1, - prod*nu2);
	    m_add_val(A,i,k+2, - prod*nu3);

	}
}

/* givens -- returns c,s parameters for Givens rotation to
		eliminate y in the vector [ x y ]' */
#ifndef ANSI_C
void	givens(x,y,c,s)
double  x,y;
double	*c,*s;
#else
void	givens(double x, double y, double *c, double *s)
#endif
{
	double	norm;

	norm = sqrt(x*x+y*y);
	if ( norm == 0.0 )
	{	*c = 1.0;	*s = 0.0;	}	/* identity */
	else
	{	*c = x/norm;	*s = y/norm;	}
}


/* schur -- computes the Schur decomposition of the matrix A in situ
	-- optionally, gives Q matrix such that Q^T.A.Q is upper triangular
	-- returns upper triangular Schur matrix */
#ifndef ANSI_C
MAT	*schur(A,Q)
MAT	*A, *Q;
#else
MAT	*schur(MAT *A, MAT *Q)
#endif
{
    int		i, j, iter, k, k_min, k_max, k_tmp, n, split;
    double	beta2, c, discrim, dummy, nu1, s, t, tmp, x, y, z;
    double	**A_me;
    double	sqrt_macheps;
    static	VEC	*diag=VNULL, *beta=VNULL;

    n = A->n;
    diag = v_resize(diag,A->n);
    beta = v_resize(beta,A->n);
    MEM_STAT_REG(diag,TYPE_VEC);
    MEM_STAT_REG(beta,TYPE_VEC);
    /* compute Hessenberg form */
    Hfactor(A,diag,beta);
    
    /* save Q if necessary */
    if ( Q )
	Q = makeHQ(A,diag,beta,Q);
    makeH(A,A);

    sqrt_macheps = sqrt(MACHEPS);

    k_min = 0;	A_me = A->me;

    while ( k_min < n )
    {
	double	a00, a01, a10, a11;
	double	scale, t, numer, denom;

	/* find k_max to suit:
	   submatrix k_min..k_max should be irreducible */
	k_max = n-1;
	for ( k = k_min; k < k_max; k++ )
	    /* if ( A_me[k+1][k] == 0.0 ) */
	    if ( m_entry(A,k+1,k) == 0.0 )
	    {	k_max = k;	break;	}

	if ( k_max <= k_min )
	{
	    k_min = k_max + 1;
	    continue;		/* outer loop */
	}

	/* check to see if we have a 2 x 2 block
	   with complex eigenvalues */
	if ( k_max == k_min + 1 )
	{
	    /* tmp = A_me[k_min][k_min] - A_me[k_max][k_max]; */
	    a00 = m_entry(A,k_min,k_min);
	    a01 = m_entry(A,k_min,k_max);
	    a10 = m_entry(A,k_max,k_min);
	    a11 = m_entry(A,k_max,k_max);
	    tmp = a00 - a11;
	    /* discrim = tmp*tmp +
		4*A_me[k_min][k_max]*A_me[k_max][k_min]; */
	    discrim = tmp*tmp + 4*a01*a10;
	    if ( discrim < 0.0 )
	    {	/* yes -- e-vals are complex
		   -- put 2 x 2 block in form [a b; c a];
		   then eigenvalues have real part a & imag part sqrt(|bc|) */
		numer = - tmp;
		denom = ( a01+a10 >= 0.0 ) ?
		    (a01+a10) + sqrt((a01+a10)*(a01+a10)+tmp*tmp) :
		    (a01+a10) - sqrt((a01+a10)*(a01+a10)+tmp*tmp);
		if ( denom != 0.0 )
		{   /* t = s/c = numer/denom */
		    t = numer/denom;
		    scale = c = 1.0/sqrt(1+t*t);
		    s = c*t;
		}
		else
		{
		    c = 1.0;
		    s = 0.0;
		}
		rot_cols(A,k_min,k_max,c,s,A);
		rot_rows(A,k_min,k_max,c,s,A);
		if ( Q != MNULL )
		    rot_cols(Q,k_min,k_max,c,s,Q);
		k_min = k_max + 1;
		continue;
	    }
	    else /* discrim >= 0; i.e. block has two real eigenvalues */
	    {	/* no -- e-vals are not complex;
		   split 2 x 2 block and continue */
		/* s/c = numer/denom */
		numer = ( tmp >= 0.0 ) ?
		    - tmp - sqrt(discrim) : - tmp + sqrt(discrim);
		denom = 2*a01;
		if ( fabs(numer) < fabs(denom) )
		{   /* t = s/c = numer/denom */
		    t = numer/denom;
		    scale = c = 1.0/sqrt(1+t*t);
		    s = c*t;
		}
		else if ( numer != 0.0 )
		{   /* t = c/s = denom/numer */
		    t = denom/numer;
		    scale = 1.0/sqrt(1+t*t);
		    c = fabs(t)*scale;
		    s = ( t >= 0.0 ) ? scale : -scale;
		}
		else /* numer == denom == 0 */
		{
		    c = 0.0;
		    s = 1.0;
		}
		rot_cols(A,k_min,k_max,c,s,A);
		rot_rows(A,k_min,k_max,c,s,A);
		/* A->me[k_max][k_min] = 0.0; */
		if ( Q != MNULL )
		    rot_cols(Q,k_min,k_max,c,s,Q);
		k_min = k_max + 1;	/* go to next block */
		continue;
	    }
	}

	/* now have r x r block with r >= 2:
	   apply Francis QR step until block splits */
	split = FALSE;		iter = 0;
	while ( ! split )
	{
	    iter++;
	    
	    /* set up Wilkinson/Francis complex shift */
	    k_tmp = k_max - 1;

	    a00 = m_entry(A,k_tmp,k_tmp);
	    a01 = m_entry(A,k_tmp,k_max);
	    a10 = m_entry(A,k_max,k_tmp);
	    a11 = m_entry(A,k_max,k_max);

	    /* treat degenerate cases differently
	       -- if there are still no splits after five iterations
	          and the bottom 2 x 2 looks degenerate, force it to
		  split */
#ifdef DEBUG
	    printf("# schur: bottom 2 x 2 = [%lg, %lg; %lg, %lg]\n",
		   a00, a01, a10, a11);
#endif
	    if ( iter >= 5 &&
		 fabs(a00-a11) < sqrt_macheps*(fabs(a00)+fabs(a11)) &&
		 (fabs(a01) < sqrt_macheps*(fabs(a00)+fabs(a11)) ||
		  fabs(a10) < sqrt_macheps*(fabs(a00)+fabs(a11))) )
	    {
	      if ( fabs(a01) < sqrt_macheps*(fabs(a00)+fabs(a11)) )
		m_set_val(A,k_tmp,k_max,0.0);
	      if ( fabs(a10) < sqrt_macheps*(fabs(a00)+fabs(a11)) )
		{
		  m_set_val(A,k_max,k_tmp,0.0);
		  split = TRUE;
		  continue;
		}
	    }

	    s = a00 + a11;
	    t = a00*a11 - a01*a10;

	    /* break loop if a 2 x 2 complex block */
	    if ( k_max == k_min + 1 && s*s < 4.0*t )
	    {
		split = TRUE;
		continue;
	    }

	    /* perturb shift if convergence is slow */
	    if ( (iter % 10) == 0 )
	    {	s += iter*0.02;		t += iter*0.02;
	    }

	    /* set up Householder transformations */
	    k_tmp = k_min + 1;
	    /********************
	    x = A_me[k_min][k_min]*A_me[k_min][k_min] +
		A_me[k_min][k_tmp]*A_me[k_tmp][k_min] -
		    s*A_me[k_min][k_min] + t;
	    y = A_me[k_tmp][k_min]*
		(A_me[k_min][k_min]+A_me[k_tmp][k_tmp]-s);
	    if ( k_min + 2 <= k_max )
		z = A_me[k_tmp][k_min]*A_me[k_min+2][k_tmp];
	    else
		z = 0.0;
	    ********************/

	    a00 = m_entry(A,k_min,k_min);
	    a01 = m_entry(A,k_min,k_tmp);
	    a10 = m_entry(A,k_tmp,k_min);
	    a11 = m_entry(A,k_tmp,k_tmp);

	    /********************
	    a00 = A->me[k_min][k_min];
	    a01 = A->me[k_min][k_tmp];
	    a10 = A->me[k_tmp][k_min];
	    a11 = A->me[k_tmp][k_tmp];
	    ********************/
	    x = a00*a00 + a01*a10 - s*a00 + t;
	    y = a10*(a00+a11-s);
	    if ( k_min + 2 <= k_max )
		z = a10* /* m_entry(A,k_min+2,k_tmp) */ A->me[k_min+2][k_tmp];
	    else
		z = 0.0;

	    for ( k = k_min; k <= k_max-1; k++ )
	    {
		if ( k < k_max - 1 )
		{
		    hhldr3(x,y,z,&nu1,&beta2,&dummy);
		    if ( Q != MNULL )
			hhldr3rows(Q,k,n-1,beta2,nu1,y,z);
		}
		else
		{
		    givens(x,y,&c,&s);
		    rot_cols(A,k,k+1,c,s,A);
		    rot_rows(A,k,k+1,c,s,A);
		    if ( Q )
			rot_cols(Q,k,k+1,c,s,Q);
		}
		/* if ( k >= 2 )
		    m_set_val(A,k,k-2,0.0); */
		/* x = A_me[k+1][k]; */
		x = m_entry(A,k+1,k);
		if ( k <= k_max - 2 )
		    /* y = A_me[k+2][k];*/
		    y = m_entry(A,k+2,k);
		else
		    y = 0.0;
		if ( k <= k_max - 3 )
		    /* z = A_me[k+3][k]; */
		    z = m_entry(A,k+3,k);
		else
		    z = 0.0;
	    }
	    /* if ( k_min > 0 )
		m_set_val(A,k_min,k_min-1,0.0);
	    if ( k_max < n - 1 )
		m_set_val(A,k_max+1,k_max,0.0); */
	    for ( k = k_min; k <= k_max-2; k++ )
	    {
		/* zero appropriate sub-diagonals */
		m_set_val(A,k+2,k,0.0);
		if ( k < k_max-2 )
		    m_set_val(A,k+3,k,0.0);
	    }

	    /* test to see if matrix should split */
	    for ( k = k_min; k < k_max; k++ )
		if ( fabs(A_me[k+1][k]) < MACHEPS*
		    (fabs(A_me[k][k])+fabs(A_me[k+1][k+1])) )
		{	A_me[k+1][k] = 0.0;	split = TRUE;	}
	}
    }
    
    /* polish up A by zeroing strictly lower triangular elements
       and small sub-diagonal elements */
    for ( i = 0; i < A->m; i++ )
	for ( j = 0; j < i-1; j++ )
	    A_me[i][j] = 0.0;
    for ( i = 0; i < A->m - 1; i++ )
	if ( fabs(A_me[i+1][i]) < MACHEPS*
	    (fabs(A_me[i][i])+fabs(A_me[i+1][i+1])) )
	    A_me[i+1][i] = 0.0;

#ifdef	THREADSAFE
    V_FREE(diag);	V_FREE(beta);
#endif

    return A;
}

/* schur_vals -- compute real & imaginary parts of eigenvalues
	-- assumes T contains a block upper triangular matrix
		as produced by schur()
	-- real parts stored in real_pt, imaginary parts in imag_pt */
#ifndef ANSI_C
void	schur_evals(T,real_pt,imag_pt)
MAT	*T;
VEC	*real_pt, *imag_pt;
#else
void	schur_evals(MAT *T, VEC *real_pt, VEC *imag_pt)
#endif
{
	int	i, n;
	double	discrim, **T_me;
	double	diff, sum, tmp;

	n = T->n;	T_me = T->me;
	real_pt = v_resize(real_pt,(unsigned int)n);
	imag_pt = v_resize(imag_pt,(unsigned int)n);

	i = 0;
	while ( i < n )
	{
		if ( i < n-1 && T_me[i+1][i] != 0.0 )
		{   /* should be a complex eigenvalue */
		    sum  = 0.5*(T_me[i][i]+T_me[i+1][i+1]);
		    diff = 0.5*(T_me[i][i]-T_me[i+1][i+1]);
		    discrim = diff*diff + T_me[i][i+1]*T_me[i+1][i];
		    if ( discrim < 0.0 )
		    {	/* yes -- complex e-vals */
			real_pt->ve[i] = real_pt->ve[i+1] = sum;
			imag_pt->ve[i] = sqrt(-discrim);
			imag_pt->ve[i+1] = - imag_pt->ve[i];
		    }
		    else
		    {	/* no -- actually both real */
			tmp = sqrt(discrim);
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

complex double *getEigvalues(MAT *A){
    	MAT *T = MNULL, *Q = MNULL;
    	VEC *evals_re = VNULL, *evals_im = VNULL;

    	complex double *z;

    	Q = m_get(A->m,A->n);
    	T = m_copy(A,MNULL);

    	/* compute Schur form: A = Q.T.Q^T */
    	schur(T,Q);
    	/* extract eigenvalues */
    	evals_re = v_get(A->m);
    	evals_im = v_get(A->m);
    	schur_evals(T,evals_re,evals_im);

    	z=malloc(evals_re->dim*sizeof(complex double));
    	for(int i=0;i<evals_re->dim;i++){
    	  	z[i]=evals_re->ve[i]+I*evals_im->ve[i];
    	}
    	return z;
    }

    double maxMagEigVal(complex double *z, size_t size){
    	double maximum=0, aux;
    	for (int c = 1; c < size; c++)
    	  {
    		aux=cabs(z[c]);
    	    if (aux > maximum)
    	    {
    	       maximum  = aux;
    	    }
    	  }
    	return (double)maximum;
    }

    double log_b(double base, double x) {
        return (double)(log(x) / log(base));
    }

/* m_zero -- zero the matrix A */
#ifndef ANSI_C
MAT	*m_zero(A)
MAT	*A;
#else
MAT	*m_zero(MAT *A)
#endif
{
	int	i, A_m, A_n;
	double	**A_me;

	A_m = A->m;	A_n = A->n;	A_me = A->me;
	for ( i = 0; i < A_m; i++ )
		__zero__(A_me[i],A_n);
		/* for ( j = 0; j < A_n; j++ )
			A_me[i][j] = 0.0; */

	return A;
}

/* mat_id -- set A to being closest to identity matrix as possible
	-- i.e. A[i][j] == 1 if i == j and 0 otherwise */
#ifndef ANSI_C
MAT	*m_ident(A)
MAT	*A;
#else
MAT	*m_ident(MAT *A)
#endif
{
	int	i, size;

	m_zero(A);
	size = min(A->m,A->n);
	for ( i = 0; i < size; i++ )
		A->me[i][i] = 1.0;

	return A;
}

/* __sub__ -- subtract arrays c.f. v_sub() */
#ifndef ANSI_C
void	__sub__(dp1,dp2,out,len)
register double	*dp1, *dp2, *out;
register int	len;
#else
void	__sub__(const double *dp1, const double *dp2, double *out, int len)
#endif
{
    register int	i;
    for ( i = 0; i < len; i++ )
	out[i] = dp1[i] - dp2[i];
}

/* m_sub -- matrix subtraction -- may be in-situ */
#ifndef ANSI_C
MAT	*m_sub(mat1,mat2,out)
MAT	*mat1,*mat2,*out;
#else
MAT	*m_sub(const MAT *mat1, const MAT *mat2, MAT *out)
#endif
{
	unsigned int	m,n,i;

	if ( out==(MAT *)NULL || out->m != mat1->m || out->n != mat1->n )
		out = m_resize(out,mat1->m,mat1->n);
	m = mat1->m;	n = mat1->n;
	for ( i=0; i<m; i++ )
	{
		__sub__(mat1->me[i],mat2->me[i],out->me[i],(int)n);
		/**************************************************
		for ( j=0; j<n; j++ )
			out->me[i][j] = mat1->me[i][j]-mat2->me[i][j];
		**************************************************/
	}

	return (out);
}

/* __add__ -- add arrays c.f. v_add() */
#ifndef ANSI_C
void	__add__(dp1,dp2,out,len)
register double	*dp1, *dp2, *out;
register int	len;
#else
void	__add__(const double *dp1, const double *dp2, double *out, int len)
#endif
{
    register int	i;
    for ( i = 0; i < len; i++ )
	out[i] = dp1[i] + dp2[i];
}

/* m_add -- matrix addition -- may be in-situ */
#ifndef ANSI_C
MAT	*m_add(mat1,mat2,out)
MAT	*mat1,*mat2,*out;
#else
MAT	*m_add(const MAT *mat1, const MAT *mat2, MAT *out)
#endif
{
	unsigned int	m,n,i;

	if ( out==(MAT *)NULL || out->m != mat1->m || out->n != mat1->n )
		out = m_resize(out,mat1->m,mat1->n);
	m = mat1->m;	n = mat1->n;
	for ( i=0; i<m; i++ )
	{
		__add__(mat1->me[i],mat2->me[i],out->me[i],(int)n);
		/**************************************************
		for ( j=0; j<n; j++ )
			out->me[i][j] = mat1->me[i][j]+mat2->me[i][j];
		**************************************************/
	}

	return (out);
}

/* m_mlt -- matrix-matrix multiplication */
#ifndef ANSI_C
MAT	*m_mlt(A,B,OUT)
MAT	*A,*B,*OUT;
#else
MAT	*m_mlt(const MAT *A, const MAT *B, MAT *OUT)
#endif
{
	unsigned int	i, /* j, */ k, m, n, p;
	double	**A_v, **B_v /*, *B_row, *OUT_row, sum, tmp */;

	m = A->m;	n = A->n;	p = B->n;
	A_v = A->me;		B_v = B->me;

	if ( OUT==(MAT *)NULL || OUT->m != A->m || OUT->n != B->n )
		OUT = m_resize(OUT,A->m,B->n);

/****************************************************************
	for ( i=0; i<m; i++ )
		for  ( j=0; j<p; j++ )
		{
			sum = 0.0;
			for ( k=0; k<n; k++ )
				sum += A_v[i][k]*B_v[k][j];
			OUT->me[i][j] = sum;
		}
****************************************************************/
	m_zero(OUT);
	for ( i=0; i<m; i++ )
		for ( k=0; k<n; k++ )
		{
		    if ( A_v[i][k] != 0.0 )
		        __mltadd__(OUT->me[i],B_v[k],A_v[i][k],(int)p);
		    /**************************************************
		    B_row = B_v[k];	OUT_row = OUT->me[i];
		    for ( j=0; j<p; j++ )
			(*OUT_row++) += tmp*(*B_row++);
		    **************************************************/
		}

	return OUT;
}

/* px_get -- gets a PERM of given 'size' by dynamic memory allocation
	-- Note: initialized to the identity permutation
	-- the permutation is on the set {0,1,2,...,size-1} */
#ifndef ANSI_C
PERM	*px_get(size)
int	size;
#else
PERM	*px_get(int size)
#endif
{
   PERM	*permute;
   int	i;

   mem_bytes(TYPE_PERM,0,sizeof(PERM));
   mem_numvar(TYPE_PERM,1);
   
   permute->size = permute->max_size = size;
   
   mem_bytes(TYPE_PERM,0,size*sizeof(unsigned int));
   
   for ( i=0; i<size; i++ )
     permute->pe[i] = i;
   
   return (permute);
}

/* px_resize -- returns the permutation px with size new_size
   -- px is set to the identity permutation */
#ifndef ANSI_C
PERM	*px_resize(px,new_size)
PERM	*px;
int	new_size;
#else
PERM	*px_resize(PERM *px, int new_size)
#endif
{
   int	i;

   if ( ! px )
     return px_get(new_size);
   
   /* nothing is changed */
   if (new_size == px->size)
     return px;

   if ( new_size > px->max_size )
   {
	 mem_bytes(TYPE_PERM,px->max_size*sizeof(unsigned int),
		      new_size*sizeof(unsigned int));
      px->pe = RENEW(px->pe,new_size,unsigned int);
      px->max_size = new_size;
   }
   if ( px->size <= new_size )
     /* extend permutation */
     for ( i = px->size; i < new_size; i++ )
       px->pe[i] = i;
   else
     for ( i = 0; i < new_size; i++ )
       px->pe[i] = i;
   
   px->size = new_size;
   
   return px;
}

/* m_inverse -- returns inverse of A, provided A is not too rank deficient
	-- uses LU factorisation */
#ifndef ANSI_C
MAT	*m_inverse(A,out)
MAT	*A, *out;
#else
MAT	*m_inverse(const MAT *A, MAT *out)
#endif
{
	int	i;
	static VEC	*tmp = VNULL, *tmp2 = VNULL;
	static MAT	*A_cp = MNULL;
	static PERM	*pivot = PNULL;

	if ( ! out || out->m < A->m || out->n < A->n )
	    out = m_resize(out,A->m,A->n);

	A_cp = m_resize(A_cp,A->m,A->n);
	A_cp = m_copy(A,A_cp);
	tmp = v_resize(tmp,A->m);
	tmp2 = v_resize(tmp2,A->m);
	pivot = px_resize(pivot,A->m);
	MEM_STAT_REG(A_cp,TYPE_MAT);
	MEM_STAT_REG(tmp, TYPE_VEC);
	MEM_STAT_REG(tmp2,TYPE_VEC);
	MEM_STAT_REG(pivot,TYPE_PERM);
	for ( i = 0; i < A->n; i++ )
	{
	    v_zero(tmp);
	    tmp->ve[i] = 1.0;
	    set_col(out,i,tmp2);
	}

#ifdef	THREADSAFE
	V_FREE(tmp);	V_FREE(tmp2);
	M_FREE(A_cp);	PX_FREE(pivot);
#endif

	return out;
}

    int k_ss(MAT *A, MAT *B, MAT *C, MAT *D, double p, double u){
    	double k_ss, y_ss, x, c_bar;
    	int k_bar;
    	MAT *AUX = MNULL, *AUX1 = MNULL, *AUX2 = MNULL, *AUX3 = MNULL, *temp = MNULL;
    	MAT *Id = MNULL;
    	double lambdaMax;
    	int n=A->max_n;
    	Id = m_get(A->m,A->n);
    	//get the expression y_ss=(C(I-A)^(-1)B+D)u
    	Id = m_ident(Id);
    	AUX1 = m_sub(Id,A,AUX1);
    	AUX2 = m_inverse(AUX1,AUX2);
    	AUX = m_mlt(C,AUX2,AUX);;
    	AUX3 = m_mlt(AUX,B,AUX3);
    	temp = m_add(AUX3,D,temp);

    	y_ss=((double)temp->me[0][0])*u;
    	printf("y_ss=%f\n",y_ss);
    	lambdaMax=maxMagEigVal(getEigvalues(A),(size_t)sizeof(getEigvalues(A)));
    	printf("LambMax=%f\n",lambdaMax);
    	if(y_ss<0){
    		y_ss=-y_ss;
    	}

    	c_bar = m_norm_inf(A);printf("c_bar=%f\n",c_bar);
    	x = (p*y_ss)/(100*n*c_bar);
    	printf("x=%f\n",x);
    	k_ss=log_b(lambdaMax,x);
    	k_ss=(int)ceil(k_ss);
		//AUX = m_pow(A);
    	return k_bar=(int)k_ss;
    }

    double my_max(double* v){
    	double x=v[0];
    	for(int i=1;i<size(v);i++){
    		if(v[i]>x){
    			x=v[i];
    		}
    	}
    	return x;
    }

/* _m_pow -- computes integer powers of a square matrix A, A^p
   -- uses tmp as temporary workspace */
#ifndef ANSI_C
MAT	*_m_pow(A, p, tmp, out)
MAT	*A, *tmp, *out;
int	p;
#else
MAT	*_m_pow(const MAT *A, int p, MAT *tmp, MAT *out)
#endif
{
   int		it_cnt, k, max_bit;
   
   /*
     File containing routines for evaluating matrix functions
     esp. the exponential function
     */

#define	Z(k)	(((k) & 1) ? tmp : out)
   
   out = m_resize(out,A->m,A->n);
   tmp = m_resize(tmp,A->m,A->n);
   
   if ( p == 0 )
     m_ident(out);
   else if ( p > 0 )
   {
      it_cnt = 1;
      for ( max_bit = 0; ; max_bit++ )
	if ( (p >> (max_bit+1)) == 0 )
	  break;
      tmp = m_copy(A,tmp);
      
      for ( k = 0; k < max_bit; k++ )
      {
	 m_mlt(Z(it_cnt),Z(it_cnt),Z(it_cnt+1));
	 it_cnt++;
	 if ( p & (1 << (max_bit-1)) )
	 {
	    m_mlt(A,Z(it_cnt),Z(it_cnt+1));
	    /* m_copy(Z(it_cnt),out); */
	    it_cnt++;
	 }
	 p <<= 1;
      }
      if (it_cnt & 1)
	out = m_copy(Z(it_cnt),out);
   }

   return out;

#undef Z   
}

/* m_free -- returns MAT & asoociated memory back to memory heap */
#ifndef ANSI_C
int	m_free(mat)
MAT	*mat;
#else
int	m_free(MAT *mat)
#endif
{
#ifdef SEGMENTED
   int	i;
#endif
   
   if ( mat==(MAT *)NULL || (int)(mat->m) < 0 ||
       (int)(mat->n) < 0 )
     /* don't trust it */
     return (-1);
   
#ifndef SEGMENTED
   if ( mat->base != (double *)NULL ) {
	 mem_bytes(TYPE_MAT,mat->max_m*mat->max_n*sizeof(double),0);
      
      free((char *)(mat->base));
   }
#else
   for ( i = 0; i < mat->max_m; i++ )
     if ( mat->me[i] != (double *)NULL ) {
	   mem_bytes(TYPE_MAT,mat->max_n*sizeof(double),0);
	
	free((char *)(mat->me[i]));
     }
#endif
   if ( mat->me != (double **)NULL ) {
	 mem_bytes(TYPE_MAT,mat->max_m*sizeof(double *),0);
      
      free((char *)(mat->me));
   }

      mem_bytes(TYPE_MAT,sizeof(MAT),0);
      mem_numvar(TYPE_MAT,-1);
   free((char *)mat);
   
   return (0);
}



/* px_free -- returns PERM & asoociated memory back to memory heap */
#ifndef ANSI_C
int	px_free(px)
PERM	*px;
#else
int	px_free(PERM *px)
#endif
{
   if ( px==(PERM *)NULL || (int)(px->size) < 0 )
     /* don't trust it */
     return (-1);
   
   if ( px->pe == (unsigned int *)NULL ) {
	 mem_bytes(TYPE_PERM,sizeof(PERM),0);
	 mem_numvar(TYPE_PERM,-1);
      
      free((char *)px);
   }
   else
   {
	 mem_bytes(TYPE_PERM,sizeof(PERM)+px->max_size*sizeof(unsigned int),0);
	 mem_numvar(TYPE_PERM,-1);
     free((char *)px->pe);
     free((char *)px);
   }
   
   return (0);
}



/* v_free -- returns VEC & asoociated memory back to memory heap */
#ifndef ANSI_C
int	v_free(vec)
VEC	*vec;
#else
int	v_free(VEC *vec)
#endif
{
   if ( vec==(VEC *)NULL || (int)(vec->dim) < 0 )
     /* don't trust it */
     return (-1);
   
   if ( vec->ve == (double *)NULL ) {
	 mem_bytes(TYPE_VEC,sizeof(VEC),0);
	 mem_numvar(TYPE_VEC,-1);
      free((char *)vec);
   }
   else
   {
	 mem_bytes(TYPE_VEC,sizeof(VEC)+vec->max_dim*sizeof(double),0);
	 mem_numvar(TYPE_VEC,-1);
      free((char *)vec->ve);
      free((char *)vec);
   }
   
   return (0);
}

/* macros that also check types and sets pointers to NULL */
#define	M_FREE(mat)	( m_free(mat),	(mat)=(MAT *)NULL )
#define V_FREE(vec)	( v_free(vec),	(vec)=(VEC *)NULL )
#define	PX_FREE(px)	( px_free(px),	(px)=(PERM *)NULL )

/* m_pow -- computes integer powers of a square matrix A, A^p */
#ifndef ANSI_C
MAT	*m_pow(A, p, out)
MAT	*A, *out;
int	p;
#else
MAT	*m_pow(const MAT *A, int p, MAT *out)
#endif
{
   static MAT	*wkspace=MNULL, *tmp=MNULL;
   
   
   wkspace = m_resize(wkspace,A->m,A->n);
   MEM_STAT_REG(wkspace,TYPE_MAT);
   if ( p < 0 )
   {
       tmp = m_resize(tmp,A->m,A->n);
       MEM_STAT_REG(tmp,TYPE_MAT);
       out = _m_pow(tmp, -p, wkspace, out);
   }
   else
       out = _m_pow(A, p, wkspace, out);

#ifdef	THREADSAFE
   M_FREE(wkspace);	M_FREE(tmp);
#endif

   return out;
}

static const char    *format = "%14.9g ";

/* m_foutput -- prints a representation of the matrix a onto file/stream fp */
#ifndef ANSI_C
void    m_foutput(fp,a)
FILE    *fp;
MAT     *a;
#else
void    m_foutput(FILE *fp, const MAT *a)
#endif
{
     unsigned int      i, j, tmp;
     
     if ( a == (MAT *)NULL )
     {  fprintf(fp,"Matrix: NULL\n");   return;         }
     fprintf(fp,"Matrix: %d by %d\n",a->m,a->n);
     if ( a->me == (double **)NULL )
     {  fprintf(fp,"NULL\n");           return;         }
     for ( i=0; i<a->m; i++ )   /* for each row... */
     {
	  fprintf(fp,"row %u: ",i);
	  for ( j=0, tmp=2; j<a->n; j++, tmp++ )
	  {             /* for each col in row... */
	       fprintf(fp,format,a->me[i][j]);
	       if ( ! (tmp % 5) )       putc('\n',fp);
	  }
	  if ( tmp % 5 != 1 )   putc('\n',fp);
     }
}

void main(){
    MAT *A = MNULL, *B = MNULL, *C = MNULL, *D = MNULL, *T = MNULL, *Q = MNULL, *X_re = MNULL, *X_im = MNULL, *Q1 = MNULL, *Q1_inv = MNULL;
    MAT *Q1_temp, *Test = MNULL;
    VEC *evals_re = VNULL, *evals_im = VNULL;
    MAT *F = MNULL, *G = MNULL, *H = MNULL;
    int k=3;
    double y, x0;
    complex double *z;
    //ZMAT *ZQ = ZMNULL, *ZQ_temp, *ZQ_inv = ZMNULL, *ZH, *ZF;

   //setting up A matrix
//    A=m_get(4,4);
//    A->me[0][0]=-0.5000;A->me[0][1]=0.6000;A->me[0][2]=0;A->me[0][3]=0;
//    A->me[1][0]=-0.6000;A->me[1][1]=-0.5000;A->me[1][2]=0;A->me[1][3]=0;
//    A->me[2][0]=0;A->me[2][1]=0;A->me[2][2]=0.2000;A->me[2][3]=0.8000;
//    A->me[3][0]=0;A->me[3][1]=0;A->me[3][2]=-0.8000;A->me[3][3]=0.2000;printf("A ");m_output(A);
    A=m_get(5,5);
    A->me[0][0]=-0.5000;A->me[0][1]=0.6000;A->me[0][2]=0;A->me[0][3]=0;A->me[0][4]=0;
    A->me[1][0]=-0.6000;A->me[1][1]=-0.5000;A->me[1][2]=0;A->me[1][3]=0;A->me[1][4]=0;
    A->me[2][0]=0;A->me[2][1]=0;A->me[2][2]=0.2000;A->me[2][3]=0.8000;A->me[2][4]=0;
    A->me[3][0]=0;A->me[3][1]=0;A->me[3][2]=-0.8000;A->me[3][3]=0.2000;A->me[3][4]=0;
    A->me[4][0]=0;A->me[4][1]=0;A->me[4][2]=0;A->me[4][3]=0;A->me[4][4]=0.6;printf("A ");m_output(A);
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
//    C->me[0][0]=0;C->me[0][1]=2.6;C->me[0][2]=0.5;C->me[0][3]=1.2;printf("C ");m_output(C);
        C=m_get(1,5);
        C->me[0][0]=0;C->me[0][1]=2.6;C->me[0][2]=0.5;C->me[0][3]=1.2;C->me[0][4]=0;printf("C ");m_output(C);
    //setting up D matrix
    D=m_get(1,1);
    D->me[0][0]=0;printf("D ");m_output(D);
    printf("-----------------------------------------------------------\n");
    printf("k_ss=%d\n",k_ss(A,B,C,D,5,1.0f));
    Test = m_pow(A,2,Test);
	m_output(Test);

//    /* read in A matrix */
//    printf("Input A matrix:\n");
//
//    A = m_input(MNULL);     /* A has whatever size is input */
//    //B = m_input(MNULL);     /* B has whatever size is input */
//
//    if ( A->m < A->n )
//    {
//        printf("Need m >= n to obtain least squares fit\n");
//        exit(0);
//    }
//    printf("# A =\n");       m_output(A);
//
//    //zm_output(zm_A_bar(A));
//
//    Q = m_get(A->m,A->n);
//    T = m_copy(A,MNULL);
//
//    /* compute Schur form: A = Q.T.Q^T */
//    schur(T,Q);
//    /* extract eigenvalues */
//    evals_re = v_get(A->m);
//    evals_im = v_get(A->m);
//    schur_evals(T,evals_re,evals_im);
//
//    z=malloc(evals_re->dim*sizeof(complex double));
//    for(int i=0;i<evals_re->dim;i++){
//    	z[i]=evals_re->ve[i]+I*evals_im->ve[i];
//    	printf("Z[%d]=%f + i%f\n", i, creal(z[i]), cimag(z[i]));
//    }
//
//    size_t size=(size_t)sizeof(z);
//    printf("Maximum:%f\n",maxMagEigVal(z,size));

}