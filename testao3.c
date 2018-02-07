#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <malloc.h>

/* miscellaneous constants */
#define	VNULL	((VEC *)NULL)
#define	MNULL	((MAT *)NULL)
#define	PNULL	((PERM *)NULL)

/* macros that also check types and sets pointers to NULL */
#define	M_FREE(mat)	( m_free(mat),	(mat)=(MAT *)NULL )
#define V_FREE(vec)	( v_free(vec),	(vec)=(VEC *)NULL )
#define	PX_FREE(px)	( px_free(px),	(px)=(PERM *)NULL )

/* available standard types */
#define TYPE_NULL              (-1)
//#define TYPE_MAT    	        0
#define TYPE_BAND               1
//#define TYPE_PERM		2
//#define TYPE_VEC		3
#define TYPE_IVEC		4

#define	m_output(mat)	m_foutput(stdout,mat)

static const char    *format = "%14.9g ";

#ifndef ANSI_C
#define ANSI_C 1
#endif
#ifndef CHAR0ISDBL0
#define CHAR0ISDBL0 1
#endif

#define SEGMENTED

#ifndef THREADSAFE	/* for use as a shared library */
#define	THREADSAFE 1
#endif

#define TYPE_MAT mem_bytes(0,0,sizeof(MAT))
#define TYPE_VEC mem_bytes(0,0,sizeof(VEC))
#define TYPE_PERM mem_bytes(0,0,sizeof(PERM))
#define MEM_CONNECT_MAX_LISTS 3

/* standard copy & zero functions */
#define	MEM_COPY(from,to,size)	memmove((to),(from),(size))
#define	MEM_ZERO(where,size)	memset((where),'\0',(size))

/* allocate one object of given type */
#define	NEW(type)	((type *)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type */
#define	NEW_A(num,type)	((type *)calloc((size_t)(num),(size_t)sizeof(type)))

 /* re-allocate arry to have num objects of the given type */
#define	RENEW(var,num,type) \
    ((var)=(type *)((var) ? \
		    realloc((char *)(var),(size_t)((num)*sizeof(type))) : \
		    calloc((size_t)(num),(size_t)sizeof(type))))

#define	MEMCOPY(from,to,n_items,type) \
 MEM_COPY((char *)(from),(char *)(to),(unsigned)(n_items)*sizeof(type))

/* type independent min and max operations */
#ifndef max
#define	max(a,b)	((a) > (b) ? (a) : (b))
#endif /* max */
#ifndef min
#define	min(a,b)	((a) > (b) ? (b) : (a))
#endif /* min */

#ifndef THREADSAFE
#define mem_stat_reg(var,type)  mem_stat_reg_list((void **)var,type,0,__FILE__,__LINE__)
#define MEM_STAT_REG(var,type)  mem_stat_reg_list((void **)&(var),type,0,__FILE__,__LINE__)
#define mem_stat_free(mark)   mem_stat_free_list(mark,0)
#else
#define mem_stat_reg(var,type)
#define MEM_STAT_REG(var,type)
#define mem_stat_free(mark)
#endif

/* matrix definition */
typedef	struct	{
		unsigned int	m, n;
		unsigned int	max_m, max_n, max_size;
		double	**me,*base;	/* base is base of alloc'd mem */
} MAT;

/* vector definition */
typedef	struct	{
		unsigned int	dim, max_dim;
		double	*ve;
} VEC;

/* permutation definition */
typedef	struct	{
		unsigned int	size, max_size, *pe;
} PERM;

/* m_add -- matrix addition -- may be in-situ */
#ifndef ANSI_C
MAT	*m_add(mat1,mat2,out)
MAT	*mat1,*mat2,*out;
#else
MAT	*m_add(const MAT *mat1, const MAT *mat2, MAT *out)
#endif
{
	unsigned int	m,n,i,j;

	/*if ( out==(MAT *)NULL || out->m != mat1->m || out->n != mat1->n )
		out = m_resize(out,mat1->m,mat1->n);*/
	m = mat1->m;	n = mat1->n;
	for ( i=0; i<m; i++ )
	{
//		__add__(mat1->me[i],mat2->me[i],out->me[i],(int)n);
		/**************************************************/
		for ( j=0; j<n; j++ )
			out->me[i][j] = mat1->me[i][j]+mat2->me[i][j];
		/**************************************************/
	}

	return (out);
}

/* m_sub -- matrix subtraction -- may be in-situ */
#ifndef ANSI_C
MAT	*m_sub(mat1,mat2,out)
MAT	*mat1,*mat2,*out;
#else
MAT	*m_sub(const MAT *mat1, const MAT *mat2, MAT *out)
#endif
{
	unsigned int	m,n,i,j;

	/*if ( out==(MAT *)NULL || out->m != mat1->m || out->n != mat1->n )
		out = m_resize(out,mat1->m,mat1->n);*/
	m = mat1->m;	n = mat1->n;
	for ( i=0; i<m; i++ )
	{
//		__sub__(mat1->me[i],mat2->me[i],out->me[i],(int)n);
		/**************************************************/
		for ( j=0; j<n; j++ )
			out->me[i][j] = mat1->me[i][j]-mat2->me[i][j];
		/**************************************************/
	}

	return (out);
}

/* m_get -- gets an mxn matrix (in MAT form) by dynamic memory allocation
	-- normally ALL matrices should be obtained this way
	-- if either m or n is negative this will raise an error
	-- note that 0 x n and m x 0 matrices can be created */
MAT	*m_get(int m, int n)
{
   MAT	*matrix = malloc(sizeof *matrix);
   int	i,j;

   if (m < 0 || n < 0)
        printf("The matrix dimensions must be positive!\n");
   if ((matrix=NEW(MAT)) == (MAT *)NULL )
	   printf("The matrix is NULL!\n");

   matrix->m = m;		matrix->n = matrix->max_n = n;
   matrix->max_m = m;	matrix->max_size = m*n;

   matrix->me = (double **)malloc(m * sizeof(double*));
   for(int i = 0; i < m; i++)
	   matrix->me[i] = (double *)malloc(n * sizeof(double));

   return (matrix);
}

/* m_resize -- returns the matrix A of size new_m x new_n; A is zeroed
   -- if A == NULL on entry then the effect is equivalent to m_get() */
MAT	*m_resize(MAT *A,int new_m, int new_n)
{
   int	i;
   int	new_max_m, new_max_n, old_m, old_n, add_rows;
   double **tmp;

   if (new_m < 0 || new_n < 0)
     printf("The size must be positive!");

   if ( ! A )
     return m_get(new_m,new_n);

   // nothing was changed
   if (new_m == A->m && new_n == A->n)
     return A;

   old_m = A->m;	old_n = A->n;	add_rows = new_m-old_m;
   if ( new_m > A->max_m )
   {	// re-allocate A->me
	   tmp = realloc( A->me, sizeof *A->me * (new_m) );
	   if ( tmp )
	   {
		 A->me = tmp;
	     for (i = 0; i < add_rows; i++ )
	     {
	    	 A->me[old_m + i] = malloc( sizeof *A->me[old_m + i] * old_n );
	     }
	   }
	   if ( new_n > A->max_n )
	   {
		   double *tmp;
		   for ( int i = 0; i < old_m; i++ )
		   {
		     tmp = realloc( A->me[i], sizeof *A->me[i] * (new_n) );
		     if ( tmp )
		     {
		    	A->me[i] = tmp;
		     }
		   }
	   }
	   else if ( new_n < A->max_n )
	   {
		   double *tmp;
		   for ( int i = 0; i < old_n; i++ )
		   {
		     tmp = realloc( A->me[i], sizeof *A->me[i] * (new_n) );
		     if ( tmp )
		    	A->me[i] = tmp;
		   }
	   }
   }
   else if ( new_m < A->max_m )
   {
	   int del_rows = old_m-new_m;
	   double *tmp;
	   for ( int i = 1; i <= del_rows; i++ ){
	     free( A->me[old_m - i] );

	     tmp = realloc( A->me, old_m - del_rows );
	     if ( tmp )
		   A->me[i] = tmp;
	   }
	   if ( new_n > A->max_n )
	   {
	     double *tmp;
	   	 for ( int i = 0; i < old_m; i++ )
	   	 {
	   	   tmp = realloc( A->me[i], sizeof *A->me[i] * (new_n) );
	   	   if ( tmp )
	   	   {
	   		 A->me[i] = tmp;
	   	   }
	   	 }
	   }
	   else if ( new_n < A->max_n )
	   {
	     double *tmp;
	   	 for ( int i = 0; i < old_n; i++ )
	   	 {
	   	   tmp = realloc( A->me[i], sizeof *A->me[i] * (new_n) );
	   	   if ( tmp )
	   	   	 A->me[i] = tmp;
	   	 }
	   }
   }

   new_max_m = max(new_m,A->max_m);
   new_max_n = max(new_n,A->max_n);

   A->max_m = new_max_m;
   A->max_n = new_max_n;
   A->max_size = A->max_m*A->max_n;
   A->m = new_m;	A->n = new_n;

   return A;
}

/* m_zero -- zero the matrix A */
#ifndef ANSI_C
MAT	*m_zero(A)
MAT	*A;
#else
MAT	*m_zero(MAT *A)
#endif
{
	int	i, j, A_m, A_n;
	double	**A_me;

	A_m = A->m;	A_n = A->n;	A_me = A->me;
	for ( i = 0; i < A_m; i++ )
//		__zero__(A_me[i],A_n);
		for ( j = 0; j < A_n; j++ )
			A_me[i][j] = 0.0;
	return A;
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

/* m_mlt -- matrix-matrix multiplication */
#ifndef ANSI_C
MAT	*m_mlt(A,B,OUT)
MAT	*A,*B,*OUT;
#else
MAT	*m_mlt(const MAT *A, const MAT *B, MAT *OUT)
#endif
{
	unsigned int	i, /* j, */ k, m, n, p;
	double	**A_v, **B_v, *B_row, *OUT_row, sum, tmp;

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

/* m_copy -- copies matrix into new area
	-- out(i0:m,j0:n) <- in(i0:m,j0:n) */
#ifndef ANSI_C
MAT	*m_copy(in,out/*,i0,j0*/)
MAT	*in,*out;
// unsigned int	i0,j0;
#else
MAT	*m_copy(const MAT *in, MAT *out)
#endif
{
	unsigned int i0 = 0, j0 = 0;
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

/* v_zero -- zero the vector x */
#ifndef ANSI_C
VEC	*v_zero(x)
VEC	*x;
#else
VEC	*v_zero(VEC *x)
#endif
{

	// __zero__(x->ve,x->dim);
	for (int i = 0; i < x->dim; i++ )
		x->ve[i] = 0.0; 

	return x;
}

/* v_get -- gets a VEC of dimension 'size'
   -- Note: initialized to zero */
#ifndef ANSI_C
VEC	*v_get(size)
int	size;
#else
VEC	*v_get(int size)
#endif
{
   VEC	*vector = malloc(sizeof *vector);
   int	i,j;

   vector->dim = vector->max_dim = size;
   if (size < 0)
        printf("The vector dimension must be positive!\n");
   if ((vector->ve=NEW_A(size,double)) == (double *)NULL )
   {
      free(vector);
   }
   else 
   {
	  vector->ve = (double *)malloc(size * sizeof(double));
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
   double *ptr;
   if ( ! x )
     return v_get(new_dim);

   /* nothing is changed */
   if (new_dim == x->dim)
     return x;

   if ( x->max_dim == 0 )	/* assume that it's from sub_vec */
     return v_get(new_dim);
   ptr = x->ve;
   if ( new_dim > x->max_dim )
   {
     ptr = realloc( ptr, new_dim * sizeof *ptr );
   }
   if ( new_dim > x->dim )
   {
     for (int i = 1; i < (new_dim - x->dim); i++ )
		x->ve[new_dim-i] = 0.0;
   }
   else if( new_dim < x->dim )
   {
	 ptr = realloc( ptr, new_dim * sizeof *ptr );
   }

   x->dim = new_dim;
   
   return x;
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
   PERM	*permute = malloc(sizeof *permute);
   int	i;
   
   permute->size = permute->max_size = size;

   permute->pe = (unsigned int *)malloc(size * sizeof(unsigned int));
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
   unsigned int *ptr;

   if ( ! px )
     return px_get(new_size);
   
   /* nothing is changed */
   if (new_size == px->size)
     return px;
   
   ptr = px->pe;
   if ( new_size > px->max_size )
   {
     ptr = realloc( ptr, new_size * sizeof *ptr );
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

/* set_col -- sets column of matrix to values given in vec (in situ)
	-- that is, mat(i0:lim,col) <- vec(i0:lim) */
#ifndef ANSI_C
MAT	*set_col(mat,col,vec)
MAT	*mat;
VEC	*vec;
unsigned int	col;
#else
MAT	*set_col(MAT *mat, unsigned int col, const VEC *vec/*, unsigned int i0*/)
#endif
{
   unsigned int	i,lim,i0;
   
   lim = min(mat->m,vec->dim);
   for ( i=i0; i<lim; i++ )
     mat->me[i][col] = vec->ve[i];
   
   return (mat);
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
     mat->base = (double *)malloc(mat->max_m*mat->max_n * sizeof(double));
     free((char *)(mat->base));
   }
#else
   for ( i = 0; i < mat->max_m; i++ )
     if ( mat->me[i] != (double *)NULL ) {
	   mat->me[i] = (double *)malloc(mat->max_n * sizeof(double));
	free((char *)(mat->me[i]));
     }
#endif
   if ( mat->me != (double **)NULL ) {
     mat->me = (double **)malloc(mat->max_m * sizeof(double*));
      free((char *)(mat->me));
   }

	  mat = malloc(sizeof *mat);
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
	 px = malloc(sizeof *px);
      
      free((char *)px);
   }
   else
   {
	 px = malloc(sizeof *px);
	 px->pe = (unsigned int *)malloc(px->max_size*sizeof(unsigned int));
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
	 vec = malloc(sizeof *vec);
      free((char *)vec);
   }
   else
   {
	 vec = malloc(sizeof *vec);
	 vec->ve = (double *)malloc(vec->max_dim*sizeof(double));
      free((char *)vec->ve);
      free((char *)vec);
   }
   
   return (0);
}

/* _v_copy -- copies vector into new area
	-- out(i0:dim) <- in(i0:dim) */
#ifndef ANSI_C
VEC	*v_copy(in,out)
VEC	*in,*out;
#else
VEC	*v_copy(const VEC *in, VEC *out/*, unsigned int i0*/)
#endif
{
	/* unsigned int	i,j; */
	unsigned int i0 = 0;

	if ( in==out )
		return (out);
	if ( out==VNULL || out->dim < in->dim )
		out = v_resize(out,in->dim);

	MEM_COPY(&(in->ve[i0]),&(out->ve[i0]),(in->dim - i0)*sizeof(double));
	/* for ( i=i0; i < in->dim; i++ )
		out->ve[i] = in->ve[i]; */

	return (out);
}

/* px_vec -- permute vector */
#ifndef ANSI_C
VEC	*px_vec(px,vector,out)
PERM	*px;
VEC	*vector,*out;
#else
VEC	*px_vec(PERM *px, const VEC *vector, VEC *out)
#endif
{
    unsigned int	old_i, i, size, start;
    double	tmp;
    
    // if ( px==PNULL || vector==VNULL )
	// error(E_NULL,"px_vec");
    // if ( px->size > vector->dim )
	// error(E_SIZES,"px_vec");
    if ( out==VNULL || out->dim < vector->dim )
	out = v_resize(out,vector->dim);
    
    size = px->size;
    if ( size == 0 )
	return v_copy(vector,out);
    if ( out != vector )
    {
	for ( i=0; i<size; i++ )
	    // if ( px->pe[i] >= size )
		// error(E_BOUNDS,"px_vec");
	    // else
		out->ve[i] = vector->ve[px->pe[i]];
    }
    else
    {	/* in situ algorithm */
	start = 0;
	while ( start < size )
	{
	    old_i = start;
	    i = px->pe[old_i];
	    if ( i >= size )
	    {
		start++;
		continue;
	    }
	    tmp = vector->ve[start];
	    while ( 1 )
	    {
		vector->ve[old_i] = vector->ve[i];
		px->pe[old_i] = i+size;
		old_i = i;
		i = px->pe[old_i];
		if ( i >= size )
		    break;
		if ( i == start )
		{
		    vector->ve[old_i] = tmp;
		    px->pe[old_i] = i+size;
		    break;
		}
	    }
	    start++;
	}

	for ( i = 0; i < size; i++ )
		px->pe[i] = px->pe[i]-size;
    }
    
    return out;
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

/* px_transp -- transpose elements of permutation
		-- Really multiplying a permutation by a transposition */
#ifndef ANSI_C
PERM	*px_transp(px,i1,i2)
PERM	*px;		/* permutation to transpose */
unsigned int	i1,i2;		/* elements to transpose */
#else
PERM	*px_transp(PERM *px, unsigned int i1, unsigned int i2)
#endif
{
	unsigned int	temp;

	if ( px==(PERM *)NULL )
		printf("PERM is NULL!\n");

	if ( i1 < px->size && i2 < px->size )
	{
		temp = px->pe[i1];
		px->pe[i1] = px->pe[i2];
		px->pe[i2] = temp;
	}

	return px;
}

/* LUfactor -- gaussian elimination with scaled partial pivoting
		-- Note: returns LU matrix which is A */
#ifndef ANSI_C
MAT	*LUfactor(A,pivot)
MAT	*A;
PERM	*pivot;
#else
MAT	*LUfactor(MAT *A, PERM *pivot)
#endif
{
	unsigned int	i, j, m, n;
	int	i_max, k, k_max;
	double	**A_v, *A_piv, *A_row;
	double	max1, temp, tiny;
	static	VEC	*scale = VNULL;

	// if ( A==(MAT *)NULL || pivot==(PERM *)NULL )
	// 	error(E_NULL,"LUfactor");
	// if ( pivot->size != A->m )
	// 	error(E_SIZES,"LUfactor");
	m = A->m;	n = A->n;
	scale = v_resize(scale,A->m);
	MEM_STAT_REG(scale,TYPE_VEC);
	A_v = A->me;

	tiny = 10.0/HUGE_VAL;

	/* initialise pivot with identity permutation */
	for ( i=0; i<m; i++ )
		pivot->pe[i] = i;

	/* set scale parameters */
	for ( i=0; i<m; i++ )
	{
		max1 = 0.0;
		for ( j=0; j<n; j++ )
		{
			temp = fabs(A_v[i][j]);
			max1 = max(max1,temp);
		}
		scale->ve[i] = max1;
	}

	/* main loop */
	k_max = min(m,n)-1;
	for ( k=0; k<k_max; k++ )
	{
	    /* find best pivot row */
	    max1 = 0.0;	i_max = -1;
	    for ( i=k; i<m; i++ )
		if ( fabs(scale->ve[i]) >= tiny*fabs(A_v[i][k]) )
		{
		    temp = fabs(A_v[i][k])/scale->ve[i];
		    if ( temp > max1 )
		    { max1 = temp;	i_max = i;	}
		}
	    
	    /* if no pivot then ignore column k... */
	    if ( i_max == -1 )
	    {
		/* set pivot entry A[k][k] exactly to zero,
		   rather than just "small" */
		A_v[k][k] = 0.0;
		continue;
	    }
	    
	    /* do we pivot ? */
	    if ( i_max != k )	/* yes we do... */
	    {
		px_transp(pivot,i_max,k);
		for ( j=0; j<n; j++ )
		{
		    temp = A_v[i_max][j];
		    A_v[i_max][j] = A_v[k][j];
		    A_v[k][j] = temp;
		}
	    }
	    
	    /* row operations */
	    for ( i=k+1; i<m; i++ )	/* for each row do... */
	    {	/* Note: divide by zero should never happen */
		temp = A_v[i][k] = A_v[i][k]/A_v[k][k];
		A_piv = &(A_v[k][k+1]);
		A_row = &(A_v[i][k+1]);
		if ( k+1 < n )
		    __mltadd__(A_row,A_piv,-temp,(int)(n-(k+1)));
		/*********************************************
		  for ( j=k+1; j<n; j++ )
		  A_v[i][j] -= temp*A_v[k][j];
		  (*A_row++) -= temp*(*A_piv++);
		  *********************************************/
	    }
	    
	}

#ifdef	THREADSAFE
	V_FREE(scale);
#endif

	return A;
}

/* Usolve -- back substitution with optional over-riding diagonal
		-- can be in-situ but doesn't need to be */
#ifndef ANSI_C
VEC	*Usolve(matrix,b,out,diag)
MAT	*matrix;
VEC	*b, *out;
double	diag;
#else
VEC	*Usolve(const MAT *matrix, const VEC *b, VEC *out, double diag)
#endif
{
	unsigned int	dim /* , j */;
	int	i, i_lim;
	double	**mat_ent, *mat_row, *b_ent, *out_ent, *out_col, sum, tiny;

	// if ( matrix==MNULL || b==VNULL )
	// 	error(E_NULL,"Usolve");
	// dim = min(matrix->m,matrix->n);
	// if ( b->dim < dim )
	// 	error(E_SIZES,"Usolve");
	if ( out==VNULL || out->dim < dim )
		out = v_resize(out,matrix->n);
	mat_ent = matrix->me;	b_ent = b->ve;	out_ent = out->ve;

	tiny = 10.0/HUGE_VAL;

	for ( i=dim-1; i>=0; i-- )
		if ( b_ent[i] != 0.0 )
		    break;
		else
		    out_ent[i] = 0.0;
	i_lim = i;

	for (    ; i>=0; i-- )
	{
		sum = b_ent[i];
		mat_row = &(mat_ent[i][i+1]);
		out_col = &(out_ent[i+1]);
		sum -= __ip__(mat_row,out_col,i_lim-i);
		/******************************************************
		for ( j=i+1; j<=i_lim; j++ )
			sum -= mat_ent[i][j]*out_ent[j];
			sum -= (*mat_row++)*(*out_col++);
		******************************************************/
		if ( diag==0.0 )
		{
			if ( fabs(mat_ent[i][i]) <= tiny*fabs(sum) )
				printf("Element too small!\n");
			else
				out_ent[i] = sum/mat_ent[i][i];
		}
		else
			out_ent[i] = sum/diag;
	}

	return (out);
}

/* Lsolve -- forward elimination with (optional) default diagonal value */
#ifndef ANSI_C
VEC	*Lsolve(matrix,b,out,diag)
MAT	*matrix;
VEC	*b,*out;
double	diag;
#else
VEC	*Lsolve(const MAT *matrix, const VEC *b, VEC *out, double diag)
#endif
{
	unsigned int	dim, i, i_lim /* , j */;
	double	**mat_ent, *mat_row, *b_ent, *out_ent, *out_col, sum, tiny;

	// if ( matrix==(MAT *)NULL || b==(VEC *)NULL )
	// 	error(E_NULL,"Lsolve");
	// dim = min(matrix->m,matrix->n);
	// if ( b->dim < dim )
	// 	error(E_SIZES,"Lsolve");
	if ( out==(VEC *)NULL || out->dim < dim )
		out = v_resize(out,matrix->n);
	mat_ent = matrix->me;	b_ent = b->ve;	out_ent = out->ve;

	for ( i=0; i<dim; i++ )
		if ( b_ent[i] != 0.0 )
		    break;
		else
		    out_ent[i] = 0.0;
	i_lim = i;

	tiny = 10.0/HUGE_VAL;

	for (    ; i<dim; i++ )
	{
		sum = b_ent[i];
		mat_row = &(mat_ent[i][i_lim]);
		out_col = &(out_ent[i_lim]);
		sum -= __ip__(mat_row,out_col,(int)(i-i_lim));
		/*****************************************************
		for ( j=i_lim; j<i; j++ )
			sum -= mat_ent[i][j]*out_ent[j];
			sum -= (*mat_row++)*(*out_col++);
		******************************************************/
		if ( diag==0.0 )
		{
			if ( fabs(mat_ent[i][i]) <= tiny*fabs(sum) )
				printf("Element too small!\n");
			else
				out_ent[i] = sum/mat_ent[i][i];
		}
		else
			out_ent[i] = sum/diag;
	}

	return (out);
}

/* LUsolve -- given an LU factorisation in A, solve Ax=b */
#ifndef ANSI_C
VEC	*LUsolve(LU,pivot,b,x)
MAT	*LU;
PERM	*pivot;
VEC	*b,*x;
#else
VEC	*LUsolve(const MAT *LU, PERM *pivot, const VEC *b, VEC *x)
#endif
{
    x = v_resize(x,b->dim);
	px_vec(pivot,b,x);	/* x := P.b */
	Lsolve(LU,x,x,1.0);	/* implicit diagonal = 1 */
	Usolve(LU,x,x,0.0);	/* explicit diagonal */

	return (x);
}

// /* m_inverse -- returns inverse of A, provided A is not too rank deficient
// 	-- uses LU factorisation */
// #ifndef ANSI_C
// MAT	*m_inverse(A,out)
// MAT	*A, *out;
// #else
// MAT	*m_inverse(const MAT *A, MAT *out)
// #endif
// {
// 	printf("here!\n");
// 	int	i;
// 	static VEC	*tmp = VNULL, *tmp2 = VNULL;
// 	static MAT	*A_cp = MNULL;
// 	static PERM	*pivot = PNULL;

// 	if ( ! out || out->m < A->m || out->n < A->n )
// 	{
// 	    out = m_resize(out,A->m,A->n);
// 		printf("hereT!\n");
// 	}
// 	A_cp = m_resize(A_cp,A->m,A->n);
// 	printf("here2!\n");
// 	A_cp = m_copy(A,A_cp);
// 	printf("here3!\n");
// 	tmp = v_resize(tmp,A->m);
// 	printf("here4!\n");
// 	tmp2 = v_resize(tmp2,A->m);
// 	printf("here5!\n");
// 	pivot = px_resize(pivot,A->m);
// 	printf("here6!\n");
// 	MEM_STAT_REG(A_cp,TYPE_MAT);
// 	printf("here7!\n");
// 	MEM_STAT_REG(tmp, TYPE_VEC);
// 	printf("here8!\n");
// 	MEM_STAT_REG(tmp2,TYPE_VEC);
// 	printf("here9!\n");
// 	MEM_STAT_REG(pivot,TYPE_PERM);
// 	printf("here10!\n");
// 	LUfactor(A_cp,pivot);
// 	printf("here11!\n");
// 	for ( i = 0; i < A->n; i++ )
// 	{
// 	    v_zero(tmp);
// 	    tmp->ve[i] = 1.0;
// 		LUsolve(A_cp,pivot,tmp,tmp2);
// 		printf("temp2=%f!\n",tmp->ve[i]);
// 		printf("here12!\n");
// 	    set_col(out,i,tmp2);
// 		printf("here13!\n");
// 	}

// #ifdef	THREADSAFE
// 	V_FREE(tmp);	V_FREE(tmp2);
// 	M_FREE(A_cp);	PX_FREE(pivot);
// #endif
// printf("here14!\n");
// 	return out;
// }

/* m_inverse -- returns inverse of A, provided A is not too rank deficient
	-- uses LU factorisation */
#ifndef ANSI_C
MAT	*m_inverse(A,out)
MAT	*A, *out;
#else
MAT	*m_inverse(const MAT *A, MAT *out)
#endif
{
  int i,j,k,matsize;
  float temp;

  matsize = A->m;
  
	for(i=0;i<matsize;i++)									//automatically initialize the unit matrix, e.g.
		for(j=0;j<matsize;j++)								//	-       -
			if(i==j)										// | 1  0  0 |
				out->me[i][j]=1;								// | 0  1  0 |
			else											// | 0  0  1 |
				out->me[i][j]=0;								//  -       -
/*---------------LoGiC starts here------------------*/		//procedure to make the matrix A to unit matrix
	for(k=0;k<matsize;k++)									//by some row operations,and the same row operations of
	{														//Unit mat. I gives the inverse of matrix A
		temp=A->me[k][k];										//'temp' stores the A[k][k] value so that A[k][k] will not change
		for(j=0;j<matsize;j++)								//during the operation A[i][j]/=A[k][k] when i=j=k
		{
			A->me[k][j]/=temp;									//it performs the following row operations to make A to unit matrix
			out->me[k][j]/=temp;									//R0=R0/A[0][0],similarly for I also R0=R0/A[0][0]
		}													//R1=R1-R0*A[1][0] similarly for I
		for(i=0;i<matsize;i++)								//R2=R2-R0*A[2][0]		,,
		{
			temp=A->me[i][k];									//R1=R1/A[1][1]
			for(j=0;j<matsize;j++)							//R0=R0-R1*A[0][1]
			{												//R2=R2-R1*A[2][1]
				if(i==k)
					break;									//R2=R2/A[2][2]
				A->me[i][j]-=A->me[k][j]*temp;						//R0=R0-R2*A[0][2]
				out->me[i][j]-=out->me[k][j]*temp;						//R1=R1-R2*A[1][2]
			}
		}
	}
/*---------------LoGiC ends here--------------------*/

	return out;
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

void main(){
	printf("testing \n");
    MAT *A = MNULL, *A2 = MNULL, *A3 = MNULL, *A4 = MNULL, *A5 = MNULL, *A6 = MNULL, *B = MNULL, *C = MNULL, *D = MNULL, *T = MNULL, *Q = MNULL, *X_re = MNULL, *X_im = MNULL, *Q1 = MNULL, *Q1_inv = MNULL;
    MAT *Q1_temp, *Test = MNULL;
    //VEC *evals_re = VNULL, *evals_im = VNULL;
    MAT *F = MNULL, *G = MNULL, *H = MNULL;
    int k=3;
    double y, x0;
    //complex double *z;
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
    A->me[4][0]=0;A->me[4][1]=0;A->me[4][2]=0;A->me[4][3]=0;A->me[4][4]=0.6;printf("A ");
    m_output(A);
    A2=m_get(5,5);
    A2 = m_add(A, A, A2);
    m_output(A2);
    A3=m_get(5,5);
    A3 = m_sub(A, A, A3);
    m_output(A3);
    A4=m_get(5,5);
    A4 = m_mlt(A, A, A4);
    m_output(A4);
	A5=m_get(5,5);
    A5 = m_inverse(A,A5);
    m_output(A5);
	A6=m_get(5,5);
    A6 = m_pow(A,3,A6);
    m_output(A6);
//    printf("testing /n");
    //setting up B matrix
//    B=m_get(4,1);
//    B->me[0][0]=0;
//    B->me[1][0]=0;
//    B->me[2][0]=2.5;
//    B->me[3][0]=1;printf("B ");m_output(B);
    /*B=m_get(5,1);
    B->me[0][0]=0;
    B->me[1][0]=0;
    B->me[2][0]=2.5;
    B->me[3][0]=1;
    B->me[4][0]=0;printf("B ");m_output(B);*/
    //setting up C matrix
//    C=m_get(1,4);
//    C->me[0][0]=0;C->me[0][1]=2.6;C->me[0][2]=0.5;C->me[0][3]=1.2;printf("C ");m_output(C);
        /*C=m_get(1,5);
        C->me[0][0]=0;C->me[0][1]=2.6;C->me[0][2]=0.5;C->me[0][3]=1.2;C->me[0][4]=0;printf("C ");m_output(C);*/
    //setting up D matrix
    /*D=m_get(1,1);
    D->me[0][0]=0;printf("D ");m_output(D);
    printf("-----------------------------------------------------------\n");
    printf("k_ss=%d\n",k_ss(A,B,C,D,5,1.0f));
    Test = m_pow(A,2,Test);
	m_output(Test);*/

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
