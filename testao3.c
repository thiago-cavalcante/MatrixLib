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

#define	m_output(mat)	m_foutput(stdout,mat)

static const char    *format = "%14.9g ";

#ifndef ANSI_C
#define ANSI_C 1
#endif
#ifndef CHAR0ISDBL0
#define CHAR0ISDBL0 1
#endif

#define SEGMENTED

#define TYPE_MAT mem_bytes(0,0,sizeof(MAT))
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


/* matrix definition */
typedef	struct	{
		unsigned int	m, n;
		unsigned int	max_m, max_n, max_size;
		double	**me,*base;	/* base is base of alloc'd mem */
} MAT;


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

void main(){
	printf("testing \n");
    MAT *A = MNULL, *A2 = MNULL, *A3 = MNULL, *A4 = MNULL, *B = MNULL, *C = MNULL, *D = MNULL, *T = MNULL, *Q = MNULL, *X_re = MNULL, *X_im = MNULL, *Q1 = MNULL, *Q1_inv = MNULL;
    MAT *Q1_temp, *Test = MNULL;
    //VEC *evals_re = VNULL, *evals_im = VNULL;
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
    printf("testing \n");
    A->me[0][0]=-0.5000;A->me[0][1]=0.6000;A->me[0][2]=0;A->me[0][3]=0;A->me[0][4]=0;
    A->me[1][0]=-0.6000;A->me[1][1]=-0.5000;A->me[1][2]=0;A->me[1][3]=0;A->me[1][4]=0;
    A->me[2][0]=0;A->me[2][1]=0;A->me[2][2]=0.2000;A->me[2][3]=0.8000;A->me[2][4]=0;
    A->me[3][0]=0;A->me[3][1]=0;A->me[3][2]=-0.8000;A->me[3][3]=0.2000;A->me[3][4]=0;
    A->me[4][0]=0;A->me[4][1]=0;A->me[4][2]=0;A->me[4][3]=0;A->me[4][4]=0.6;printf("A ");
    printf("testing %f \n",A->me[0][0]);
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