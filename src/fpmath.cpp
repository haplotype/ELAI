#include "fpmath.h"
#include "string.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#if defined (MPI_ENABLED)
#include "mpi.h"
#endif

using namespace std;

double sum(double * ta, int n)
{
	double res = 0.0; 
	for (int i = 0; i < n; i++)
		res += ta[i]; 
	return res; 
}

int sample_pmass(int n, double *pmass)
{
	int res = n-1; 
	double prob = gsl_rng_uniform(gsl_r);
	double cur_prob = 0;
	for (int j = 0; j < n; j++)
	{
		cur_prob += pmass[j];
		if(prob < cur_prob) 
		{
			res = j; 
			break;
		}
	}       	 
	return (res); 
}

double sumlog(double * ta, int n)
{
	double dmax = ta[0]; 
	for (int i = 1; i < n; i++)
	{
        if(dmax < ta[i]) 
			dmax = ta[i]; 
	}

	double dmul = 0; 
	for (int i = 0; i < n; i++)
	{
		double diff = ta[i] - dmax; 
        dmul += exp(diff); 
	}
	return(dmax + log(dmul)); 
}

int geometric(int len, int meanrank)
{
	if(meanrank > len) 
		meanrank = len;
	double tmin = exp(-(double) len / (double) meanrank); 
	double x = gsl_rng_uniform(gsl_r) * (1.0 - tmin) + tmin; 
	return ((int) floor(-log(x) * meanrank));  
}

double pgeometric(int pick, int len, int meanrank)
{
	if(meanrank > len) 
		meanrank = len;   
	double p = 1.0 - exp(-1.0 / (double) meanrank);
	double q = 1.0 - p;

	double temp = len * log(q); 
	if (temp < -100) 
		temp = 0; 
	else 
		temp = exp(temp); 
	double logp = pick * log(q) + log(p) - log(1.0 - temp); 
	return(logp); 
}

string getline(streambuf * pbuf)
{
    char ch;
	string str;
	size_t pos; 
	while((ch = pbuf->sgetc()) != EOF)
	{
		if(ch != '\n' && ch != '\r')
		{
			str.push_back(ch);
			ch = pbuf->snextc();
		}
		else {
			pbuf->sbumpc();  //chomp;
			pos = str.find_first_not_of(";, \t", 0); 
			if(str.empty() || pos  == string::npos || str.at(pos) == '#')
			{
				str.clear();
				continue;
			}
			else
				break;
		}
	}
	return str;
}   //this getline use ;, \t as delimit, and ignore the lines either full of delimit or starting with #. 

int imap(int nK, int x, int y)
{
	int tmax = x; 
	int tmin = y; 
	if(x < y) {
		tmax = y; tmin = x; 
	}
	int res = ((tmax*(tmax+1)) >> 1) + tmin;  
//	int res1 = (tmin) * nK - (tmin * tmin - tmin)/2; 
//	res1 += (tmax - tmin);  
//	if(res1 != res) {cout << "worng in imap " << x << " " << y << " || " <<  res1 << " " << res << endl; } 
	return (res); 
}

void safe_exit() 
{
#if defined (MPI_ENABLED)
	MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Finalize(); 
#endif
	exit(0); 
}

int compare_pair(const void * a, const void * b)
{
	Pair * ta = (Pair *) a;
	Pair * tb = (Pair *) b; 
	if((tb->bf) > (ta->bf)) return 1;
	else if((tb->bf) == (ta->bf)) return 0;
	else return -1;
}

int compare(const void * a, const void * b)
{
	int *ta = (int *) a;
	int *tb = (int *) b;
	if ((*ta) > (*tb)) return 1;
	else if ((*ta) == (*tb)) return 0;
	else return -1;
}

void * Allocate1D(size_t us, int dim)
{
	size_t size = dim * us; 
	char * m = (char *) malloc(size);
	memset(m, 0, size); 
	return (void *) m; 
}

void Free1D(void * m)
{
	if(m == NULL) return; 
	free(m); 
	m = NULL; 
}

void ** Allocate2D(size_t us,  int dim1, int dim2)
{
	char ** m;
	m = (char **) malloc((size_t)(dim1 * us));
	m[0] = (char *) malloc((size_t)(dim1 *dim2 * us));
   	memset(m[0], 0, (size_t) (dim1 * dim2 * us));  
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D double matrix. \n");
		safe_exit();
	}
	for(int i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + (size_t) (dim2 * us);
	}
	return ((void **) (m));
}

void Free2D(void ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}

double ** Allocate2DMatrix(int dim1, int dim2)
{
	int i;
	double ** m;
	
	m = (double **) malloc((size_t)((dim1)*sizeof(double*)));
	m[0] = (double *) malloc((size_t)((dim1*dim2)*sizeof(double)));
	memset(m[0], 0, (dim1*dim2)*sizeof(double)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D double matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}

//
//void *** Allocate3DMatrix(size_t us, int dim1, int dim2, int dim3)
//{
//	void *** m; 
//	m = (void ***) malloc((size_t)(dim1 * us));
//	m[0] = (void **) malloc((size_t)(dim1 * dim2 * us));
//	m[0][0] = (void *) malloc((size_t)(dim1 * dim2 * dim3 * us));
//	if (!(m && m[0] &&  m[0][0]))
//	{
//		printf("Error: Problem allocating a 3D double matrix. \n");
//		safe_exit();
//	}
//
//	for (int j = 1; j < dim2; j++)
//	{
//		m[0][j] = m[0][j-1] + dim3;
//	}
//	for (int i = 1; i < dim1; i++)
//	{
//		m[i] = m[i-1] + dim2;
//		m[i][0] = m[i-1][dim2 - 1] + dim3;
//		for(int j = 1; j < dim2; j++)
//		{
//			m[i][j] = m[i][j-1] + dim3;
//		}
//	}
//	return (m);
//}
//
//void Free3DMatrix(void *** m)
//{
//	if(m == NULL) return;
//	free(m[0][0]);
//	free(m[0]);
//	free(m);
//	m = NULL; 
//}

double *** Allocate3DMatrix(int dim1, int dim2, int dim3)
{
	int i, j;
	double *** m; 

	m = (double ***) malloc((size_t)((dim1)*sizeof(double**)));
	m[0] = (double **) malloc((size_t)((dim1) * (dim2) * sizeof(double *)));
	m[0][0] = (double *) malloc((size_t)((dim1) * (dim2) * (dim3) * sizeof(double)));
	memset(m[0][0], 0, (dim1) * (dim2) * (dim3) * sizeof(double));  
	if (!(m && m[0] &&  m[0][0]))
	{
		printf("Error: Problem allocating a 3D double matrix. \n");
		safe_exit();
	}

	for (j = 1; j < dim2; j++)
	{
		m[0][j] = m[0][j-1] + dim3;
	}
	for (i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
		m[i][0] = m[i-1][dim2 - 1] + dim3;
		for(j = 1; j < dim2; j++)
		{
			m[i][j] = m[i][j-1] + dim3;
		}
	}
	return (m);
}

void Free2DMatrix(double ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}

void Free3DMatrix(double *** m)
{
	if(m == NULL) return;
	free(m[0][0]);
	free(m[0]);
	free(m);
	m = NULL; 
}
 

int ** Allocate2DIntMatrix(int dim1, int dim2)
{
	int i;
	int ** m;
	
	m = (int **) malloc((size_t)((dim1)*sizeof(int*)));
	m[0] = (int *) malloc((size_t)((dim1*dim2)*sizeof(int)));
	memset(m[0], 0, (dim1*dim2)*sizeof(int)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D int matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}

short int ** Allocate2DShortMatrix(int dim1, int dim2)
{
	int i;
	short int ** m;
	
	m = (short int **) malloc((size_t)((dim1)*sizeof(short int*)));
	m[0] = (short int *) malloc((size_t)((dim1*dim2)*sizeof(short int)));
	memset(m[0], 0, (dim1*dim2)*sizeof(short)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D int matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}    

void Free2DShortMatrix(short int ** m)
{
	if(m == NULL) return;
	free(m[0]);
	free(m);
	m = NULL; 
}

void Free2DIntMatrix(int ** m)
{                        
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}  

char  ** Allocate2DCharMatrix(int dim1, int dim2)
{
	int i;
	char ** m;
	
	m = (char **) malloc((size_t)((dim1)*sizeof(char*)));
	m[0] = (char *) malloc((size_t)((dim1*dim2)*sizeof(char)));
	memset(m[0], 0, (dim1*dim2)*sizeof(char)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D char matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}    


void Free2DCharMatrix(char ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}  

void DistinctIntArray(int low, int high, int dim, int * A)
{
    if (dim >= high - low)
		return;  
	
	int i, j, r; 
	int bingle; 
	for (i=0; i<dim; i++)
		A[i] = -1;
	
	int howmany = 0;
	for (i=high - dim; i<high; i++)
	{
		bingle = 0;
		r = low + gsl_rng_uniform_int(gsl_r, i-low);
		for (j = 0; j < howmany; j++)
		{
			if (r == A[j])
			{
				bingle = 1; 
				break;
			}
		}

		if (bingle) A[howmany] = i;
		else A[howmany] = r; 
		howmany++;
	}
}

//void lu_decomp(double ** a, int n, int *indx, int *d)
//{
//	int i, j, k;
//	int imax = -1;
//	double big, dum, sum, temp;
//	double *vv = new double[n];
//	*d = 1;
//
//	for (i = 0; i < n; i++)
//	{
//		big = 0.0; 
//		for (j = 0; j < n; j++)
//			if((temp = fabs(a[i][j])) > big) big = temp;
//		if (big == 0.0)
//		{
//			cout << "singular matrix in routine ludcmp" << endl;
//			exit(0); 
//		}
//		vv[i] = 1.0 / big; 
//	}
//
//	for (j = 0; j < n; j++)
//	{
//		for (i = 0; i < j; i++)
//		{
//			sum = a[i][j];
//			for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
//			a[i][j] = sum;
//		}
//
//		big = 0.0;
//		for (i = j; i < n; i++)
//		{
//			sum = a[i][j];
//			for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
//			a[i][j] = sum;
//			if ((dum = vv[i] * fabs(sum)) >= big)
//			{
//				big = dum; 
//				imax = i; 
//			}
//		}
//
//		if ( j != imax)
//		{
//			for (k = 0; k < n; k++)
//			{
//				dum = a[imax][k];
//				a[imax][k] = a[j][k];
//				a[j][k] = dum; 
//			}
//			*d = -(*d); 
//			vv[imax] = vv[j];
//		}
//		indx[j] = imax;
//		if(a[j][j] == 0.0) a[j][j] = TINY;
//		if (j != n)
//		{
//			dum = 1.0 / a[j][j];
//			for (i=j+1; i < n; i++) a[i][j] *= dum; 
//		}
//	}
//	delete[] vv; 
//}
//
//void lu_back_sub(double ** a, int n, int *indx, double b[])
//{
//	int i, ii = -1, ip, j; 
//	double sum;
//
//	for (i = 0; i < n; i++)
//	{
//		ip = indx[i];
//		sum = b[ip];
//		b[ip] = b[i];
//		if(ii>=0)
//			for (j=ii;j <i;j++) sum -= a[i][j] * b[j];
//		else if(sum) ii = i; 
//		b[i] = sum; 
//	}
//
//	for (i = n-1; i >= 0; i--) 
//	{
//		sum = b[i];
//		for (j=i+1; j < n; j++) sum -= a[i][j] * b[j];
//		b[i] = sum/a[i][i];
//	}
//} 

void center(double * xx, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += xx[i];
	sum /= (double) n;
	for (int i = 0; i < n; i++)
		xx[i] -= sum;
}

void center_col(double ** xx, int nrow, int ncol)
{
    for (int j = 0; j< ncol; j++)
    {
        double sum = 0;
        for (int i = 0; i < nrow; i++)
            sum += xx[i][j];
        sum /= (double) nrow;
        for (int i = 0; i < nrow; i++)
            xx[i][j] -= sum;
    }

}

int nchoosek(int n, int k) 
{
  	if (k <= 0 || k >= n)
  		return 1;
  	else
  		return nchoosek(n-1, k-1) + nchoosek(n-1, k);
}

int normalize(double * v, int n)
{
	double total = 0; 
	int neg = 0; 
	for (int i = 0; i < n; i++)
	{
		if(v[i] < 0) { neg = 1; break;}
		total += v[i]; 
	}
	if(neg || total < 1e-100) 
	{
		for (int i = 0; i < n; i++)
			v[i] = 0; 
	}
	for (int i = 0; i < n; i++)
		v[i] /= total; 
	return 1; 
}

int find_max_index(double *v, int n)
{
	int index = 0; 
	double mval = -1e100; 
	for (int i = 0; i < n; i++)
	{
        if(v[i] > mval) 
		{
			mval = v[i];
			index = i; 
		}
	}
	return index; 
}

double find_max_value(double *v, int n)
{
	double mval = -1e100; 
	for (int i = 0; i < n; i++)
	{
        if(v[i] > mval) 
			mval = v[i];
	}
	return mval; 
}

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

/******************************************************************************/
double local_min ( double a, double b, double eps, double t, 
  double fd(double x, void * par) , double * x , void * par)


/******************************************************************************/
/*
  Purpose:

    LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].

  Discussion:

    The method used is a combination of golden section search and
    successive parabolic interpolation.  Convergence is never much slower
    than that for a Fibonacci search.  If F has a continuous second
    derivative which is positive at the minimum (which is not at A or
    B), then convergence is superlinear, and usually of the order of
    about 1.324....

    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
    F is never evaluated at two points closer than TOL.  

    If F is a unimodal function and the computed values of F are always
    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
    LOCAL_MIN approximates the abscissa of the global minimum of F on the 
    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.  

    If F is not unimodal, then LOCAL_MIN may approximate a local, but 
    perhaps non-global, minimum to the same accuracy.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    Orignal FORTRAN77 version by Richard Brent.
    C version by John Burkardt.

  Reference:

    Richard Brent,
    Algorithms for Minimization Without Derivatives,
    Dover, 2002,
    ISBN: 0-486-41998-3,
    LC: QA402.5.B74.

  Parameters:

    Input, double A, B, the endpoints of the interval.

    Input, double EPS, a positive relative error tolerance.
    EPS should be no smaller than twice the relative machine precision,
    and preferably not much less than the square root of the relative
    machine precision.

    Input, double T, a positive absolute error tolerance.

    Input, double F ( double x ), a user-supplied
    function whose local minimum is being sought.

    Output, double *X, the estimated value of an abscissa
    for which F attains a local minimum value in [A,B].

    Output, double LOCAL_MIN, the value F(X).
*/
{
  double c;
  double d = 0;
  double e;
  double fu;
  double fv;
  double fw;
  double fx;
  double m;
  double p;
  double q;
  double r;
  double sa;
  double sb;
  double t2;
  double tol;
  double u;
  double v;
  double w;
  int iter = 0; 
/*
  C is the square of the inverse of the golden ratio.
*/
  c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

  sa = a;
  sb = b;
  *x = sa + c * ( b - a );
  w = *x;
  v = w;
  e = 0.0;
  fx = (*fd) ( *x, par);
  fw = fx;
  fv = fw;

  for ( ; ; )
  { 
	  iter++; 
	  if(iter > 200) 
		  break; 

    m = 0.5 * ( sa + sb ) ;
    tol = eps * fabs ( *x ) + t;
    t2 = 2.0 * tol;
/*
  Check the stopping criterion.
*/
    if ( fabs ( *x - m ) <= t2 - 0.5 * ( sb - sa ) )
    {
      break;
    }
/*
  Fit a parabola.
*/
    r = 0.0;
    q = r;
    p = q;

    if ( tol < fabs ( e ) )
    {
      r = ( *x - w ) * ( fx - fv );
      q = ( *x - v ) * ( fx - fw );
      p = ( *x - v ) * q - ( *x - w ) * r;
      q = 2.0 * ( q - r );
      if ( 0.0 < q )
      {
        p = - p;
      }
      q = fabs ( q );
      r = e;
      e = d;
    }

    if ( fabs ( p ) < fabs ( 0.5 * q * r ) && 
         q * ( sa - *x ) < p && 
         p < q * ( sb - *x ) )
    {
/*
  Take the parabolic interpolation step.
*/
      d = p / q;
      u = *x + d;
/*
  F must not be evaluated too close to A or B.
*/
      if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
      {
        if ( *x < m )
        {
          d = tol;
        }
        else
        {
          d = - tol;
        }
      }
    }
/*
  A golden-section step.
*/
    else
    {
      if ( *x < m )
      {
        e = sb - *x;
      }
      else
      {
        e = a - *x;
      }
      d = c * e;
    }
/*
  F must not be evaluated too close to X.
*/
    if ( tol <= fabs ( d ) )
    {
      u = *x + d;
    }
    else if ( 0.0 < d )
    {
      u = *x + tol;
    }
    else
    {
      u = *x - tol;
    }

    fu = (*fd) ( u, par );
/*
  Update A, B, V, W, and X.
*/
    if ( fu <= fx )
    {
      if ( u < *x )
      {
        sb = *x;
      }
      else
      {
        sa = *x;
      }
      v = w;
      fv = fw;
      w = *x;
      fw = fx;
      *x = u;
      fx = fu;
    }
    else
    {
      if ( u < *x )
      {
        sa = u;
      }
      else
      {
        sb = u;
      }

      if ( fu <= fw || w == *x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == *x || v== w )
      {
        v = u;
        fv = fu;
      }
    }
  }
  return fx;
}
