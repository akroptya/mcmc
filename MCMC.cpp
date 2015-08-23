// MCMC.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>

#define nummod 26
#define STEP_INIT 10

double* FSample = new double[844*3]; 

double* GumbelML( );

bool SteppingOut( double (*f)( double* x ), int n, double* x_0, double y, double w, int m, double* L, double* R );

bool VariablesInTurn( double (*f)( double* x ), int n, double* x_0, double y, double w, int m, double* L, double* R );

double* MultivariateSliceOpt( double (*f)( double* x ), int n );
double derivative( double (*F)( int n, double* variable, double* par ), int n, double* variable, double* par );
double argmaximize( double (*F)(  double variable, int np, double* parameters ), int np, double* parameters );

double Fextr( double threshold, double x, double theta, double beta, int Nu, int n  );
double fextr( double threshold, double x, double theta, double beta, int Nu, int n  );

double Fnorm( double x, double mu, double D );
double fnorm( double x, double mu, double D );

double F( double cross, double threshold, double x, double mu, double D, double theta, double beta, int Nu, int n  );
double f( double cross, double threshold, double x, double mu, double D, double theta, double beta, int Nu, int n  );

double CrossFound( double mu, double D, double threshold, double theta, double beta, int Nu, int n );
double CrossFunction( double variable, int np, double* par );

double FindLeftPoint( double start, double mu, double D, double value );

double* USample( int numobs, double* Sample );

double GumbelLikelihood( double* par );
double GumbelPDF( int n, double* u, double* par );
double GumbelCDF( int n, double* u, double* par );

void FnFill( double* Fn, double mu, double D );
void MultivariateSliceGen( double (*f)( double* x ), int n );
double PDF3( double* x );


double* Fn1 = new double[2000002];
double* Fn2 = new double[2000002]; 
double* Fn3 = new double[2000002]; 

	int n;
	int Nu; 
	double* mu1 =  new double[nummod];
	double* mu2 =  new double[nummod];
	double* mu3 =  new double[nummod];
	double* D1  =  new double[nummod];
	double* D2  =  new double[nummod];
	double* D3  =  new double[nummod];
	double* threshold1  =  new double[nummod];
	double* theta1  =  new double[nummod]; 
	double* beta1   =  new double[nummod];
	double* threshold2  =  new double[nummod];
	double* theta2  =  new double[nummod]; 
	double* beta2   =  new double[nummod];
	double* threshold3  =  new double[nummod];
	double* theta3  =  new double[nummod]; 
	double* beta3   =  new double[nummod];	
	double cross1;
	double cross2;
	double cross3;

	double leftpoint1;
	double leftpoint2;
	double leftpoint3;

	int model;
	
	double* par = new double[2]; // параметри копули

void SetPar( );

void BackTesting()
{
	SetPar( );	

	FILE* fp = fopen("d:\\data\\UAH25\\u.txt","r");
	float x1, x2, x3;
	double* sample = new double[772*3];
	double* curs = sample;
	int j;
	for( j = 0; j < 772; j++ )
	{
		fscanf(fp, "%f \t %f \t %f", &x1, &x2, &x3 ); 
		*(curs++) = x1;
		*(curs++) = x2;
		*(curs++) = x3;
	}
	fclose(fp);

	srand( (unsigned)1 );

	for(model = 15; model < nummod; model++)
	{
		FnFill( Fn1, mu1[model], D1[model] );
		FnFill( Fn2, mu2[model], D2[model] );
		FnFill( Fn3, mu3[model], D3[model] );

		n = 522 + model*10;
		Nu = (int)n*0.03;

		cross1 = CrossFound( mu1[model], D1[model], threshold1[model], theta1[model], beta1[model], Nu, n );
		cross2 = CrossFound( mu2[model], D2[model], threshold2[model], theta2[model], beta2[model], Nu, n );
		cross3 = CrossFound( mu3[model], D3[model], threshold3[model], theta3[model], beta3[model], Nu, n );
		
		double value1 = Fnorm( cross1, mu1[model], D1[model] ) - Fextr( threshold1[model], cross1, theta1[model], beta1[model], Nu, n  );
		double value2 = Fnorm( cross2, mu2[model], D2[model] ) - Fextr( threshold2[model], cross2, theta2[model], beta2[model], Nu, n  ); 
		double value3 = Fnorm( cross3, mu3[model], D3[model] ) - Fextr( threshold3[model], cross3, theta3[model], beta3[model], Nu, n  );

		leftpoint1 = FindLeftPoint( mu1[model], mu1[model], D1[model], value1 );
		leftpoint2 = FindLeftPoint( mu2[model], mu2[model], D2[model], value2 );
		leftpoint3 = FindLeftPoint( mu3[model], mu3[model], D3[model], value3 );

		curs = FSample;
		double* obs = sample;
		for( int i = 0; i < n; i++ )
		{	
			if( *obs < leftpoint1 ){ *obs = leftpoint1; }		
			*curs = F( cross1, threshold1[model], *obs, mu1[model], D1[model], theta1[model], beta1[model], Nu, n  );
			if(*curs < 0.02){ *curs = 0.5; }
			curs++; obs++;

			if( *obs < leftpoint2 ){ *obs = leftpoint2; }
			*curs = F( cross2, threshold2[model], *obs, mu2[model], D2[model], theta2[model], beta2[model], Nu, n  );
			if(*curs < 0.02){ *curs = 0.5; }
			curs++; obs++;

			if( *obs < leftpoint3 ){ *obs = leftpoint3; }
			*curs = F( cross3, threshold3[model], *obs, mu3[model], D3[model], theta3[model], beta3[model], Nu, n  );
			if(*curs < 0.02){ *curs = 0.5; }
			curs++; obs++;
		}		
		
		delete[] par;
		par = GumbelML();

		MultivariateSliceGen( PDF3, 3 );

	}


	delete[] sample;
	delete[] mu1;
	delete[] D1;
	delete[] threshold1;
	delete[] theta1;
	delete[] beta1;
	delete[] mu2;
	delete[] D2;
	delete[] threshold2;
	delete[] theta2;
	delete[] beta2;
	delete[] mu3;
	delete[] D3;
	delete[] threshold3;
	delete[] theta3;
	delete[] beta3;
}



double fnorm1( double x )
{
	return exp( -(x - 0.053066)*(x - 0.053066)/(5208.404) ) / 127.916693;
}
double fnorm2( double x, double mu, double D )
{
	return exp( -(x + .29)*(x + .29)/(89.274) ) / 16.747;
}
double fnorm3( double x, double mu, double D )
{
	return exp( -(x + 0.00141844)*(x + 0.00141844)/(1.178174) ) / 1.92388721;
}



double PDF3( double* x )
{
	double* u = new double[3];
	u[0] = F( cross1, threshold1[model], x[0], mu1[model], D1[model], theta1[model], beta1[model], Nu, n  );	
	u[1] = F( cross2, threshold2[model], x[1], mu2[model], D2[model], theta2[model], beta2[model], Nu, n  );
	u[2] = F( cross3, threshold3[model], x[2], mu3[model], D3[model], theta3[model], beta3[model], Nu, n  );			
	
	double* fm = new double[3];
	fm[0] = f( cross1, threshold1[model], x[0], mu1[model], D1[model], theta1[model], beta1[model], Nu, n  );	
	fm[1] = f( cross2, threshold2[model], x[1], mu2[model], D2[model], theta2[model], beta2[model], Nu, n  );
	fm[2] = f( cross3, threshold3[model], x[2], mu3[model], D3[model], theta3[model], beta3[model], Nu, n  );			

	double res = GumbelPDF( 3, u, par )*fm[0]*fm[1]*fm[2];

	delete[] u;
	delete[] fm;

	return res;
}



double PDF2( double* x )
{
	double* u = new double[3];
	u[0] = F( 133.196648,	64,		x[0], 0.053066,		2604.202,	-0.1574473, 97.12882, 42, 844 );
	u[1] = F( 14.588010,	10.9,	x[1], -0.29,		44.637,		-0.2032381, 5.147675, 42, 844 );
	u[2] = F( 2.574822,		0.2,	x[2], -0.00141844,	0.589087,	 0.8424977, 0.3019789,42, 844 );			
	
	double* fm = new double[3];
	fm[0] = f( 133.196648,	64,		x[0], 0.053066,		2604.202,	-0.1574473, 97.12882, 42, 844 );
	fm[1] = f( 14.588010,	10.9,	x[1], -0.29,		44.637,		-0.2032381, 5.147675, 42, 844 );
	fm[2] = f( 2.574822,	0.2,	x[2], -0.00141844,	0.589087,	 0.8424977, 0.3019789,42, 844 );			

	double* par = new double[2]; 
	par[0] = 1.1966612750633259; 
	par[1] = 1.1971495712149418;

	double res = GumbelPDF( 3, u, par )*fm[0]*fm[1]*fm[2];

	delete[] u;
	delete[] par;
	delete[] fm;

	return res;
}


// multidimensional slice-sampling methods
void MultivariateSliceGen( double (*f)( double* x ), int n )
{	
	int j;
	//srand( (unsigned)time( NULL ) );
	//srand( (unsigned)1 );

	double* x_i = (double*)new double[n];
/*	x_i[0] = 400;
	x_i[1] = 20;
	x_i[2] = 10;*/
	for( j = 0; j < n; j++ )
	{
		x_i[j] = 0;
	}
	
	double f_i = f(x_i);
	double fmax = 0;	
	double y;
	double* L = (double*)new double[n];
	L[0] = leftpoint1;
	L[1] = leftpoint2;
	L[2] = leftpoint3; 

	double* R = (double*)new double[n];
	R[0] = 410;
	R[1] = 25;
	R[2] = 10; 
	double w = 10;
	int count = 0;


	double test;

//	time_t ltime, ltinit;
	//time( &ltinit );
	//double t1 = Fnorm( 5000, 0.053066, 2604.202); 

	FILE* fp2 = fopen("d:\\data\\vyborka.txt","w");

	for( int i=0; i < 10000; i++ )
	{
		y = ((double)rand()/(double)RAND_MAX)*f_i;
		
		// нахождение интервала содержащего slice
		// VariablesInTurn( f, n, x_i, y, w, 40, L, R );
		// SteppingOut( f, n, x_i, y, w, 10, L, R );
		//L[0] = x_i[0]-200;
		//L[1] = x_i[1]-10;
		//L[2] = x_i[2]-5; 
/*
		if( L[0] < -107.735300 )
		{L[0] = -107.735300;}
		if( L[1] < -16.091800 )
		{L[1] = -16.091800;}
		if( L[2] < -2.612000 )
		{L[2] = -2.612000; }
*/
	/*	R[0] = x_i[0]+400;
		R[1] = x_i[1]+20;
		R[2] = x_i[2]+10; 
	*/	
		if( R[0] > 420)
		{R[0] = 420;}
		if( R[1] > 25)
		{R[1] = 25;}
		if( R[2] > 10)
		{R[2] = 10;} 

		bool flag = true;
		do
		{
			// sampling from interval
			for( j=0; j < n; j++ )
			{
				x_i[j] = ((double)rand()/(double)RAND_MAX)*(R[j] - L[j]) + L[j];
			}			
			f_i = f(x_i);
			
			/*if( y < f_i )
			{
				// shrinkage
				for( j=0; j < n; j++ )
				{
					if( x_i[j] < x_0[j] )
					{LR[j*2] = x_1[j];	}
					else
					{LR[j*2+1] = x_1[j]; }
				}
			}*/

			if(_isnan(f_i))
			{
				test = x_i[0];
			}
			if(x_i[0] > 200)
			{
				test = x_i[0];
			}

		}
		while( (y > f_i) || (_isnan(f_i)) || (f_i > 1) );
		
		for( j = 0; j < n; j++ )
		{			
			fprintf(fp2, "%f\t", x_i[j] );
			printf("%f\t", x_i[j] );
		}		
		fprintf(fp2, "\n" );
		//time( &ltime );
		//printf("\t %ld\n", ltime-ltinit );
		printf("%d\n",count++ );
	}  
	fclose(fp2);

	printf("%f\n",test);
	delete[] L;
	delete[] R;
	delete[] x_i;
}


-- multidimensional optimization
double* MultivariateSliceOpt( double (*f)( double* x ), int n )
{	
	int j;
	//srand( (unsigned)time( NULL ) );
//	srand( (unsigned)1 );


	double* x_i = (double*)new double[n];
	double* argmax = (double*)new double[n];
	for( j = 0; j < n; j++ )
	{
		x_i[j] = 1;
		argmax[j] = -100000; 
	}
	
	double f_i = f(x_i);
	double fmax = 0;	
	double y;
	double* L = (double*)new double[n];
	double* R = (double*)new double[n];
	double w = 1;

	for( int i=0; i < 100; i++ )
	{
		y = ((double)rand()/(double)RAND_MAX)*f_i;
		
		// нахождение интервала содержащего slice
		// VariablesInTurn( f, n, x_i, y, w, 40, L, R );
		SteppingOut( f, n, x_i, y, w, 10, L, R );
		
		bool flag = true;
		do
		{
			// sampling from interval
			for( j=0; j < n; j++ )
			{
				x_i[j] = ((double)rand()/(double)RAND_MAX)*(R[j] - L[j]) + L[j];
			}			
			f_i = f(x_i);
			/*if( y < f_i )
			{
				// shrinkage
				for( j=0; j < n; j++ )
				{
					if( x_i[j] < x_0[j] )
					{LR[j*2] = x_1[j];	}
					else
					{LR[j*2+1] = x_1[j]; }
				}
			}*/
		}
		while( y > f_i );
		
		if( f_i > fmax )
		{
			for( j = 0; j < n; j++ )
			{
				argmax[j] = x_i[j]; 
			}
			fmax = f_i;
		}		
	}  

	return argmax;
}


double* GumbelMultivariateSliceOpt( double (*f)( double* x ), int n )
{	
	int j;
	//srand( (unsigned)time( NULL ) );
	srand( (unsigned)1 );

	double* L = (double*)new double[n];
	double* R = (double*)new double[n];

	double* x_i = (double*)new double[n];
	double* argmax = (double*)new double[n];

	for( j = 0; j < n; j++ )
	{
		x_i[j] = 1;
		argmax[j] = -100000; 
		L[j] = 1;
		R[j] = 3;
	}
	
	double f_i = f(x_i);
	double fmax = 0;	
	double y;
	double w = 1;


	for( int i=0; i < 10000; i++ )
	{
		y = ((double)rand()/(double)RAND_MAX)*f_i;
		
		bool flag = true;
		do
		{
			do
			{
				for( j=0; j < n; j++ )
				{
					x_i[j] = ((double)rand()/(double)RAND_MAX)*(R[j] - L[j]) + L[j];
				}			
			}
			while( x_i[0] > x_i[1]);

			f_i = f(x_i);
		}
		while( (y > f_i) || (_isnan(f_i)) );
		
		if( f_i > fmax )
		{
			for( j = 0; j < n; j++ )
			{
				argmax[j] = x_i[j]; 
			}
			fmax = f_i;
		}		
	}  

	delete[] L;
	delete[] R;
	delete[] x_i;

	return argmax;
}


bool VariablesInTurn( double (*f)( double* x ), int n, double* x_0, double y, double w, int m, double* L, double* R )
{
	double* J = (double*)new double[n];
	double* K = (double*)new double[n];

//	srand( (unsigned)time( NULL ) );
//	srand( (unsigned)2 );
	
	int j;
	for( j=0; j < n; j++ )
	{
		L[j] = x_0[j] - w*((double)rand()/(double)RAND_MAX);
		R[j] = L[j] + w;

		J[j] = floor(m*((double)rand()/(double)RAND_MAX));
		K[j] = (m-1) - J[j];

		while( (J[j] > 0) && (f(L) > y) )
		{
			L[j] -= w;
			J[j]--;
		}
		
		while( (K[j] > 0) && (f(R) > y) )
		{
			R[j] += w;
			K[j]--;
		}			
	}

	delete[] J;
	delete[] K;

	return true;
}

bool SteppingOut( double (*f)( double* x ), int n, double* x_0, double y, double w, int m, double* L, double* R )
{
	double* J = (double*)new double[n];
	double* K = (double*)new double[n];
	double* x_1 = (double*)new double[n];
	//srand( (unsigned)time( NULL ) );
	//srand( (unsigned)2 );
	
	int j;
	for( j=0; j < n; j++ )
	{
		L[j] = x_0[j] - w*((double)rand()/(double)RAND_MAX);
		R[j] = L[j] + w;

		J[j] = floor(m*((double)rand()/(double)RAND_MAX));
		K[j] = (m-1) - J[j];
	}

	bool flag;
	while( flag )
	{
		flag = false;
		for( j=0; j < n; j++ )
		{
			if( J[j] != 0 )
			{
				L[j] -= w;
				J[j]--;
				flag = true;
			}
			
			if( K[j] != 0 )
			{
				R[j] += w;
				K[j]--;
				flag = true;
			}
		}

		// Сгенерировать точку
		for( j=0; j < n; j++ )
		{
			x_1[j] = ((double)rand()/(double)RAND_MAX)*(R[j] - L[j]) + L[j];
		}
		if( y > f(x_1))
		{		
			for( j=0; j < n; j++ )
			{
				if( x_1[j] < x_0[j] )
				{J[j] = 0; }
				else
				{K[j] = 0; }
			}		
		}
	}
	
	delete[] J;
	delete[] K;
	delete[] x_1;

	return true;
}


double Copula(double* u) 
{
//	return u[0]*u[1];
/*	if(u[0]+u[1]-1 > 0) 
		{return u[0]+u[1]-1;}
	else 
		{return 0;}*/
	if(u[0]> u[1]) 
		{return u[1];}
	else 
		{return u[0];}
}


int main(int argc, char* argv[])
{  
	int n = 2;
	double* x_i = (double*)new double[n];
	x_i[0]=0; x_i[1]=0; 

	double f_i = Copula(x_i);
	double y;

	FILE* fp2 = fopen("d:\\temp\\vyborka.txt","w");

	srand( (unsigned)2 );
		
	for( int i=0; i < 10000; i++ )
	{
		y = ((double)rand()/(double)RAND_MAX)*f_i;
		fprintf(fp2,"%f,%f\n",x_i[0],x_i[1]);
		do
		{
			for(int j=0; j < n; j++ )
			{
				x_i[j] = ((double)rand()/(double)RAND_MAX);
			}			
			f_i = Copula(x_i);			
		}
		while( y > f_i );		
	}
	fclose(fp2);
	return 0;
}



double GumbelCDF( int n, double* u, double* par )
{
	return exp( -pow( pow( pow(log((1.0/u[0])),par[1]) + pow(log((1.0/u[1])),par[1]) , (par[0]/par[1])) + pow(log((1.0/u[2])),par[0]) , (1.0/par[0]) ) );
}

double GumbelPDF( int n, double* u, double* par )
{
	return derivative( GumbelCDF, n, u, par );
}


double* GumbelML( )
{
	double* theta = new double[2];
	theta = GumbelMultivariateSliceOpt( GumbelLikelihood, 2 );

	return theta;
}



double GumbelLikelihood( double* par ) 
{ 
	double gl = 0;
	double* obs = FSample;
	
	for(int i = 0; i < n; i++ )
	{			
		gl += log(GumbelPDF( 3, obs, par ));
		if( _isnan(gl))
		{
			printf("Hello world!");
		}
		obs += 3;
	}
	
	return gl+1000;
}




// numerical differention of multivariable function
// n - number of variables 
// variable - array of varibales
double derivative( double (*F)( int n, double* variable, double* par ), int n, double* variable, double* par )
{
	double res = 0;
	double* x_h = (double*)new double[n];
	
	_int32 signs = 0;
	_int32 sign = 1;
	unsigned _int32 members = 1;
	double h_n = 1;
	members = (members << n);

	int up = 0; 
	int down = 0;
	
	int k;
	double step = 10;
	for( k = 0; k < n; k++ )
	{
		if( variable[k] < step )
		{
			step = variable[k];
		}
		if( (1-variable[k]) < step )
		{
			step = 1-variable[k];
		}

	}

	step *= 0.1;

	h_n = 1;
	for( k = 0; k < n; k++ )
	{
		h_n = h_n *step; 
	}
	double res1;
	for(int i = 0; i < members; i++)
	{
		sign = 1;
		for(int j = 0; j < n; j++)
		{	
			signs = 1;
			signs = 1 - 2*(((signs << j) & i) >> j);

			x_h[j] = variable[j] + signs*step; 			
			sign *= signs;
		}
		res1 = F(n, x_h, par );
		res += sign*F(n, x_h, par );		
	}

	res = fabs(res / (members * h_n));

	delete[] x_h;

	return res; 
}


double Fextr( double threshold, double x, double theta, double beta, int Nu, int n  )
{
	double Fextr_x;
	double Fu;
	if( theta != 0 )
	{	Fu = 1 - pow( 1 + theta*(x - threshold)/beta, -1/theta ); }
	else
	{	Fu = 1 - exp(-(x - threshold)/beta);	}
		
	Fextr_x = ((double)Nu/(double)n)*Fu + ((double)(n - Nu))/(double)n;
	
	return Fextr_x;
}


double fextr( double threshold, double x, double theta, double beta, int Nu, int n  )
{
	double fextr_x;
	double fu;
	if( theta != 0 )
	{	fu = (1/beta)*pow( 1 + theta*(x - threshold)/beta, -(1 + 1/theta) ); }
	else
	{	fu = (1/beta)*exp(-(x - threshold)/beta);	}
		
	fextr_x = ((double)Nu/(double)n)*fu;

	return fextr_x;
}

double fnorm( double x, double mu, double D )
{
	return exp( -(x - mu)*(x - mu)/(2*D) ) / sqrt(2*3.141592*D);
}

double Fnorm( double x, double mu, double D )
{
	double* Fn = Fn1;
	if( mu == Fn2[0] ){
		Fn = Fn2;		
	}
	else { if(mu == Fn3[0]){ Fn = Fn3;} }
	
	double res = Fn[(int)((x+10*D-mu)*1000000/fabs(10*D))];
	return res;	
}

void FnFill( double* Fn, double mu, double D )
{
	Fn[0] = mu;
	Fn[1] = D;
	Fn[2] = 0;

	double Fnorm_x = 0;
	double t = -10*D + mu;
	double h = fabs(10*D) / (double)1000000;
    for( int i = 3; t < 10*D + mu; i++ )
	{
        Fn[i] = Fn[i-1] + fnorm( t, mu, D )*h;
        t += h;
	}
}



double F( double cross, double threshold, double x, double mu, double D, double theta, double beta, int Nu, int n  )
{
	double Fx = 0;
	if( x >= cross )
	{
		Fx = Fextr( threshold, x, theta, beta, Nu, n );
	}
	else
	{
		Fx = Fnorm( x, mu, D );
	}

	return Fx;
}

double f( double cross, double threshold, double x, double mu, double D, double theta, double beta, int Nu, int n  )
{
	double fx = 0;
	if( x >= cross )
	{
		fx = fextr( threshold, x, theta, beta, Nu, n );
	}
	else
	{
		fx = fnorm( x, mu, D );
	}

	return fx;
}


// maximization single-variable function by MCMC
// most probable values
double argmaximize( double (*F)( double variable, int np, double* parameters ), int np, double* parameters )
{
	srand( (unsigned)time( NULL ) );

	double x_i = 0;
	double f_i = F(x_i, np, parameters );
	double fmax = 0;
	double argmax = 0;
	double y;
	double w_s, w_f, w;
	w = 1;

	for( int i=0; i < 100000; i++ )
	{
		y = ((double)rand()/(double)RAND_MAX)*f_i;
		
		// нахождение интервала содержащего slice
		// stepping out procedure
		w_s = x_i-(w/2);
		w_f = x_i+(w/2);

		while( F(w_s, np, parameters ) > y )
		{
			w_s -= w;
		}
		
		while( F(w_f, np, parameters ) > y )
		{
			w_f += w;
		}		
		
		do
		{
			x_i = ((double)rand()/(double)RAND_MAX)*(w_f-w_s) + w_s;
			f_i = F(x_i, np, parameters );
		}
		while( f_i < y );
		
		if( f_i > fmax)
		{
			argmax = x_i;
			fmax = f_i;
		}		
	}  

	return argmax;
}

double CrossFound( double mu, double D, double threshold, double theta, double beta, int Nu, int n )
{
	double* par = new double[7];
	par[0] = mu;	
	par[1] = D;
	par[2] = threshold; 
	par[3] = theta; 
	par[4] = beta; 
	par[5] = Nu; 
	par[6] = n;

	double x = threshold;
	double h = 0.000001;

	double f_p = CrossFunction( x, 7, par ); 
	x += h;
	double f_n = CrossFunction( x, 7, par );
	
	while( f_n > 0 )
	{
		f_n = CrossFunction( x, 7, par );
		x += h;
		f_p = f_n;
	}

	return x-h;
}


double CrossFunction( double variable, int np, double* par )
{
	return fnorm(variable, par[0], par[1]) - fextr( par[2], variable, par[3], par[4], par[5], par[6] );
}


double FindLeftPoint( double start, double mu, double D, double value )
{
	double x = start;	
	double h = 0.0001;
	
	double f_n = Fnorm( x, mu, D );
	x -= h;

	while( f_n > value )
	{
		f_n = Fnorm( x, mu, D );
		x -= h;
	}

	return x-h;
}

double* USample( int numobs, double* Sample )
{
	double* US = new double[numobs*3];
	double* curs = US;
	double* obs = Sample;
	
	for( int i =0; i < numobs; i++ )
	{		
		*curs = F( 133.196648, 64, *obs, 0.053066, 2604.202, -0.1574473, 97.12882, 42, 844 );
		curs++; obs++;
		*curs = F( 14.588010, 10.9, *obs, -0.29, 44.637, -0.2032381, 5.147675, 42, 844 );
		curs++; obs++;
		*curs = F( 2.574822, 0.2, *obs, -0.00141844, 0.589087, 0.8424977, 0.3019789, 42, 844 );
		curs++; obs++;
	}

	return US;
}



// main
	/*
	FILE* fp2 = fopen("d:\\data\\gl2.txt","w");
	double* theta;
	theta = GumbelML();
	fprintf(fp2, "%f\n", theta[0]);
	fprintf(fp2, "%f\n", theta[1]);

	*/
/*
	FILE* fp2 = fopen("d:\\data\\gl2.txt","w");
	double* par = new double[2];
	par[0] = 1.01;
	par[1] = 2;
	

	double* obs = FSample;
	double g, gcdf;
	double L = 0;
	int j;

/*	for( j = 0; j < 828; j++ )
	{
		obs += 3;
		}*/
/*	for( j = 0; j < 843; j++ )
	{
		int  decimal, sign;
		char *buffer;

		gcdf = GumbelPDF( 3, obs, par ); 
		//g = GumbelPDF( 3, obs, par ); 
		
		
		//buffer = _fcvt( g, 16, &decimal, &sign );
		//fprintf(fp2, "%s \n", buffer );
		obs += 3;
		buffer1 = _fcvt( gcdf, 16, &decimal, &sign );
		fprintf(fp2, "%f \n", gcdf );	
		//fprintf(fp2, ".%s \n", buffer1 );			
	}

	for( j = 0; j < 1; j++ )
	{		
		while( par[0] < par[1] )
		{			
			L = 0;
			obs = FSample;
			for( i = 0; i < 843; i++ )
			{
				g = GumbelPDF( 3, obs, par ); 
				L += log(g);
				fprintf(fp2, "%f\n", g);
				fprintf(fp2, "%f\n", L);
				obs += 3;
			}
			fprintf(fp2, "%f\n", L);	
			par[0] += 0.01;	
		}
		par[1] += 0.1;
	}
	fprintf(fp2, "\n");
	
	
	
	/*
	
	double* par = new double[2];
	par[0] = 1;
	par[1] = 1;

	for( int j = 0; j < 100000; j++ )
	{
		par[1] += 0.1;
		while( par[0] < par[1] )
		{
			par[0] += 0.01;	
			double g = GumbelLikelihood( par ); 
			if( g > 0 )
			{
				fprintf(fp2, "%f\n", g);	
			}
		}
	}
	*/
/*	double* obs = FSample;
	for( int j = 0; j < 848; j++ )
	{
		fprintf(fp2, "%f\n", GumbelPDF( 3, obs, par ));
		obs += 3;
	}
fclose(fp2);
	*/
	
	


	/*
	FILE* fp = fopen("d:\\data\\UEC.txt","r");

	float x1, x2, x3;
	double* sample = new double[848*3];
	double* curs = sample;
	int i;
	for( i =0; i<848; i++ )
	{
		fscanf(fp, "%f \t %f \t %f", &x1, &x2, &x3 ); 
		*(curs++) = x1;
		*(curs++) = x2;
		*(curs++) = x3;
	}

	fclose(fp);

	curs = USample( 848, sample );

	fp = fopen("d:\\data\\out.txt","w");
	for( i =0; i < 848; i++ )
	{
		fprintf(fp, "%f \t %f \t %f\n", *curs, *(curs+1), *(curs+2) ); 
		curs += 3;
	}
	fclose(fp);

	delete[] sample;
*/
	/*
	printf("%f\n", 1-Fextr( 64, 133.196648, -0.1574473, 97.12882, 42, 848 ));
	printf("%f\n", 1-Fextr( 10.9, 14.588010, -0.2032381, 5.147675, 42, 848 ));
	printf("%f\n", 1-Fextr( 0.2, 2.574822, 0.8424977, 0.3019789, 42, 848 ));
*/
/*
	printf("%f\n", FindLeftPoint( -107, 0.053066, 2604.202, 0.0178783413 ));
	printf("%f\n", FindLeftPoint( -16, -0.29, 44.637, 0.0090264438 ));
	printf("%f\n", FindLeftPoint( -2, -0.00141844, 0.589087, 0.000335319 ));
*/
	/*
	printf("UAH\n");
	printf("%f\n", 1-Fnorm( 91.8, 0.053066, 2604.202 ));
	printf("%f\n", 1-Fextr( 64, 91.8, -0.1574473, 97.12882, 42, 848 ));
*/
//	printf("%f\n",CrossFound( 0.053066, 2604.202, 64, -0.1574473, 97.12882, 42, 848 ));
//	printf("%f\n",CrossFound( -0.29, 44.637, 10.9, -0.2032381, 5.147675, 42, 848 ));
//	printf("%f\n",CrossFound( -0.00141844, 0.589087, 0.2, 0.8424977, 0.3019789, 42, 848 ));
/*
	printf("%f\n", fnorm( 133.196700, 0.053066, 2604.202 ));
	printf("%f\n", fextr( 64, 133.196700, -0.1574473, 97.12882, 42, 848 ));

	printf("EUR\n");
	printf("%f\n", 1-Fnorm( 20, -0.29, 44.637 ));
	
	printf("CNY\n");
	printf("%f\n", 1-Fnorm( 2, -0.00141844, 0.589087 ));
	*/

/*
	printf("%f\n", argmaximize());
	
	printf("%f\n", argmax[0]);
	printf("%f", argmax[1]);
*/
	/*
	int n = 2;
	double* variable = (double*)new double[n];
	variable[0] = 5;
	variable[1] = 5;	

	printf("%f", derivative( n, variable ));
*/
	/*
	FILE* fp = fopen("d:\\norm.txt","w");

	srand( (unsigned)time( NULL ) );

	double x_i = 0;
	double y;
	double w_s, w_f, w;
	w = 0.1;

	for( int i=0; i < 101000; i++ )
	{
		y = ((double)rand()/(double)RAND_MAX)*PDF1(x_i);
		
		// нахождение интервала содержащего slice
		w_s = x_i-(w/2);
		w_f = x_i+(w/2);
		while( PDF1(w_s) > y )
		{
			w_s -= w;
		}
		
		while( PDF1(w_f) > y )
		{
			w_f += w;
		}		
		
		do
		{
			x_i = ((double)rand()/(double)RAND_MAX)*(w_f-w_s) + w_s;
		}
		while( PDF1(x_i) < y );
		
		if( i > 100000 )
			fprintf(fp,"%6f\n", x_i );	    
	}  

	fclose(fp);*/

double PDF1( double x )
{
	// -1*x*x+4
/*	if( abs(x) > 2 )
	{
		return 0;
	}
	else
	{
		return -1*x*x+4;
	}*/

	//return (1/sqrt(2*3.14))*exp(-(x*x)/2);

/*
	if( abs(x) > 15 )
	{
		return 0;
	}
	else
	{
		return x*sin(x)*cos(x)+10;
	}*/
	return (1/sqrt(2*3.14))*exp(-(x*x)/2);
}

/*
double F( int n, double* variable )
{
	double res = 2*variable[0]*variable[0]*variable[1] + 3*variable[1]*variable[1]*variable[0];
	//double res = variable[0]*variable[1]/100;
	//double res = 3*variable[0];

	return res;
}
*/



void SetPar( )
{
//u0
mu1[0] = 0.132183908;	
mu2[0] = -0.209578544;	
mu3[0] = -0.002298851;

D1[0] = 2515.933895;	
D2[0] = 49.89533994;	
D3[0] = 0.832739427;

//u1
mu1[1] = 0.132894737;	
mu2[1] = -0.214473684;	
mu3[1] = 0.003195489;

D1[1] = 2478.080215;	
D2[1] = 49.598792;	
D3[1] = 0.836505777;

//u2
mu1[2] = 0.131549815;	
mu2[2] = -0.253874539;	
mu3[2] = -0.00295203;

D1[2] = 2434.310297;	
D2[2] = 48.92943967;	
D3[2] = 0.837994966;

//u3
mu1[3] = 0.134057971;	
mu2[3] = -0.235507246;	
mu3[3] = -0.001449275;

D1[3] = 2400.254773;	
D2[3] = 48.41300192;	
D3[3] = 0.891867224;

//u4
mu1[4] = 0.122241993;	
mu2[4] = -0.236476868;	
mu3[4] = -0.001423488;

D1[4] = 2365.772178;
D2[4] = 48.09946921;	
D3[4] = 0.87596945;

//u5
mu1[5] = 0.111538462;	
mu2[5] = -0.272552448;	
mu3[5] = -0.002272727;

D1[5] = 2328.453317;	
D2[5] = 47.54048875;	
D3[5] = 0.860782917;

//u6
mu1[6] = 0.208762887;	
mu2[6] = -0.21185567;	
mu3[6] = -0.002061856;

D1[6] = 2305.411885	;
D2[6] = 47.64032392	;
D3[6] = 0.846191955;

//u7
mu1[7] = 0.23597973	;
mu2[7] = -0.216047297;	
mu3[7] = -0.002368866;

D1[7] = 2278.551207;	
D2[7] = 47.46493663;	
D3[7] = 0.833824887;

//u8
mu1[8] = 0.117607973;	
mu2[8] = -0.249335548;	
mu3[8] = -0.002159468;

D1[8] = 2250.355064;
D2[8] = 47.06406778;	
D3[8] = 0.81904691;

//u9
mu1[9] = 0.089379085;	
mu2[9] = -0.305228758;	
mu3[9] = -0.002124183;

D1[9] = 2317.772309;	
D2[9] = 46.78141288;
D3[9] = 0.805805628;

//u10
mu1[10] = 0.088424437;	mu2[10] = -0.33392283;	mu3[10] = -0.002411576;
D1[10] = 2636.708143;	D2[10] = 46.84160101;	D3[10] = 0.793085962;

//u11
mu1[11] = 0.089082278;	mu2[11] = -0.33164557;	mu3[11] = -0.00221519;
D1[11] = 2954.975237;	D2[11] = 46.85154848;	D3[11] = 0.781104435;

//u12
mu1[12] = 0.073208723;	mu2[12] = -0.332087227;	mu3[12] = -0.002336449;
D1[12] = 2909.035849;	D2[12] = 46.83778316;	D3[12] = 0.769401709;

//u13
mu1[13] = 0.079601227;	mu2[13] = -0.301993865;	mu3[13] = -0.001226994;
D1[13] = 2927.296373;	D2[13] = 46.6180759;	D3[13] = 0.757847955;

//u14
mu1[14] = 0.096827795;	mu2[14] = -0.252416918;	mu3[14] = -0.001661631;
D1[14] = 2889.53002;	D2[14] = 46.3693058;	D3[14] = 0.746759716;

//u15
mu1[15] = 0.080654762;	mu2[15] = -0.279464286;	mu3[15] = -0.001041667;
D1[15] = 2851.700519;	D2[15] = 46.04819166;	D3[15] = 0.736139003;

//u16
mu1[16] = 0.08372434;	mu2[16] = -0.257624633;	mu3[16] = -0.001906158;
D1[16] = 2810.214228;	D2[16] = 45.826792;	D3[16] = 0.725943498;

//u17
mu1[17] = 0.246820809;	mu2[17] = -0.205202312;	mu3[17] = -0.00216763;
D1[17] = 2786.747182;	D2[17] = 45.71869938;	D3[17] = 0.715668232;

//u18
mu1[18] = 0.088746439;	mu2[18] = -0.244301994;	mu3[18] = -0.001851852;
D1[18] = 2815.493867;	D2[18] = 45.99211442;	D3[18] = 0.705831088;

//u19
mu1[19] = 0.08511236;	mu2[19] = -0.264325843;	mu3[19] = -0.001404494;
D1[19] = 2803.869539;	D2[19] = 45.87653147;	D3[19] = 0.696088039;

//u20
mu1[20] = 0.076731302;	mu2[20] = -0.293905817;	mu3[20] = -0.001800554;
D1[20] = 2765.235574;	D2[20] = 45.7603789;	D3[20] = 0.686779;

//u21
mu1[21] = 0.076775956;	mu2[21] = -0.278278689;	mu3[21] = -0.002459016;
D1[21] = 2727.601457;	D2[21] = 45.42785859;	D3[21] = 0.677504205;

//u22
mu1[22] = 0.078167116;	mu2[22] = -0.258086253;	mu3[22] = -0.002156334;
D1[22] = 2690.979685;	D2[22] = 45.23752562;	D3[22] = 0.668497368;

//u23
mu1[23] = 0.19481383;	mu2[23] = -0.299202128;	mu3[23] = -0.00212766;
D1[23] = 2755.423049;	D2[23] = 45.20492613;	D3[23] = 0.659755787;

//u24
mu1[24] = 0.2;	mu2[24] = -0.306824147;	mu3[24] = -0.00144357;
D1[24] = 2752.530565;	D2[24] = 45.00828451;	D3[24] = 0.651338255;

//u25
mu1[25] = 0.057901554;	mu2[25] = -0.324222798;	mu3[25] = -0.001554404;
D1[25] = 2730.375203;	D2[25] = 44.55003506;	D3[25] = 0.6430326;



theta1[0] = -0.4427414;
beta1[0] = 154.9238;
theta2[0] = -0.1228940;
beta2[0] = 3.357390;
theta3[0] = 0.7284898;
beta3[0] = 0.5598324;

theta1[1] = -0.4427414;
beta1[1] = 154.9238;
theta2[1] = -0.1228940;
beta2[1] = 3.357390;
theta3[1] = 0.8481006;
beta3[1] = 0.5188907;

theta1[2] = -0.4256245;
beta1[2] = 154.5698;
theta2[2] = -0.2675161;
beta2[2] = 4.219264;
theta3[2] = 0.7177012;
beta3[2] = 0.6333643;

theta1[3] = -0.4256245;
beta1[3] = 154.5698;
theta2[3] = -0.2675161;
beta2[3] = 4.219264;
theta3[3] = 1.001071;
beta3[3] = 0.5081818;

theta1[4] = -0.4256245;
beta1[4] = 154.5698;
theta2[4] = -0.2675161;
beta2[4] = 4.219264;
theta3[4] = 1.001071;
beta3[4] = 0.5081818;

theta1[5] = -0.3512079;
beta1[5] = 139.2945;
theta2[5] = -0.4557341;
beta2[5] = 5.776231;
theta3[5] = 0.7818423;
beta3[5] = 0.674986;

theta1[6] = -0.3512079;
beta1[6] = 139.2945;
theta2[6] = -0.2333952;
beta2[6] = 3.975059;
theta3[6] = 0.7818423;
beta3[6] = 0.674986;

theta1[7] = -0.3512079; 
beta1[7] = 139.2945;
theta2[7] = -0.2333952;
beta2[7] = 3.975059;
theta3[7] = 0.7818423;
beta3[7] = 0.674986;

theta1[8] = -0.3207777;
beta1[8] = 134.5586;
theta2[8] = -0.4245923;
beta2[8] = 5.508326;
theta3[8] = 0.6479136;
beta3[8] = 0.8061088;

theta1[9] = -0.3302577;
beta1[9] = 133.3287;
theta2[9] = -0.4245923;
beta2[9] = 5.508326;
theta3[9] = 0.6479136;
beta3[9] = 0.8061088;

theta1[10] = -0.4328279;
beta1[10] = 148.1962;
theta2[10] = -0.4245923;
beta2[10] = 5.508326;
theta3[10] = 0.6479136;
beta3[10] = 0.8061088;

theta1[11] = -0.2624649;
beta1[11] = 102.2657;
theta2[11] = -0.4245923;
beta2[11] = 5.508326;
theta3[11] = 0.6479136;
beta3[11] = 0.8061088;

theta1[12] = -0.4063643;
beta1[12] = 133.1133;
theta2[12] = -0.3581708;
beta2[12] = 5.124195;
theta3[12] = 0.7528904;
beta3[12] = 0.5086749;

theta1[13] = -0.4063643;
beta1[13] = 133.1133;
theta2[13] = -0.3581708;
beta2[13] = 5.124195;
theta3[13] = 0.7528904;
beta3[13] = 0.5086749;

theta1[14] = -0.4063643;
beta1[14] = 133.1133;
theta2[14] = -0.3581708;
beta2[14] = 5.124195;
theta3[14] = 0.7528904;
beta3[14] = 0.5086749;

theta1[15] = -0.6518546; 
beta1[15] = 194.8907;
theta2[15] = -0.3581708;
beta2[15] = 5.124195;
theta3[15] = 0.7528904;
beta3[15] = 0.5086749;

theta1[16] = -0.6518546;
beta1[16] = 194.8907;
theta2[16] = -0.3581708;
beta2[16] = 5.124195;
theta3[16] = 0.7528904;
beta3[16] = 0.5086749;

theta1[17] = -0.3992525;
beta1[17] = 134.0689;
theta2[17] = -0.3375199;
beta2[17] = 4.939455;
theta3[17] = 0.7528904;
beta3[17] = 0.5086749;

theta1[18] = -0.4464179;
beta1[18] = 148.6064;
theta2[18] = -0.3375199;
beta2[18] = 4.939455;
theta3[18] = 0.7528904;
beta3[18] = 0.5086749;

theta1[19] = -0.4229385;
beta1[19] = 142.4821; 
theta2[19] = -0.3375199; 
beta2[19] = 4.939455; 
theta3[19] = 0.7528904;
beta3[19] = 0.5086749; 

theta1[20] = -0.4229385;
beta1[20] = 142.4821; 
theta2[20] = -0.3375199;
beta2[20] = 4.939455; 
theta3[20] = 0.7528904;
beta3[20] = 0.5086749;

theta1[21] = -0.4229385;
beta1[21] = 142.4821;
theta2[21] = -0.3375199;
beta2[21] = 4.939455;
theta3[21] = 0.7528904;
beta3[21] = 0.5086749;

theta1[22] = -0.3897278;
beta1[22] = 136.6802;
theta2[22] = -0.3574722;
beta2[22] = 5.191864;
theta3[22] = 0.7528904;
beta3[22] = 0.5086749;

theta1[23] = -0.2376162;
beta1[23] = 102.5908;
theta2[23] = -0.3574722;
beta2[23] = 5.191864;
theta3[23] = 0.7528904;
beta3[23] = 0.5086749; 

theta1[24] = -0.1715222;
beta1[24] = 91.29588;
theta2[24] = -0.3574722;
beta2[24] = 5.191864;
theta3[24] = 0.7528904;
beta3[24] = 0.5086749;

theta1[25] = -0.1855790;
beta1[25] = 94.2092;
theta2[25] = -0.6208403;
beta2[25] = 7.914634; 
theta3[25] = 0.7528904;
beta3[25] = 0.5086749;


threshold1[0] = 96.7;	threshold2[0] = 15.7;	threshold3[0] = 0.6;
threshold1[1] = 96.7;	threshold2[1] = 15.7;	threshold3[1] = 0.7;
threshold1[2] = 92.4;	threshold2[2] = 15.4;	threshold3[2] = 0.6;
threshold1[3] = 92.4;	threshold2[3] = 15.4;	threshold3[3] = 0.8;
threshold1[4] = 92.4;	threshold2[4] = 15.4;	threshold3[4] = 0.8;
threshold1[5] = 84.6;	threshold2[5] = 14.9;	threshold3[5] = 0.7;
threshold1[6] = 84.6;	threshold2[6] = 15.4;	threshold3[6] = 0.7;
threshold1[7] = 84.6;	threshold2[7] = 15.4;	threshold3[7] = 0.7;
threshold1[8] = 83.5;	threshold2[8] = 14.9;	threshold3[8] = 0.6;
threshold1[9] = 92.4;	threshold2[9] = 14.9;	threshold3[9] = 0.6;
threshold1[10] = 114.3;	threshold2[10] = 14.9;	threshold3[10] = 0.6;
threshold1[11] = 140.9;	threshold2[11] = 14.9;	threshold3[11] = 0.6;
threshold1[12] = 132.2;	threshold2[12] = 14.1;	threshold3[12] = 0.5;
threshold1[13] = 132.2;	threshold2[13] = 14.1;	threshold3[13] = 0.5;
threshold1[14] = 132.2;	threshold2[14] = 14.1;	threshold3[14] = 0.5;
threshold1[15] = 114.3;	threshold2[15] = 14.1;	threshold3[15] = 0.5;
threshold1[16] = 114.3;	threshold2[16] = 14.1;	threshold3[16] = 0.5;
threshold1[17] = 114.3;	threshold2[17] = 14.1;	threshold3[17] = 0.5;
threshold1[18] = 108.4;	threshold2[18] = 14.1;	threshold3[18] = 0.5;
threshold1[19] = 108.4;	threshold2[19] = 14.1;	threshold3[19] = 0.5;
threshold1[20] = 108.4;	threshold2[20] = 14.1;	threshold3[20] = 0.5;
threshold1[21] = 108.4;	threshold2[21] = 14.1;	threshold3[21] = 0.5;
threshold1[22] = 99.4;	threshold2[22] = 13.8;	threshold3[22] = 0.5;
threshold1[23] = 119.2;	threshold2[23] = 13.8;	threshold3[23] = 0.5;
threshold1[24] = 119.4;	threshold2[24] = 13.8;	threshold3[24] = 0.5;
threshold1[25] = 119.2;	threshold2[25] = 13.5;	threshold3[25] = 0.5;

}