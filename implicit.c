#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grRDF.h"

/// utilitzant la funcion exemplar del campus




// u= x^2+xy^2-2t-2tx
/*double f(double t,double x, double y)
{
    return -4.-4.*x;
}

double g(double t,double x, double y)
{
    return x*x + x*y*y -2.*t -2.*t*x;
}

double h(double x,double y)
{
    return x*x +x*y*y;
}*/



double f(double t,double x, double y)
{
    if ((x-0.5)*(x-0.5)*(y-0.5)*(y-0.5)<0.2*0.2)
    {
        return 100./0.13;
    }
    return 0.;
}

double g(double t,double x, double y)
{
    return 0.;
}
double h(double x,double y)
{
    return 0.;
}


//P=t(x^3+y^3+x^2y+xy^2+x+y+1)
/*double f(double t, double x, double y)
{
    double u_t = x*x*x + y*y*y + x*x*y + x*y*y + x + y + 1 ;
    double laplace_u = t * (8*x + 8*y);
    return u_t - laplace_u;
}

double g(double t, double x, double y)
{
    return t*(x*x*x+y*y*y+x*x*y+x*y*y+x+y+1.);
}

double h(double x, double y)
{
    return 0.;
}*/















int main()
{
    grRDF *gr=malloc(sizeof(grRDF));

    double Lx=1.; // on acaba x
    double Ly=1.; // on acaba y
    double T=1.;  // on acaba t
    int nx=100;
    int ny=100;
    int nt=100;
    double tol=1e-12;
    int maxit=1000;
    double w=1.7;

    double dx,dy,dt; // magnitud en cada pas
    dx=Lx/nx;
    dy=Ly/ny;
    dt=T/nt;

    double sol[(nx+1)*(ny+1)];
    for(int k=0;k<(nx+1)*(ny+1);k++)
        sol[k]=0.;


    FILE *fp=fopen("dades.txt","w");

    grRDF_init(gr,dx,dy,dt,nx,ny,h);

    for(int k=0;k<nt;k++)
    {

        printf("\n\n");
        printf("in the iteration %i:\n",k);
        //printf("our numerical solution is:\n");
        //grRDF_show(gr);

        grRDF_escriure(gr,fp);
        grRDF_pasCalCN(gr,w,tol,maxit,f,g);
    }
    grRDF_allib(gr);
    fclose(fp);
    return 0;
}
