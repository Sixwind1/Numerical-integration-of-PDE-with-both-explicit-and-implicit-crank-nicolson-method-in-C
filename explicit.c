#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grRDF.h"

/// utilitzant la funcion exemplar del campus




// u= x^2+xy^2-2t-2tx

double f(double t,double x, double y)
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
}




// problema de la difusio de calor

/*double f(double t,double x, double y)
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
}*/






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
    int nx=5;
    int ny=5;
    int nt=120;


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
        printf("our numerical solution is:\n");
        grRDF_show(gr);
        printf("\n the true solution is:\n");

        for(int j=gr->ny;j>=0;j--)
        {
            for(int i=0;i<gr->nx+1;i++)
            {
                sol[i*(ny+1)+j]=g(k*gr->dt,i*gr->dx,j*gr->dy);
                printf("%.3G ",sol[i*(ny+1)+j]);
            }
            printf("\n");
        }
        printf("\n the absolute error at every position is:\n");
        for(int j=gr->ny;j>=0;j--)
        {
            for(int i=0;i<gr->nx+1;i++)
            {
                printf("%.3G ",fabs(sol[i*(ny+1)+j]-gr->u[i*(gr->ny+1)+j]));
            }
            printf("\n");
        }

        grRDF_escriure(gr,fp);
        grRDF_pasCalExpl(gr,f,g);
    }
    grRDF_allib(gr);
    fclose(fp);
    return 0;
}
