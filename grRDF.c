#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define IDX(i,j,colsize) (i*colsize+j)

typedef struct {
   double dx;
   double dy;
   double mux;
   double muy;
   int nx;
   int ny;
   double t;
   double dt;
   double *u;
} grRDF;

void grRDF_show(grRDF *gr)
{
    printf("\n dx =%.3G    dy =%.3G    mux =%.3G    muy =%.3G \n",gr->dx,gr->dy,gr->mux,gr->muy);
    printf("\n nx =%i    ny =%i   t =%.3G      dt =%.3G \n",gr->nx,gr->ny,gr->t,gr->dt);
    printf("\n the layer u is shown below in the analogous Euclidean space ordering: \n");

    int colsize=gr->ny+1;
    for(int j=gr->ny;j>=0;j--)
    {
        for(int i=0;i<gr->nx+1;i++)
        {
            printf("%.3G ",gr->u[IDX(i,j, colsize)]);
        }
        printf("\n");
    }
}

void grRDF_init (grRDF *gr, double dx, double dy, double dt,
                 int nx, int ny, double (*h)(double x, double y))
{
    /// we assume that gr already points to an adequate memory block
    /// this function fills all the constant data of gr and sets u to the initial conditions (time 0)
    gr->dx=dx;
    gr->dy=dy;
    gr->t=0.;
    gr->dt=dt;
    gr->mux= gr->dt/ ((gr->dx)*(gr->dx));
    gr->muy= gr->dt/ ((gr->dy)*(gr->dy));
    gr->nx=nx; //  # 0,1,...,nx for a total of nx+1 nodes
    gr->ny=ny;
    gr->u=malloc(sizeof(double)*(nx+1)*(ny+1));
    int colsize=ny+1;
    for(int j=gr->ny;j>=0;j--)
    {
        for(int i=0;i<gr->nx+1;i++)
        {
            gr->u[ IDX(i,j,colsize)]=h(i*dx,j*dy);
        }
    }
    printf("in grRDF init: successfully initialized\n");
}

void grRDF_pasCalExpl (grRDF *gr,
      double (*f)(double t, double x, double y),
      double (*g)(double t, double x, double y)) // run one step and update gr
{
    int rowsize=gr->nx+1;
    int colsize=gr->ny+1;
    double tpdt=gr->t + gr->dt; // tpdt = t + dt
    double u1[rowsize * colsize]; // auxiliary grid to store the next iteration

    // fill boundary conditions
    for(int j=0;j<colsize;j++)
    {
        u1[IDX(0,j,colsize)]=g(tpdt,  0.,   j*gr->dy);
        u1[IDX(gr->nx,j,colsize)]=g(tpdt,  (gr->nx)*gr->dx,  j*gr->dy);
    }
    for(int i=0;i<rowsize;i++)
    {
        u1[IDX(i,0,colsize)]= g(tpdt,  i*gr->dx,  0.);
        u1[IDX(i,gr->ny,colsize)]= g(tpdt,  i*gr->dx,  (gr->ny)*gr->dy);
    }

    // apply algorithm 1.1 to fill the interior of the grid
    for(int i=1;i<rowsize-1;i++)
    {
        for(int j=1;j<colsize-1;j++)
        {
            u1[IDX(i,j,colsize)] = (1.-2.*(gr->mux+gr->muy))*gr->u[IDX(i,j,colsize)];
            u1[IDX(i,j,colsize)]+= gr->mux*(gr->u[IDX((i+1),j,colsize)]+gr->u[IDX((i-1),j,colsize)]);
            u1[IDX(i,j,colsize)]+= gr->muy*(gr->u[IDX(i,(j+1),colsize)]+gr->u[IDX(i,(j-1),colsize)]);
            u1[IDX(i,j,colsize)]+= gr->dt * f(gr->dt,gr->dx*i,gr->dy*j);
        }
    }
    memcpy(gr->u,u1,sizeof(double)*colsize*rowsize); // copy matrix u1 into u
    gr->t += gr->dt; // update time
}

void grRDF_escriure (grRDF *gr, FILE *fp) // copy the data of gr and write it to fp
{
    // NOTE: we assume that fp IS ALREADY OPENED IN WRITE MODE!!!

    /* to draw a single frame we ONLY need the set of points and not their positions (which are already implicit in the data)
    # t time1       the data are distributed by columns, starting from the leftmost column
    x1,y1,u11
    x1,y2,u12
    x1,y2,u12
    .
    .
    .
    x1,ym,u1m
    x2,y1,u21
    x2,y2,u22
    .
    .
    .
    x2,ym,u2m
    .
    .
    .
    .
    .
    xn,ym,unm


    # t time2
    etc...*/

    int colsize=gr->ny+1;
    fprintf(fp,"# t %.3G\n",gr->t);
    for(int i=0;i<gr->nx+1;i++) //
    {
        for(int j=0;j<gr->ny+1;j++)
        {
            fprintf(fp,"%.3G %.3G %.3G \n",i*gr->dx,j*gr->dy,gr->u[IDX(i,j,colsize)]);
        }

    } // when exiting we are on a new line
    fprintf(fp,"\n \n"); // two extra line breaks to start writing again
    //printf("in grRDF write: data successfully written to fp\n");
}

void grRDF_allib (grRDF *gr) // completely free a gr
{
    free(gr->u);
    free(gr);
}

void grRDF_pasCalCN (grRDF *gr, double w, double tol, int maxit,
      double (*f)(double t, double x, double y),
      double (*g)(double t, double x, double y)) // for part 2 of Crank–Nicolson
{
    int rowsize=gr->nx+1;
    int colsize=gr->ny+1;
    double tpdt=gr->t + gr->dt; // tpdt = t + dt
    double UM[rowsize][colsize];
    double UM1[rowsize][colsize];
    double mx=gr->mux;
    double my=gr->muy;

    // fill boundary conditions in UM (and copy them into UM1)
    for (int j=0;j<colsize;j++)
    {
        UM[0][j]=g(tpdt,0.,j*gr->dy);
        UM[gr->nx][j] = g(tpdt,gr->nx*gr->dx,j*gr->dy);
        UM1[0][j]=UM[0][j];
        UM1[gr->nx][j]=UM[gr->nx][j];
    }

    for (int i=0; i<rowsize;i++)
    {
        UM[i][0]=g(tpdt,i*gr->dx,0.);
        UM[i][gr->ny]=g(tpdt,i*gr->dx,gr->ny*gr->dy);
        UM1[i][0]=UM[i][0];
        UM1[i][gr->ny]=UM[i][gr->ny];
    }
    for(int i=1;i<rowsize-1;i++) // interior of UM and UM1 = interior of gr->u
    {
        for(int j=1;j<colsize-1;j++)
        {
            UM[i][j]=gr->u[IDX(i,j,colsize)];
            UM1[i][j]=gr->u[IDX(i,j,colsize)];
        }
    }

    // apply iteration (6) of the pdf
    int it=0;
    double err=0.;
    double aux=0.; // the part inside brackets of the algorithm
    while (it<maxit)
    {
        it++;
        //printf("in grRDF_pasCalCN: while loop, starting iteration %i\n",it);
        err=0.;
        for(int i=1;i<rowsize-1;i++)
        {
            for(int j=1;j<colsize-1;j++)
            {
                aux=(1.-mx-my)*gr->u[IDX(i,j,colsize)];
                aux+= mx/2. * (UM[i+1][j]+UM1[i-1][j]+gr->u[IDX((i+1),j,colsize)]+gr->u[IDX((i-1),j,colsize)]);
                aux+= my/2. * (UM[i][j+1]+UM1[i][j-1]+gr->u[IDX(i,(j+1),colsize)]+gr->u[IDX(i,(j-1),colsize)]);
                aux+= gr->dt/2. * (f(tpdt,i*gr->dx,j*gr->dy) + f(gr->t,i*gr->dx,j*gr->dy));
                UM1[i][j]= (1.-w)*UM[i][j]+ w*aux/ (1.+mx+my);

                if(err<fabs(UM[i][j]-UM1[i][j])) /// err = infinity norm (UM - UM1)
                    err=fabs(UM[i][j]-UM1[i][j]);
            }
        }

        /*printf("in grRDF_pasCalCN: while loop, the difference matrix UM1-UM is:\n");
        for(int j=colsize-1;j>=0;j--)
        {
            for(int i=0;i<rowsize;i++)
            {
                printf("%.3G ",UM1[i][j]-UM[i][j]);
            }
            printf("\n");
        }*/

        for (int i=1;i<rowsize-1;i++) // UM = UM1  (indices: UM and UM1 only differ in the interior)
        {
            for (int j=1;j<colsize-1;j++)
            {
                UM[i][j] = UM1[i][j];
            }
        }

        if (err<tol)
            break;
    }
    printf(" \n in pasCalCN, goal achieved in %i SOR iterations and the maximal norm error in this last iteration is %.8G \n",it,err);
    if(it>=maxit)
    {
        printf("\n in pasCalCN: error, reached maxit\n");
        return;
    }

    for (int i=0;i<rowsize;i++) // gr->u = UM1
    {
        for (int j=0;j<colsize;j++)
        {
            gr->u[IDX(i,j,colsize)] = UM1[i][j];
        }
    }
    gr->t=tpdt;
}
