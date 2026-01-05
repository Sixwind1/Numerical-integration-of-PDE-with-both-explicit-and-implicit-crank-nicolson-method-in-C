
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

void grRDF_init (grRDF *gr, double dx, double dy, double dt,
                 int nx, int ny, double (*h)(double x, double y));

void grRDF_pasCalExpl (grRDF *gr,
      double (*f)(double t, double x, double y),
      double (*g)(double t, double x, double y));

void grRDF_escriure (grRDF *gr, FILE *fp);

void grRDF_allib (grRDF *gr);

int grRDF_pasCalCN (grRDF *gr, double w, double tol, int maxit,
      double (*f)(double t, double x, double y),
      double (*g)(double t, double x, double y));

