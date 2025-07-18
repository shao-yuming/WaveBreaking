#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  
#include "view.h"
#include "adapt_wavelet_leave_interface.h"
#include "tag.h"
#include "output_surfaces.h"
#include "output_vtu_foreach.h"
#include "liutex_omega.h"

// Simulation parameters
// ak: wave steepness, BO: Bond number, RE: Reynolds number
double ak = 0.55;
double BO = 1000.;
double RE = 40000.;
int LEVEL = dimension == 2 ? 9 : 9;
double uemax = 0.001;// Max velocity error for AMR

// Fluid properties (density and viscosity ratios)
#define RATIO (1./998.)
#define MURATIO (17.4e-6/8.9e-4)
int DIRAC = 0;

// Physical constants
#define k_  (2.*pi)      // Wave number
#define h_   0.5         // Water depth
#define g_   1           // Gravitational acceleration
#define sig (1./(BO*sq(k_)))  // Surface tension coefficient

// Main entry point
int main (int argc, char * argv[])
{
  // Read command-line parameters
  if (argc > 1) LEVEL = atoi (argv[1]);
  if (argc > 2) ak = atof(argv[2]);
  if (argc > 3) BO = atof(argv[3]);
  if (argc > 4) RE = atof(argv[4]);    
  if (argc > 5) DIRAC = atof(argv[5]);

  // Set simulation domain and periodicity
  origin (-L0/2, -L0/2, -L0/2);
  periodic (right);
#if dimension > 2
  periodic (front);
#endif

  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; 
  mu2 = 1.0/RE*MURATIO;
  f.sigma = 1./(BO*sq(k_));
  G.y = -g_;

#if TREE  
  N = 32;
#else
  N = 1 << LEVEL;
#endif
  run();
}

// Initial free surface profile (third-order Stokes wave)
double wave (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
}

double eta (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3;
}

// Gaussian weighting function for vorticity initialization
double gaus (double y, double yc, double T){
  double deltaw = sqrt(2.0/RE)/k_;
  double deltaa = sqrt(2.0/RE*MURATIO/RATIO)/k_;
  double r = y - yc;
  return 2.0/(sqrt(2.0*pi*sq(deltaa)) + sqrt(2.0*pi*sq(deltaw))) *
    (T*exp(-sq(r)/(2.0*sq(deltaw))) + (1.0 - T)*exp(-sq(r)/(2.0*sq(deltaa))));
}

// Initialization event: sets free surface and velocity field
event init (i = 0)
{
  if (!restore ("restart")) {
    do {
      fraction (f, wave(x,y));

      scalar Phi[];
      foreach() {
	double alpa = 1./tanh(k_*h_);
	double a_ = ak/k_;
	double sgma = sqrt(g_*k_*tanh(k_*h_)*
			   (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					      (sq(alpa) - 1.) + sq(alpa))));
	double A_ = a_*g_/sgma;
	double phi1 = A_*cosh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x);
	double phi2 = 3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
	  cosh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_);
	double phi3 = 1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
	  (9.*sq(alpa) - 13.)*
	  cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
	Phi[] = phi1 + ak*phi2 + ak*ak*phi3;
      }
      boundary ({Phi});
     
     // Velocity field initialization (Dirac-based or potential-based)
      if (DIRAC) {
	scalar vort2[];
	scalar psi[];
	foreach() {
	  vort2[] = -2.0*gaus(y,wave(x,y)+y,f[])*(Phi[1,0]-Phi[-1,0])/(2.*Delta);
	  psi[] = 0.0;
	}
	boundary ({vort2,psi});
	psi[top] = dirichlet(0.);
	psi[bottom] = dirichlet(0.);
	poisson (psi, vort2);
	foreach() {
	  u.x[] = (psi[0,1] - psi[0,-1])/(2.*Delta);
	  u.y[] = -(psi[1] - psi[-1])/(2.*Delta);
	}
      }
      else {
	foreach()
	  foreach_dimension()
	    u.x[] = (Phi[1] - Phi[-1])/(2.0*Delta) * f[];
      }
      boundary ((scalar *){u});
    }
#if TREE  

   while (adapt_wavelet_leave_interface ((scalar *){u},{f}, (double[]){uemax,uemax,uemax}, LEVEL, 6, 2).nf);

#else
    while (0);

#endif
  }
}

// Save checkpoint fields at regular time intervals
event dumpforwave (t = 0; t <= 10.5; t += 0.1) {
   char dname[100];
   sprintf (dname, "field/dump-%g", t);
   dump (file = dname);
 }

// Adaptive mesh refinement event
event adapt (i++) {
  adapt_wavelet_leave_interface ((scalar *){u},{f}, (double[]){uemax,uemax,uemax}, LEVEL, 6, 2);
  }

// Export free-surface field for ParaView (VTP)
void save_surface(int nf, int j)
{
    char name[80];
    FILE *fp;
    nf > 0 ? sprintf(name, "surface/surface_%6.6d_n%3.3d.vtp", j, pid()) : sprintf(name, "surface/surface_%6.6d.vtp", j);
    fp = fopen(name, "w");
    output_vtp(f, fp, (vector *){u});
    fclose(fp);

    char subname[80];
    sprintf(subname, "surface/surface_%6.6d.pvtp", j);
    char base_name[80];
    sprintf(base_name, "surface_%6.6d", j);
    fp = fopen(subname, "w");
    output_pvtp((vector *){u}, fp, base_name);
    fclose(fp);
}

// Export full-field variables (VTU) for ParaView, including Liutex vortex field
void save_vtu(int nf, int j)
{
  vector liutex[];
  scalar omegar[];
  omega_liutex (u,liutex,omegar);

    FILE *fp;
    char name[111];
    nf > 0 ? sprintf(name, "vtu/all_%6.6d_n%3.3d.vtu", j, pid()) : sprintf(name, "vtu/all_%6.6d.vtu", j);
    fp = fopen(name, "w");
    output_vtu_pid((scalar *){f, p,omegar}, (vector *){u,liutex}, t, fp, false);
    fclose(fp);

    char subname[112];
    sprintf(subname, "vtu/all_%6.6d.pvtu", j);
    char base_name[112];
    sprintf(base_name, "all_%6.6d", j);
    fp = fopen(subname, "w");
   output_pvtu((scalar *){f, p,omegar}, (vector *){u,liutex}, t, fp, base_name);
    fclose(fp);
}

// Regularly save surface and volume data
event savedata(t = 0; t <= 10.5; t +=0.001)
{
    static int j = 0;
    save_surface(1, j);
    save_vtu(1, j);
    j++;
    
}

// Compute viscous dissipation rates in water and air
int dissipation_rate (double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			     sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    rateWater += mu1/rho[]*f[]*sqterm; //water
    rateAir   += mu2/rho[]*(1. - f[])*sqterm; //air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}

// Output energy budget (kinetic, potential, surface energy, and dissipation)
event graphs (t = 0; t <= 11; t +=0.01) {
  static FILE * fpwater = fopen("budgetWater.dat", "w");
  static FILE * fpair = fopen("budgetAir.dat", "w");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  double se=0;
  foreach(reduction(+:ke) reduction(+:gpe) 
	  reduction(+:keAir) reduction(+:gpeAir) reduction(+:se)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += rho[]*norm2*f[]*dv();
    keAir += rho[]*norm2*(1.0-f[])*dv();
    gpe += rho1*g_*y*f[]*dv();
    gpeAir += rho2*g_*y*(1.0-f[])*dv();

    for (scalar s in interfaces)
      if (s[] > 1e-7 && s[] < 1 - 1e-7) {
        coord  normal = mycs(point, s), parea;
        double alpha = plane_alpha(s[], normal);
        double dS = pow(Delta, dimension - 1) * plane_area_center(normal, alpha, &parea);
        se += sig*dS;
      }

  }
  double rates[2];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
    if (i == 0) {
    fprintf (fpwater, "t ke gpe se dissipation\n");
    fprintf (fpair, "t ke gpe dissipation\n");
    }
  fprintf (fpwater, "%g %g %g %g %g\n",
	   t/(k_/sqrt(g_*k_)), ke/2., gpe + 0.125, se, dissWater);
fflush(fpwater);
  fprintf (fpair, "%g %g %g %g\n",
	   t/(k_/sqrt(g_*k_)), keAir/2., gpeAir + 0.125, dissAir);
fflush(fpair);

}

// Identify and log droplet and bubble statistics (volume, centroid, velocity)
event countDropsBubble(t = 0; t <= 10.5; t +=0.001)
{
  scalar m1[]; //droplets
  scalar m2[]; //bubbles
  foreach(){
    m1[] = f[] > 0.5; //i.e. set m true if f[] is close to unity (droplets)
    m2[] = f[] < 0.5; //m true if f[] close to zero (bubbles)
  }
  int n1 = tag(m1);
  int n2 = tag(m2);

  double v1[n1]; //droplet
  coord b1[n1];  //droplet
  double ux1[n1];
  double uy1[n1];
  double uz1[n1];
  double v2[n2]; //bubble
  coord b2[n2];  //bubble
  double ux2[n2];
  double uy2[n2];
  double uz2[n2];

  for (int j=0; j<n1; j++)
     // v1[j] = b1[j].x = b1[j].y = b1[j].z = 0.0;
      v1[j] = b1[j].x = b1[j].y = b1[j].z =ux1[j]=uy1[j]=uz1[j]= 0.0;
  for (int j=0; j<n2; j++)
      //v2[j] = b2[j].x = b2[j].y = b2[j].z = 0.0;
      v2[j] = b2[j].x = b2[j].y = b2[j].z =ux2[j]=uy2[j]=uz2[j]= 0.0;

  foreach_leaf() {
    // droplets
    if (m1[] > 0) {
      int j = m1[] - 1;
      v1[j] += dv()*f[]; //increment the volume of the droplet
      coord p = {x,y,z};

      ux1[j] += dv()*u.x[]*f[];
      uy1[j]+= dv()*u.y[]*f[];
      uz1[j] += dv()*u.z[]*f[];

      foreach_dimension()
	      b1[j].x += dv()*f[]*p.x;

    }
    // bubbles
    if (m2[] > 0) {
      int j = m2[] - 1;
      v2[j] += dv()*(1.0-f[]);
      coord p = {x,y,z};

       ux2[j] += dv()*u.x[]*(1.0-f[]);
  	   uy2[j] += dv()*u.y[]*(1.0-f[]);
       uz2[j] += dv()*u.z[]*(1.0-f[]);

      foreach_dimension()
	     b2[j].x += dv()*(1.0-f[])*p.x;
       
    }
  }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b1, 3*n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (MPI_IN_PLACE, ux1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uy1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uz1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (MPI_IN_PLACE, v2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (MPI_IN_PLACE, ux2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uy2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uz2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif

  static FILE * fdrop = fopen("droplets.dat","w");
  static FILE * fbubb = fopen("bubbles.dat","w");
  for (int j=0; j<n1; j++)
    fprintf (fdrop, "%d %g %d %g %g %g %g %g %g %g\n", i, t,
	     j, v1[j], b1[j].x/v1[j], b1[j].y/v1[j], b1[j].z/v1[j], ux1[j]/v1[j], uy1[j]/v1[j], uz1[j]/v1[j]);
  for (int j=0; j<n2; j++)
    fprintf (fbubb, "%d %g %d %g %g %g %g %g %g %g\n", i, t,
	     j, v2[j], b2[j].x/v2[j], b2[j].y/v2[j], b2[j].z/v2[j], ux2[j]/v2[j], uy2[j]/v2[j], uz2[j]/v2[j]);
fflush(fdrop);
fflush(fbubb);
}

