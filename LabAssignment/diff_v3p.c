#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef REAL
#define REAL float
#endif

#ifndef M_PI
#define M_PI (3.1415926535897932384626)
#endif

// v1: change the order of memory storage to speed up the Cache acces
#define F1(x,y,z) F1[((ny)*(nx)*(z))+((nx)*(y))+(x)]
#define F2(x,y,z) F2[((ny)*(nx)*(z))+((nx)*(y))+(x)]

void diffusion (REAL * F1, REAL * F2, 
                int nx, int ny, int nz,
                REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                REAL cb, REAL cc, int time)
{
  int t, x, y, z;
  int west, east, north, south, top, down;
  for (t=0; t < time; ++t)
  {
    // Change x, y and z for constants pointing to the right points for vertices, sides and faces
    F2(0,0,0) =
       F1(0,0,0) *cc + F1(0,0,0) *cw +
       F1(1,0,0) *ce + F1(0,0,0) *cn +
       F1(0,1,0) *cs + F1(0,0,0) *cb +
       F1(0,0,1) *ct;
    F2(nx-1,ny-1,nz-1) =
       F1(nx-1,ny-1,nz-1) *cc + F1(nx-2,ny-1,nz-1) *cw +
       F1(nx-1,ny-1,nz-1) *ce + F1(nx-1,ny-2,nz-1) *cn +
       F1(nx-1,ny-1,nz-1) *cs + F1(nx-1,ny-1,nz-2) *cb +
       F1(nx-1,ny-1,nz-1) *ct;
    F2(nx-1,0,0) =
       F1(nx-1,0,0) *cc + F1(nx-2,0,0) *cw +
       F1(nx-1,0,0) *ce + F1(nx-1,0,0) *cn +
       F1(nx-1,1,0) *cs + F1(nx-1,0,0) *cb +
       F1(nx-1,0,1) *ct;
    F2(0,ny-1,0) =
       F1(0,ny-1,0) *cc + F1(0,ny-1,0) *cw +
       F1(1,ny-1,0) *ce + F1(0,ny-2,0) *cn +
       F1(0,ny-1,0) *cs + F1(0,ny-1,0) *cb +
       F1(0,ny-1,1) *ct;
    F2(nx-1,ny-1,0) =
       F1(nx-1,ny-1,0) *cc + F1(nx-2,ny-1,0) *cw +
       F1(nx-1,ny-1,0) *ce + F1(nx-1,ny-2,0) *cn +
       F1(nx-1,ny-1,0) *cs + F1(nx-1,ny-1,0) *cb +
       F1(nx-1,ny-1,1) *ct;
    F2(0,0,nz-1) =
       F1(0,0,nz-1) *cc + F1(0,0,nz-1) *cw +
       F1(1,0,nz-1) *ce + F1(0,0,nz-1) *cn +
       F1(0,1,nz-1) *cs + F1(0,0,nz-2) *cb +
       F1(0,0,nz-1) *ct;
    F2(nx-1,0,nz-1) =
       F1(nx-1,0,nz-1) *cc + F1(nx-2,0,nz-1) *cw +
       F1(nx-1,0,nz-1) *ce + F1(nx-1,0,nz-1) *cn +
       F1(nx-1,1,nz-1) *cs + F1(nx-1,0,nz-2) *cb +
       F1(nx-1,0,nz-1) *ct;
    F2(0,ny-1,nz-1) =
       F1(0,ny-1,nz-1) *cc + F1(0,ny-1,nz-1) *cw +
       F1(1,ny-1,nz-1) *ce + F1(0,ny-2,nz-1) *cn +
       F1(0,ny-1,nz-1) *cs + F1(0,ny-1,nz-2) *cb +
       F1(0,ny-1,nz-1) *ct;
    for (x=1; x < nx-1; ++x)
    {
      F2(x,0,0) =
         F1(x,0,0)   *cc + F1(x-1,0,0) *cw +
         F1(x+1,0,0) *ce + F1(x,0,0)   *cn +
         F1(x,1,0)   *cs + F1(x,0,0)   *cb +
         F1(x,0,1)   *ct;
      F2(x,ny-1,nz-1) =
         F1(x,ny-1,nz-1)   *cc + F1(x-1,ny-1,nz-1) *cw +
         F1(x+1,ny-1,nz-1) *ce + F1(x,ny-2,nz-1)   *cn +
         F1(x,ny-1,nz-1)   *cs + F1(x,ny-1,nz-2)   *cb +
         F1(x,ny-1,nz-1)   *ct;
      F2(x,0,nz-1) =
         F1(x,0,nz-1)   *cc + F1(x-1,0,nz-1) *cw +
         F1(x+1,0,nz-1) *ce + F1(x,0,nz-1)   *cn +
         F1(x,1,nz-1)   *cs + F1(x,0,nz-2)   *cb +
         F1(x,0,nz-1)   *ct;
      F2(x,ny-1,0) =
         F1(x,ny-1,0)   *cc + F1(x-1,ny-1,0) *cw +
         F1(x+1,ny-1,0) *ce + F1(x,ny-2,0)   *cn +
         F1(x,ny-1,0)   *cs + F1(x,ny-1,0)   *cb +
         F1(x,ny-1,1)   *ct;
    }
    for (y=1; y < ny-1; ++y)
    {
      F2(0,y,0) =
         F1(0,y,0)   *cc + F1(0,y,0)   *cw +
         F1(1,y,0)   *ce + F1(0,y-1,0) *cn +
         F1(0,y+1,0) *cs + F1(0,y,0)   *cb +
         F1(0,y,1)   *ct;
      F2(nx-1,y,nz-1) =
         F1(nx-1,y,nz-1)   *cc + F1(nx-2,y,nz-1)   *cw +
         F1(nx-1,y,nz-1)   *ce + F1(nx-1,y-1,nz-1) *cn +
         F1(nx-1,y+1,nz-1) *cs + F1(nx-1,y,nz-2)   *cb +
         F1(nx-1,y,nz-1)   *ct;
      F2(0,y,nz-1) =
         F1(0,y,nz-1)   *cc + F1(0,y,nz-1)   *cw +
         F1(1,y,nz-1)   *ce + F1(0,y-1,nz-1) *cn +
         F1(0,y+1,nz-1) *cs + F1(0,y,nz-2)   *cb +
         F1(0,y,nz-1)   *ct;
      F2(nx-1,y,0) =
         F1(nx-1,y,0)   *cc + F1(nx-2,y,0)   *cw +
         F1(nx-1,y,0)   *ce + F1(nx-1,y-1,0) *cn +
         F1(nx-1,y+1,0) *cs + F1(nx-1,y,0)   *cb +
         F1(nx-1,y,1)   *ct;
    }
    for (z=1; z < nz-1; ++z)
    {
      F2(0,0,z) =
         F1(0,0,z)   *cc + F1(0,0,z)   *cw +
         F1(1,0,z)   *ce + F1(0,0,z)   *cn +
         F1(0,1,z)   *cs + F1(0,0,z-1) *cb +
         F1(0,0,z+1) *ct;
      F2(nx-1,ny-1,z) =
         F1(nx-1,ny-1,z)   *cc + F1(nx-2,ny-1,z)   *cw +
         F1(nx-1,ny-1,z)   *ce + F1(nx-1,ny-2,z)   *cn +
         F1(nx-1,ny-1,z)   *cs + F1(nx-1,ny-1,z-1) *cb +
         F1(nx-1,ny-1,z+1) *ct;
      F2(0,ny-1,z) =
         F1(0,ny-1,z)   *cc + F1(0,ny-1,z)   *cw +
         F1(1,ny-1,z)   *ce + F1(0,ny-2,z)   *cn +
         F1(0,ny-1,z)   *cs + F1(0,ny-1,z-1) *cb +
         F1(0,ny-1,z+1) *ct;
      //!!!!
      F2(nx-1,0,z) = 
         F1(nx-1,0,z)   *cc + F1(nx-2,0,z)   *cw +
         F1(nx-1,0,z)   *ce + F1(nx-1,0,z)   *cn +
         F1(nx-1,1,z)   *cs + F1(nx-1,0,z-1) *cb +
         F1(nx-1,0,z+1) *ct;
    }
    for (y=1; y < ny-1; ++y)
    {
      for (x=1; x < nx-1; ++x)
      {
        F2(x,y,0) =
             F1(x,y,0)   *cc + F1(x-1,y,0) *cw +
             F1(x+1,y,0) *ce + F1(x,y-1,0) *cn +
             F1(x,y+1,0) *cs + F1(x,y,0)   *cb +
             F1(x,y,1) *ct;
        F2(x,y,nz-1) =
             F1(x,y,nz-1)   *cc + F1(x-1,y,nz-1) *cw +
             F1(x+1,y,nz-1) *ce + F1(x,y-1,nz-1) *cn +
             F1(x,y+1,nz-1) *cs + F1(x,y,nz-2) *cb +
             F1(x,y,nz-1)   *ct;
      }
    }
    for (z=1; z < nz-1; ++z)
    {
      for (x=1; x < nx-1; ++x)
      {
        F2(x,0,z) =
             F1(x,0,z)   *cc + F1(x-1,0,z) *cw +
             F1(x+1,0,z) *ce + F1(x,0,z)   *cn +
             F1(x,1,z) *cs + F1(x,0,z-1) *cb +
             F1(x,0,z+1) *ct;
        F2(x,ny-1,z) =
             F1(x,ny-1,z)   *cc + F1(x-1,ny-1,z) *cw +
             F1(x+1,ny-1,z) *ce + F1(x,ny-2,z) *cn +
             F1(x,ny-1,z)   *cs + F1(x,ny-1,z-1) *cb +
             F1(x,ny-1,z+1) *ct;
      }
      for (y=1; y < ny-1; ++y)
      {
        F2(0,y,z) =
             F1(0,y,z)   *cc + F1(0,y,z)   *cw +
             F1(1,y,z)   *ce + F1(0,y-1,z) *cn +
             F1(0,y+1,z) *cs + F1(0,y,z-1) *cb +
             F1(0,y,z+1) *ct;
        F2(nx-1,y,z) =
             F1(nx-1,y,z)   *cc + F1(nx-2,y,z)   *cw +
             F1(nx-1,y,z)   *ce + F1(nx-1,y-1,z) *cn +
             F1(nx-1,y+1,z) *cs + F1(nx-1,y,z-1) *cb +
             F1(nx-1,y,z+1) *ct;
        for (x = 1; x < nx-1; ++x)
        {
          F2(x,y,z) =
             F1(x,y,z)   *cc + F1(x-1,y,z) *cw +
             F1(x+1,y,z) *ce + F1(x,y-1,z) *cn +
             F1(x,y+1,z) *cs + F1(x,y,z-1) *cb +
             F1(x,y,z+1) *ct;
        }
      }
    }
    REAL *tt = F1;   F1 = F2;  F2 = tt;  // swap matrices
  }
}

void init (REAL *F1, const int nx, const int ny, const int nz,
           const REAL kx, const REAL ky, const REAL kz,
           const REAL dx, const REAL dy, const REAL dz,
           const REAL kappa, const REAL time)
{
  REAL ax, ay, az;
  int jz, jy, jx;
  ax = exp(-kappa*time*(kx*kx));
  ay = exp(-kappa*time*(ky*ky));
  az = exp(-kappa*time*(kz*kz));
  for (jz = 0; jz < nz; jz++) 
    for (jy = 0; jy < ny; jy++) 
      for (jx = 0; jx < nx; jx++)
      {
        REAL x = dx*((REAL)(jx + 0.5));
        REAL y = dy*((REAL)(jy + 0.5));
        REAL z = dz*((REAL)(jz + 0.5));
        REAL f0 = (REAL)0.125
          *(1.0 - ax*cos(kx*x))
          *(1.0 - ay*cos(ky*y))
          *(1.0 - az*cos(kz*z));
        F1(jx,jy,jz) = f0;
      }
}

REAL sum_values (REAL *F1, const int nx, const int ny, const int nz)
{
  REAL sum=0.0;
  int jz, jy, jx;
  for (jz = 0; jz < nz; jz++)
    for (jy = 0; jy < ny; jy++)
      for (jx = 0; jx < nx; jx++)
        sum += F1(jx,jy,jz);
  return sum;
}


void printmatrix(REAL *F1, const int nx, const int ny, const int nz)
{
  int i, j, k;
  for (k = 0; k < nz; ++k)
  {
    printf("\nk = %d\n\n", k);
    for (j = 0; j < ny; ++j)
    {
      for (i = 0; i < nx; ++i)
      {
        printf("%f\t", F1(i,j,k));
      }
      printf("\n");
    }
  }
}

int main(int argc, char *argv[]) 
{ 
  int  jz, jy, jx;
  int  NX=11, NY=11, NZ=11;

  if (argc>1) { NX= atoi(argv[1]); } // get  first command line parameter
  if (argc>2) { NY= atoi(argv[2]); } // get second command line parameter
  if (argc>3) { NZ= atoi(argv[3]); } // get  third command line parameter
  if (NX < 1 || NY < 1 || NZ < 1)
  {
    printf("arguments: NX NY NZ\n");
    return 1;
  }

  REAL *f1 = (REAL *) malloc(sizeof(REAL)*NX*NY*NZ);
  REAL *f2 = (REAL *) malloc(sizeof(REAL)*NX*NY*NZ);

  REAL *f_final = NULL;

  REAL  time  = 0.0;
  int   count = 0;  

  REAL l, dx, dy, dz, kx, ky, kz, kappa, dt;
  REAL ce, cw, cn, cs, ct, cb, cc;

  l = 1.0;
  kappa = 0.1;
  dx = l / NX;  dy = l / NY;  dz = l / NZ;
  kx = ky = kz = 2.0 * M_PI;
  dt = dx*dy*dz / kappa;
  count = 0.01 / dt;
  f_final = (count % 2)? f2 : f1;

  init(f1, NX, NY, NZ, kx, ky, kz, dx, dy, dz, kappa, time);

  REAL err = sum_values(f1, NX, NY, NZ);

  ce = cw = kappa*dt/(dx*dx);
  cn = cs = kappa*dt/(dy*dy);
  ct = cb = kappa*dt/(dz*dz);
  cc = 1.0 - (ce + cw + cn + cs + ct + cb);

  printf("Running diffusion kernel with NX=%d, NY=%d, NZ=%d, %d times\n", 
         NX, NY, NZ, count);

  diffusion(f1, f2, NX, NY, NZ, ce, cw, cn, cs, ct, cb, cc, count);

  err = err - sum_values(f_final, NX, NY, NZ);
  fprintf(stderr, "Accuracy     : %E\n", err);
  printmatrix(f_final, NX, NY, NZ);
  
  free(f1); free(f2);
  return 0;
}
