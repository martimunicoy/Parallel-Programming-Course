/*
 * Copyright (c) 2016, NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef FP32
#define FP64
#endif

#ifdef FP64
typedef double real;
const real SOFTENING_SQUARED = 0.01;
#define RSQRT(x) 1.0 / sqrt((x))
#define POW(x,y) pow((x),(y))
#else
typedef float real;
const real SOFTENING_SQUARED = 0.01f;
#define RSQRT(x) 1.0f / sqrtf((x))
#define POW(x,y) powf((x),(y))
#endif

typedef struct { real x, y, z; }    real3;
typedef struct { real x, y, z, w; } real4;

real3 bodyBodyInteraction(real4 iPos, real4 jPos)
{
  real rx, ry, rz;

  rx = jPos.x - iPos.x;
  ry = jPos.y - iPos.y;
  rz = jPos.z - iPos.z;

  real distSqr = rx*rx+ry*ry+rz*rz;
  distSqr += SOFTENING_SQUARED;

  real s = jPos.w / POW(distSqr,3.0/2.0);

  real3 f;
  f.x = rx * s;
  f.y = ry * s;
  f.z = rz * s;

  return f;
}

void integrate(real4 * out, real4 * in,
               real3 * vel, real3 * force,
               real    dt,  int n)
{
  int i, j;
  for (i = 0; i < n; i++)
  {
    real fx=0, fy=0, fz=0;

    for (j = 0; j < n; j++)
    {
      real3 ff = bodyBodyInteraction(in[i], in[j]);
      fx += ff.x;  fy += ff.y; fz += ff.z;
    }

    force[i].x = fx;  force[i].y = fy;  force[i].z = fz;
  }

  for (i = 0; i < n; i++)
  {
    real fx = force[i].x, fy = force[i].y, fz = force[i].z;
    real px = in[i].x,    py = in[i].y,    pz = in[i].z,    invMass = in[i].w;
    real vx = vel[i].x,   vy = vel[i].y,   vz = vel[i].z;

    // acceleration = force / mass; 
    // new velocity = old velocity + acceleration * deltaTime
    vx += (fx * invMass) * dt;
    vy += (fy * invMass) * dt;
    vz += (fz * invMass) * dt;

    // new position = old position + velocity * deltaTime
    px += vx * dt;
    py += vy * dt;
    pz += vz * dt;

    out[i].x = px;
    out[i].y = py;
    out[i].z = pz;
    out[i].w = invMass;

    vel[i].x = vx;
    vel[i].y = vy;
    vel[i].z = vz;
  }
}

real dot(real v0[3], real v1[3])
{
  return v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2];
}

real normalize(real vector[3])
{
  float dist = sqrt(dot(vector, vector));
  if (dist > 1e-6)
  {
    vector[0] /= dist;
    vector[1] /= dist;
    vector[2] /= dist;
  }
  return dist;
}

void cross(real out[3], real v0[3], real v1[3])
{
  out[0] = v0[1]*v1[2]-v0[2]*v1[1];
  out[1] = v0[2]*v1[0]-v0[0]*v1[2];
  out[2] = v0[0]*v1[1]-v0[1]*v1[0];
}

void randomizeBodies(real4* pos, 
                     real3* vel, 
                     float clusterScale, 
                     float velocityScale, 
                     int   n)
{
  srand(42);
  float scale = clusterScale;
  float vscale = scale * velocityScale;
  float inner = 2.5f * scale;
  float outer = 4.0f * scale;

  int p = 0, v=0;
  int i = 0;
  while (i < n)
  {
    real x, y, z;
    x = rand() / (float) RAND_MAX * 2 - 1;
    y = rand() / (float) RAND_MAX * 2 - 1;
    z = rand() / (float) RAND_MAX * 2 - 1;

    real point[3] = {x, y, z};
    real len = normalize(point);
    if (len > 1)
      continue;

    pos[i].x= point[0] * (inner + (outer - inner)*rand() / (real) RAND_MAX);
    pos[i].y= point[1] * (inner + (outer - inner)*rand() / (real) RAND_MAX);
    pos[i].z= point[2] * (inner + (outer - inner)*rand() / (real) RAND_MAX);
    pos[i].w= 1.0f;

    x = 0.0f; 
    y = 0.0f; 
    z = 1.0f; 
    real axis[3] = {x, y, z};
    normalize(axis);

    if (1 - dot(point, axis) < 1e-6)
    {
      axis[0] = point[1];
      axis[1] = point[0];
      normalize(axis);
    }
    real vv[3] = {(real)pos[i].x, (real)pos[i].y, (real)pos[i].z};
    cross(vv, vv, axis);
    vel[i].x = vv[0] * vscale;
    vel[i].y = vv[1] * vscale;
    vel[i].z = vv[2] * vscale;

    i++;
  }
}

real3 average(real4 * p, int n)
{
  int i;
  real3 av= {0.0, 0.0, 0.0};
  for (i = 0; i < n; i++)
  {
    av.x += p[i].x;
    av.y += p[i].y;
    av.z += p[i].z;
  }
  av.x /= n;
  av.y /= n;
  av.z /= n;
  return av;
}

int main(int argc, char** argv)
{
  int i, n = 20000;
  int iterations = 10;
  real dt = 0.01667;

  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) iterations = atoi(argv[2]);

  real4 *pin  = (real4*)malloc(n * sizeof(real4));
  real4 *pout = (real4*)malloc(n * sizeof(real4));
  real3 *v    = (real3*)malloc(n * sizeof(real3));
  real3 *f    = (real3*)malloc(n * sizeof(real3));

  randomizeBodies(pin, v,  1.54f, 8.0f, n);

  printf("n=%d bodies for %d iterations:\n", n, iterations);

  for (i = 0; i < iterations; i++)
  {
    integrate (pout, pin, v, f, dt, n);
    real4 *t = pout;
    pout = pin; 
    pin = t;
  }

  real3 p_av= average( pout, n);
  printf("Average position: (%f,%f,%f)\n", p_av.x, p_av.y, p_av.z);

  free(pin);  free(pout);  free(v);  free(f);

  return 0;
}

