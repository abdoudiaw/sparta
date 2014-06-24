/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "compute_surf.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NUM,PRESS,XPRESS,YPRESS,ZPRESS,XSHEAR,YSHEAR,ZSHEAR,KE,EROT,EVIB,ETOT};

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputeSurf::ComputeSurf(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute surf command");

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Compute surf mixture ID does not exist");

  nvalue = narg - 3;
  which = new int[nvalue];

  nvalue = 0;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) which[nvalue++] = NUM;
    else if (strcmp(arg[iarg],"press") == 0) which[nvalue++] = PRESS;
    else if (strcmp(arg[iarg],"px") == 0) which[nvalue++] = XPRESS;
    else if (strcmp(arg[iarg],"py") == 0) which[nvalue++] = YPRESS;
    else if (strcmp(arg[iarg],"pz") == 0) which[nvalue++] = ZPRESS;
    else if (strcmp(arg[iarg],"shx") == 0) which[nvalue++] = XSHEAR;
    else if (strcmp(arg[iarg],"shy") == 0) which[nvalue++] = YSHEAR;
    else if (strcmp(arg[iarg],"shz") == 0) which[nvalue++] = ZSHEAR;
    else if (strcmp(arg[iarg],"ke") == 0) which[nvalue++] = KE;
    else if (strcmp(arg[iarg],"erot") == 0) which[nvalue++] = EROT;
    else if (strcmp(arg[iarg],"evib") == 0) which[nvalue++] = EVIB;
    else if (strcmp(arg[iarg],"etot") == 0) which[nvalue++] = ETOT;
    else error->all(FLERR,"Illegal compute surf command");
    iarg++;
  }

  surf_tally_flag = 1;
  timeflag = 1;
  per_surf_flag = 1;
  ngroup = particle->mixture[imix]->ngroup;
  ntotal = ngroup*nvalue;
  size_per_surf_cols = ntotal;

  nlocal = maxlocal = 0;
  glob2loc = loc2glob = NULL;
  array_surf = NULL;
  normflux = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSurf::~ComputeSurf()
{
  delete [] which;
  memory->destroy(glob2loc);
  memory->destroy(loc2glob);
  memory->destroy(array_surf);
  memory->destroy(normflux);
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute surf mixture has changed");

  // local copies

  dimension = domain->dimension;
  lines = surf->lines;
  tris = surf->tris;

  // allocate and initialize glob2loc indices

  if (dimension == 2) nsurf = surf->nline;
  else nsurf = surf->ntri;

  memory->destroy(glob2loc);
  memory->create(glob2loc,nsurf,"surf:glob2loc");
  for (int i = 0; i < nsurf; i++) glob2loc[i] = -1;

  // normflux for each surface element I own, based on area and timestep size
  // one-time only initialization
  // NOTE: eventually will be only for surfs I know about (local + ghost)

  if (normflux == NULL)  {
    memory->create(normflux,nsurf,"surf:normflux");

    double nfactor = update->dt/update->fnum;
    int dimension = domain->dimension;
    double tmp;

    for (int i = 0; i < nsurf; i++) {
      if (dimension == 2) normflux[i] = surf->line_size(i);
      else normflux[i] = surf->tri_size(i,tmp);
      normflux[i] *= nfactor;
    }
  }

  // initialize tally array in case accessed before a tally timestep

  clear();
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::compute_per_surf()
{
  invoked_per_surf = update->ntimestep;

  // normalize all local tallies (except NUM) by normflux

  int iglobal,m;
  double norm;

  for (int i = 0; i < nlocal; i++) {
    iglobal = loc2glob[i];
    norm = normflux[iglobal];
    for (int j = 0; j < ntotal; j++) {
      m = j % nvalue;
      if (which[m] != NUM) array_surf[i][j] /= norm;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::clear()
{
  int i,j;
  double *norm;

  // reset all set glob2loc values to -1

  for (i = 0; i < nlocal; i++)
    glob2loc[loc2glob[i]] = -1;
  nlocal = 0;
}

/* ----------------------------------------------------------------------
   tally values for a single particle colliding with Isurf
   eold = vector of velocity/energy values before collision
        = 0,1,2 = velocity components, 3 = rotational eng, 4 = vibrational eng
   p = particle ptr with velocity/energy values after collision
------------------------------------------------------------------------- */

void ComputeSurf::surf_tally(int isurf, double *eold, Particle::OnePart *p)
{
  // skip species not in mixture group

  int ispecies = p->ispecies;
  int igroup = particle->mixture[imix]->species2group[ispecies];
  if (igroup < 0) return;

  // ilocal = local index of global isurf
  // if 1st particle hitting isurf, add isurf to local list
  // grow local list if needed

  double *vec;

  int ilocal = glob2loc[isurf];
  if (ilocal < 0) {
    if (nlocal == maxlocal) grow();
    ilocal = nlocal++;
    loc2glob[ilocal] = isurf;
    glob2loc[isurf] = ilocal;
    vec = array_surf[ilocal];
    for (int i = 0; i < ntotal; i++) vec[i] = 0.0;
  }

  // tally all values associated with group into array
  // set nflag and tflag after normal and tangent computation is done once

  double pre,post,vsqpre,vsqpost;
  double vnorm[3],vdelta[3],vtang[3];
  double *norm;

  if (dimension == 2) norm = lines[isurf].norm;
  else norm = tris[isurf].norm;

  double mass = particle->species[ispecies].mass;
  double mvv2e = update->mvv2e;

  vec = array_surf[ilocal];
  int k = igroup*nvalue;
  int nflag = 0;
  int tflag = 0;

  for (int m = 0; m < nvalue; m++) {
    switch (which[m]) {
    case NUM:
      vec[k++] += 1.0;
      break;
    case PRESS:
      pre = MathExtra::dot3(eold,norm);
      post = MathExtra::dot3(p->v,norm);
      vec[k++] += mass * (post-pre);
      break;
    case XPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtra::sub3(p->v,eold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
      }
      vec[k++] -= mass * vnorm[0];
      break;
    case YPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtra::sub3(p->v,eold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
      }
      vec[k++] -= mass * vnorm[1];
      break;
    case ZPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtra::sub3(p->v,eold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
      }
      vec[k++] -= mass * vnorm[2];
      break;
    case XSHEAR:
      if (!tflag) {
        tflag = 1;
        MathExtra::sub3(p->v,eold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
        MathExtra::sub3(vdelta,vnorm,vtang);
      }
      vec[k++] -= mass * vtang[0];
      break;
    case YSHEAR:
      if (!tflag) {
        tflag = 1;
        MathExtra::sub3(p->v,eold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
        MathExtra::sub3(vdelta,vnorm,vtang);
      }
      vec[k++] -= mass * vtang[1];
      break;
    case ZSHEAR:
      if (!tflag) {
        tflag = 1;
        MathExtra::sub3(p->v,eold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
        MathExtra::sub3(vdelta,vnorm,vtang);
      }
      vec[k++] -= mass * vtang[2];
      break;
    case KE:
      vsqpre = MathExtra::lensq3(eold);
      vsqpost = MathExtra::lensq3(p->v);
      vec[k++] -= 0.5*mvv2e*mass * (vsqpost-vsqpre);
      break;
    case EROT:
      vec[k++] -= p->erot - eold[3];
      break;
    case EVIB:
      vec[k++] -= p->evib - eold[4];
      break;
    case ETOT:
      vsqpre = MathExtra::lensq3(eold);
      vsqpost = MathExtra::lensq3(p->v);
      vec[k++] -= 0.5*mvv2e*mass*(vsqpost-vsqpre) + 
        (p->erot-eold[3]) + (p->evib-eold[4]);
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
   input N is value from 1 to Ncols
------------------------------------------------------------------------- */

double *ComputeSurf::normptr(int n)
{
  return NULL;
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
------------------------------------------------------------------------- */

int ComputeSurf::surfinfo(int *&locptr)
{
  locptr = loc2glob;
  return nlocal;
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::grow()
{
  maxlocal += DELTA;
  memory->grow(loc2glob,maxlocal,"surf:loc2glob");
  memory->grow(array_surf,maxlocal,ntotal,"surf:array_surf");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

bigint ComputeSurf::memory_usage()
{
  bigint bytes = 0;
  bytes += ntotal*maxlocal * sizeof(double);
  bytes += maxlocal * sizeof(int);
  bytes += nsurf * sizeof(int);
  return bytes;
}
