/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */
#include <chrono>
#include "spatype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "update.h"
#include "math_const.h"
#include "particle.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "domain.h"
#include "comm.h"
#include "collide.h"
#include "grid.h"
#include "surf.h"
#include "surf_collide.h"
#include "surf_react.h"
#include "input.h"
#include "output.h"
#include "geometry.h"
#include "random_mars.h"
#include "timer.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "random_knuth.h"
#include <string>  
#include <cmath>
#include <iostream>
#include <string>
#include <map>

#include "/opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense"
// #include "/usr/include/eigen3/Eigen/Dense"

#include <random>
#include "math_const.h"
#include "react_bird.h"
#include "surfaceDataInterpolator.h"


static const double eV2Kelvin = 11604.505;
static const std::string filename = "lim.txt";
static const double kB = 1.38064852e-23;
static const double protonCharge = 1.60217662e-19;

using namespace SPARTA_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};      // several files
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid
enum{TALLYAUTO,TALLYREDUCE,TALLYRVOUS};         // same as Surf
enum{PERAUTO,PERCELL,PERSURF};                  // several files
enum{NOFIELD,CFIELD,PFIELD,GFIELD};             // several files



const double OXYGEN_MASS = 16.0;
const double TUNGSTEN_MASS = 184.0;
const double NEUTRAL_CHARGE = 0.0;


#define MAXSTUCK 20
#define EPSPARAM 1.0e-7

// either set ID or PROC/INDEX, set other to -1

//#define MOVE_DEBUG 1              // un-comment to debug one particle
#define MOVE_DEBUG_ID 308143534  // particle ID
#define MOVE_DEBUG_PROC -1        // owning proc
#define MOVE_DEBUG_INDEX -1   // particle index on owning proc
#define MOVE_DEBUG_STEP 4107    // timestep

int ionization_mapping[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
int recombination_mapping[] = {1, 2, 3, 4, 5, 6, 7, 8, 0};
double electron_charge = 1.60217662e-19;
double proton_mass = 1.6726219e-27;
// / Set background parameters
const double background_Z = 1.0; // Deuterium
const double background_mass = 2.0; // Deuterium
const double Q = 1.60217662e-19;
const double EPS0 = 8.854187e-12;
const double MI = 1.6737236e-27;
const double ME = 9.10938356e-31;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);
std::normal_distribution<double> distribution(0.0,1.0);
const double TWO_PI = 2.0 * M_PI;


/* ---------------------------------------------------------------------- */

Update::Update(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  ntimestep = 0;
  firststep = laststep = 0;
  beginstep = endstep = 0;
  runflag = 0;

  unit_style = NULL;
  set_units("si");

  fnum = 1.0;
  nrho = 1.0;
  vstream[0] = vstream[1] = vstream[2] = 0.0;
  temp_thermal = 273.15;
  optmove_flag = 0;
  fstyle = NOFIELD;
  fieldID = NULL;

  maxmigrate = 0;
  mlist = NULL;

  nslist_compute = nblist_compute = 0;
  slist_compute = blist_compute = NULL;
  slist_active = blist_active = NULL;

  nulist_surfcollide  = 0;
  ulist_surfcollide = NULL;

  ranmaster = new RanMars(sparta);

  reorder_period = 0;
  global_mem_limit = 0;
  mem_limit_grid_flag = 0;

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
  if (copymode) return;

  delete [] unit_style;
  delete [] fieldID;
  memory->destroy(mlist);
  delete [] slist_compute;
  delete [] blist_compute;
  delete [] slist_active;
  delete [] blist_active;
  delete [] ulist_surfcollide;
  delete ranmaster;
}

/* ---------------------------------------------------------------------- */

void Update::set_units(const char *style)
{
  // physical constants from:
  // http://physics.nist.gov/cuu/Constants/Table/allascii.txt

  if (strcmp(style,"cgs") == 0) {
    boltz = 1.380649e-16;
    mvv2e = 1.0;
    dt = 1.0;

  } else if (strcmp(style,"si") == 0) {
    boltz = 1.380649e-23;
    mvv2e = 1.0;
    dt = 1.0;

  } else error->all(FLERR,"Illegal units command");

  delete [] unit_style;
  int n = strlen(style) + 1;
  unit_style = new char[n];
  strcpy(unit_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
  // init the Update class if performing a run, else just return
  // only set first_update if a run is being performed

  if (runflag == 0) return;
  first_update = 1;

  if (optmove_flag) {
    if (!grid->uniform)
      error->all(FLERR,"Cannot use optimized move with non-uniform grid");
    else if (surf->exist)
      error->all(FLERR,"Cannot use optimized move when surfaces are defined");
    else {
      for (int ifix = 0; ifix < modify->nfix; ifix++) {
        if (strstr(modify->fix[ifix]->style,"adapt") != NULL)
          error->all(FLERR,"Cannot use optimized move with fix adapt");
      }
    }
  }

  // choose the appropriate move method

  if (domain->dimension == 3) {
    if (surf->exist)
      moveptr = &Update::move<3,1,0>;
    else {
      if (optmove_flag) moveptr = &Update::move<3,0,1>;
      else moveptr = &Update::move<3,0,0>;
    }
  } else if (domain->axisymmetric) {
    if (surf->exist)
      moveptr = &Update::move<1,1,0>;
    else {
      if (optmove_flag) moveptr = &Update::move<1,0,1>;
      else moveptr = &Update::move<1,0,0>;
    }
  } else if (domain->dimension == 2) {
    if (surf->exist)
      moveptr = &Update::move<2,1,0>;
    else {
      if (optmove_flag) moveptr = &Update::move<2,0,1>;
      else moveptr = &Update::move<2,0,0>;
    }
  }

  // checks on external field options

  if (fstyle == CFIELD) {
    if (domain->dimension == 2 && field[2] != 0.0)
      error->all(FLERR,"External field in z not allowed for 2d");
    if (domain->axisymmetric && field[1] != 0.0)
      error->all(FLERR,
                 "External field in y not allowed for axisymmetric model");
  } else if (fstyle == PFIELD) {
    ifieldfix = modify->find_fix(fieldID);
    if (ifieldfix < 0) error->all(FLERR,"External field fix ID not found");
    if (!modify->fix[ifieldfix]->per_particle_field)
      error->all(FLERR,"External field fix does not compute necessary field");
  } else if (fstyle == GFIELD) {
    ifieldfix = modify->find_fix(fieldID);
    if (ifieldfix < 0) error->all(FLERR,"External field fix ID not found");
    if (!modify->fix[ifieldfix]->per_grid_field)
      error->all(FLERR,"External field fix does not compute necessary field");
  }

  // moveperturb method is set if external field perturbs particle motion

  moveperturb = NULL;

  if (fstyle == CFIELD) {
    if (domain->dimension == 2) moveperturb = &Update::field2d;
    if (domain->dimension == 3) moveperturb = &Update::field3d;
  } else if (fstyle == PFIELD) {
    moveperturb = &Update::field_per_particle;
    field_active = modify->fix[ifieldfix]->field_active;
  } else if (fstyle == GFIELD) {
    moveperturb = &Update::field_per_grid;
    field_active = modify->fix[ifieldfix]->field_active;
  }

  if (moveperturb) perturbflag = 1;
  else perturbflag = 0;
}

/* ---------------------------------------------------------------------- */

void Update::setup()
{
  // initialize counters in case stats outputs them
  // initialize running stats before each run

  ntouch_one = ncomm_one = 0;
  nboundary_one = nexit_one = 0;
  nscheck_one = nscollide_one = 0;
  surf->nreact_one = 0;

  first_running_step = update->ntimestep;
  niterate_running = 0;
  nmove_running = ntouch_running = ncomm_running = 0;
  nboundary_running = nexit_running = 0;
  nscheck_running = nscollide_running = 0;
  surf->nreact_running = 0;
  nstuck = naxibad = 0;

  collide_react = collide_react_setup();
  bounce_tally = bounce_setup();

  dynamic = 0;
  dynamic_setup();

  modify->setup();
  output->setup(1);
}

/* ---------------------------------------------------------------------- */

void Update::run(int nsteps)
{
  int n_start_of_step = modify->n_start_of_step;
  int n_end_of_step = modify->n_end_of_step;

  // external per grid cell field
  // only evaluate once at beginning of run b/c time-independent
  // fix calculates field acting at center point of all grid cells

  if (fstyle == GFIELD && fieldfreq == 0)
    modify->fix[ifieldfix]->compute_field();

  // cellweightflag = 1 if grid-based particle weighting is ON

  int cellweightflag = 0;
  if (grid->cellweightflag) cellweightflag = 1;

  // loop over timesteps

  for (int i = 0; i < nsteps; i++) {

    ntimestep++;

    if (collide_react) collide_react_reset();
    if (bounce_tally) bounce_set(ntimestep);

    timer->stamp();

    // dynamic parameter updates

    if (dynamic) dynamic_update();

    // start of step fixes

    if (n_start_of_step) {
      modify->start_of_step();
      timer->stamp(TIME_MODIFY);
    }

    // move particles

    if (cellweightflag) particle->pre_weight();
    (this->*moveptr)();
    timer->stamp(TIME_MOVE);

    // communicate particles

    comm->migrate_particles(nmigrate,mlist);
    if (cellweightflag) particle->post_weight();
    timer->stamp(TIME_COMM);

    if (collide) {
      particle->sort();
      timer->stamp(TIME_SORT);

      collide->collisions();
      timer->stamp(TIME_COLLIDE);
    }

    if (collide_react) collide_react_update();

    // diagnostic fixes

    if (n_end_of_step) {
      modify->end_of_step();
      timer->stamp(TIME_MODIFY);
    }

    // all output

    if (ntimestep == output->next) {
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
}

/* ----------------------------------------------------------------------
   advect particles thru grid
   DIM = 2/3 for 2d/3d, 1 for 2d axisymmetric
   SURF = 0/1 for no surfs or surfs
   use multiple iterations of move/comm if necessary
------------------------------------------------------------------------- */

template < int DIM, int SURF, int OPT > void Update::move()
{
  bool hitflag;
  int m,icell,icell_original,nmask,outface,bflag,nflag,pflag,itmp;
  int side,minside,minsurf,nsurf,cflag,isurf,exclude,stuck_iterate;
  int pstart,pstop,entryexit,any_entryexit,reaction;
  surfint *csurfs;
  cellint *neigh;
  double dtremain,frac,newfrac,param,minparam,rnew,dtsurf,tc,tmp;
  double xnew[3],xhold[3],xc[3],vc[3],minxc[3],minvc[3];
  double *x,*v,*lo,*hi;
  double Lx,Ly,Lz,dx,dy,dz;
  double *boxlo, *boxhi;
  Grid::ParentCell *pcell;
  Surf::Tri *tri;
  Surf::Line *line;
  Particle::OnePart iorig;
  Particle::OnePart *particles;
  Particle::OnePart *ipart,*jpart;

  if (OPT) {
    boxlo = domain->boxlo;
    boxhi = domain->boxhi;
    Lx = boxhi[0] - boxlo[0];
    Ly = boxhi[1] - boxlo[1];
    Lz = boxhi[2] - boxlo[2];
    dx = Lx/grid->unx;
    dy = Ly/grid->uny;
    dz = Lz/grid->unz;
  }

  // for 2d and axisymmetry only
  // xnew,xc passed to geometry routines which use or set z component

  if (DIM < 3) xnew[2] = xc[2] = 0.0;

  // extend migration list if necessary

  int nlocal = particle->nlocal;
  int maxlocal = particle->maxlocal;

  if (nlocal > maxmigrate) {
    maxmigrate = maxlocal;
    memory->destroy(mlist);
    memory->create(mlist,maxmigrate,"particle:mlist");
  }

  // counters

  niterate = 0;
  ntouch_one = ncomm_one = 0;
  nboundary_one = nexit_one = 0;
  nscheck_one = nscollide_one = 0;
  surf->nreact_one = 0;

  // move/migrate iterations

  Grid::ChildCell *cells = grid->cells;
  Grid::ParentCell *pcells = grid->pcells;
  Surf::Tri *tris = surf->tris;
  Surf::Line *lines = surf->lines;
  double dt = update->dt;


    // ci->volume = MY_PI * (hi[1]*hi[1]-lo[1]*lo[1]) * (hi[0]-lo[0]);
    // printf(" ass ci->volume = %f\n", ci->volume);

        Grid::ChildInfo *cinfo = grid->cinfo;

  // external per particle field
  // fix calculates field acting on all owned particles

  if (fstyle == PFIELD) modify->fix[ifieldfix]->compute_field();

  // external per grid cell field
  // evaluate once every fieldfreq steps b/c time-dependent
  // fix calculates field acting at center point of all grid cells

  if (fstyle == GFIELD && fieldfreq && ((ntimestep-1) % fieldfreq == 0))
    modify->fix[ifieldfix]->compute_field();

  // one or more loops over particles
  // first iteration = all my particles
  // subsequent iterations = received particles

  while (1) {

    niterate++;
    particles = particle->particles;
    nmigrate = 0;
    entryexit = 0;

    if (niterate == 1) {
      pstart = 0;
      pstop = nlocal;
    }

    for (int i = pstart; i < pstop; i++) {
      pflag = particles[i].flag;

      // received from another proc and move is done
      // if first iteration, PDONE is from a previous step,
      //   set pflag to PKEEP so move the particle on this step
      // else do nothing

      if (pflag == PDONE) {
        pflag = particles[i].flag = PKEEP;
        if (niterate > 1) continue;
      }

      x = particles[i].x;
      v = particles[i].v;

      // get particle physical properties

      ipart = &particles[i];
      Particle::Species *species = particle->species;

      int isp = ipart->ispecies;
      double mass = species[isp].mass;
      double charge = species[isp].charge * electron_charge;
      // GET CELL
      icell = ipart->icell;
      // printf("Particle %d is in cell %d\n", i, icell);
      // get cell coordinates
      lo = cells[icell].lo;
      hi = cells[icell].hi;


      // cell xc is center of cell
      xc[0] = 0.5*(lo[0] + hi[0]);
      xc[1] = 0.5*(lo[1] + hi[1]);
      if (DIM == 3) xc[2] = 0.5*(lo[2] + hi[2]);

      
    // get electric and magnetic fields
    // printf("xc %g %g\n",x[0], xc[1]);
    PlasmaParams params = interpolatePlasmaData( xc[0], xc[1]);

    double Bfield[3];
    double E[3];
    Bfield[0] = params.b_r ;
    Bfield[1] = params.b_z;
    Bfield[2] = params.b_phi;
    double te = params.temp_e;
    double ne = params.dens_e;
    double flow = params.parr_flow;
        // printf(" bfields te ne %g %g %g %g %g %g\n", Bfield[0], Bfield[1], Bfield[2], te, ne, flow);
        // exit(0);
      exclude = -1;

      // apply moveperturb() to PKEEP and PINSERT since are computing xnew
      // not to PENTRY,PEXIT since are just re-computing xnew of sender
      // set xnew[2] to linear move for axisymmetry, will be remapped later
      // let pflag = PEXIT persist to check during axisymmetric cell crossing

      if (pflag == PKEEP) {
        // printf("Advecting pflag == PKEEP particle %d\n", i);
        dtremain = dt;
        // check whether charge is zero
        if (charge == 0.0) {
          // printf("Charge is zero\n");
          xnew[0] = x[0] + dtremain*v[0];
          xnew[1] = x[1] + dtremain*v[1];
          if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
          if (perturbflag)
            (this->*moveperturb)(i,particles[i].icell,dtremain,xnew,v);
        }
        else{
          // Apply collisional diffusion
          // get plasma background data
//      backgroundCollisions(xc, v, dtremain, species[isp].molwt, species[isp].charge, params);
      pusher_boris(xc, x, v, xnew, charge, mass,  dtremain , params);    
//      crossFieldDiffusion(xc, xnew, dtremain, Bfield );
        }

      // printf("B field at particle %d: %g %g %g\n", i, B[0], B[1], B[2]);
      } else if (pflag == PINSERT) {
        // printf("Advecting pflag == PINSERT particle %d\n", i);
        dtremain = particles[i].dtremain;
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
        if (perturbflag)
          (this->*moveperturb)(i,particles[i].icell,dtremain,xnew,v);
      } else if (pflag == PENTRY) {
//        printf("Advecting pflag == PENTRY particle %d\n", i);
        icell = particles[i].icell;
        if (cells[icell].nsplit > 1) {
          if (DIM == 3 && SURF) icell = split3d(icell,x);
          if (DIM < 3 && SURF) icell = split2d(icell,x);
          particles[i].icell = icell;
        }
        dtremain = particles[i].dtremain;
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
      } else if (pflag == PEXIT) {
        dtremain = particles[i].dtremain;
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
      } else if (pflag >= PSURF) {
        // get particle icell
        icell = particles[i].icell;
       
//// *** ///
    // // get temperature of surface
    double Eb = 11.1; // eV        
    double converteV2K = 11604.505;
    double wall_mass = 184.0 * 1.6726219e-27;
    double twall = Eb * converteV2K/2;
    double vrm = sqrt(2.0*update->boltz * twall / wall_mass);
    double vnorm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    // Normalize v
    v[0] /= vnorm;
    v[1] /= vnorm;
    v[2] /= vnorm;

    // Scale the normalized v by vrm
    v[0] *= vrm;
    v[1] *= vrm;
    v[2] *= vrm;

        // printf("Advecting pflag == PSURF particle %d\n", i);
        dtremain = particles[i].dtremain;
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
        if (pflag > PSURF) exclude = pflag - PSURF - 1;
      }

      // optimized move

      if (OPT) {
        int optmove = 1;

        if (xnew[0] < boxlo[0] || xnew[0] > boxhi[0])
          optmove = 0;

        if (xnew[1] < boxlo[1] || xnew[1] > boxhi[1])
          optmove = 0;

        if (DIM == 3) {
          if (xnew[2] < boxlo[2] || xnew[2] > boxhi[2])
            optmove = 0;
        }

        if (optmove) {
          const int ip = static_cast<int>((xnew[0] - boxlo[0])/dx);
          const int jp = static_cast<int>((xnew[1] - boxlo[1])/dy);
          int kp = 0;
          if (DIM == 3) kp = static_cast<int>((xnew[2] - boxlo[2])/dz);

          int cellIdx = (kp*grid->uny + jp)*grid->unx + ip + 1;

          // particle outside ghost grid halo must use standard move

          if (grid->hash->find(cellIdx) != grid->hash->end()) {
          
            int icell = (*(grid->hash))[cellIdx];

            // reset particle cell and coordinates

            particles[i].icell = icell;
            particles[i].flag = PKEEP;
            x[0] = xnew[0];
            x[1] = xnew[1];
            x[2] = xnew[2];

            if (cells[icell].proc != me) {
              mlist[nmigrate++] = i;
              particles[i].flag = PDONE;
              ncomm_one++;
            }

            continue;
          }
        }
      }

      particles[i].flag = PKEEP;
      icell = particles[i].icell;
      lo = cells[icell].lo;
      hi = cells[icell].hi;
      neigh = cells[icell].neigh;
      nmask = cells[icell].nmask;
      stuck_iterate = 0;
      ntouch_one++;

      // advect one particle from cell to cell and thru surf collides til done

      //int iterate = 0;

      while (1) {

#ifdef MOVE_DEBUG
        if (DIM == 3) {
          if (ntimestep == MOVE_DEBUG_STEP &&
              (MOVE_DEBUG_ID == particles[i].id ||
               (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
            printf("PARTICLE %d %ld: %d %d: %d: x %g %g %g: xnew %g %g %g: %d "
                   CELLINT_FORMAT ": lo %g %g %g: hi %g %g %g: DTR %g\n",
                   me,update->ntimestep,i,particles[i].id,
                   cells[icell].nsurf,
                   x[0],x[1],x[2],xnew[0],xnew[1],xnew[2],
                   icell,cells[icell].id,
                   lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],dtremain);
        }
        if (DIM == 2) {
          if (ntimestep == MOVE_DEBUG_STEP &&
              (MOVE_DEBUG_ID == particles[i].id ||
               (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
            printf("PARTICLE %d %ld: %d %d: %d: x %g %g: xnew %g %g: %d "
                   CELLINT_FORMAT ": lo %g %g: hi %g %g: DTR: %g\n",
                   me,update->ntimestep,i,particles[i].id,
                   cells[icell].nsurf,
                   x[0],x[1],xnew[0],xnew[1],
                   icell,cells[icell].id,
                   lo[0],lo[1],hi[0],hi[1],dtremain);
        }
        if (DIM == 1) {
          if (ntimestep == MOVE_DEBUG_STEP &&
              (MOVE_DEBUG_ID == particles[i].id ||
               (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
            printf("PARTICLE %d %ld: %d %d: %d: x %g %g: xnew %g %g: %d "
                   CELLINT_FORMAT ": lo %g %g: hi %g %g: DTR: %g\n",
                   me,update->ntimestep,i,particles[i].id,
                   cells[icell].nsurf,
                   x[0],x[1],xnew[0],sqrt(xnew[1]*xnew[1]+xnew[2]*xnew[2]),
                   icell,cells[icell].id,
                   lo[0],lo[1],hi[0],hi[1],dtremain);
        }
#endif

        // check if particle crosses any cell face
        // frac = fraction of move completed before hitting cell face
        // this section should be as efficient as possible,
        //   since most particles won't do anything else
        // axisymmetric y cell face crossings:
        //   these faces are curved cylindrical shells
        //   axi_horizontal_line() checks for intersection of
        //     straight-line y,z move with circle in y,z
        //   always check move against lower y face
        //     except when particle starts on face and
        //     PEXIT is set (just received) or particle is moving downward in y
        //   only check move against upper y face
        //     if remapped final y position (rnew) is within cell,
        //     or except when particle starts on face and
        //     PEXIT is set (just received) or particle is moving upward in y
        //   unset pflag so not checked again for this particle

        outface = INTERIOR;
        frac = 1.0;

        if (xnew[0] < lo[0]) {
          frac = (lo[0]-x[0]) / (xnew[0]-x[0]);
          outface = XLO;
        } else if (xnew[0] >= hi[0]) {
          frac = (hi[0]-x[0]) / (xnew[0]-x[0]);
          outface = XHI;
        }

        if (DIM != 1) {
          if (xnew[1] < lo[1]) {
            newfrac = (lo[1]-x[1]) / (xnew[1]-x[1]);
            if (newfrac < frac) {
              frac = newfrac;
              outface = YLO;
            }
          } else if (xnew[1] >= hi[1]) {
            newfrac = (hi[1]-x[1]) / (xnew[1]-x[1]);
            if (newfrac < frac) {
              frac = newfrac;
              outface = YHI;
            }
          }
        }

        if (DIM == 1) {
          if (x[1] == lo[1] && (pflag == PEXIT || v[1] < 0.0)) {
            frac = 0.0;
            outface = YLO;
          } else if (Geometry::
                     axi_horizontal_line(dtremain,x,v,lo[1],itmp,tc,tmp)) {
            newfrac = tc/dtremain;
            if (newfrac < frac) {
              frac = newfrac;
              outface = YLO;
            }
          }

          if (x[1] == hi[1] && (pflag == PEXIT || v[1] > 0.0)) {
            frac = 0.0;
            outface = YHI;
          } else {
            rnew = sqrt(xnew[1]*xnew[1] + xnew[2]*xnew[2]);
            if (rnew >= hi[1]) {
              if (Geometry::
                  axi_horizontal_line(dtremain,x,v,hi[1],itmp,tc,tmp)) {
                newfrac = tc/dtremain;
                if (newfrac < frac) {
                  frac = newfrac;
                  outface = YHI;
                }
              }
            }
          }

          pflag = 0;
        }

        if (DIM == 3) {
          if (xnew[2] < lo[2]) {
            newfrac = (lo[2]-x[2]) / (xnew[2]-x[2]);
            if (newfrac < frac) {
              frac = newfrac;
              outface = ZLO;
            }
          } else if (xnew[2] >= hi[2]) {
            newfrac = (hi[2]-x[2]) / (xnew[2]-x[2]);
            if (newfrac < frac) {
              frac = newfrac;
              outface = ZHI;
            }
          }
        }

        //if (iterate == 10) exit(1);
        //iterate++;

#ifdef MOVE_DEBUG
        if (ntimestep == MOVE_DEBUG_STEP &&
            (MOVE_DEBUG_ID == particles[i].id ||
             (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX))) {
          if (outface != INTERIOR)
            printf("  OUTFACE %d out: %d %d, frac %g\n",
                   outface,grid->neigh_decode(nmask,outface),
                   neigh[outface],frac);
          else
            printf("  INTERIOR %d %d\n",outface,INTERIOR);
        }
#endif

        // START of code specific to surfaces

        if (SURF) {

          // skip surf checks if particle flagged as EXITing this cell
          // then unset pflag so not checked again for this particle

          nsurf = cells[icell].nsurf;
          if (pflag == PEXIT) {
            nsurf = 0;
            pflag = 0;
          }
          nscheck_one += nsurf;

          if (nsurf) {

            // particle crosses cell face, reset xnew exactly on face of cell
            // so surface check occurs only for particle path within grid cell
            // xhold = saved xnew so can restore below if no surf collision

            if (outface != INTERIOR) {
              xhold[0] = xnew[0];
              xhold[1] = xnew[1];
              if (DIM != 2) xhold[2] = xnew[2];

              xnew[0] = x[0] + frac*(xnew[0]-x[0]);
              xnew[1] = x[1] + frac*(xnew[1]-x[1]);
              if (DIM != 2) xnew[2] = x[2] + frac*(xnew[2]-x[2]);

              if (outface == XLO) xnew[0] = lo[0];
              else if (outface == XHI) xnew[0] = hi[0];
              else if (outface == YLO) xnew[1] = lo[1];
              else if (outface == YHI) xnew[1] = hi[1];
              else if (outface == ZLO) xnew[2] = lo[2];
              else if (outface == ZHI) xnew[2] = hi[2];
            }

            // for axisymmetric, dtsurf = time that particle stays in cell
            // used as arg to axi_line_intersect()

            if (DIM == 1) {
              if (outface == INTERIOR) dtsurf = dtremain;
              else dtsurf = dtremain * frac;
            }

            // check for collisions with triangles or lines in cell
            // find 1st surface hit via minparam
            // skip collisions with previous surf, but not for axisymmetric
            // not considered collision if 2 params are tied and one INSIDE surf
            // if collision occurs, perform collision with surface model
            // reset x,v,xnew,dtremain and continue single particle trajectory

            cflag = 0;
            minparam = 2.0;
            csurfs = cells[icell].csurfs;

            for (m = 0; m < nsurf; m++) {
              isurf = csurfs[m];

              if (DIM > 1) {
                if (isurf == exclude) continue;
              }
              if (DIM == 3) {
                tri = &tris[isurf];
                hitflag = Geometry::
                  line_tri_intersect(x,xnew,tri->p1,tri->p2,tri->p3,
                                     tri->norm,xc,param,side);
              }
              if (DIM == 2) {
                line = &lines[isurf];
                hitflag = Geometry::
                  line_line_intersect(x,xnew,line->p1,line->p2,
                                      line->norm,xc,param,side);
              }
              if (DIM == 1) {
                line = &lines[isurf];
                hitflag = Geometry::
                  axi_line_intersect(dtsurf,x,v,outface,lo,hi,line->p1,line->p2,
                                     line->norm,exclude == isurf,
                                     xc,vc,param,side);
              }

#ifdef MOVE_DEBUG
              if (DIM == 3) {
                if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("SURF COLLIDE: %d %d %d %d: "
                         "P1 %g %g %g: P2 %g %g %g: "
                         "T1 %g %g %g: T2 %g %g %g: T3 %g %g %g: "
                         "TN %g %g %g: XC %g %g %g: "
                         "Param %g: Side %d\n",
                         MOVE_DEBUG_INDEX,icell,nsurf,isurf,
                         x[0],x[1],x[2],xnew[0],xnew[1],xnew[2],
                         tri->p1[0],tri->p1[1],tri->p1[2],
                         tri->p2[0],tri->p2[1],tri->p2[2],
                         tri->p3[0],tri->p3[1],tri->p3[2],
                         tri->norm[0],tri->norm[1],tri->norm[2],
                         xc[0],xc[1],xc[2],param,side);
              }
              if (DIM == 2) {
                if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("SURF COLLIDE: %d %d %d %d: P1 %g %g: P2 %g %g: "
                         "L1 %g %g: L2 %g %g: LN %g %g: XC %g %g: "
                         "Param %g: Side %d\n",
                         MOVE_DEBUG_INDEX,icell,nsurf,isurf,
                         x[0],x[1],xnew[0],xnew[1],
                         line->p1[0],line->p1[1],line->p2[0],line->p2[1],
                         line->norm[0],line->norm[1],
                         xc[0],xc[1],param,side);
              }
              if (DIM == 1) {
                if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("SURF COLLIDE %d %ld: %d %d %d %d: P1 %g %g: P2 %g %g: "
                         "L1 %g %g: L2 %g %g: LN %g %g: XC %g %g: "
                         "VC %g %g %g: Param %g: Side %d\n",
                         hitflag,ntimestep,MOVE_DEBUG_INDEX,icell,nsurf,isurf,
                         x[0],x[1],
                         xnew[0],sqrt(xnew[1]*xnew[1]+xnew[2]*xnew[2]),
                         line->p1[0],line->p1[1],line->p2[0],line->p2[1],
                         line->norm[0],line->norm[1],
                         xc[0],xc[1],vc[0],vc[1],vc[2],param,side);
                double edge1[3],edge2[3],xfinal[3],cross[3];
                MathExtra::sub3(line->p2,line->p1,edge1);
                MathExtra::sub3(x,line->p1,edge2);
                MathExtra::cross3(edge2,edge1,cross);
                if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                    MOVE_DEBUG_ID == particles[i].id)
                  printf("CROSSSTART %g %g %g\n",cross[0],cross[1],cross[2]);
                xfinal[0] = xnew[0];
                xfinal[1] = sqrt(xnew[1]*xnew[1]+xnew[2]*xnew[2]);
                xfinal[2] = 0.0;
                MathExtra::sub3(xfinal,line->p1,edge2);
                MathExtra::cross3(edge2,edge1,cross);
                if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                    MOVE_DEBUG_ID == particles[i].id)
                  printf("CROSSFINAL %g %g %g\n",cross[0],cross[1],cross[2]);
              }
#endif

              if (hitflag && param < minparam && side == OUTSIDE) {
                cflag = 1;
                minparam = param;
                minside = side;
                minsurf = isurf;
                minxc[0] = xc[0];
                minxc[1] = xc[1];
                if (DIM == 3) minxc[2] = xc[2];
                if (DIM == 1) {
                  minvc[1] = vc[1];
                  minvc[2] = vc[2];
                }
              }

            } // END of for loop over surfs

            // tri/line = surf that particle hit first

            if (cflag) {
              if (DIM == 3) tri = &tris[minsurf];
              if (DIM != 3) line = &lines[minsurf];

              // set x to collision point
              // if axisymmetric, set v to remapped velocity at collision pt

              x[0] = minxc[0];
              x[1] = minxc[1];
              if (DIM == 3) x[2] = minxc[2];
              if (DIM == 1) {
                v[1] = minvc[1];
                v[2] = minvc[2];
              }

              // perform surface collision using surface collision model
              // surface chemistry may destroy particle or create new one
              // must update particle's icell to current icell so that
              //   if jpart is created, it will be added to correct cell
              // if jpart, add new particle to this iteration via pstop++
              // tally surface statistics if requested using iorig

              ipart = &particles[i];
              ipart->icell = icell;
              dtremain *= 1.0 - minparam*frac;

              if (nsurf_tally)
                memcpy(&iorig,&particles[i],sizeof(Particle::OnePart));

              if (DIM == 3)
                jpart = surf->sc[tri->isc]->
                  collide(ipart,dtremain,minsurf,tri->norm,tri->isr,reaction);
              if (DIM != 3)
                jpart = surf->sc[line->isc]->
                  collide(ipart,dtremain,minsurf,line->norm,line->isr,reaction);

              if (jpart) {
                particles = particle->particles;
                x = particles[i].x;
                v = particles[i].v;
                jpart->flag = PSURF + 1 + minsurf;
                jpart->dtremain = dtremain;
                jpart->weight = particles[i].weight;
                pstop++;
              }

              if (nsurf_tally)
                for (m = 0; m < nsurf_tally; m++)
                  slist_active[m]->surf_tally(minsurf,icell,reaction,
                                              &iorig,ipart,jpart);

              // stuck_iterate = consecutive iterations particle is immobile

              if (minparam == 0.0) stuck_iterate++;
              else stuck_iterate = 0;

              // reset post-bounce xnew

              xnew[0] = x[0] + dtremain*v[0];
              xnew[1] = x[1] + dtremain*v[1];
              if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];

              exclude = minsurf;
              nscollide_one++;

#ifdef MOVE_DEBUG
              if (DIM == 3) {
                if (ntimestep == MOVE_DEBUG_STEP &&
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("POST COLLISION %d: %g %g %g: %g %g %g: %g %g %g\n",
                         MOVE_DEBUG_INDEX,
                         x[0],x[1],x[2],xnew[0],xnew[1],xnew[2],
                         minparam,frac,dtremain);
              }
              if (DIM == 2) {
                if (ntimestep == MOVE_DEBUG_STEP &&
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("POST COLLISION %d: %g %g: %g %g: %g %g %g\n",
                         MOVE_DEBUG_INDEX,
                         x[0],x[1],xnew[0],xnew[1],
                         minparam,frac,dtremain);
              }
              if (DIM == 1) {
                if (ntimestep == MOVE_DEBUG_STEP &&
                    (MOVE_DEBUG_ID == particles[i].id ||
                     (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
                  printf("POST COLLISION %d: %g %g: %g %g: vel %g %g %g: %g %g %g\n",
                         MOVE_DEBUG_INDEX,
                         x[0],x[1],
                         xnew[0],sqrt(xnew[1]*xnew[1]+xnew[2]*xnew[2]),
                         v[0],v[1],v[2],
                         minparam,frac,dtremain);
              }
#endif

              // if ipart = NULL, particle discarded due to surface chem
              // else if particle not stuck, continue advection while loop
              // if stuck, mark for DISCARD, and drop out of SURF code

              if (ipart == NULL) particles[i].flag = PDISCARD;
              else if (stuck_iterate < MAXSTUCK) continue;
              else {
                particles[i].flag = PDISCARD;
                nstuck++;
              }

            } // END of cflag if section that performed collision

            // no collision, so restore saved xnew if changed it above

            if (outface != INTERIOR) {
              xnew[0] = xhold[0];
              xnew[1] = xhold[1];
              if (DIM != 2) xnew[2] = xhold[2];
            }

          } // END of if test for any surfs in this cell
        } // END of code specific to surfaces

        // break from advection loop if discarding particle

        if (particles[i].flag == PDISCARD) break;

        // no cell crossing and no surface collision
        // set final particle position to xnew, then break from advection loop
        // for axisymmetry, must first remap linear xnew and v
        // for axisymmetry, check if final particle position is within cell
        //   can be rare epsilon round-off cases where particle ends up outside
        //     of final cell curved surf when move logic thinks it is inside
        //   example is when Geom::axi_horizontal_line() says no crossing of cell edge
        //     but axi_remap() puts particle outside the cell
        //   in this case, just DISCARD particle and tally it to naxibad
        // if migrating to another proc,
        //   flag as PDONE so new proc won't move it more on this step

        if (outface == INTERIOR) {
          if (DIM == 1) axi_remap(xnew,v);
          x[0] = xnew[0];
          x[1] = xnew[1];
          if (DIM == 3) x[2] = xnew[2];
          if (DIM == 1) {
            if (x[1] < lo[1] || x[1] > hi[1]) {
              particles[i].flag = PDISCARD;
              naxibad++;
              break;
            }
          }
          if (cells[icell].proc != me) particles[i].flag = PDONE;
          break;
        }

        // particle crosses cell face
        // decrement dtremain in case particle is passed to another proc
        // for axisymmetry, must then remap linear x and v
        // reset particle x to be exactly on cell face
        // for axisymmetry, must reset xnew for next iteration since v changed

        dtremain *= 1.0-frac;
        exclude = -1;

        x[0] += frac * (xnew[0]-x[0]);
        x[1] += frac * (xnew[1]-x[1]);
        if (DIM != 2) x[2] += frac * (xnew[2]-x[2]);
        if (DIM == 1) axi_remap(x,v);

        if (outface == XLO) x[0] = lo[0];
        else if (outface == XHI) x[0] = hi[0];
        else if (outface == YLO) x[1] = lo[1];
        else if (outface == YHI) x[1] = hi[1];
        else if (outface == ZLO) x[2] = lo[2];
        else if (outface == ZHI) x[2] = hi[2];

        if (DIM == 1) {
          xnew[0] = x[0] + dtremain*v[0];
          xnew[1] = x[1] + dtremain*v[1];
          xnew[2] = x[2] + dtremain*v[2];
        }

        // nflag = type of neighbor cell: child, parent, unknown, boundary
        // if parent, use id_find_child to identify child cell
        //   result can be -1 for unknown cell, occurs when:
        //   (a) particle hits face of ghost child cell
        //   (b) the ghost cell extends beyond ghost halo
        //   (c) cell on other side of face is a parent
        //   (d) its child, which the particle is in, is entirely beyond my halo
        // if new cell is child and surfs exist, check if a split cell

        nflag = grid->neigh_decode(nmask,outface);
        icell_original = icell;

        if (nflag == NCHILD) {
          icell = neigh[outface];
          if (DIM == 3 && SURF) {
            if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
              icell = split3d(icell,x);
          }
          if (DIM < 3 && SURF) {
            if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
              icell = split2d(icell,x);
          }
        } else if (nflag == NPARENT) {
          pcell = &pcells[neigh[outface]];
          icell = grid->id_find_child(pcell->id,cells[icell].level,
                                      pcell->lo,pcell->hi,x);
          if (icell >= 0) {
            if (DIM == 3 && SURF) {
              if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                icell = split3d(icell,x);
            }
            if (DIM < 3 && SURF) {
              if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                icell = split2d(icell,x);
            }
          }
        } else if (nflag == NUNKNOWN) icell = -1;

        // neighbor cell is global boundary
        // tally boundary stats if requested using iorig
        // collide() updates x,v,xnew as needed due to boundary interaction
        //   may also update dtremain (piston BC)
        // for axisymmetric, must recalculate xnew since v may have changed
        // surface chemistry may destroy particle or create new one
        // if jpart, add new particle to this iteration via pstop++
        // OUTFLOW: exit with particle flag = PDISCARD
        // PERIODIC: new cell via same logic as above for child/parent/unknown
        // OTHER: reflected particle stays in same grid cell

        else {
          ipart = &particles[i];

          if (nboundary_tally)
            memcpy(&iorig,&particles[i],sizeof(Particle::OnePart));

          bflag = domain->collide(ipart,outface,icell,xnew,dtremain,
                                  jpart,reaction);

          if (jpart) {
            particles = particle->particles;
            x = particles[i].x;
            v = particles[i].v;
          }

          if (nboundary_tally)
            for (m = 0; m < nboundary_tally; m++)
              blist_active[m]->
                boundary_tally(outface,bflag,reaction,&iorig,ipart,jpart);

          if (DIM == 1) {
            xnew[0] = x[0] + dtremain*v[0];
            xnew[1] = x[1] + dtremain*v[1];
            xnew[2] = x[2] + dtremain*v[2];
          }

          if (bflag == OUTFLOW) {
            particles[i].flag = PDISCARD;
            nexit_one++;
            break;

          } else if (bflag == PERIODIC) {
            if (nflag == NPBCHILD) {
              icell = neigh[outface];
              if (DIM == 3 && SURF) {
                if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                  icell = split3d(icell,x);
              }
              if (DIM < 3 && SURF) {
                if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                  icell = split2d(icell,x);
              }
            } else if (nflag == NPBPARENT) {
              pcell = &pcells[neigh[outface]];
              icell = grid->id_find_child(pcell->id,cells[icell].level,
                                          pcell->lo,pcell->hi,x);
              if (icell >= 0) {
                if (DIM == 3 && SURF) {
                  if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                    icell = split3d(icell,x);
                }
                if (DIM < 3 && SURF) {
                  if (cells[icell].nsplit > 1 && cells[icell].nsurf >= 0)
                    icell = split2d(icell,x);
                }
              } else domain->uncollide(outface,x);
            } else if (nflag == NPBUNKNOWN) {
              icell = -1;
              domain->uncollide(outface,x);
            }

          } else if (bflag == SURFACE) {
            if (ipart == NULL) {
              particles[i].flag = PDISCARD;
              break;
            } else if (jpart) {
              jpart->flag = PSURF;
              jpart->dtremain = dtremain;
              jpart->weight = particles[i].weight;
              pstop++;
            }
            nboundary_one++;
            ntouch_one--;    // decrement here since will increment below

          } else {
            nboundary_one++;
            ntouch_one--;    // decrement here since will increment below
          }
        }

        // neighbor cell is unknown
        // reset icell to original icell which must be a ghost cell
        // exit with particle flag = PEXIT, so receiver can identify neighbor

        if (icell < 0) {
          icell = icell_original;
          particles[i].flag = PEXIT;
          particles[i].dtremain = dtremain;
          entryexit = 1;
          break;
        }

        // if nsurf < 0, new cell is EMPTY ghost
        // exit with particle flag = PENTRY, so receiver can continue move

        if (cells[icell].nsurf < 0) {
          particles[i].flag = PENTRY;
          particles[i].dtremain = dtremain;
          entryexit = 1;
          break;
        }

        // move particle into new grid cell for next stage of move

        lo = cells[icell].lo;
        hi = cells[icell].hi;
        neigh = cells[icell].neigh;
        nmask = cells[icell].nmask;
        ntouch_one++;
      }

      // END of while loop over advection of single particle

#ifdef MOVE_DEBUG
      if (ntimestep == MOVE_DEBUG_STEP &&
          (MOVE_DEBUG_ID == particles[i].id ||
           (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
        printf("MOVE DONE %d %d %d: %g %g %g: DTR %g\n",
               MOVE_DEBUG_INDEX,particles[i].flag,icell,
               x[0],x[1],x[2],dtremain);
#endif

      // move is complete, or as much as can be done on this proc
      // update particle's grid cell
      // if particle flag set, add particle to migrate list
      // if discarding, migration will delete particle

      particles[i].icell = icell;

      if (particles[i].flag != PKEEP) {
        mlist[nmigrate++] = i;
        if (particles[i].flag != PDISCARD) {
          if (cells[icell].proc == me) {
            char str[128];
            sprintf(str,
                    "Particle %d on proc %d being sent to self "
                    "on step " BIGINT_FORMAT,
                    i,me,update->ntimestep);
            error->one(FLERR,str);
          }
          ncomm_one++;
        }
      }
    }

    // END of pstart/pstop loop advecting all particles

    // if gridcut >= 0.0, check if another iteration of move is required
    // only the case if some particle flag = PENTRY/PEXIT
    //   in which case perform particle migration
    // if not, move is done and final particle comm will occur in run()
    // if iterating, reset pstart/pstop and extend migration list if necessary

    if (grid->cutoff < 0.0) break;

    timer->stamp(TIME_MOVE);
    MPI_Allreduce(&entryexit,&any_entryexit,1,MPI_INT,MPI_MAX,world);
    timer->stamp();

    if (any_entryexit) {
      timer->stamp(TIME_MOVE);
      pstart = comm->migrate_particles(nmigrate,mlist);
      timer->stamp(TIME_COMM);
      pstop = particle->nlocal;
      if (pstop-pstart > maxmigrate) {
        maxmigrate = pstop-pstart;
        memory->destroy(mlist);
        memory->create(mlist,maxmigrate,"particle:mlist");
      }
    } else break;

    // END of single move/migrate iteration

  }

  // END of all move/migrate iterations

  particle->sorted = 0;

  // accumulate running totals

  niterate_running += niterate;
  nmove_running += particle->nlocal;
  ntouch_running += ntouch_one;
  ncomm_running += ncomm_one;
  nboundary_running += nboundary_one;
  nexit_running += nexit_one;
  nscheck_running += nscheck_one;
  nscollide_running += nscollide_one;
  surf->nreact_running += surf->nreact_one;
}

/* ----------------------------------------------------------------------
   calculate motion perturbation for a single particle I
     due to external per particle field
   array in fix[ifieldfix] stores per particle perturbations for x and v
------------------------------------------------------------------------- */

void Update::field_per_particle(int i, int icell, double dt, double *x, double *v)
{
  double dtsq = 0.5*dt*dt;
  double **array = modify->fix[ifieldfix]->array_particle;

  int icol = 0;
  if (field_active[0]) {
    x[0] += dtsq*array[i][icol];
    v[0] += dt*array[i][icol];
    icol++;
  }
  if (field_active[1]) {
    x[1] += dtsq*array[i][icol];
    v[1] += dt*array[i][icol];
    icol++;
  }
  if (field_active[2]) {
    x[2] += dtsq*array[i][icol];
    v[2] += dt*array[i][icol];
    icol++;
  }
};

/* ----------------------------------------------------------------------
   calculate motion perturbation for a single particle I in grid cell Icell
     due to external per grid cell field
   array in fix[ifieldfix] stores per grid cell perturbations for x and v
------------------------------------------------------------------------- */

void Update::field_per_grid(int i, int icell, double dt, double *x, double *v)
{
  double dtsq = 0.5*dt*dt;
  double **array = modify->fix[ifieldfix]->array_grid;

  int icol = 0;
  if (field_active[0]) {
    x[0] += dtsq*array[icell][icol];
    v[0] += dt*array[icell][icol];
    icol++;
  }
  if (field_active[1]) {
    x[1] += dtsq*array[icell][icol];
    v[1] += dt*array[icell][icol];
    icol++;
  }
  if (field_active[2]) {
    x[2] += dtsq*array[icell][icol];
    v[2] += dt*array[icell][icol];
    icol++;
  }
};

/* ----------------------------------------------------------------------
   particle is entering split parent icell at x
   determine which split child cell it is in
   return index of sub-cell in ChildCell
------------------------------------------------------------------------- */

int Update::split3d(int icell, double *x)
{
  int m,cflag,isurf,hitflag,side,minsurfindex;
  double param,minparam;
  double xc[3];
  Surf::Tri *tri;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  Surf::Tri *tris = surf->tris;

  // check for collisions with lines in cell
  // find 1st surface hit via minparam
  // only consider tris that are mapped via csplits to a split cell
  //   unmapped tris only touch cell surf at xnew
  //   another mapped tri should include same xnew
  // NOTE: these next 2 lines do not seem correct compared to code
  // not considered a collision if particles starts on surf, moving out
  // not considered a collision if 2 params are tied and one is INSIDE surf

  int nsurf = cells[icell].nsurf;
  surfint *csurfs = cells[icell].csurfs;
  int isplit = cells[icell].isplit;
  int *csplits = sinfo[isplit].csplits;
  double *xnew = sinfo[isplit].xsplit;

  cflag = 0;
  minparam = 2.0;

  for (m = 0; m < nsurf; m++) {
    if (csplits[m] < 0) continue;
    isurf = csurfs[m];
    tri = &tris[isurf];
    hitflag = Geometry::
      line_tri_intersect(x,xnew,tri->p1,tri->p2,tri->p3,
                         tri->norm,xc,param,side);

    if (hitflag && side != INSIDE && param < minparam) {
      cflag = 1;
      minparam = param;
      minsurfindex = m;
    }
  }

  if (!cflag) return sinfo[isplit].csubs[sinfo[isplit].xsub];
  int index = csplits[minsurfindex];
  return sinfo[isplit].csubs[index];
}

/* ----------------------------------------------------------------------
   particle is entering split ICELL at X
   determine which split sub-cell it is in
   return index of sub-cell in ChildCell
------------------------------------------------------------------------- */

int Update::split2d(int icell, double *x)
{
  int m,cflag,isurf,hitflag,side,minsurfindex;
  double param,minparam;
  double xc[3];
  Surf::Line *line;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  Surf::Line *lines = surf->lines;

  // check for collisions with lines in cell
  // find 1st surface hit via minparam
  // only consider lines that are mapped via csplits to a split cell
  //   unmapped lines only touch cell surf at xnew
  //   another mapped line should include same xnew
  // NOTE: these next 2 lines do not seem correct compared to code
  // not considered a collision if particle starts on surf, moving out
  // not considered a collision if 2 params are tied and one is INSIDE surf

  int nsurf = cells[icell].nsurf;
  surfint *csurfs = cells[icell].csurfs;
  int isplit = cells[icell].isplit;
  int *csplits = sinfo[isplit].csplits;
  double *xnew = sinfo[isplit].xsplit;

  cflag = 0;
  minparam = 2.0;
  for (m = 0; m < nsurf; m++) {
    if (csplits[m] < 0) continue;
    isurf = csurfs[m];
    line = &lines[isurf];
    hitflag = Geometry::
      line_line_intersect(x,xnew,line->p1,line->p2,line->norm,xc,param,side);

    if (hitflag && side != INSIDE && param < minparam) {
      cflag = 1;
      minparam = param;
      minsurfindex = m;
    }
  }

  if (!cflag) return sinfo[isplit].csubs[sinfo[isplit].xsub];
  int index = csplits[minsurfindex];
  return sinfo[isplit].csubs[index];
}

/* ----------------------------------------------------------------------
   check if any surface collision or reaction models are defined
   return 1 if there are any, 0 if not
------------------------------------------------------------------------- */

int Update::collide_react_setup()
{
  nsc = surf->nsc;
  sc = surf->sc;
  nsr = surf->nsr;
  sr = surf->sr;

  if (nsc || nsr) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   zero counters for tallying surface collisions/reactions
   done at start of each timestep
   done within individual SurfCollide and SurfReact instances
------------------------------------------------------------------------- */

void Update::collide_react_reset()
{
  for (int i = 0; i < nsc; i++) sc[i]->tally_reset();
  for (int i = 0; i < nsr; i++) sr[i]->tally_reset();
}

/* ----------------------------------------------------------------------
   update cummulative counters for tallying surface collisions/reactions
   done at end of each timestep
   this is done within individual SurfCollide and SurfReact instances
------------------------------------------------------------------------- */

void Update::collide_react_update()
{
  for (int i = 0; i < nsc; i++) sc[i]->tally_update();
  for (int i = 0; i < nsr; i++) sr[i]->tally_update();
}

/* ----------------------------------------------------------------------
   setup lists of all computes that tally surface and boundary bounce info
   return 1 if there are any, 0 if not
------------------------------------------------------------------------- */

int Update::bounce_setup()
{
  delete [] slist_compute;
  delete [] blist_compute;
  delete [] slist_active;
  delete [] blist_active;
  slist_compute = blist_compute = NULL;

  nslist_compute = nblist_compute = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->surf_tally_flag) nslist_compute++;
    if (modify->compute[i]->boundary_tally_flag) nblist_compute++;
  }

  if (nslist_compute) slist_compute = new Compute*[nslist_compute];
  if (nblist_compute) blist_compute = new Compute*[nblist_compute];
  if (nslist_compute) slist_active = new Compute*[nslist_compute];
  if (nblist_compute) blist_active = new Compute*[nblist_compute];

  nslist_compute = nblist_compute = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->surf_tally_flag)
      slist_compute[nslist_compute++] = modify->compute[i];
    if (modify->compute[i]->boundary_tally_flag)
      blist_compute[nblist_compute++] = modify->compute[i];
  }

  if (nslist_compute || nblist_compute) return 1;
  nsurf_tally = nboundary_tally = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   set bounce tally flags for current timestep
   nsurf_tally = # of surface computes needing bounce info on this step
   nboundary_tally = # of boundary computes needing bounce info on this step
   clear accumulators in computes that will be invoked this step
------------------------------------------------------------------------- */

void Update::bounce_set(bigint ntimestep)
{
  int i;

  nsurf_tally = 0;
  if (nslist_compute) {
    for (i = 0; i < nslist_compute; i++)
      if (slist_compute[i]->matchstep(ntimestep)) {
        slist_active[nsurf_tally++] = slist_compute[i];
        slist_compute[i]->clear();
      }
  }

  nboundary_tally = 0;
  if (nblist_compute) {
    for (i = 0; i < nblist_compute; i++)
      if (blist_compute[i]->matchstep(ntimestep)) {
        blist_active[nboundary_tally++] = blist_compute[i];
        blist_compute[i]->clear();
      }
  }
}

/* ----------------------------------------------------------------------
   make list of classes that reset dynamic parameters
   currently only surf collision models
------------------------------------------------------------------------- */

void Update::dynamic_setup()
{
  delete [] ulist_surfcollide;
  ulist_surfcollide = NULL;

  nulist_surfcollide = 0;
  for (int i = 0; i < surf->nsc; i++)
    if (surf->sc[i]->dynamicflag) nulist_surfcollide++;

  if (nulist_surfcollide)
    ulist_surfcollide = new SurfCollide*[nulist_surfcollide];

  nulist_surfcollide = 0;
  for (int i = 0; i < surf->nsc; i++)
    if (surf->sc[i]->dynamicflag)
      ulist_surfcollide[nulist_surfcollide++] = surf->sc[i];

  if (nulist_surfcollide) dynamic = 1;
}

/* ----------------------------------------------------------------------
   invoke class methods that reset dynamic parameters
------------------------------------------------------------------------- */

void Update::dynamic_update()
{
  if (nulist_surfcollide) {
    for (int i = 0; i < nulist_surfcollide; i++)
      ulist_surfcollide[i]->dynamic();
  }
}

/* ----------------------------------------------------------------------
   set global properites via global command in input script
------------------------------------------------------------------------- */

void Update::global(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal global command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"fnum") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      fnum = input->numeric(FLERR,arg[iarg+1]);
      if (fnum <= 0.0) error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"optmove") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"yes") == 0) optmove_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) optmove_flag = 0;
      else error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nrho") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      nrho = input->numeric(FLERR,arg[iarg+1]);
      if (nrho <= 0.0) error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"vstream") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal global command");
      vstream[0] = input->numeric(FLERR,arg[iarg+1]);
      vstream[1] = input->numeric(FLERR,arg[iarg+2]);
      vstream[2] = input->numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      temp_thermal = input->numeric(FLERR,arg[iarg+1]);
      if (temp_thermal <= 0.0) error->all(FLERR,"Illegal global command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"field") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"none") == 0) {
        fstyle = NOFIELD;
        iarg += 2;
      } else if (strcmp(arg[iarg+1],"constant") == 0) {
        if (iarg+6 > narg) error->all(FLERR,"Illegal global field command");
        fstyle = CFIELD;
        double fmag = input->numeric(FLERR,arg[iarg+2]);
        field[0] = input->numeric(FLERR,arg[iarg+3]);
        field[1] = input->numeric(FLERR,arg[iarg+4]);
        field[2] = input->numeric(FLERR,arg[iarg+5]);
        if (fmag <= 0.0) error->all(FLERR,"Illegal global field command");
        if (field[0] == 0.0 && field[1] == 0.0 && field[2] == 0.0)
          error->all(FLERR,"Illegal global field command");
        MathExtra::snorm3(fmag,field);
        iarg += 6;
      } else if (strcmp(arg[iarg+1],"particle") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal global field command");
        delete [] fieldID;
        fstyle = PFIELD;
        int n = strlen(arg[iarg+2]) + 1;
        fieldID = new char[n];
        strcpy(fieldID,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"grid") == 0) {
        if (iarg+4 > narg) error->all(FLERR,"Illegal global field command");
        delete [] fieldID;
        fstyle = GFIELD;
        int n = strlen(arg[iarg+2]) + 1;
        fieldID = new char[n];
        strcpy(fieldID,arg[iarg+2]);
        fieldfreq = input->inumeric(FLERR,arg[iarg+3]);
        if (fieldfreq < 0) error->all(FLERR,"Illegal global field command");
        iarg += 4;
      } else error->all(FLERR,"Illegal global field command");

    } else if (strcmp(arg[iarg],"surfs") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      surf->global(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"surfgrid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (surf->exist)
        error->all(FLERR,
                   "Cannot set global surfgrid when surfaces already exist");
      if (strcmp(arg[iarg+1],"auto") == 0) grid->surfgrid_algorithm = PERAUTO;
      else if (strcmp(arg[iarg+1],"percell") == 0)
        grid->surfgrid_algorithm = PERCELL;
      else if (strcmp(arg[iarg+1],"persurf") == 0)
        grid->surfgrid_algorithm = PERSURF;
      else error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"surfmax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (surf->exist)
        error->all(FLERR,
                   "Cannot set global surfmax when surfaces already exist");
      grid->maxsurfpercell = atoi(arg[iarg+1]);
      if (grid->maxsurfpercell <= 0) error->all(FLERR,"Illegal global command");
      // reallocate paged data structs for variable-length surf info
      grid->allocate_surf_arrays();
      iarg += 2;
    } else if (strcmp(arg[iarg],"splitmax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (surf->exist)
        error->all(FLERR,
                   "Cannot set global splitmax when surfaces already exist");
      grid->maxsplitpercell = atoi(arg[iarg+1]);
      if (grid->maxsplitpercell <= 0) error->all(FLERR,"Illegal global command");
      // reallocate paged data structs for variable-length cell info
      grid->allocate_surf_arrays();
      iarg += 2;
    } else if (strcmp(arg[iarg],"gridcut") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      grid->cutoff = input->numeric(FLERR,arg[iarg+1]);
      if (grid->cutoff < 0.0 && grid->cutoff != -1.0)
        error->all(FLERR,"Illegal global command");
      // force ghost info to be regenerated with new cutoff
      grid->remove_ghosts();
      iarg += 2;
    } else if (strcmp(arg[iarg],"weight") == 0) {
      // for now assume just one arg after "cell"
      // may need to generalize later
      if (iarg+3 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"cell") == 0) grid->weight(1,&arg[iarg+2]);
      else error->all(FLERR,"Illegal weight command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"comm/sort") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"yes") == 0) comm->commsortflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) comm->commsortflag = 0;
      else error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm/style") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"neigh") == 0) comm->commpartstyle = 1;
      else if (strcmp(arg[iarg+1],"all") == 0) comm->commpartstyle = 0;
      else error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"surftally") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"auto") == 0) surf->tally_comm = TALLYAUTO;
      else if (strcmp(arg[iarg+1],"reduce") == 0) surf->tally_comm = TALLYREDUCE;
      else if (strcmp(arg[iarg+1],"rvous") == 0) surf->tally_comm = TALLYRVOUS;
      else error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"particle/reorder") == 0) {
      reorder_period = input->inumeric(FLERR,arg[iarg+1]);
      if (reorder_period < 0) error->all(FLERR,"Illegal global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mem/limit") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal global command");
      if (strcmp(arg[iarg+1],"grid") == 0) mem_limit_grid_flag = 1;
      else {
        double factor = input->numeric(FLERR,arg[iarg+1]);
        bigint global_mem_limit_big = static_cast<bigint> (factor * 1024*1024);
        if (global_mem_limit_big < 0) error->all(FLERR,"Illegal global command");
        if (global_mem_limit_big > MAXSMALLINT)
          error->all(FLERR,"Global mem/limit setting cannot exceed 2GB");
        global_mem_limit = global_mem_limit_big;
      }
      iarg += 2;
    } else error->all(FLERR,"Illegal global command");
  }
}

/* ----------------------------------------------------------------------
   reset timestep as called from input script
------------------------------------------------------------------------- */

void Update::reset_timestep(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal reset_timestep command");
  bigint newstep = ATOBIGINT(arg[0]);
  reset_timestep(newstep);
}

/* ----------------------------------------------------------------------
   reset timestep
   set atimestep to new timestep, so future update_time() calls will be correct
   trigger reset of timestep for output and for fixes that require it
   do not allow any timestep-dependent fixes to be defined
   reset eflag/vflag global so nothing will think eng/virial are current
   reset invoked flags of computes,
     so nothing will think they are current between runs
   clear timestep list of computes that store future invocation times
   called from rerun command and input script (indirectly)
------------------------------------------------------------------------- */

void Update::reset_timestep(bigint newstep)
{
  ntimestep = newstep;
  if (ntimestep < 0) error->all(FLERR,"Timestep must be >= 0");
  if (ntimestep > MAXBIGINT) error->all(FLERR,"Too big a timestep");

  output->reset_timestep(ntimestep);

  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i]->time_depend)
      error->all(FLERR,
                 "Cannot reset timestep with a time-dependent fix defined");
    //modify->fix[i]->reset_timestep(ntimestep);
  }

  for (int i = 0; i < modify->ncompute; i++) {
    modify->compute[i]->invoked_scalar = -1;
    modify->compute[i]->invoked_vector = -1;
    modify->compute[i]->invoked_array = -1;
    modify->compute[i]->invoked_per_particle = -1;
    modify->compute[i]->invoked_per_grid = -1;
    modify->compute[i]->invoked_per_surf = -1;
  }

  for (int i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->timeflag) modify->compute[i]->clearstep();
}

/* ----------------------------------------------------------------------
   get mem/limit based on grid memory
------------------------------------------------------------------------- */

void Update::set_mem_limit_grid(int gnlocal)
{
  if (gnlocal == 0) gnlocal = grid->nlocal;

  bigint global_mem_limit_big = static_cast<bigint> (gnlocal*sizeof(Grid::ChildCell));

  if (global_mem_limit_big > MAXSMALLINT)
    error->one(FLERR,"Global mem/limit setting cannot exceed 2GB");

  global_mem_limit = global_mem_limit_big;
}

/* ----------------------------------------------------------------------
   get mem/limit based on grid memory
------------------------------------------------------------------------- */

int Update::have_mem_limit()
{
  if (mem_limit_grid_flag)
    set_mem_limit_grid();

  int mem_limit_flag = 0;

  if (global_mem_limit > 0 || (mem_limit_grid_flag && !grid->nlocal))
    mem_limit_flag = 1;

  return mem_limit_flag;
}


/* rate */
vector<DataPoint2> ReadAndCacheData(const string& filename) {
    vector<DataPoint2> data;
    ifstream file(filename);
    
    if (!file.is_open()) {
        throw runtime_error("Unable to open the file: " + filename);
    }
    
    string line;
    while (getline(file, line)) {
        DataPoint2 point;
        istringstream iss(line);
        if (!(iss >> point.center_rr >> point.center_zz >> point.ne >> point.te >> point.ti >> point.vflow >> point.br >> point.bz >> point.bt)) {
            cerr << "Error parsing line: " << line << endl;
            continue;
        }
        data.push_back(point);
    }
    
    file.close();
    
    return data;
}

// Function to perform linear interpolation
DataPoint2 LinearInterpolate(const DataPoint2& p1, const DataPoint2& p2, double t) {
    DataPoint2 result;
    double t1 = 1.0 - t;
    result.center_rr = p1.center_rr * t1 + p2.center_rr * t;
    result.center_zz = p1.center_zz * t1 + p2.center_zz * t;
    result.ne = p1.ne * t1 + p2.ne * t;
    result.te = p1.te * t1 + p2.te * t;
    result.ti = p1.ti * t1 + p2.ti * t;
    result.vflow = p1.vflow * t1 + p2.vflow * t;
    result.br = p1.br * t1 + p2.br * t;
    result.bz = p1.bz * t1 + p2.bz * t;
    result.bt = p1.bt * t1 + p2.bt * t;
    return result;
}

/* rate */
PlasmaData Update::readPlasmaData(const std::string& filePath) {

    auto it = plasmaDataCache.find(filePath);
    if (it != plasmaDataCache.end()) {
        return it->second;  // Return cached content
    }

    H5::H5File file(filePath, H5F_ACC_RDONLY);

    auto rawDataTo2DVector = [](const std::vector<float>& rawData, const std::vector<hsize_t>& dims) {
        std::vector<std::vector<float>> data(dims[0], std::vector<float>(dims[1]));
        for (hsize_t i = 0; i < dims[0]; ++i)
            for (hsize_t j = 0; j < dims[1]; ++j)
                data[i][j] = rawData[i * dims[1] + j];
        return data;
    };

    auto read2DDataSet = [&file, &rawDataTo2DVector](const std::string& datasetPath) {
        H5::DataSet ds = file.openDataSet(datasetPath);
        H5::DataSpace space = ds.getSpace();
        std::vector<hsize_t> dims(2);
        space.getSimpleExtentDims(dims.data(), NULL);

        std::vector<float> rawData(dims[0] * dims[1]);
        ds.read(rawData.data(), H5::PredType::NATIVE_FLOAT);

        return rawDataTo2DVector(rawData, dims);
    };

    auto read2DSliceFrom3DDataSet = [&file, &rawDataTo2DVector](const std::string& datasetPath) {
        H5::DataSet ds = file.openDataSet(datasetPath);
        H5::DataSpace space = ds.getSpace();
        std::vector<hsize_t> dims(3);
        space.getSimpleExtentDims(dims.data(), NULL);

        hsize_t offset[3] = {0, 0, dims[2] - 1};
        hsize_t count[3] = {dims[0], dims[1], 1};
        space.selectHyperslab(H5S_SELECT_SET, count, offset);
    
        std::vector<float> rawData(dims[0] * dims[1]);
        H5::DataSpace mem_space(2, count);
        ds.read(rawData.data(), H5::PredType::NATIVE_FLOAT, mem_space, space);

        return rawDataTo2DVector(rawData, {dims[0], dims[1]});
    };

    PlasmaData data;
    data.dens_e = read2DDataSet("n_e/dens");
    data.temp_e = read2DDataSet("n_e/temp");
    data.dens_i = read2DDataSet("n_i/dens");
    data.temp_i = read2DDataSet("n_i/temp");
    data.parr_flow = read2DDataSet("n_i/parr_flow");
    data.b_r = read2DDataSet("bfield/b_r");
    data.b_z = read2DDataSet("bfield/b_z");
    data.b_phi = read2DDataSet("bfield/b_phi");
    data.r = read2DSliceFrom3DDataSet("solps_like/r");
    data.z = read2DSliceFrom3DDataSet("solps_like/z");

    file.close();

    // Cache the content
    plasmaDataCache[filePath] = data;

    return data;
}


// Helper function for bilinear interpolation
float bilinearInterpolation(float alpha, float beta,
                            float q11, float q12, float q21, float q22) {
    return (1 - alpha) * (1 - beta) * q11 +
           alpha * (1 - beta) * q21 +
           (1 - alpha) * beta * q12 +
           alpha * beta * q22;
}

PlasmaParams interpolate(const PlasmaData2& lower, const PlasmaData2& upper, double alpha) {
    PlasmaParams params;
    params.dens_e = lower.ne + alpha * (upper.ne - lower.ne);
    params.temp_e = lower.te + alpha * (upper.te - lower.te);
    params.temp_i = lower.ti + alpha * (upper.ti - lower.ti);
    params.parr_flow = lower.vflow + alpha * (upper.vflow - lower.vflow);
      params.b_r = lower.br + alpha * (upper.br - lower.br);
    params.b_z = lower.bz + alpha * (upper.bz - lower.bz);
    params.b_phi = lower.bt + alpha * (upper.bt - lower.bt);
    return params;
}

PlasmaParams Update::interpolatePlasmaData(float r_val, float z_val) {
    std::pair<float, float> key = {r_val, z_val};
    if (plasma_cache.find(key) != plasma_cache.end()) {
        return plasma_cache[key];
    }

//    std::ifstream dataFile("data.txt");

//     // Skip the header line
//     std::string header;
//     std::getline(dataFile, header);

//     std::vector<PlasmaData2> data2;

//     PlasmaData2 tempData;
//     while (dataFile >> tempData.center_rr >> tempData.center_zz >> tempData.ne >> tempData.te >> tempData.ti
//                    >> tempData.vflow >> tempData.br >> tempData.bz >> tempData.bt) {
//         data2.push_back(tempData);
//     }

//     dataFile.close();

//     PlasmaParams params;

//     // plasma_cache[key] = params;
//     int i_r = -1, j_r = -1, i_z = -1, j_z = -1;

//   // Find i_r based on r_val
//   for (int i = 0; i < data2.size() - 1; i++) {
//       if (data2[i].center_rr <= r_val && r_val < data2[i + 1].center_rr) {
//           i_r = i;
//           break;
//       }
//   }

//   // Find j_z based on z_val
//   for (int j = 0; j < data2.size() - 1; j++) {
//       if (data2[j].center_zz <= z_val && z_val < data2[j + 1].center_zz) {
//           j_z = j;
//           break;
//       }
//   }

// // ...

// // Check if either coordinate was not found
// if (i_r == -1 || j_z == -1) {
//     params = {0, 0, 0, 0, 0, 0, 0, 0};  // Assign values to the existing 'params' variable
// }

// // Now you can perform interpolation using i_r, j_r, i_z, and j_z
// double alpha_r = (r_val - data2[i_r].center_rr) / (data2[i_r + 1].center_rr - data2[i_r].center_rr);
// double alpha_z = (z_val - data2[j_z].center_zz) / (data2[j_z + 1].center_zz - data2[j_z].center_zz);

// // PlasmaParams params;
// params.dens_e = bilinearInterpolation(alpha_r, alpha_z, data2[i_r].ne, data2[i_r + 1].ne, data2[j_z].ne, data2[j_z + 1].ne);
// params.temp_e = bilinearInterpolation(alpha_r, alpha_z, data2[i_r].te, data2[i_r + 1].te, data2[j_z].te, data2[j_z + 1].te);
// params.dens_i = bilinearInterpolation(alpha_r, alpha_z, data2[i_r].ne, data2[i_r + 1].ne, data2[j_z].ne, data2[j_z + 1].ti);
// params.temp_i = bilinearInterpolation(alpha_r, alpha_z, data2[i_r].ti, data2[i_r + 1].ti, data2[j_z].ti, data2[j_z + 1].ti);
// params.b_phi = bilinearInterpolation(alpha_r, alpha_z, data2[i_r].bt, data2[i_r + 1].bt, data2[j_z].bt, data2[j_z + 1].bt);
// params.b_r = bilinearInterpolation(alpha_r, alpha_z, data2[i_r].br, data2[i_r + 1].br, data2[j_z].br, data2[j_z + 1].br);
// params.b_z = bilinearInterpolation(alpha_r, alpha_z, data2[i_r].bz, data2[i_r + 1].bz, data2[j_z].bz, data2[j_z + 1].bz);
// params.parr_flow = bilinearInterpolation(alpha_r, alpha_z, data2[i_r].vflow, data2[i_r + 1].vflow, data2[j_z].vflow, data2[j_z + 1].vflow);

// set manually
// temp_e (float r_val, float z_val) 
//  r_val < 2.9 -.4 < and z < 0.4 temp_e = 40
// else temp_e  decreases exponentially with r and z. starts

// magnetic field 
// br 
// for z> 0 and br increases linearly with z 0.001 to -.2
// for z< 0 and br increases linearly with z 0.001 to .2
// bz is 0.3 at z=0 and decreases linearly with z to 0.1 
// bphi is constant 3
PlasmaParams params;
// if ()
// br
 if (z_val > 0) {
        // For z > 0, br increases linearly with z from 0.001 to -0.2
        params.b_r = -0.2 - 0.001 * z_val;
    } else {
        // For z < 0, br increases linearly with z from 0.001 to 0.2
        params.b_r = 0.2 - 0.001 * z_val;
    }
// bz
if (z_val > 0) {
        // For z > 0, bz decreases linearly with z from 0.3 to 0.1
        params.b_z = 0.3 - 0.2 * z_val;
    } else {
        // For z < 0, bz increases linearly with z from 0.3 to 0.1
        params.b_z = 0.3 + 0.2 * z_val;
    }
// bphi
params.b_phi = 3.0;
if (r_val < 2.9 && r_val > -0.4 && z_val < 0.4) {
        params.temp_e =  40.0f;
        params.dens_e = 0.8e19;
        params.parr_flow = 500.0f;

    } else {
        // Calculate temp_e decreasing exponentially with r and z
        float reference_r = 2.9f;
        float reference_z = 0.4f;
        float initial_temp_e = 40.0f;
        float intial_dens_e = 0.8e19;
        float initial_parr_flow = 500.0f;
        
        
        // Calculate the distance from the reference point (r, z)
        float distance = sqrt((r_val - reference_r) * (r_val - reference_r) + (z_val - reference_z) * (z_val - reference_z));
        
        // Calculate the decrease rate (you can adjust the rate as needed)
        float decrease_rate = 0.1f; // Adjust this rate
        
        // Calculate the new temp_e based on the exponential decrease
        params.temp_e= initial_temp_e * exp(-decrease_rate * distance);
        params.dens_e = intial_dens_e * exp(-decrease_rate * distance);

        // parflow increaaSE
      float  increase_rate  = 0.01f;
      params.parr_flow = -initial_parr_flow * exp(increase_rate * distance);
        
    }
    // params.temp_e = temp_e;
    // params.temp_i = temp_e;
    // params.dens_e = dens_e;
    // params.dens_i = dens_e;
    // params.parr_flow = parr_flow;
    // params.b_r = br;
    // params.b_z = bz;
    // params.b_phi = bphi;


// Cache the interpolated result
plasma_cache[key] = params;

// Return the interpolated result
return params;


}

#include <random>

// Assuming this is inside a class, these could be member variables
std::default_random_engine generator;  // you might want to seed it
//std::normal_distribution<double> distribution(0.0,1.0);

double rand_normal() {
    return distribution(generator);
}

void Update::pusher_boris(double *xc, double *x, double *v, double *xnew, double charge, double mass, double dt, PlasmaParams params)
{
    // Directly compute coefficients that will be reused
    const double q_over_2m = charge * dt / (2.0 * mass);

    // Get Electric Field
    double E[3];
//    auto [eField, angle] = this->getElectricPotential(xc[0], xc[1], charge, mass, params);
    E[1] = 0;
    double E0= 300.73307680233894;
//    printf("E0 %g\n", E0);
    double Debye=1.6619560612438563e-3;
    E[0] = E0 * charge/1.602e-19; //*std::exp(-0.5 * xc[0]/Debye); //eField;
    E[2] = 0;
//    printf("charge E %g %g\n", charge/1.602e-19, E[0]);
    
//    double B[3] = {params.b_r, params.b_z, params.b_phi};
    double B[3] = {0,0,0};
    B[2]= 0;
    
    
    // Half electric field update
    for (int i = 0; i < 3; i++) {
        v[i] += E[i] * q_over_2m;
    }

    // Compute the t-vector for magnetic field rotation
    double t[3];
    double qom = 0.5 * dt * charge / mass; // charge-to-mass ratio times dt/2
    for (int i = 0; i < 3; i++) {
        t[i] = qom * B[i];
    }

    double tsq = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
    
    // Calculate v-prime (v' = v + v x t)
    double vprime[3] = {
        v[1] * t[2] - v[2] * t[1] + v[0],
        v[2] * t[0] - v[0] * t[2] + v[1],
        v[0] * t[1] - v[1] * t[0] + v[2]
    };

    // Calculate s-vector (s = 2t / (1 + t^2))
    double s[3];
    double denom = 1.0 / (1.0 + tsq);
    for (int i = 0; i < 3; i++) {
        s[i] = 2.0 * t[i] * denom;
    }

    // Update velocity using s (v_new = v' + v' x s)
    for (int i = 0; i < 3; i++) {
        v[i] = vprime[1] * s[2] - vprime[2] * s[1] + vprime[0];
        v[i] += E[i] * q_over_2m;  // Combine second half electric field update here
    }

    
//    // Update position
//    for (int i = 0; i < 3; i++) {
//        xnew[i] = x[i] + v[i] * dt;
//    }
//
    
    const double D_perp = domain->t_D_perp;
    double diffusion_term = std::sqrt(2 * D_perp * dt) * (rand_normal() - 0.5);  // Assuming rand_normal() returns a number between 0 and 1
//    printf("diffusion_term %g\n", diffusion_term);
    xnew[0] = x[0] + v[0] * dt + diffusion_term;
    xnew[1] = x[1] + v[1] * dt;
    xnew[2] = x[2] + v[2] * dt;
//    // Update position
//    for (int i = 0; i < 3; i++) {
//        xnew[i] = x[i] + v[i] * dt;
//    }
}


// void Update::pusher_boris(double *xc, double *x, double *v, double *xnew, double charge, double mass, double dt, 
//  PlasmaParams params)
// {
//   double E[3];
//   auto [eField, angle] = this->getElectricPotential(xc[0], xc[1], charge, mass, params);
//   E[1] = eField;
//   E[0] = 0 ;
//   E[2] = 0;

//   double B[3] = {params.b_r, params.b_phi, params.b_z};

//     // Half electric field update
//     for (int i = 0; i < 3; i++) {
//         v[i] += charge * E[i] * dt / (2.0 * mass);
//     }

//     // Magnetic field rotation
//     double qom = 0.5 * dt * charge / mass; // charge-to-mass ratio times dt/2

//     // Compute the t-vector for magnetic field rotation
//     double t[3] = {
//         qom * B[0],
//         qom * B[1],
//         qom * B[2]
//     };

//     // Compute the magnitude of t squared
//     double tsq = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];

//     // Calculate v-prime (v' = v + v x t)
//     double vprime[3] = {
//         v[1] * t[2] - v[2] * t[1] + v[0],
//         v[2] * t[0] - v[0] * t[2] + v[1],
//         v[0] * t[1] - v[1] * t[0] + v[2]
//     };

//     // Calculate s-vector (s = 2t / (1 + t^2))
//     double s[3] = {
//         2.0 * t[0] / (1.0 + tsq),
//         2.0 * t[1] / (1.0 + tsq),
//         2.0 * t[2] / (1.0 + tsq)
//     };

//     // Update velocity using s (v_new = v' + v x s)
//     v[0] = vprime[1] * s[2] - vprime[2] * s[1] + vprime[0];
//     v[1] = vprime[2] * s[0] - vprime[0] * s[2] + vprime[1];
//     v[2] = vprime[0] * s[1] - vprime[1] * s[0] + vprime[2];

//     // Half electric field update
//     for (int i = 0; i < 3; i++) {
//         v[i] += charge * E[i] * dt / (2.0 * mass);
//     }

//     // Update position
//     for (int i = 0; i < 3; i++) {
//         xnew[i]  = x[i] + v[i] * dt;
//     }
// }

/* ----------------------------------------------------------------------
   get electric potential for particle I from x,y,z fields
------------------------------------------------------------------------- */

void  Update::getNeighboringPoints(const Point& closestPoint, Point& previousPoint, Point& nextPoint, const std::vector<Point>& points) const {
    for (size_t i = 0; i < points.size(); i++) {
        if (points[i].id == closestPoint.id) {
            // Check bounds before accessing neighboring points
            if (i > 0) {
                previousPoint = points[i - 1];
            }
            if (i + 1 < points.size()) {
                nextPoint = points[i + 1];
            }
            break;
        }
    }
}

std::pair<double, double> Update::getElectricPotential(double x, double y, double charge, double mass, PlasmaParams params) const {
        cachePointsFromFile();

        if (cachedPoints.empty()) {
            std::cerr << "No points were read from the file." << std::endl;
            return std::make_pair(0.0, 0.0);  // Or any other default/error value you'd like to use
        }

        // Get the closest point
        Point target = {0, x, y};
        Result closest;
        closest.distance = std::numeric_limits<double>::max();
        for (const Point& point : cachedPoints) {
            double distSq = squaredDistance(target, point);
            if (distSq < closest.distance) {
                closest.distance = distSq;
                closest.point = point;
            }
        }

        Point previousPoint, nextPoint;
        getNeighboringPoints(closest.point, previousPoint, nextPoint, cachedPoints); 
        // // Compute the tangent (deltaX, deltaY)
        double deltaX = nextPoint.x - previousPoint.x;
        double deltaY = nextPoint.y - previousPoint.y;

        double normal[2] = {-deltaY, deltaX};
        // Normalize the normal vector
        double normalMagnitude = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
        normal[0] /= normalMagnitude;
        normal[1] /= normalMagnitude;

        // Compute plasma parameters
        double B[3] = {params.b_r, params.b_z, params.b_phi};
        double Bmagnitude = std::sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);

        // Compute the dot product between B and the normal
        double dotProduct = B[0] * normal[0] + B[1] * normal[1];

        // Compute the angle
        double angleRadians = std::acos(dotProduct / Bmagnitude);
        if (angleRadians > M_PI) {
            angleRadians =  M_PI;
        }

        closest.distance = std::sqrt(closest.distance);  // Compute actual distance once

        double te = params.temp_e;
        double ti = params.temp_i;
        double debyeLength = 7.43 * std::sqrt(te / (params.dens_e * 1e-6));
        double  cs = std::sqrt((te + ti) * kB * eV2Kelvin/ (mass ));
        double omega_ci =  charge * Bmagnitude / (mass );
        double larmorRadius = cs / omega_ci;
        const double t_width = domain->t_sheath;
        double lmps = std::abs(0.5 * larmorRadius * std::cos(angleRadians) * t_width ); 
        double potential_wall = 0.5 * std::log( 2 * M_PI * ME/ (proton_mass) * ( 1.0 + ti/te));
        double potential_vmps = std::log(std::sin(angleRadians)); // log(sin(angleRadians))
        double potential_ds = potential_wall  - potential_vmps;

        double Efield = std::abs(potential_brooks(closest.distance, debyeLength, lmps, potential_ds, potential_vmps, Bmagnitude ) * kB * te * eV2Kelvin / (protonCharge ));
        if (Efield < 1e-10) {
            Efield = 0.0;
        }
        
        electricPotentialCache[{x, y}] = std::abs(Efield);


        // return Efield;
        return std::make_pair(Efield, angleRadians);

    };

double Update::potential_brooks(double z, double lambda_ds, double lmps, double v_ds, double v_wall, double Bmagnitude) const {

    if (Bmagnitude < 1e-10 || z/lambda_ds < 100 ) {
        return v_ds * std::exp(- 0.5 * z/lambda_ds)/lambda_ds;
    } 

    if (z/lambda_ds > 50 ) {
        return v_wall * std::exp(-z/lmps) / lmps;
    } 

    return v_ds * std::exp(- 0.5 * z/lambda_ds)/lambda_ds  + v_wall * std::exp(-z/lmps) / lmps;
}


double Update::potential_PIC(double alpha, double minDistance) const {

    const double t_width = domain->t_sheath;
    const double C1_alpha = -0.00281407f * alpha - 2.31655435f;
    const double C2_alpha = 0.00640402f * alpha + 0.01023915f;
    return std::abs(C1_alpha * C2_alpha) * std::exp(-C2_alpha * std::abs(minDistance * t_width));
}                                                                                                                 

void Update::crossFieldDiffusion(double *xc, double *xnew, double dt, double *Bfield ) {
  //
    const double D_perp = domain->t_D_perp;

    double B[3];
    B[0] = Bfield[0];
    B[1] = Bfield[1];
    B[2] = Bfield[2];

    // Normalize B-field
    double Bnorm = std::sqrt(Bfield[0]*Bfield[0] + Bfield[1]*Bfield[1] + Bfield[2]*Bfield[2]);
    double Bx = Bfield[0] / Bnorm;
    double By = Bfield[1] / Bnorm;
    double Bz = Bfield[2] / Bnorm;
double P1x, P1y, P1z, P2x, P2y, P2z;
    if (std::abs(Bz) < 1e-10) {
        P1x = 0;
        P1y = -Bz;
        P1z = By;
    } else {
        P1x = By;
        P1y = -Bx;
        P1z = 0;
    }
    double P1norm = std::sqrt(P1x*P1x + P1y*P1y + P1z*P1z);
    P1x /= P1norm;
    P1y /= P1norm;
    P1z /= P1norm;

    P2x = By * P1z - Bz * P1y;
    P2y = Bz * P1x - Bx * P1z;
    P2z = Bx * P1y - By * P1x;

    std::normal_distribution<double> dist(0.0, std::sqrt(2 * D_perp * dt));
    double delta_x = dist(gen);
    double delta_y = dist(gen);

    xnew[0] += P1x * delta_x + P2x * delta_y;
    xnew[1] += P1y * delta_x + P2y * delta_y;
    xnew[2] += P1z * delta_x + P2z * delta_y;
}

/* define background collisions */
void Update::getSlowDownFrequencies(double& nu_friction, double& nu_deflection, 
    double& nu_parallel, double& nu_energy, 
    double *v, double charge, double amu, double te_eV, double ti_eV, double density) {

  const double vx = v[0];
  const double vy = v[1];
  const double vz = v[2];

  const double flowVelocity[3]= {0.0, 0.0, 0.0};
  const double relativeVelocity[3] = {vx - flowVelocity[0], vy - flowVelocity[1], vz - flowVelocity[2]};
  const double velocityNorm = std::sqrt(relativeVelocity[0]*relativeVelocity[0] + 
                                         relativeVelocity[1]*relativeVelocity[1] + 
                                         relativeVelocity[2]*relativeVelocity[2]);

  const double lam_d = std::sqrt(EPS0 * te_eV / (density * background_Z * background_Z * Q));
  const double lam = 12.0 * M_PI * density * lam_d * lam_d * lam_d / charge;
  double lam_gitr = 12.0*M_PI*density*std::pow(lam_d,3)/charge;


  double gam_electron_background = 0.238762895 * charge * charge * std::log(lam) / (amu * amu);
  double gam_ion_background = 0.238762895 * charge * charge * background_Z * background_Z * std::log(lam) / (amu * amu);

  if(gam_electron_background < 0.0) gam_electron_background=0.0;
  if(gam_ion_background < 0.0) gam_ion_background=0.0;

  const double a_ion = background_mass * MI / (2 * ti_eV * Q);
  const double a_electron = ME / (2 * te_eV * Q);

  const double xx = velocityNorm * velocityNorm * a_ion;
  const double psi_prime = 2.0 * std::sqrt(xx/M_PI) * std::exp(-xx);
  const double psi_psiprime = std::erf(std::sqrt(xx));
  const double psi = psi_psiprime - psi_prime;

  const double xx_e = velocityNorm * velocityNorm * a_electron;
  const double psi_prime_e = 1.128379 * std::sqrt(xx_e);
  const double psi_e = 0.75225278 * std::pow(xx_e,1.5);
  const double psi_psiprime_e = psi_e + psi_prime_e;

  const double nu_0_i = gam_electron_background * density / std::pow(velocityNorm, 3);
  const double nu_0_e = gam_ion_background * density / std::pow(velocityNorm, 3);


  nu_friction = (1 + amu / background_mass) * psi * nu_0_i;
  nu_deflection = 2 * (psi_psiprime - psi / (2 * xx)) * nu_0_i;
  nu_parallel = psi / xx * nu_0_i;
  nu_energy = 2 * (amu / background_mass * psi - psi_prime) * nu_0_i;

  if(te_eV <= 0.0 || density <= 0.0) {
    nu_friction = 0.0;
    nu_deflection = 0.0;
    nu_parallel = 0.0;
    nu_energy = 0.0;
  }
}

void Update::getSlowDownDirections2(double parallel_direction[], double perp_direction1[], 
    double perp_direction2[], double vx, double vy, double vz) {

  double v = std::sqrt(vx * vx + vy * vy + vz * vz);
  
  // Handle the case where v == 0.0
  if (v == 0.0) {
    v = 1.0;
    vz = 1.0;
    vx = 0.0;
    vy = 0.0;
  }

  double ez1 = vx / v;
  double ez2 = vy / v;
  double ez3 = vz / v;
    
  // Compute the first perpendicular direction
  double ex1 = ez2;
  double ex2 = -ez1;
  double ex3 = 0.0;
  
  // Handle edge case for particles moving purely in z-direction
  double exnorm = std::sqrt(ex1 * ex1 + ex2 * ex2);
  if (std::abs(exnorm) < 1.0e-12) {
    ex1 = -ez3;
    ex2 = 0.0;
    ex3 = ez1;
  }

  // Normalize ex
  exnorm = std::sqrt(ex1 * ex1 + ex2 * ex2 + ex3 * ex3);
  ex1 /= exnorm;
  ex2 /= exnorm;
  ex3 /= exnorm;
  
  // Compute the second perpendicular direction
  double ey1 = ez2 * ex3 - ez3 * ex2;
  double ey2 = ez3 * ex1 - ez1 * ex3;
  double ey3 = ez1 * ex2 - ez2 * ex1;

  // Populate output arrays
  parallel_direction[0] = ez1; 
  parallel_direction[1] = ez2;
  parallel_direction[2] = ez3;
  
  perp_direction1[0] = ex1; 
  perp_direction1[1] = ex2;
  perp_direction1[2] = ex3;
  
  perp_direction2[0] = ey1; 
  perp_direction2[1] = ey2;
  perp_direction2[2] = ey3;
}

void Update::backgroundCollisions(double *xc, double *v, double dt, double mass, double charge,  PlasmaParams params) {
    // Get plasma data
  //  auto start = std::chrono::high_resolution_clock::now(); 
      double n1 = distribution(gen);
    double n2 = distribution(gen);
    double xsi = dist(gen);
    
    double flowVelocity[3] = {0.0, 0.0, 0.0};
    double t_flow_scale = domain->t_flow;
    double B[3];
    B[0] = params.b_r ;
    B[1] = params.b_z;
    B[2] = params.b_phi;
    double density = params.dens_i;
    double ti_eV = params.temp_i;
    double te_eV = params.temp_e;
    double vflow_parr = params.parr_flow * t_flow_scale;
    
    double Bnorm = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

    for(int i=0; i<3; i++) {
        flowVelocity[i] = vflow_parr * (B[i] / Bnorm);
     
    }
    
    double relativeVelocity[3];
    for(int i=0; i<3; i++) {
        relativeVelocity[i] = v[i] - flowVelocity[i];
    }
    double velocityRelativeNorm = std::sqrt(std::inner_product(relativeVelocity, relativeVelocity+3, relativeVelocity, 0.0));


    // getSlowDownFrequencies
    double nu_friction, nu_deflection, nu_parallel, nu_energy;
    getSlowDownFrequencies(nu_friction, nu_deflection, nu_parallel, nu_energy, v, charge, mass, te_eV, ti_eV,  density);

    // getSlowDownDirections2
    double parallel_direction[3], perp_direction1[3], perp_direction2[3];
    getSlowDownDirections2(parallel_direction, perp_direction1, perp_direction2, v[0], v[1], v[2]);
    
    double coeff_par = n1 * std::sqrt(2.0 * nu_parallel * dt);
    double cosXsi = std::min(std::max(cos(2.0 * M_PI * xsi), -1.0), 1.0);
    double coeff_perp1 = cosXsi * std::sqrt(nu_deflection * dt * 0.5);
    double coeff_perp2 = sin(2.0 * M_PI * xsi) * std::sqrt(nu_deflection * dt * 0.5);
    double nuEdt = std::min(std::max(nu_energy * dt, -1.0), 1.0);

    for(int i=0; i<3; i++) {
        double V_rel = velocityRelativeNorm * (1.0 - 0.5 * nuEdt * dt) * 
                      (coeff_par * parallel_direction[i] + 
                       std::abs(n2) * coeff_perp1 * perp_direction1[i] +
                       std::abs(n2) * coeff_perp2 * perp_direction2[i]) -
                      velocityRelativeNorm * dt * nu_friction * parallel_direction[i];
        v[i] = V_rel + flowVelocity[i];
    }

}

