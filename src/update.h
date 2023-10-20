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

#ifndef SPARTA_UPDATE_H
#define SPARTA_UPDATE_H

#include "math.h"
#include "pointers.h"
#include <vector>
#include <string>
#include <H5Cpp.h>
#include <map>
#include <tuple>
#include "particle.h"
#include <fstream> // for std::ifstream
#include <iostream> // for std::cerr and std::endl

using namespace std;

namespace SPARTA_NS {

struct PointRZ {
    float x;
    float z;

        PointRZ(float r_val, float z_val) : x(r_val), z(z_val) {}
};

struct DataPoint2 {
    double center_rr, center_zz, ne, te, ti, vflow, br, bz, bt;
};

// Define a simple Polygon structure with a containment check
struct Polygon {
    std::vector<PointRZ> vertices;

    bool contains(const PointRZ& PointRZ) const {
        bool inside = false;
        int j = vertices.size() - 1;

        for (size_t i = 0; i < vertices.size(); ++i) {
            if ((vertices[i].z > PointRZ.z) != (vertices[j].z > PointRZ.z) &&
                (PointRZ.x < (vertices[j].x - vertices[i].x) * (PointRZ.z - vertices[i].z) / (vertices[j].z - vertices[i].z) + vertices[i].x)) {
                inside = !inside;
            }
            j = i;
        }

        return inside;
    }
};

struct PlasmaData2 {
    double center_rr;
    double center_zz;
    double ne;
    double te;
    double ti;
    double vflow;
    double br;
    double bz;
    double bt;
};

struct PlasmaData {
    std::vector<std::vector<float>> dens_e, temp_e, dens_i, temp_i, parr_flow, b_r, b_z, b_phi;
    // r, z are 1D arrays
    std::vector<float> r, z;
};


struct PlasmaParams {
    double dens_e, temp_e, dens_i, temp_i, parr_flow, b_r, b_z, b_phi;
};

struct Point {
    int id;           // surface id
    double x;
    double y;
};

struct Result {
    double distance;
    Point point;
};

struct DataPoint {    // Added for representing the magnetic field data
    double r, z, y, br, bz, bt;
};


struct DataPointPlasma {
    double r, z, vflow, ti, te, ni, ne;

    double getValue(int index) const {
        switch(index) {
            case 0: return vflow;
            case 1: return ti;
            case 2: return te;
            case 3: return ni;
            case 4: return ne;
            default: return 0;
        }
    }
};
struct DataPointRate {    // Added for representing the magnetic field data
    double ne, te, rate;
};

struct DataPointReflectionSputtering {    
    double angle, energy, rpyld, spyld;
};

class Update : protected Pointers {
 public:
  bigint ntimestep;               // current timestep
  int nsteps;                     // # of steps to run
  int runflag;                    // 0 for unset, 1 for run
  bigint firststep,laststep;      // 1st & last step of this run
  bigint beginstep,endstep;       // 1st and last step of multiple runs
  int first_update;               // 0 before initial update, 1 after
  double dt;                      // timestep size

  char *unit_style;      // style of units used throughout simulation
  double boltz;          // Boltzmann constant (eng/degree K)
  double mvv2e;          // conversion of mv^2 to energy

  double fnum;           // ratio of real particles to simulation particles
  double nrho;           // number density of background gas
  double vstream[3];     // streaming velocity of background gas
  double temp_thermal;   // thermal temperature of background gas
  int optmove_flag;      // global optmove option set

  int fstyle;            // external field: NOFIELD, CFIELD, PFIELD, GFIELD
  double field[3];       // constant external field
  char *fieldID;         // fix ID for PFIELD or GFIELD
  int ifieldfix;         // index of external field fix
  int *field_active;     // ptr to field_active flags in fix
  int fieldfreq;         // update GFIELD every this many timsteps

  int nmigrate;          // # of particles to migrate to new procs
  int *mlist;            // indices of particles to migrate

                         // current step counters
  int niterate;          // iterations of move/comm
  int ntouch_one;        // particle-cell touches
  int ncomm_one;         // particles migrating to new procs
  int nboundary_one;     // particles colliding with global boundary
  int nexit_one;         // particles exiting outflow boundary
  int nscheck_one;       // surface elements checked for collisions
  int nscollide_one;     // particle/surface collisions

  bigint first_running_step; // timestep running counts start on
  int niterate_running;      // running count of move/comm interations
  bigint nmove_running;      // running count of total particle moves
  bigint ntouch_running;     // running count of current step counters
  bigint ncomm_running;
  bigint nboundary_running;
  bigint nexit_running;
  bigint nscheck_running;
  bigint nscollide_running;

  int nstuck;                // # of particles stuck on surfs and deleted
  int naxibad;               // # of particles where axisymm move was bad
                             // in this case, bad means particle ended up
                             // outside of final cell curved surf by epsilon
                             // when move logic thinks it is inside cell

  int reorder_period;        // # of timesteps between particle reordering
  int global_mem_limit;      // max # of bytes in arrays for rebalance and reordering
  int mem_limit_grid_flag;   // 1 if using size of grid as memory limit
  void set_mem_limit_grid(int gnlocal = 0);
  int have_mem_limit();      // 1 if have memory limit

  int copymode;          // 1 if copy of class (prevents deallocation of
                         //  base class when child copy is destroyed)

  class RanMars *ranmaster;   // master random number generator

  double rcblo[3],rcbhi[3];    // debug info from RCB for dump image

  // this info accessed by SurfReactAdsorb to do on-surface reaction tallying

  int nsurf_tally;         // # of Cmp tallying surf bounce info this step
  int nboundary_tally;     // # of Cmp tallying boundary bounce info this step
  class Compute **slist_active;   // list of active surf Computes this step
  class Compute **blist_active;   // list of active boundary Computes this step


  // Caching mechanism additions:
  mutable std::vector<Point> CachedSurfData; // Marked mutable since we modify it in a const function.
  mutable std::string cachedSurfDataFilename; // Same reason as above.
  mutable bool isSurfDataInitialized = false; // Same reason as above.

  mutable std::vector<Point> cachedPoints;
  mutable bool pointsCached = false;


  struct pair_hash {
    std::size_t operator() (const std::pair<double, double>& p) const {
        auto h1 = std::hash<double>{}(p.first); 
        auto h2 = std::hash<double>{}(p.second); 
        return h1 ^ h2;
    }
};
  struct VectorHash {
        size_t operator()(const std::array<double, 3>& vec) const {
            // Hash the values in xc to produce a unique key
            return std::hash<double>()(vec[0]) ^ std::hash<double>()(vec[1]) ^ std::hash<double>()(vec[2]);
        }
    };

struct ArrayHash {
    std::size_t operator()(const std::array<double, 3>& arr) const {
        std::hash<double> hashFn;
        std::size_t h0 = hashFn(arr[0]);
        std::size_t h1 = hashFn(arr[1]);
        std::size_t h2 = hashFn(arr[2]);
        return h0 ^ (h1 << 1) ^ (h2 << 2); // Just one possible combination
    }
};


    mutable std::unordered_map<std::array<double, 3>, double, VectorHash> potentialCache;
// std::unordered_map<std::array<double, 3>, std::array<double, 3>> EfieldCache;
// std::unordered_map<std::array<double, 3>, std::array<double, 3>, ArrayHash> EfieldCache;
mutable std::unordered_map<std::array<double, 3>, std::array<double, 3>, ArrayHash> EfieldCache;


    mutable std::unordered_map<std::pair<double, double>, double, pair_hash> electricPotentialCache;
    std::map<std::pair<float, float>, PlasmaParams> plasma_cache;
    std::map<std::string, PlasmaData> plasmaDataCache;
vector<DataPoint2> ReadAndCacheData(const string& filename);

    PlasmaData readPlasmaData(const std::string& filePath);
    
    PlasmaParams interpolatePlasmaData(double r_val, double z_val);

DataPoint2 LinearInterpolate2(const DataPoint2& p1, const DataPoint2& p2, double t);

PlasmaParams interpolate(const PlasmaData2& lower, const PlasmaData2& upper, double alpha);

    const std::vector<DataPointReflectionSputtering>& getCachedDataReflectionSputtering(int, int);

    // double interp2dCombined(double x, double y, double z, int nx, int nz,
    // float* gridx, float* gridz, float* data);

    double interp2dCombined(double x, double y, double z, int nx, int nz,
    const float* gridx, const float* gridz, const float* data);



// void loadPlasmaData(const std::string &filename) {
//     plasma_data = readPlasmaData(filename);
    
//     flattened_temp_e = flatten(plasma_data.temp_e);
//     flattened_temp_i = flatten(plasma_data.temp_i);
//     flattened_dens_e = flatten(plasma_data.dens_e);
//     flattened_dens_i = flatten(plasma_data.dens_i);
//     flattened_parr_flow = flatten(plasma_data.parr_flow);
//     flattened_b_r = flatten(plasma_data.b_r);
//     flattened_b_z = flatten(plasma_data.b_z);
//     flattened_b_phi = flatten(plasma_data.b_phi);

//     r_double.assign(plasma_data.r.begin(), plasma_data.r.end());
//     z_double.assign(plasma_data.z.begin(), plasma_data.z.end());
// }

  int getMaxChargeNumber(double );
  bool validateSpeciesChange(double , int , int );

  void process_particle(Particle::OnePart *p, Particle::Species *species, int sp,
  double te, double ne); //, RateData &rateData);


void getNeighboringPoints(const Point& closestPoint, Point& previousPoint, Point& nextPoint, const std::vector<Point>& points) const;

    // double potential_brooks(double , double , double , double , double, double ) const ;
    // double potential_brooks(double minDistance, double alpha, PlasmaParams params) const;
    // double potential_brooks(double *xc, double minDistance, double alpha, const PlasmaParams params) const ;
// double potential_brooks(double *xc, double minDistance, double alpha, const PlasmaParams params) const;
// std::array<double, 3> potential_brooks(double *xc, const PlasmaParams params);
std::array<double, 3> potential_brooks(double *xc, const PlasmaParams params) const;

    double potential_PIC(double , double  ) const;
    // double getElectricPotential(double  , double , double, double, PlasmaParams ) const;
    std::pair<double, double> getElectricPotential(double x, double y, double charge, double mass, PlasmaParams params) const;

    void crossFieldDiffusion( double *, double *, double , double * );
    void getSlowDownFrequencies(double& , double& , double& , double& , double *, double , double, double , double , double);
    void getSlowDownDirections2(double [], double [], double [], double , double , double );
    void backgroundCollisions(double *, double *, double, double, double, PlasmaParams );
//   double getSoledgeData(float, float);
//   PlasmaParams getSoledgeData(double r_val, double z_val) const ;

  Update(class SPARTA *);
  ~Update();
  void set_units(const char *);
  virtual void init();
  virtual void setup();
  virtual void run(int);
  void global(int, char **);
  void reset_timestep(int, char **);

  int split3d(int, double *);
  int split2d(int, double *);


 protected:


  void cachePointsFromFile() const {
          if (!pointsCached) {
              static const std::string filename = "lim.txt";
              std::ifstream file(filename);
              if (!file) {
                  std::cerr << "Failed to open the file: " << filename << std::endl;
                  return;
              }
              Point point;
              while (file >> point.id >> point.x >> point.y) {
                  cachedPoints.push_back(point);
              }
              pointsCached = true;
          }
      }

      double squaredDistance(const Point& p1, const Point& p2) const {
          return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
      }

  int me,nprocs;
  int maxmigrate;            // max # of particles in mlist
  class RanKnuth *random;     // RNG for particle timestep moves

  int collide_react;         // 1 if any SurfCollide or React classes defined
  int nsc,nsr;               // copy of Collide/React data in Surf class
  class SurfCollide **sc;
  class SurfReact **sr;

  int bounce_tally;               // 1 if any bounces are ever tallied
  int nslist_compute;             // # of computes that tally surf bounces
  int nblist_compute;             // # of computes that tally boundary bounces
  class Compute **slist_compute;  // list of all surf bounce Computes
  class Compute **blist_compute;  // list of all boundary bounce Computes

  int surf_pre_tally;       // 1 to log particle stats before surf collide
  int boundary_pre_tally;   // 1 to log particle stats before boundary collide

  int collide_react_setup();
  void collide_react_reset();
  void collide_react_update();

  int bounce_setup();
  virtual void bounce_set(bigint);

  int nulist_surfcollide;
  SurfCollide **ulist_surfcollide;

  int dynamic;              // 1 if any classes do dynamic updates of params
  void dynamic_setup();
  void dynamic_update();

  void reset_timestep(bigint);

  //int axi_vertical_line(double, double *, double *, double, double, double,
  //                     double &);

  // remap x and v components into axisymmetric plane
  // input x at end of linear move (x = xold + dt*v)
  // change x[1] = sqrt(x[1]^2 + x[2]^2), x[2] = 0.0
  // change vy,vz by rotation into axisymmetric plane

  inline void axi_remap(double *x, double *v) {
    double ynew = x[1];
    double znew = x[2];
    x[1] = sqrt(ynew*ynew + znew*znew);
    x[2] = 0.0;
    double rn = ynew / x[1];
    double wn = znew / x[1];
    double vy = v[1];
    double vz = v[2];
    v[1] = vy*rn + vz*wn;
    v[2] = -vy*wn + vz*rn;
  };

  typedef void (Update::*FnPtr)();
  FnPtr moveptr;             // ptr to move method
  template < int, int, int > void move();

  int perturbflag;
  typedef void (Update::*FnPtr2)(int, int, double, double *, double *);
  FnPtr2 moveperturb;        // ptr to moveperturb method

  // variants of moveperturb method
  // adjust end-of-move x,v due to perturbation on straight-line advection

  inline void field2d(int i, int icell, double dt, double *x, double *v) {
    double dtsq = 0.5*dt*dt;
    x[0] += dtsq*field[0];
    x[1] += dtsq*field[1];
    v[0] += dt*field[0];
    v[1] += dt*field[1];
  };

  inline void field3d(int i, int icell, double dt, double *x, double *v) {
    double dtsq = 0.5*dt*dt;
    x[0] += dtsq*field[0];
    x[1] += dtsq*field[1];
    x[2] += dtsq*field[2];
    v[0] += dt*field[0];
    v[1] += dt*field[1];
    v[2] += dt*field[2];
  };

  // NOTE: cannot be inline b/c ref to modify->fix[] is not supported
  //       unless possibly include modify.h and fix.h in this file
  void field_per_particle(int, int, double, double *, double *);
  void field_per_grid(int, int, double, double *, double *);

  //

//   void pusher_boris3D(double *, double *, double *, double *, double , double , double , PlasmaParams);
//   void pusher_boris3D(double *xc, double *x, double *v, double *xnew, double charge, double mass, double dt, PlasmaParams params,  double *E, double *B);
    void pusher_boris2D(double *, double *, double *, double *, double , double , double, double *, double *);

   std::vector<DataPointReflectionSputtering> loadDataSurfaceData(const std::string& filename); // Helper function to load the data from the file

   std::string cachedFilename;  // <-- Add this line
   std::string cachedIonFilename;
   std::string cachedRecomFilename;
   std::string cachedSurfaceDataFilename;
  //  std::string cachedSurfDataFilename;
};
}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Gravity in z not allowed for 2d

Self-explanatory.

E: Gravity in y not allowed for axi-symmetric model

Self-explanatory.

E: Particle %d on proc %d hit inside of surf %d on step %ld

This error should not happen if particles start outside of physical
objects.  Please report the issue to the SPARTA developers.

E: Sending particle to self

This error should not occur.  Please report the issue to the SPARTA
developers.

E: Cannot set global surfmax when surfaces already exist

This setting must be made before any surfac elements are
read via the read_surf command.

E: Global mem/limit setting cannot exceed 2GB

Self-expanatory, prevents 32-bit interger overflow

E: Timestep must be >= 0

Reset_timestep cannot be used to set a negative timestep.

E: Too big a timestep

Reset_timestep timestep value must fit in a SPARTA big integer, as
specified by the -DSPARTA_SMALL, -DSPARTA_BIG, or -DSPARTA_BIGBIG
options in the low-level Makefile used to build SPARTA.  See
Section 2.2 of the manual for details.

E: Cannot reset timestep with a time-dependent fix defined

The timestep cannot be reset when a fix that keeps track of elapsed
time is in place.

*/
