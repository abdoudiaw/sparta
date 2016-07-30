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

#ifdef FIX_CLASS

FixStyle(ave/grid,FixAveGrid)

#else

#ifndef SPARTA_FIX_AVE_GRID_H
#define SPARTA_FIX_AVE_GRID_H

#include "fix.h"

namespace SPARTA_NS {

class FixAveGrid : public Fix {
 public:
  int ave;                   // public so FixAdapt can check it

  FixAveGrid(class SPARTA *, int, char **);
  ~FixAveGrid();
  int setmask();
  void init();
  void setup();
  void end_of_step();

  void add_grid_one(int, int);
  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void compress_grid();
  double memory_usage();

 private:
  int nvalues,maxvalues;
  int nrepeat,irepeat,nsample;
  bigint nvalid;

  char **ids;                // ID/name of compute,fix,variable to access   
  int *which;                // COMPUTE or FIX or VARIABLE
  int *argindex;             // which column from compute or fix to access
  int *value2index;          // index of compute,fix,variable
  int *post_process;         // 1 if need compute->post_process() on value

  int ntotal;                // total # of columns in tally array
  double **tally;            // array of tally quantities, cells by ntotal
                             // can be multiple tally quants per value

                             // used when normalizing tallies
  int *nmap;                 // # of tally quantities for each value
                             //   these may not be unique
  int **map;                 // indices of non-unique tally quantities
                             //   in tally, for each value

                             // used when accumulating tallies
  int *numap;                // # of unique tally quantities for each value
  int **umap;                // indices of tally quants in tally for each value
  int **uomap;               // indices of corresponding quantities (0 to N-1)
                             //   in compute/fix tally array, for each value

  int nglocal;               // # of owned grid cells
  int nglocalmax;            // max size of per-cell vectors/arrays

  int pack_one(int, char *, int);
  int unpack_one(char *, int);
  void options(int, char **);
  void grow();
  bigint nextvalid();
  void grow_percell(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute ID for fix ave/grid does not exist

Self-explanatory.

E: Fix ID for fix ave/grid does not exist

Self-explanatory.

E: Fix ave/grid compute does not calculate per-grid values

Self-explanatory.

E: Fix ave/grid compute does not calculate a per-grid vector

Self-explanatory.

E: Fix ave/grid compute does not calculate a per-grid array

Self-explanatory.

E: Fix ave/grid compute array is accessed out-of-range

Self-explanatory.

E: Fix ave/grid fix does not calculate per-grid values

Self-explanatory.

E: Fix ave/grid fix does not calculate a per-grid vector

Self-explanatory.

E: Fix ave/grid fix does not calculate a per-grid array

Self-explanatory.

E: Fix ave/grid fix array is accessed out-of-range

Self-explanatory.

E: Fix for fix ave/grid not computed at compatible time

Fixes generate values on specific timesteps.  Fix ave/grid is
requesting a value on a non-allowed timestep.

E: Variable name for fix ave/grid does not exist

Self-explanatory.

E: Fix ave/grid variable is not grid-style variable

Self-explanatory.

*/
