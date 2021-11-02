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

FixStyle(temp/rescale,FixTempRescale)

#else

#ifndef SPARTA_FIX_TEMP_RESCALE_H
#define SPARTA_FIX_TEMP_RESCALE_H

#include "stdio.h"
#include "fix.h"

namespace SPARTA_NS {

class FixTempRescale : public Fix {
 public:
  FixTempRescale(class SPARTA *, int, char **);
  ~FixTempRescale() {}
  int setmask();
  void init();
  void end_of_step();

 private:
  double tstart,tstop;
  double tprefactor;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot open fix print file %s

The output file generated by the fix print command cannot be opened

*/
