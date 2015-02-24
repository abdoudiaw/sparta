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

FixStyle(inflow,FixInflow)

#else

#ifndef SPARTA_FIX_INFLOW_H
#define SPARTA_FIX_INFLOW_H

#include "fix.h"

namespace SPARTA_NS {

class FixInflow : public Fix {
 public:
  FixInflow(class SPARTA *, int, char **);
  ~FixInflow();
  int setmask();
  void init();
  void start_of_step();
  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void compress_grid();
  void post_compress_grid();
  double compute_vector(int);

 private:
  int nevery,imix;
  int faces[6];
  int np,perspecies;
  int npercell,nthresh;
  int nsingle,ntotal;

  struct CellFace {
    double lo[3];               // lower-left corner of face
    double hi[3];               // upper-right corner of face
    double normal[3];           // inward normal from external boundary
    double ntarget;             // # of mols to insert for all species
    double *ntargetsp;          // # of mols to insert for each species
    cellint id;                 // ID of cell with insertion face
    int pcell;                  // associated cell index for particles
                                // unsplit or sub cell (not split cell)
    int icell;                  // associated cell index, unsplit or split cell
    int iface;                  // which face of unsplit or split cell
    int ndim;                   // dim (0,1,2) normal to face
    int pdim,qdim;              // 2 dims (0,1,2) parallel to face
  };

  CellFace *cellface;           // cell/face pairs to insert particles on
  int ncf;                      // # of cell/face pairs
  int ncfmax;                   // max # of cell/face pairs allocated

  int **c2f;                    // I,J = cellface index of Jth face of cell I
                                // only for unsplit and split cells, not sub
                                // -1 if no insertions on that face
  int nglocal;                  // # of owned grid cells
  int nglocalmax;               // max size of per-cell vectors/arrays

  class RanPark *random;

  double mol_inflow(double, double, double);
  int split(int, int);
  void grow_percell(int);
  void grow_cellface(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Fix inflow mixture ID does not exist

Self-explanatory.

E: Cannot use fix inflow in z dimension for 2d simulation

Self-explanatory.

E: Cannot use fix inflow in y dimension for axisymmetric

This is because the y dimension boundaries cannot be
inflow boundaries for an axisymmetric model.

E: Cannot use fix inflow n > 0 with perspecies yes

This is because the perspecies option calculates the
number of particles to insert itself.

E: Cannot use fix inflow on periodic boundary

Self-explanatory.

E: Fix inflow used on outflow boundary

Self-explanatory.

*/
