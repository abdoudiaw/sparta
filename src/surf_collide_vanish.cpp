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
#include <fstream> 
#include "math.h"
#include "surf_collide_vanish.h"
#include "error.h"
#include <iostream> // Include header for standard I/O operations

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfCollideVanish::SurfCollideVanish(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal surf_collide vanish command");

  allowreact = 0;
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   simply return ip = NULL to delete particle
   return reaction = 0 = no reaction took place
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideVanish::
collide(Particle::OnePart *&ip, double &, int, double *, int, int &)
{
        Particle::Species *species = particle->species;

  // exit(0);
  nsingle++;

  // save particle info and surface normal

  double x[3],v[3];
  x[0] = ip->x[0];
  x[1] = ip->x[1];
  x[2] = ip->x[2];
  v[0] = ip->v[0];
  v[1] = ip->v[1];
  v[2] = ip->v[2];
  int isp = ip->ispecies;
  double mass = species[isp].molwt; 
  double charge = species[isp].charge;


    // Open a file in append mode
    std::ofstream file("particle_data.txt", std::ios::app);

    // Check if the file is open and write data
    if (file.is_open()) {
        // Writing headers if the file is empty
        file.seekp(0, std::ios::end); // Go to the end of file
        if (file.tellp() == 0) { // Check if file is empty
            file << "x,y,z,vx,vy,vz,charge,mass" << std::endl; // Write headers
        }

        // Write data in columns
        file << x[0] << ',' << x[1] << ',' << x[2] << ','
             << v[0] << ',' << v[1] << ',' << v[2] << ','
             << charge << ',' << mass << std::endl;

        // Close the file
        file.close();
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }

// save 
  ip = NULL;
  return NULL;
}
