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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "react_tce.h"
#include "particle.h"
#include "collide.h"
#include "random_knuth.h"
#include "error.h"
#include "grid.h"
#include "update.h"
#include <cmath>
#include <iostream>
using namespace SPARTA_NS;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files

/* ---------------------------------------------------------------------- */

ReactTCE::ReactTCE(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactTCE::init()
{
  if (!collide || strcmp(collide->style,"vss") != 0)
    error->all(FLERR,"React tce can only be used with collide vss");

  ReactBird::init();
}

/* ---------------------------------------------------------------------- */

int ReactTCE::attempt(Particle::OnePart *ip, Particle::OnePart *jp,
                      double pre_etrans, double pre_erot, double pre_evib,
                      double &post_etotal, int &kspecies)
{
  double pre_etotal,ecc,e_excess;
  OneReaction *r;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int icell = jp->icell;
    // loop over cells I own
  Grid::ChildInfo *cinfo = grid->cinfo;


  std::unordered_map<int, std::pair<double, double>> cellCenters;

     // 3. Check the cellCenters cache
        std::pair<double, double> xc;
        if (cellCenters.find(icell) == cellCenters.end()) {
            double *lo = grid->cells[icell].lo;
            double *hi = grid->cells[icell].hi;
            xc = {0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])};
            cellCenters[icell] = xc;
        } else {
            xc = cellCenters[icell];
        }

          // 4. Use the xc values to get plasma parameters
          PlasmaParams params = update->interpolatePlasmaData(xc.first, xc.second);

          // 5. Continue with your existing code
          double te = params.temp_e;
          double ne = params.dens_e;
          // get log10 of electron temperature and density
          double log_te = log10(te);
          double log_ne = log10(ne);



  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.0;


  int n = reactions[isp][jsp].n;
  // print species isp and jsp
  printf(" isp = %d\n",isp);
  printf(" jsp = %d\n",jsp);
  printf(" reactiions number n = %d\n",n);


  // if (n == 0) return 0;


  int *list = reactions[isp][jsp].list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform();

  // loop over possible reactions for these 2 species


  for (int i = 0; i < n; i++) {
    r = &rlist[list[i]];

    // // ignore energetically impossible reactions

    // pre_etotal = pre_etrans + pre_erot + pre_evib;

    // ecc = pre_etrans;
    // if (pre_ave_rotdof > 0.1) ecc += pre_erot*r->coeff[0]/pre_ave_rotdof;

    // e_excess = ecc - r->coeff[1];
    // printf(" coeff[1] = %g\n",r->coeff[1]);

    // if (e_excess <= 0.0) 
    // { 
    //       printf("ignoring reaction because e_excess = %f\n",e_excess);
      
    //   continue;
    // }

    // compute probability of reaction

    // get rate for ionization:
    // double rate =
    // get reaction type


    switch (r->type) {
    // case DISSOCIATION:
    case IONIZATION:
    // case EXCHANGE:
      {
        // get coefficient for reaction Ionization(double x, double a, double x0, double sigma, double alpha, double y0) 
        // IonizationRate(double x, double a, double x0, double sigma, double alpha, double y0)

        // I A 613.579 -2.976 169.458 0.861 -1239.458 0
        double a = r->coeff[0];
        double x0 = r->coeff[1];
        double sigma = r->coeff[2];
        double alpha = r->coeff[3];
        double y0 = r->coeff[4];
        double rate = pow(10, IonizationRate(log_te, a, x0, sigma, alpha, y0));
        // P = 1.0 - exp(-rate * dt * n)  
        // Taylor expansion of exp(-rate * dt * n) = 1 - rate * dt * n + rate^2 * dt^2 * n^2 / 2
        // P = rate * dt * n - rate^2 * dt^2 * n^2 / 2
        react_prob = rate * update->dt * ne -  rate * rate * update->dt * update->dt * ne * ne / 2.0;

        // react_prob += r->coeff[2] *
        //   pow(ecc-r->coeff[1],r->coeff[3]) *
        //   pow(1.0-r->coeff[1]/ecc,r->coeff[5]);
        //   // printf("react_prob = %f\n",react_prob);
        break;
      }

    case RECOMBINATION:
      {
        // skip if no 3rd particle chosen by Collide::collisions()
        //   this includes effect of boost factor to skip recomb reactions
        // check if this recomb reaction is the same one
        //   that the 3rd particle species maps to, else skip it
        // this effectively skips all recombinations reactions
        //   if selected a 3rd particle species that matches none of them
        // scale probability by boost factor to restore correct stats

        // printf("recomb reaction\n");

         double a1 = r->coeff[0];
        double x0 = r->coeff[1];
        double sigma = r->coeff[2];
        double a2 = r->coeff[3];
        double k = r->coeff[4];
        double b = r->coeff[5];
        // printf("r->coeff[2] = %f\n", r->coeff[2]);

        // // print coefficients
        // printf("a1 = %f\n",a1);
        // printf("x0 = %f\n",x0);
        // printf("sigma = %g\n",sigma);
        // printf("a2 = %f\n",a2);
        // printf("k = %f\n",k);
        // printf("b = %f\n",b);
        // // print ne and te
        // printf("ne = %g\n",ne);
        // printf("te = %g\n",te);

        // // compute rate
        double dt = update->dt;
        double rate = pow(10, RecombinationRate(log_te, a1, x0, sigma, a2, k, b));
        // P = 1.0 - exp(-rate * dt * n) 
        // Taylor expansion of exp(-rate * dt * n) = 1 - rate * dt * n + rate^2 * dt^2 * n^2 / 2
        // P = rate * dt * n - rate^2 * dt^2 * n^2 / 2

        react_prob = rate * dt * ne -  rate * rate * dt * dt * ne * ne / 2.0;

        // exit(0);
        // printf("IONIZATION reaction\n");
        // exit(0);
        // if (recomb_species < 0) continue;
        // int *sp2recomb = reactions[isp][jsp].sp2recomb;
        // if (sp2recomb[recomb_species] != list[i]) continue;

        // react_prob += recomb_boost * recomb_density * r->coeff[2] *
        //   pow(ecc,r->coeff[3]) *
        //   pow(1.0-r->coeff[1]/ecc,r->coeff[5]);

        // react_prob *= recomb_boost_inverse;

        
        break;
      }

    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }

    if (react_prob > random_prob) {
      tally_reactions[list[i]]++;
      ip->ispecies = r->products[0];

      switch (r->type) {
      case IONIZATION:
        {
          jp->ispecies = r->products[1];
          break;
        }
      case RECOMBINATION:
        {

          jp->ispecies = -1;
          break;
        }
      }

      if (r->nproduct > 2) kspecies = r->products[2];
      else kspecies = -1;

      post_etotal = pre_etotal + r->coeff[4];

      return 1;
    }
  }

  return 0;
}


double ReactTCE::IonizationRate(double x, double a, double x0, double sigma, double alpha, double y0) {
    // Gaussian function part
    double gauss = exp(-pow(x - x0, 2) / (2 * pow(sigma, 2)));
    // Error function to introduce skew
    double erf = 1 + tanh(alpha * (x - x0));
    return y0 + a * gauss * erf;
}

// for recombination
double ReactTCE::RecombinationRate(double x, double a1, double x0, double sigma, double a2, double k, double b) {
    return a1 * exp(-pow(x - x0, 2) / (2 * pow(sigma, 2))) + a2 * exp(k * x) + b;
}
