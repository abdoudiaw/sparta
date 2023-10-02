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
#include "react_adas.h"
#include "particle.h"
#include "collide.h"
#include "random_knuth.h"
#include "error.h"
#include <vector>
#include <string>
#include "/opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense"
#include <cmath> // Include this header
#include <stdexcept> // For std::runtime_error
#include "update.h"
using namespace SPARTA_NS;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files

/* ---------------------------------------------------------------------- */

ReactADAS::ReactADAS(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactADAS::init()
{
  ReactBird::init();
}

/* ---------------------------------------------------------------------- */

int ReactADAS::attempt(Particle::OnePart *ip, Particle::OnePart *jp,
                      double pre_etrans, double pre_erot, double pre_evib,
                      double &post_etotal, int &kspecies)
{
  double pre_etotal,ecc,e_excess;
  OneReaction *r;
  double *x;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  x = ip->x;
  double imass = species[isp].molwt ;
  double jmass = species[jsp].molwt;
  double icharge = species[isp].charge;
  double jcharge = species[jsp].charge;

  double dt = update->dt;
  double react_prob_ioniziation =  update->get_ionization_rates(x, imass, icharge, dt );
  double react_prob_recombination = update->get_recombination_rates(x, imass, icharge, dt );


  // probablity to compare to reaction probability
  double current_prob=0;
  double random_prob = random->uniform();

  int n = reactions[isp][jsp].n;
  if (n == 0) return 0;
  int *list = reactions[isp][jsp].list;
  // loop over possible reactions for these 2 species
  for (int i = 0; i < n; i++) {
    r = &rlist[list[i]];
    printf(" looping over possible reactions for these 2 species %i %s\n", i, r->type);
    switch (r->type) {
    case IONIZATION:
      {
        current_prob = react_prob_ioniziation;
        break;
      }
    case RECOMBINATION:
      {
        current_prob = react_prob_recombination;
        break;
      }
    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }
    if (current_prob > random_prob) {
      tally_reactions[list[i]]++;
      ip->ispecies = r->products[0];
    printf("ionization %g %g\n", react_prob_ioniziation, react_prob_recombination);
      switch (r->type) {
      case IONIZATION:
        {
        printf("ionization %g %g\n", react_prob_ioniziation, react_prob_recombination);
          jp->ispecies = r->products[0];
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

      return 1;
    }
  }
  return 0;
}



// /* ----------------------------------------------------------------------
// // read density and temperature profiles:
// ------------------------------------------------------------------------- */

// std::vector<DataPointPlasma> ReactADAS::loadDataPlasma(const std::string& filename) {
//     std::ifstream file(filename);
//     std::vector<DataPointPlasma> data;
//     std::string line;
//     // Skip the header
//     std::getline(file, line);
//     DataPointPlasma point;
//     while (file >> point.r >> point.z >> point.ne >> point.te) {
//         data.push_back(point);
//     }
//     return data;
// }


// std::vector<DataPointRate> ReactADAS::loadDataRate(const std::string& filename) {
//     std::ifstream file(filename);
//     std::vector<DataPointRate> data;
//     std::string line;
//     // Skip the header
//     std::getline(file, line);
//     DataPointRate point;
//     while (file >> point.te >> point.ne >> point.rate) {
//         data.push_back(point);
//     }
//     return data;
// }

// Eigen::VectorXd interpolatePlasma(const std::vector<DataPointPlasma>& data, double x, double z) {
//     Eigen::VectorXd result(2);
//     result << 0, 0;

//     for (size_t i = 1; i < data.size(); ++i) {
//         if (x >= data[i-1].r && x <= data[i].r) {
//             for (size_t j = 1; j < data.size(); ++j) {
//                 if (z >= data[j-1].z && z <= data[j].z) {
//                     // Bilinear interpolation
                    
//                     double alpha = (x - data[i-1].r) / (data[i].r - data[i-1].r);
//                     double beta = (z - data[j-1].z) / (data[j].z - data[j-1].z);

//                     double  ne1 = data[i-1].ne + alpha * (data[i].ne - data[i-1].ne);
//                     double  ne2 = data[j-1].ne + beta * (data[j].ne - data[j-1].ne);
                    
//                     double te1 = data[i-1].te + alpha * (data[i].te - data[i-1].te);
//                     double te2 = data[j-1].te + beta * (data[j].te - data[j-1].te);
                    
//                     result(0) =  ne1 + beta * ( ne2 -  ne1);
//                     result(1) = te1 + beta * (te2 - te1);
                    
//                     return result;
//                 }
//             }
//         }
//     }
//     return result;
// }

// const std::vector<DataPointPlasma>& ReactADAS::getCachedPlasmaData() {
//     if (cachedDataPlasma.empty()) {
//         cachedDataPlasma = loadDataPlasma("plasma.txt");
//     }
//     return cachedDataPlasma;
// }
//     std::vector<double> ReactADAS::get_density_temperature(double *x) {
//     std::vector<DataPointPlasma> data = getCachedPlasmaData();
//     auto interpolatedValues = interpolatePlasma(data, x[0], x[1]);

//     std::vector<double> plasma = { interpolatedValues(0), interpolatedValues(1), 0 };
//     return plasma;
// }


// const std::vector<DataPointRate>& ReactADAS::getCachedIonizationRates(int mass, int charge) {
//     std::string material;

//     if (mass == 184) {
//         material = "tungsten";
//     } else if (mass == 16) {
//         material = "oxygen";
//     } else {
//         throw std::runtime_error("Material does not exist in the DB!");
//     }

//     std::string filename = "rates/" + material + "_ionization." + std::to_string(charge) + ".txt";

//     // If the cache is empty or if the filename has changed (i.e., new mass or charge), reload data.
//     if (cachedDataIonizationRates.empty() || filename != cachedIonFilename) {
//         // std::cout << "Reading file: " << filename << std::endl;
//         cachedDataIonizationRates = loadDataRate(filename);
//         cachedIonFilename = filename;  // Store the current filename to check against future calls.
//     }

//     return cachedDataIonizationRates;
// }

// const std::vector<DataPointRate>& ReactADAS::getCachedRecombRates(int mass, int charge) {
//     std::string material;

//     if (mass == 184) {
//         material = "tungsten";
//     } else if (mass == 16) {
//         material = "oxygen";
//     } else {
//         throw std::runtime_error("Material does not exist in the DB!");
//     }

//     std::string filename = "rates/" + material + "_recombination." + std::to_string(charge) + ".txt";

//     // If the cache is empty or if the filename has changed (i.e., new mass or charge), reload data.
//     if (cachedDataRecombRates.empty() || filename != cachedRecomFilename) {
//         // std::cout << "Reading file: " << filename << std::endl;
//         cachedDataRecombRates = loadDataRate(filename);
//         cachedRecomFilename = filename;  // Store the current filename to check against future calls.
//     }

//     return cachedDataRecombRates;
// }


// double ReactADAS::get_ionization_rates(double *x, int mass, int charge) {
//     // Get temperature and density:
//     std::vector<double> plasmaData = get_density_temperature(x);
//     std::vector<DataPointRate> data = getCachedIonizationRates(mass, charge);

//     double neLog = log10(plasmaData[0]);
//     double teLog = log10(plasmaData[1]);

//     double rateResult = 0.0;



//     for (size_t i = 1; i < data.size(); ++i) {
//         if (neLog >= data[i-1].ne && neLog <= data[i].ne) {
//             for (size_t j = 1; j < data.size(); ++j) {
//                 if (teLog >= data[j-1].te && teLog <= data[j].te) {
//                     // Bilinear interpolation
//                     double alpha = (neLog - data[i-1].ne) / (data[i].ne - data[i-1].ne);
//                     double beta = (teLog - data[j-1].te) / (data[j].te - data[j-1].te);

//                     double rate1 = data[i-1].rate + alpha * (data[i].rate - data[i-1].rate);
//                     double rate2 = data[j-1].rate + beta * (data[j].rate - data[j-1].rate);

//                     rateResult = rate1 + beta * (rate2 - rate1);
//                     return pow(10.0, rateResult) * plasmaData[0];
//                 }
//             }
//         }
//     }
//     return pow(10.0, rateResult) * plasmaData[0];
// }



// double ReactADAS::get_recombination_rates(double *x, int mass, int charge) {
//     // Get temperature and density:
//     std::vector<double> plasmaData = get_density_temperature(x);
//     std::vector<DataPointRate> data = getCachedRecombRates(mass, charge);

//     double neLog = log10(plasmaData[0]);
//     double teLog = log10(plasmaData[1]);
//     double rateResult = 0.0;

//     for (size_t i = 1; i < data.size(); ++i) {
//         if (neLog >= data[i-1].ne && neLog <= data[i].ne) {
//             for (size_t j = 1; j < data.size(); ++j) {
//                 if (teLog >= data[j-1].te && teLog <= data[j].te) {
//                     // Bilinear interpolation
//                     double alpha = (neLog - data[i-1].ne) / (data[i].ne - data[i-1].ne);
//                     double beta = (teLog - data[j-1].te) / (data[j].te - data[j-1].te);

//                     double rate1 = data[i-1].rate + alpha * (data[i].rate - data[i-1].rate);
//                     double rate2 = data[j-1].rate + beta * (data[j].rate - data[j-1].rate);

//                     rateResult = rate1 + beta * (rate2 - rate1);
//                     return pow(10.0, rateResult) * plasmaData[0];
//                 }
//             }
//         }
//     }
//     return pow(10.0, rateResult) * plasmaData[0];
// }


