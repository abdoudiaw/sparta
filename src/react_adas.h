// /* ----------------------------------------------------------------------
//    SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
//    http://sparta.sandia.gov
//    Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
//    Sandia National Laboratories

//    Copyright (2014) Sandia Corporation.  Under the terms of Contract
//    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
//    certain rights in this software.  This software is distributed under
//    the GNU General Public License.

//    See the README file in the top-level SPARTA directory.
// ------------------------------------------------------------------------- */

// #ifdef REACT_CLASS

// ReactStyle(adas,ReactADAS)

// #else

// #ifndef SPARTA_REACT_ADAS_H
// #define SPARTA_REACT_ADS_H

// #include "react_bird.h"
// #include "particle.h"
// #include <vector>
// #include <string>
// #include <vector>
// #include <iostream>
// #include <fstream>

// namespace SPARTA_NS {

// struct DataPointPlasma {    // Added for representing the magnetic field data
//     double r, z, ne, te;
// };

// struct DataPointRate {    // Added for representing the magnetic field data
//     double ne, te, rate;
// };

// struct PlasmaData {
//     double density;
//     double temperature;
// };


// class ReactADAS : public ReactBird {
//  public:
//   ReactADAS(class SPARTA *, int, char **);
//   void init();
//   int attempt(Particle::OnePart *, Particle::OnePart *,
//               double, double, double, double &, int &);
//    // void get_density_temperature( double *, double * );
//     std::vector<double> get_density_temperature(double *) ;
//     std::vector<DataPointPlasma> cachedDataPlasma;       // To store the cached magnetic field data
//      std::vector<DataPointRate> cachedDataIonizationRates;
//      std::vector<DataPointRate> cachedDataRecombRates;
//     const std::vector<DataPointPlasma>& getCachedPlasmaData();
//     const std::vector<DataPointRate>& getCachedIonizationRates(int, int);
//     const std::vector<DataPointRate>& getCachedRecombRates(int , int );
//     double  get_ionization_rates(double *, int , int );
//     double get_recombination_rates(double *, int , int );
// protected:

//    std::vector<DataPointPlasma> loadDataPlasma(const std::string& filename); // Helper function to load the data from the file

//    std::vector<DataPointRate> loadDataRate(const std::string& filename); // Helper function to load the data from the file

//    std::string cachedFilename;  // <-- Add this line
//    std::string cachedIonFilename;
//    std::string cachedRecomFilename;

// };

// }

// #endif
// #endif
