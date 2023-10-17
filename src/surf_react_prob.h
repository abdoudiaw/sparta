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

#ifdef SURF_REACT_CLASS

SurfReactStyle(prob,SurfReactProb)

#else

#ifndef SPARTA_SURF_REACT_PROB_H
#define SPARTA_SURF_REACT_PROB_H

#include "surf_react.h"
#include <H5Cpp.h>
#include <map>
#include <tuple>

#include <fstream> // for std::ifstream
#include <iostream> // for std::cerr and std::endl



namespace SPARTA_NS {

struct SurfaceData {
    std::vector<float> energy;
    std::vector<float> angle;
    std::vector<std::vector<float>> rfyld;
    std::vector<std::vector<float>> spyld;
};

struct SurfaceDataParams{
    double rfyld, spyld;
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first); 
        auto h2 = std::hash<T2>{}(p.second); 
        return h1 ^ h2;  
    }
};
struct tuple_hash {
    template <class T1, class T2, class T3>
    std::size_t operator() (const std::tuple<T1, T2, T3>& t) const {
        size_t hash1 = std::hash<T1>{}(std::get<0>(t));
        size_t hash2 = std::hash<T2>{}(std::get<1>(t));
        size_t hash3 = std::hash<T3>{}(std::get<2>(t));

        return hash1 ^ hash2 ^ (hash3 << 1);  // Just an example of combining the hashes
    }
};


class SurfReactProb : public SurfReact {
 public:
  SurfReactProb(class SPARTA *, int, char **);
  SurfReactProb(class SPARTA *sparta) : SurfReact(sparta) {}
  virtual ~SurfReactProb();
  virtual void init();
  int react(Particle::OnePart *&, int, double *, Particle::OnePart *&, int &);
  char *reactionID(int);
  int match_reactant(char *, int);
  int match_product(char *, int);

  // reaction info, as read from file

  struct OneReaction {
    int active;                    // 1 if reaction is active
    int type;                      // reaction type = DISSOCIATION, etc
    int style;                     // reaction style = ARRHENIUS, etc
    int ncoeff;                    // # of numerical coeffs
    int nreactant,nproduct;        // # of reactants and products
    char **id_reactants,**id_products;  // species IDs of reactants/products
    int *reactants,*products;      // species indices of reactants/products
    double *coeff;                 // numerical coeffs for reaction
    char *id;                      // reaction ID (formula)
  };

 protected:

std::unordered_map<std::string, SurfaceData> surfaceDataCache; // For reading data from the file.
// std::unordered_map<std::pair<double, double, double >, SurfaceDataParams, pair_hash> interpolatedDataCache; // For interpolated data.
std::unordered_map<std::tuple<double, double, double>, SurfaceDataParams, tuple_hash> interpolatedDataCache;


  class RanKnuth *random;     // RNG for reaction probabilities

  OneReaction *rlist;              // list of all reactions read from file
  int nlist_prob;                  // # of reactions read from file
  int maxlist_prob;                // max # of reactions in rlist

  // possible reactions a reactant species is part of

  struct ReactionI {
    int *list;           // list of indices into rlist, ptr into indices
    int n;               // # of reactions in list
  };

  ReactionI *reactions;       // reactions for all species
  int *indices;               // master list of indices

  virtual void init_reactions();
  void readfile(char *);
  int readone(char *, char *, int &, int &);
    double wierzbicki_biersack(double , double , double , double , double );
    double getElectricPotential(double , double , double , double );
double yamamura(double , double , double , double , double );
    //  std::vector<DataPointPlasma> loadDataPlasma(const std::string& filename); // Helper function to load the data from the file
  SurfaceData readSurfaceData(const std::string& filePath);
  SurfaceDataParams interpolateSurfaceData(double energy_val, double angle_val, double mass);

};

}
namespace std {
    template <>
    struct hash<std::pair<double, double>> {
        size_t operator()(const std::pair<double, double>& p) const {
            auto h1 = hash<double>{}(p.first); 
            auto h2 = hash<double>{}(p.second); 
            return h1 ^ h2; 
        }
    };
}


#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
