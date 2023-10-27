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
#include "surf_react_prob.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include "random_knuth.h"
#include <string>  
#include <algorithm>
#include "/opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense"
// #include "/usr/include/eigen3/Eigen/Dense"
#include <random>
#include <iostream>
#include "surf.h"
#include "grid.h"
// #include "surfaceDataInterpolator.h"
#include "geometry_tools.h"


#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

using namespace SPARTA_NS;

enum{DISSOCIATION,EXCHANGE,RECOMBINATION,SPUTTERING};
enum{SIMPLE};

#define MAXREACTANT 1
#define MAXPRODUCT 2
#define MAXCOEFF 2

#define MAXLINE 1024
#define DELTALIST 16

double joule2ev = 6.24150934e18;
// target info
double mass_target = 184.0;
double charge_target = 74.0;
double mproton = 1.67262158e-27;

/* ---------------------------------------------------------------------- */

SurfReactProb::SurfReactProb(SPARTA *sparta, int narg, char **arg) :
  SurfReact(sparta, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal surf_react prob command");

  // initialize reaction data structs

  nlist_prob = maxlist_prob = 0;
  rlist = NULL;
  reactions = NULL;
  indices = NULL;

  // read reaction file

  readfile(arg[2]);

  // setup the reaction tallies

  nsingle = ntotal = 0;

  nlist = nlist_prob;
  tally_single = new int[nlist];
  tally_total = new int[nlist];
  tally_single_all = new int[nlist];
  tally_total_all = new int[nlist];

  size_vector = 2 + 2*nlist;

  // initialize RNG

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfReactProb::~SurfReactProb()
{
  if (copy) return;

  delete random;

  if (rlist) {
    for (int i = 0; i < maxlist_prob; i++) {
      for (int j = 0; j < rlist[i].nreactant; j++)
        delete [] rlist[i].id_reactants[j];
      for (int j = 0; j < rlist[i].nproduct; j++)
        delete [] rlist[i].id_products[j];
      delete [] rlist[i].id_reactants;
      delete [] rlist[i].id_products;
      delete [] rlist[i].reactants;
      delete [] rlist[i].products;
      delete [] rlist[i].coeff;
      delete [] rlist[i].id;
    }
    memory->destroy(rlist);
  }

  memory->destroy(reactions);
  memory->destroy(indices);
}

/* ---------------------------------------------------------------------- */

void SurfReactProb::init()
{
  SurfReact::init();
  init_reactions();
}

/* ----------------------------------------------------------------------
   select surface reaction to perform for particle with ptr IP on surface
   return which reaction 1 to N, 0 = no reaction
   if dissociation, add particle and return ptr JP
------------------------------------------------------------------------- */


int SurfReactProb::react(Particle::OnePart *&ip, int, double *,
                         Particle::OnePart *&jp, int &)
{
  int n = reactions[ip->ispecies].n;
  if (n == 0) return 0;
  int *list = reactions[ip->ispecies].list;
  // probablity to compare to reaction probability
  double v3[3];
  double x3[3];
  double react_prob = 0.0;
  OneReaction *r;
  for (int i = 0; i < n; i++) {
    r = &rlist[list[i]];
    r->reactants[i] = particle->find_species(r->id_reactants[i]);
    // get velocities
    int isp = r->reactants[i];
    double mass_incident = particle->species[isp].mass;
    double charge_incident = particle->species[isp].charge;
    for (int j = 0; j < 3; j++) {
      v3[j] = ip->v[j];
      x3[j] = ip->x[j];
    }
    double species_mass_amu = particle->species[isp].molwt;

    // PlasmaData plasmaData = update->readPlasmaData("data/plasmadata.h5");
    PlasmaParams params = update->interpolatePlasmaData( x3[0], x3[1]);
    double B[3];
    B[0] = params.b_r ;
    B[1] = params.b_z;
    B[2] = params.b_phi;
    double te = params.temp_e;
    // printf("te = %g\n",te);
    double ti = params.temp_i;
    // printf("te = %g\n",te);
    double ne = params.dens_e;
double random_prob = random->uniform();
    // double sheathEnergy = 3.0 * te * charge_incident;
    // printf("velocities of incident particle = %g %g %g\n",v3[0],v3[1],v3[2]);
    double normVel  = std::sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
    // printf("normVel = %g\n",normVel);
    // printf("mass_incident = %g\n",mass_incident);
    //     double energy_incident_J = 0.5 * mass_incident * MathExtra::lensq3(v3)  ; 
    //     printf("energy_incident_J = %g\n",energy_incident_J);
  //  printf("mass %g\n", mass_incident);

    double energy_incident_eV = 0.5 * mass_incident * normVel * joule2ev ; 
    double energy_incident_eV_2 =  3.0 * te * charge_incident + 2 * charge_incident * te ;
    energy_incident_eV = energy_incident_eV_2;
    // printf("energy_incident_eV energy_incident_eV_2 = %g %g\n",energy_incident_eV, energy_incident_eV_2);

    // energy_incident = energy_incident_eV; 
    int icell = ip->icell;
    Grid::ChildCell *cells = grid->cells;
    double *lo, *hi, xc[3];
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    // cell xc is center of cell
    xc[0] = 0.5*(lo[0] + hi[0]);
    xc[1] = 0.5*(lo[1] + hi[1]);

    // auto [eField, angle] = update->getElectricPotential(xc[0], xc[1], charge_incident, mass_incident, params);

    GeometryTools::PointWall P = {0, xc[0], xc[1]};
    K::Vector_3 velocity(v3[0], v3[1], v3[2]);
    K::Vector_3 magneticField(B[0], B[1], B[2]);
    auto [angle_velocity_normal, angle_magneticField_normal, minDist, normal_vector] = computeAnglesWithNormal(P, velocity, magneticField);

   // convert angle to radians
    double angle_velocity_normal_rad = angle_velocity_normal * M_PI / 180.0;
    // printf("energy_incident_eV = %g\n",energy_incident_eV);

    if (angle_velocity_normal > 90.0) {
      angle_velocity_normal = std::abs(180.0 - angle_velocity_normal);
    }
    if (angle_velocity_normal < .0) {
      angle_velocity_normal =  0.0;
    }
    // printf("angle_velocity_normal = %g\n",angle_velocity_normal);

    // SurfaceDataParams  params2 = interpolateSurfaceData(energy_incident_eV, angle_velocity_normal, species_mass_amu);
      SurfaceDataParams  params2 = interpolateSurfaceData(energy_incident_eV_2, angle_velocity_normal, species_mass_amu);

  double target_Z1 = 74.;
  double target_M1 = 184.;
  double ion_Z = charge_incident;
  double ion_M = species_mass_amu;
  // printf("mass = %g\n",ion_M);
  double converteV2K = 11604.505;
  double U_S = 2.0; //8.79;
  // printf("tem = %g\n",te);
  double vrm = sqrt(0.5 * update->boltz * U_S * converteV2K / target_M1 / mproton);
  // printf("energy_incident_eV = %g\n",energy_incident_eV);
  // get reflection coefficient and sputtering coefficient
  double refl = 0;
  double sput = 0;
    if (energy_incident_eV < 10.0) {
    refl =   wierzbicki_biersack(target_Z1, target_M1, ion_Z, ion_M, energy_incident_eV);
    sput=  yamamura(target_Z1, target_M1, ion_Z, ion_M, energy_incident_eV);
    }else{
        refl = params2.rfyld;
      sput = params2.spyld;
    }
    
    

    double reflection_coefficient = refl; //
    double sputtering_coefficient = sput; 
    double react_prob_reflection = 0.0;
    double react_prob_sputtering = 0.0;

    double total_coefficient = reflection_coefficient + sputtering_coefficient;
// printf("total_coefficient = %g\n",total_coefficient);
// printf("reflection_coefficient = %g\n",reflection_coefficient);
// printf("sputtering_coefficient = %g\n",sputtering_coefficient);

if (total_coefficient == 0.0) {
    // Assuming you want to set both to zero if total_coefficient is zero
    react_prob_reflection = 0.0;
    react_prob_sputtering = 0.0;
} else {
    if (reflection_coefficient == 0.0) {
        react_prob_reflection = 0.0;
    } else {
        react_prob_reflection = reflection_coefficient / total_coefficient;
    }

    if (sputtering_coefficient == 0.0) {
        react_prob_sputtering = 0.0;
    } else {
        react_prob_sputtering = sputtering_coefficient / total_coefficient;
    }
}

    if (total_coefficient  <=0){
      // printf("total_coefficient = %g\n",total_coefficient);


        // for (int j = 0; j < 3; j++) {
        //   v3[j] = 0.0;
        // }
        // memcpy(ip->v,v3,3*sizeof(double));
        // return 0;

        

        // remove particle
        ip = NULL;
        // printf("removing particle\n");
        return 0;

    }
    else{    
      nsingle++;
      tally_single[list[i]]++;
     if (react_prob_sputtering > random_prob) { 
      switch (r->type) {
       case DISSOCIATION:
        {
          // printf("sputtering/dissociation\n");
        
          // // incident species is dissociated into 1 product species
          // printf("incident species type = %d\n",isp);
          // printf("energy of incident particle = %g\n",energy_incident_eV);
          // printf("angle of incident particle = %g\n",angle_velocity_normal);
          // printf("mass of incident particle = %g\n",species_mass_amu);
          // // print charge
          // printf("charge of incident particle = %g\n",charge_incident);
          // // print te
          // printf("te = %g\n",te);
          double x[3],v[3];
          ip->ispecies = r->products[0];
          // printf("product species type = %d\n",ip->ispecies);
          int id = MAXSMALLINT*random->uniform();
          memcpy(x,ip->x,3*sizeof(double));
          memcpy(v,ip->v,3*sizeof(double));
          // print velocity
          // printf("velocity of dissociated particle = %g %g %g\n",v[0],v[1],v[2]);
          Particle::OnePart *particles = particle->particles;;
          int reallocflag =
            particle->add_particle(id,r->products[1],ip->icell,x,v,0.0,0.0);
          if (reallocflag) ip = particle->particles + (ip - particles);
          jp = &particle->particles[particle->nlocal-1];
          return (list[i] + 1);
        }
      }
    }
    else{
      // get species mass
      // printf("reflection\n");
      // incident species mass and charge
      // printf("incident species id = %d\n",ip->ispecies);
      double mass = particle->species[ip->ispecies].molwt;
      // printf("incident mass = %g\n",mass );
      double charge_i = particle->species[ip->ispecies].charge;
      // printf("incident charge = %g\n",charge_i);
      // if mass = 16 then neutralize to oxgyen species 9 and if mass 184 tungsten species 14
      if (mass == 16.0)
      { 
        // printf("neutralizing to oxygen\n");
        ip->ispecies = 8;
      }
      if (mass == 184.0) ip->ispecies = 14;
      //  printf("reflected species id = %d\n",ip->ispecies);
      double charge2 = particle->species[ip->ispecies].charge;
      double mass2 = particle->species[ip->ispecies].molwt  ;

      double theta = 2.0 * random->uniform() * M_PI;
      double vtangent = vrm * sqrt(-log(random->uniform()));
      double vx = vtangent * sin(theta);
      double vy = vtangent * cos(theta);
      double vz = 0.0;
      double v[3];
      v[0] = vx;
      v[1] = -std::abs(vy);
      v[2] = vz;
            for (int j = 0; j < 3; j++) {
        v[j] = ip->v[j];
      }
      memcpy(ip->v,v,3*sizeof(double));

        return (list[i] + 1);
  }
  }
    }
  return 0;
}

/* ---------------------------------------------------------------------- */

char *SurfReactProb::reactionID(int m)
{
  return rlist[m].id;
}

/* ---------------------------------------------------------------------- */

int SurfReactProb::match_reactant(char *species, int m)
{
  for (int i = 0; i < rlist[m].nreactant; i++)
    if (strcmp(species,rlist[m].id_reactants[i]) == 0) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

int SurfReactProb::match_product(char *species, int m)
{
  for (int i = 0; i < rlist[m].nproduct; i++)
    if (strcmp(species,rlist[m].id_products[i]) == 0) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

void SurfReactProb::init_reactions()
{
  // convert species IDs to species indices
  // flag reactions as active/inactive depending on whether all species exist

  for (int m = 0; m < nlist_prob; m++) {
    OneReaction *r = &rlist[m];
    r->active = 1;
    for (int i = 0; i < r->nreactant; i++) {
      r->reactants[i] = particle->find_species(r->id_reactants[i]);
      if (r->reactants[i] < 0) {
        r->active = 0;
        break;
      }
    }
    for (int i = 0; i < r->nproduct; i++) {
      r->products[i] = particle->find_species(r->id_products[i]);
      if (r->products[i] < 0) {
        r->active = 0;
        break;
      }
    }
  }

  // count possible reactions for each species

  memory->destroy(reactions);
  int nspecies = particle->nspecies;
  reactions = memory->create(reactions,nspecies,
                             "surf_react:reactions");

  for (int i = 0; i < nspecies; i++) reactions[i].n = 0;

  int n = 0;
  for (int m = 0; m < nlist_prob; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active) continue;
    int i = r->reactants[0];
    reactions[i].n++;
    n++;
  }

  // allocate indices = entire list of reactions for all I species

  memory->destroy(indices);
  memory->create(indices,n,"surf_react:indices");

  // reactions[i].list = offset into full indices vector

  int offset = 0;
  for (int i = 0; i < nspecies; i++) {
    reactions[i].list = &indices[offset];
    offset += reactions[i].n;
  }

  // reactions[i].list = indices of possible reactions for each species

  for (int i = 0; i < nspecies; i++) reactions[i].n = 0;

  for (int m = 0; m < nlist_prob; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active) continue;
    int i = r->reactants[0];
    reactions[i].list[reactions[i].n++] = m;
  }

  // check that summed reaction probabilities for each species <= 1.0

  double sum;
  for (int i = 0; i < nspecies; i++) {
    sum = 0.0;
    for (int j = 0; j < reactions[i].n; j++)
      sum += rlist[reactions[i].list[j]].coeff[0];
    if (sum > 1.0)
      error->all(FLERR,"Surface reaction probability for a species > 1.0");
  }
}

/* ---------------------------------------------------------------------- */

void SurfReactProb::readfile(char *fname)
{
  int n,n1,n2,eof;
  char line1[MAXLINE],line2[MAXLINE];
  char *word;
  OneReaction *r;

  // proc 0 opens file

  if (comm->me == 0) {
    fp = fopen(fname,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open reaction file %s",fname);
      error->one(FLERR,str);
    }
  }

  // read reactions one at a time and store their info in rlist

  while (1) {
    if (comm->me == 0) eof = readone(line1,line2,n1,n2);
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;

    MPI_Bcast(&n1,1,MPI_INT,0,world);
    MPI_Bcast(&n2,1,MPI_INT,0,world);
    MPI_Bcast(line1,n1,MPI_CHAR,0,world);
    MPI_Bcast(line2,n2,MPI_CHAR,0,world);

    if (nlist_prob == maxlist_prob) {
      maxlist_prob += DELTALIST;
      rlist = (OneReaction *)
        memory->srealloc(rlist,maxlist_prob*sizeof(OneReaction),
                         "surf_react:rlist");
      for (int i = nlist_prob; i < maxlist_prob; i++) {
        r = &rlist[i];
        r->nreactant = r->nproduct = 0;
        r->id_reactants = new char*[MAXREACTANT];
        r->id_products = new char*[MAXPRODUCT];
        r->reactants = new int[MAXREACTANT];
        r->products = new int[MAXPRODUCT];
        r->coeff = new double[MAXCOEFF];
        r->id = NULL;
      }
    }

    r = &rlist[nlist_prob];

    int side = 0;
    int species = 1;

    n = strlen(line1) - 1;
    r->id = new char[n+1];
    strncpy(r->id,line1,n);
    r->id[n] = '\0';

    word = strtok(line1," \t\n");

    while (1) {
      if (!word) {
        if (side == 0) error->all(FLERR,"Invalid reaction formula in file");
        if (species) error->all(FLERR,"Invalid reaction formula in file");
        break;
      }
      if (species) {
        species = 0;
        if (side == 0) {
          if (r->nreactant == MAXREACTANT)
            error->all(FLERR,"Too many reactants in a reaction formula");
          n = strlen(word) + 1;
          r->id_reactants[r->nreactant] = new char[n];
          strcpy(r->id_reactants[r->nreactant],word);
          r->nreactant++;
        } else {
          if (r->nreactant == MAXPRODUCT)
            error->all(FLERR,"Too many products in a reaction formula");
          n = strlen(word) + 1;
          r->id_products[r->nproduct] = new char[n];
          strcpy(r->id_products[r->nproduct],word);
          r->nproduct++;
        }
      } else {
        species = 1;
        if (strcmp(word,"+") == 0) {
          word = strtok(NULL," \t\n");
          continue;
        }
        if (strcmp(word,"-->") != 0)
          error->all(FLERR,"Invalid reaction formula in file");
        side = 1;
      }
      word = strtok(NULL," \t\n");
    }

    // replace single NULL product with no products

    if (r->nproduct == 1 && strcmp(r->id_products[0],"NULL") == 0) {
      delete [] r->id_products[0];
      r->id_products[0] = NULL;
      r->nproduct = 0;
    }

    word = strtok(line2," \t\n");
    if (!word) error->all(FLERR,"Invalid reaction type in file");
    if (word[0] == 'D' || word[0] == 'd') r->type = DISSOCIATION;
    else if (word[0] == 'E' || word[0] == 'e') r->type = EXCHANGE;
    else if (word[0] == 'R' || word[0] == 'r') r->type = RECOMBINATION;
    else if (word[0] == 'P' || word[0] == 'p') r->type = SPUTTERING;
    else error->all(FLERR,"Invalid reaction type in file");

    // check that reactant/product counts are consistent with type

    if (r->type == DISSOCIATION) {
      if (r->nreactant != 1 || r->nproduct != 1)
        error->all(FLERR,"Invalid dissociation reaction");
    } else if (r->type == EXCHANGE) {
      if (r->nreactant != 1 || r->nproduct != 1)
        error->all(FLERR,"Invalid exchange reaction");
    } else if (r->type == SPUTTERING) {
      if (r->nreactant != 1 || r->nproduct != 1)
        error->all(FLERR,"Invalid sputtering reaction");
    }     else if (r->type == RECOMBINATION) {
      if (r->nreactant != 1 || r->nproduct != 0)
        error->all(FLERR,"Invalid recombination reaction");
    }

    word = strtok(NULL," \t\n");
    if (!word) error->all(FLERR,"Invalid reaction style in file");
    if (word[0] == 'S' || word[0] == 's') r->style = SIMPLE;
    else error->all(FLERR,"Invalid reaction style in file");

    if (r->style == SIMPLE) r->ncoeff = 1;

    for (int i = 0; i < r->ncoeff; i++) {
      word = strtok(NULL," \t\n");
      if (!word) error->all(FLERR,"Invalid reaction coefficients in file");
      r->coeff[i] = input->numeric(FLERR,word);
    }

    word = strtok(NULL," \t\n");
    if (word) error->all(FLERR,"Too many coefficients in a reaction formula");

    nlist_prob++;
  }

  if (comm->me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   read one reaction from file
   reaction = 2 lines
   return 1 if end-of-file, else return 0
------------------------------------------------------------------------- */

int SurfReactProb::readone(char *line1, char *line2, int &n1, int &n2)
{
  char *eof;
  while ((eof = fgets(line1,MAXLINE,fp))) {
    size_t pre = strspn(line1," \t\n");
    if (pre == strlen(line1) || line1[pre] == '#') continue;
    eof = fgets(line2,MAXLINE,fp);
    if (!eof) break;
    n1 = strlen(line1) + 1;
    n2 = strlen(line2) + 1;
    return 0;
  }

  return 1;
}


// std::vector<double> analytical_func(const std::vector<double>& E) {
//     double A = 0.3;
//     double B = 0.06;
//     double C = 0.05;
//     double epsilon = 1e-10;
//     std::vector<double> result(E.size(), 0.0);  // Initialize a vector for the result

//     for (size_t i = 0; i < E.size(); ++i) {
//         if (E[i] < 28) {
//             result[i] = 0;
//         } else if (E[i] >= 28 && E[i] < 300) {
//             result[i] = A * log(B * E[i] + C + epsilon);
//         } else {
//             result[i] = 1;
//         }
//     }

//     return result;
// }

double SurfReactProb::wierzbicki_biersack(double target_Z1, double target_M1, double ion_Z, double ion_M, double energy_eV) {
    
    // printf("energy_eV = %g\n",energy_eV);
    // printf("target_Z1 = %g\n",target_Z1);
    // printf("target_M1 = %g\n",target_M1);
    // printf("ion_Z = %g\n",ion_Z);
    // printf("ion_M = %g\n",ion_M);

    double Z1 = ion_Z;
    double Z2 = target_Z1;
    double M1 = ion_M;
    double M2 = target_M1;
    double energy_keV = energy_eV / 1E3;
 
    double reduced_energy = 32.55 * energy_keV * M2 / ((M1 + M2) * Z1 * Z2 * (pow(Z1, 0.23) + pow(Z2, 0.23)));
    double mu = M2 / M1;
    // printf("mu = %g\n",mu);
 
    double a1 = 5.9638;
    double b1 = 0.0646;
    double c1 = 52.211;
 
    double a2 = 1.26E-3;
    double b2 = -0.9305;
    double c2 = 1.235;
 
    double RN_mu = std::exp(a1 * std::sqrt(1.0 - b1 * std::pow(std::log(mu / c1), 2.0)));
    double RN_e = a2 * std::exp(b2 * std::pow(std::log(reduced_energy + 1.0), c2));

   double R = RN_mu * RN_e;
  //  printf("R = %g\n",R);
   // check if R is nan
    if (std::isnan(R)) {
        R = 0.0;
    }
    // printf("setting R  to zero = %g\n",R);
    return R;
}
 


double  SurfReactProb::getElectricPotential(double minDistance, double larmorRadius, double te, double ne)  {
        double m_e = 9.10938356e-31;
        double m_p = 1.6726219e-27;
        double debyeLength = 7.43 * std::sqrt(te / ne);
        double potential_wall = 0.5 * log ( 2 * M_PI * (m_e / m_p) * 2) ;
        double thresholdForSheathModel =  20; // 200 debye lengths
        double Lmps = 2 * larmorRadius;

        double result = potential_wall  *  std::exp(- abs(minDistance  / Lmps )) * te  ;
        return abs(result);
    }



double SurfReactProb::yamamura(double target_Z1, double target_M1, double ion_Z, double ion_M, double energy_eV) {
    double z1 = ion_Z;
    double z2 = target_Z1;
    double m1 = ion_M;
    double m2 = target_M1;
    double Us = 8.79;
    double Q = 1.10;

    double reduced_mass_2 = m2 / (m1 + m2);
    double reduced_mass_1 = m1 / (m1 + m2);

    // Lindhard's reduced energy
    double reduced_energy = 0.03255 / (z1 * z2 * std::sqrt(z1 * z1 * z1 + z2 * z2 * z2)) * reduced_mass_2 * energy_eV;

    // Yamamura empirical constants
    double K = 8.478 * z1 * z2 / std::sqrt(z1 * z1 * z1 + z2 * z2 * z2) * reduced_mass_1;
    double a_star = 0.08 + 0.164 * std::pow(m2 / m1, 0.4) + 0.0145 * std::pow(m2 / m1, 1.29);

    // Sputtering threshold energy
    double Eth = (1.9 + 3.8 * (m1 / m2) + 0.134 * std::pow(m2 / m1, 1.24)) * Us;

    // Lindhard-Scharff-Schiott nuclear cross section
    double sn = 3.441 * std::sqrt(reduced_energy) * std::log(reduced_energy + 2.718) / (1. + 6.355 * std::sqrt(reduced_energy) + reduced_energy * (-1.708 + 6.882 * std::sqrt(reduced_energy)));

    // Lindhard-Scharff electronic cross section
    double k = 0.079 * std::pow(m1 + m2, 1.5) / (std::pow(m1, 1.5) * std::sqrt(m2)) * std::pow(z1, 2.0/3.0) * std::sqrt(z2) / std::pow(z1 * z1 * z1 + z2 * z2 * z2, 3.0/4.0);
    double se = k * std::sqrt(reduced_energy);
   double Y = 0.42 * a_star * Q * K * sn / Us / (1. + 0.35 * Us * se) * std::pow(1. - std::sqrt(Eth / energy_eV), 2.8);
    // check if Y is nan
    if (std::isnan(Y)) {
        Y = 0.0;
    }
    return Y;
    // return 0.42 * a_star * Q * K * sn / Us / (1. + 0.35 * Us * se) * std::pow(1. - std::sqrt(Eth / energy_eV), 2.8);
}


SurfaceData SurfReactProb::readSurfaceData(const std::string& filePath) {

    auto it = surfaceDataCache.find(filePath);
    if (it != surfaceDataCache.end()) {
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

    auto read1DDataSet = [&file](const std::string& datasetPath) {
    H5::DataSet ds = file.openDataSet(datasetPath);
    H5::DataSpace space = ds.getSpace();
    std::vector<hsize_t> dims(1);
    space.getSimpleExtentDims(dims.data(), NULL);

    std::vector<float> rawData(dims[0]);
    ds.read(rawData.data(), H5::PredType::NATIVE_FLOAT);

    return rawData;
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
    SurfaceData data;
     data.energy = read1DDataSet("E");
     data.angle = read1DDataSet("A");
     data.rfyld = read2DDataSet("rfyld");
     data.spyld = read2DDataSet("spyld");
    file.close();

    // Cache the content
    surfaceDataCache[filePath] = data;

    return data;
}

std::string getFilenameFromMass(double mass) {
    const double mass_O = 15.999; 
    const double mass_W = 184.0;
    // printf("mass = %g\n",mass);
    if (std::fabs(mass - mass_O) < 1e-2) {
        return "O_on_W.h5";
    } else if (std::fabs(mass - mass_W) < 1e-2) {
        return "W_on_W.h5";
    } else {
        printf("Material symbol not found for mass and charge = %g\n",mass);
        exit(0);
    }
}

double interpolate(const std::vector<std::vector<double>>& data_2d,
                   const std::vector<double>& energy_values,
                   const std::vector<double>& angle_values,
                   double energy_new, double angle_new) {

    // Ensure the new energy and angle values are within range
    if (energy_new <= energy_values.front() || energy_new >= energy_values.back() ||
        angle_new <= angle_values.front() || angle_new >= angle_values.back()) {
        return 0.0;
    }

    const gsl_interp2d_type * T = gsl_interp2d_bilinear;
    gsl_interp2d * interp = gsl_interp2d_alloc(T, energy_values.size(), angle_values.size());

    // Flatten the 2D data for gsl_interp2d_init
    std::vector<double> flattened_data;
    for (const auto& row : data_2d) {
        for (double value : row) {
            flattened_data.push_back(value);
        }
    }
    gsl_interp2d_init(interp, energy_values.data(), angle_values.data(), flattened_data.data(), energy_values.size(), angle_values.size());
    double interpolated_value = gsl_interp2d_eval(interp, energy_values.data(), angle_values.data(), flattened_data.data(), energy_new, angle_new, NULL, NULL);
    gsl_interp2d_free(interp); // Don't forget to free the allocated memory!

    return interpolated_value;
}

SurfaceDataParams SurfReactProb::interpolateSurfaceData(double energy_val, double angle_val, double mass) {
    // Check the interpolated data cache for existing data
    std::tuple<double, double, double> key = {energy_val, angle_val, mass};
    auto cache_it = interpolatedDataCache.find(key);
    if (cache_it != interpolatedDataCache.end()) {
        return cache_it->second;
    }

    // Read plasma data and avoid copying vectors r and z
    std::string filename = getFilenameFromMass(mass);
    const SurfaceData& data = readSurfaceData("data/" + filename);

    // Convert the energy and angle data from float to double
    std::vector<double> energy_double(data.energy.begin(), data.energy.end());
    std::vector<double> angle_double(data.angle.begin(), data.angle.end());

    // Create a SurfaceDataParams object to store the interpolated data
    SurfaceDataParams interpolated_data;
    std::vector<std::vector<double>> rfyld_double(data.rfyld.size());
    for (size_t i = 0; i < data.rfyld.size(); ++i) {
        rfyld_double[i].assign(data.rfyld[i].begin(), data.rfyld[i].end());
    }
    std::vector<std::vector<double>> spyld_double(data.spyld.size());
    for (size_t i = 0; i < data.spyld.size(); ++i) {
        spyld_double[i].assign(data.spyld[i].begin(), data.spyld[i].end());
    }

    interpolated_data.rfyld = interpolate(rfyld_double, energy_double, angle_double, energy_val, angle_val);
    interpolated_data.spyld = interpolate(spyld_double, energy_double, angle_double, energy_val, angle_val);

    // Cache the interpolated data
    interpolatedDataCache[key] = interpolated_data;

    return interpolated_data;
}
