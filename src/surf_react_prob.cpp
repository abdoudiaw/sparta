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
#include "/opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense"
#include <random>
#include <iostream>
#include "surf.h"
#include "grid.h"
using namespace SPARTA_NS;

enum{DISSOCIATION,EXCHANGE,RECOMBINATION,SPUTTERING};
enum{SIMPLE};

#define MAXREACTANT 1
#define MAXPRODUCT 2
#define MAXCOEFF 2

#define MAXLINE 1024
#define DELTALIST 16

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
  double random_prob = random->uniform();
  double joule2ev = 6.24150934e18;
  // target info
  double mass_target = 184.0;
  double charge_target = 74.0;
  double mproton = 1.67262158e-27;
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
    std::vector<double> plasmaData = update->get_density_temperature(x3);
    double B[3];
    update->get_magnetic_field( x3, B);
    double ne = plasmaData[0];
    double te = plasmaData[1];
    double sheathEnergy = 3.0 * te * charge_incident;
    // get reflection coefficient
    double energy_incident = 0.5 * mass_incident * MathExtra::lensq3(v3) * joule2ev + sheathEnergy;
    // get sputtering coefficient
    double angle_max = 2.0;
    double angle_min = 0.0;
    double angle_degrees = angle_min + (random->uniform() * (angle_max - angle_min));

    double reflection_coefficient = update->get_reflection_coefficient(energy_incident,angle_degrees, mass_incident / mproton, charge_incident );
    double sputtering_coefficient = update->get_sputtering_coefficient(energy_incident,angle_degrees, mass_incident / mproton, charge_incident );

    double react_prob_reflection = reflection_coefficient / (reflection_coefficient + sputtering_coefficient);
    double react_prob_sputtering = sputtering_coefficient / (reflection_coefficient + sputtering_coefficient);
    

     if (react_prob_reflection > random_prob) { 
      // printf("reflection\n");
      nsingle++;
      tally_single[list[i]]++;
      switch (r->type) {
       case DISSOCIATION:
        {
          double x[3],v[3];
          ip->ispecies = r->products[0];
          printf("dissociation\n");
          printf("ip->ispecies = %d\n",ip->ispecies);
          int id = MAXSMALLINT*random->uniform();
          memcpy(x,ip->x,3*sizeof(double));
          memcpy(v,ip->v,3*sizeof(double));
          Particle::OnePart *particles = particle->particles;
          int reallocflag =
            particle->add_particle(id,r->products[1],ip->icell,x,v,0.0,0.0);
          if (reallocflag) ip = particle->particles + (ip - particles);
          jp = &particle->particles[particle->nlocal-1];
          return (list[i] + 1);
        }
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


double SurfReactProb::wierzbicki_biersack(double target_Z1, double target_M1, double ion_Z, double ion_M, double energy_eV) {
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

 
    return RN_mu * RN_e;
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