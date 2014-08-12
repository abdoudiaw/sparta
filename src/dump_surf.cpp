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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "dump_surf.h"
#include "update.h"
#include "domain.h"
#include "surf.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// customize by adding keyword

enum{ID,V1X,V1Y,V1Z,V2X,V2Y,V2Z,V3X,V3Y,V3Z,
     COMPUTE,FIX,VARIABLE};
enum{INT,DOUBLE,CELLINT,STRING};    // many files

enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain

#define INVOKED_PER_SURF 32
#define CHUNK 8

/* ---------------------------------------------------------------------- */

DumpSurf::DumpSurf(SPARTA *sparta, int narg, char **arg) :
  Dump(sparta, narg, arg)
{
  if (narg == 4) error->all(FLERR,"No dump surf attributes specified");

  clearstep = 1;
  buffer_allow = 1;
  buffer_flag = 1;

  dimension = domain->dimension;

  nevery = atoi(arg[2]);

  // scan values and set size_one
  // array entries may expand into multiple fields
  // could use ioptional to add on keyword/arg pairs

  ioptional = parse_fields(narg,arg);
  if (ioptional < narg)
    error->all(FLERR,"Invalid attribute in dump surf command");
  size_one = nfield;

  // setup format strings

  vformat = new char*[nfield];

  format_default = new char[3*nfield+1];
  format_default[0] = '\0';

  for (int i = 0; i < nfield; i++) {
    if (vtype[i] == INT) strcat(format_default,"%d ");
    else if (vtype[i] == DOUBLE) strcat(format_default,"%g ");
    vformat[i] = NULL;
  }

  // setup column string using field2arg
  // add column subscripts to array args that were expanded

  int n = 0;
  for (int i = 0; i < nfield; i++) n += strlen(arg[field2arg[i]]) + 6;
  columns = new char[n];
  columns[0] = '\0';
  char subscript[6];
  for (int i = 0; i < nfield; i++) {
    strcat(columns,arg[field2arg[i]]);
    if (argindex[i] > 0 && columns[strlen(columns)-1] != ']') {
      sprintf(subscript,"[%d]",argindex[i]);
      strcat(columns,subscript);
    }
    strcat(columns," ");
  }
}

/* ---------------------------------------------------------------------- */

DumpSurf::~DumpSurf()
{
  memory->sfree(pack_choice);
  memory->destroy(vtype);
  memory->destroy(field2arg);
  memory->destroy(field2index);
  memory->destroy(argindex);

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  memory->sfree(id_compute);
  delete [] compute;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  memory->sfree(id_fix);
  delete [] fix;

  for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
  memory->sfree(id_variable);
  delete [] variable;
  for (int i = 0; i < nvariable; i++) memory->destroy(vbuf[i]);
  delete [] vbuf;

  for (int i = 0; i < nfield; i++) delete [] vformat[i];
  delete [] vformat;

  delete [] columns;
}

/* ---------------------------------------------------------------------- */

void DumpSurf::init_style()
{
  delete [] format;
  char *str;
  if (format_user) str = format_user;
  else str = format_default;

  int n = strlen(str) + 1;
  format = new char[n];
  strcpy(format,str);

  // tokenize the format string and add space at end of each format element

  char *ptr;
  for (int i = 0; i < nfield; i++) {
    if (i == 0) ptr = strtok(format," \0");
    else ptr = strtok(NULL," \0");
    delete [] vformat[i];
    vformat[i] = new char[strlen(ptr) + 2];
    strcpy(vformat[i],ptr);
    vformat[i] = strcat(vformat[i]," ");
  }

  // setup boundary string

  int m = 0;
  for (int idim = 0; idim < 3; idim++) {
    for (int iside = 0; iside < 2; iside++) {
      if (domain->bflag[idim*2+iside] == OUTFLOW) boundstr[m++] = 'o';
      else if (domain->bflag[idim*2+iside] == PERIODIC) boundstr[m++] = 'p';
      else if (domain->bflag[idim*2+iside] == REFLECT) boundstr[m++] = 'r';
      else if (domain->bflag[idim*2+iside] == AXISYM) boundstr[m++] = 'a';
      else if (domain->bflag[idim*2+iside] == SURFACE) boundstr[m++] = 's';
    }
    boundstr[m++] = ' ';
  }
  boundstr[8] = '\0';

  // setup function ptrs

  if (binary) header_choice = &DumpSurf::header_binary;
  else header_choice = &DumpSurf::header_item;

  if (binary) write_choice = &DumpSurf::write_binary;
  else if (buffer_flag == 1) write_choice = &DumpSurf::write_string;
  else write_choice = &DumpSurf::write_text;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) 
      error->all(FLERR,"Could not find dump surf compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find dump surf fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->per_surf_freq)
      error->all(FLERR,"Dump surf and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0) 
      error->all(FLERR,"Could not find dump surf variable name");
    variable[i] = ivariable;
  }

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::header_binary(bigint ndump)
{
  fwrite(&update->ntimestep,sizeof(bigint),1,fp);
  fwrite(&ndump,sizeof(bigint),1,fp);
  fwrite(domain->bflag,6*sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&nfield,sizeof(int),1,fp);
  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::header_item(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF SURFS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: SURFS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

int DumpSurf::count()
{
  // invoke Computes for per-surf quantities

  if (ncompute) {
    for (int i = 0; i < ncompute; i++)
      if (!(compute[i]->invoked_flag & INVOKED_PER_SURF)) {
	compute[i]->compute_per_surf();
	compute[i]->invoked_flag |= INVOKED_PER_SURF;
      }
  }

  // evaluate surf-style Variables for per-surf quantities

  if (nvariable)
    for (int i = 0; i < nvariable; i++)
      input->variable->compute_surf(variable[i],vbuf[i],1,0);

  return surf->nlocal;
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack()
{
  nslocal = surf->nlocal;
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpSurf::write_text(int n, double *mybuf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT) fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == DOUBLE) fprintf(fp,vformat[j],mybuf[m]);
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpSurf::parse_fields(int narg, char **arg)
{
  // initialize per-field lists

  pack_choice = NULL;
  vtype = NULL;
  field2arg = NULL;
  field2index = NULL;
  argindex = NULL;

  int maxfield = narg - 4;
  allocate_values(maxfield);
  nfield = 0;

  ncompute = 0;
  id_compute = NULL;
  compute = NULL;

  nfix = 0;
  id_fix = NULL;
  fix = NULL;

  nvariable = 0;
  id_variable = NULL;
  variable = NULL;
  vbuf = NULL;

  // customize by adding to if statement

  for (int iarg = 4; iarg < narg; iarg++) {
    if (nfield == maxfield) {
      maxfield += CHUNK;
      allocate_values(maxfield);
    }

    argindex[nfield] = -1;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[nfield] = &DumpSurf::pack_id;
      vtype[nfield] = INT;
      field2arg[nfield] = iarg;
      nfield++;

    } else if (strcmp(arg[iarg],"v1x") == 0) {
      pack_choice[nfield] = &DumpSurf::pack_v1x;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;
      nfield++;
    } else if (strcmp(arg[iarg],"v1y") == 0) {
      pack_choice[nfield] = &DumpSurf::pack_v1y;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;
      nfield++;
    } else if (strcmp(arg[iarg],"v1z") == 0) {
      if (dimension == 2) 
	error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[nfield] = &DumpSurf::pack_v1z;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;
      nfield++;
    } else if (strcmp(arg[iarg],"v2x") == 0) {
      pack_choice[nfield] = &DumpSurf::pack_v2x;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;
      nfield++;
    } else if (strcmp(arg[iarg],"v2y") == 0) {
      pack_choice[nfield] = &DumpSurf::pack_v2y;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;
      nfield++;
    } else if (strcmp(arg[iarg],"v2z") == 0) {
      if (dimension == 2) 
	error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[nfield] = &DumpSurf::pack_v2z;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;
      nfield++;
    } else if (strcmp(arg[iarg],"v3x") == 0) {
      if (dimension == 2) 
	error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[nfield] = &DumpSurf::pack_v3x;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;
      nfield++;
    } else if (strcmp(arg[iarg],"v3y") == 0) {
      if (dimension == 2) 
	error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[nfield] = &DumpSurf::pack_v3y;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;
      nfield++;
    } else if (strcmp(arg[iarg],"v3z") == 0) {
      if (dimension == 2) 
	error->all(FLERR,"Invalid dump surf field for 2d simulation");
      pack_choice[nfield] = &DumpSurf::pack_v3z;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;
      nfield++;

    // compute value = c_ID
    // if no trailing [], then index = 0, else index = int between []
    // if index = 0 and compute stores array, expand to one value per column

    } else if (strncmp(arg[iarg],"c_",2) == 0) {
      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      int index;
      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Invalid attribute in dump surf command");
	index = atoi(ptr+1);
	*ptr = '\0';
      } else index = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump surf compute ID");
      if (modify->compute[n]->per_surf_flag == 0)
	error->all(FLERR,"Dump surf compute does not compute per-surf info");
      if (index > 0 && modify->compute[n]->size_per_surf_cols == 0)
	error->all(FLERR,
		   "Dump surf compute does not calculate per-surf array");
      if (index > 0 && index > modify->compute[n]->size_per_surf_cols)
	error->all(FLERR,"Dump surf compute array is accessed out-of-range");

      if (index == 0 && modify->compute[n]->size_per_surf_cols > 0) {
	int ncol = modify->compute[n]->size_per_surf_cols;
	for (int i = 0; i < ncol; i++) {
	  if (nfield == maxfield) {
	    maxfield += CHUNK;
	    allocate_values(maxfield);
	  }
	  pack_choice[nfield] = &DumpSurf::pack_compute;
	  vtype[nfield] = DOUBLE;
	  argindex[nfield] = i+1;
	  field2arg[nfield] = iarg;
	  field2index[nfield] = add_compute(suffix);
	  nfield++;
	} 
      } else {
	pack_choice[nfield] = &DumpSurf::pack_compute;
	vtype[nfield] = DOUBLE;
	argindex[nfield] = index;
	field2arg[nfield] = iarg;
	field2index[nfield] = add_compute(suffix);
	nfield++;
      }

      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then index = 0, else index = int between []
    // if index = 0 and compute stores array, expand to one value per column

    } else if (strncmp(arg[iarg],"f_",2) == 0) {
      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      int index;
      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Invalid attribute in dump surf command");
	index = atoi(ptr+1);
	*ptr = '\0';
      } else index = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump surf fix ID");
      if (modify->fix[n]->per_surf_flag == 0)
	error->all(FLERR,"Dump surf fix does not compute per-surf info");
      if (index > 0 && modify->fix[n]->size_per_surf_cols == 0)
	error->all(FLERR,"Dump surf fix does not compute per-surf array");
      if (index > 0 && index > modify->fix[n]->size_per_surf_cols)
	error->all(FLERR,"Dump surf fix array is accessed out-of-range");

      if (index == 0 && modify->fix[n]->size_per_surf_cols > 0) {
	int ncol = modify->fix[n]->size_per_surf_cols;
	for (int i = 0; i < ncol; i++) {
	  if (nfield == maxfield) {
	    maxfield += CHUNK;
	    allocate_values(maxfield);
	  }
	  pack_choice[nfield] = &DumpSurf::pack_fix;
	  vtype[nfield] = DOUBLE;
	  argindex[nfield] = i+1;
	  field2arg[nfield] = iarg;
	  field2index[nfield] = add_fix(suffix);
	  nfield++;
	} 
      } else {
	pack_choice[nfield] = &DumpSurf::pack_fix;
	vtype[nfield] = DOUBLE;
	argindex[nfield] = index;
	field2arg[nfield] = iarg;
	field2index[nfield] = add_fix(suffix);
	nfield++;
      }

      delete [] suffix;

    // variable value = v_name

    } else if (strncmp(arg[iarg],"v_",2) == 0) {
      pack_choice[nfield] = &DumpSurf::pack_variable;
      vtype[nfield] = DOUBLE;
      field2arg[nfield] = iarg;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      argindex[nfield] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump surf variable name");
      if (input->variable->surf_style(n) == 0)
	error->all(FLERR,"Dump surf variable is not surf-style variable");

      field2index[nfield] = add_variable(suffix);
      delete [] suffix;
      nfield++;

    } else return iarg;
  }

  return narg;
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpSurf::add_compute(char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) return icompute;
  
  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  delete [] compute;
  compute = new Compute*[ncompute+1];

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpSurf::add_fix(char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,id_fix[ifix]) == 0) break;
  if (ifix < nfix) return ifix;
  
  id_fix = (char **)
    memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  delete [] fix;
  fix = new Fix*[nfix+1];

  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add Variable to list of Variables used by dump
   return index of where this Variable is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpSurf::add_variable(char *id)
{
  int ivariable;
  for (ivariable = 0; ivariable < nvariable; ivariable++)
    if (strcmp(id,id_variable[ivariable]) == 0) break;
  if (ivariable < nvariable) return ivariable;
  
  id_variable = (char **)
    memory->srealloc(id_variable,(nvariable+1)*sizeof(char *),
		     "dump:id_variable");
  delete [] variable;
  variable = new int[nvariable+1];
  delete [] vbuf;
  vbuf = new double*[nvariable+1];
  for (int i = 0; i <= nvariable; i++) vbuf[i] = NULL;

  int n = strlen(id) + 1;
  id_variable[nvariable] = new char[n];
  strcpy(id_variable[nvariable],id);
  nvariable++;
  return nvariable-1;
}

/* ----------------------------------------------------------------------
   reallocate vectors for each input value, of length N
------------------------------------------------------------------------- */

void DumpSurf::allocate_values(int n)
{
  pack_choice = (FnPtrPack *) 
    memory->srealloc(pack_choice,n*sizeof(FnPtrPack),"dump:pack_choice");
  memory->grow(vtype,n,"dump:vtype");
  memory->grow(field2arg,n,"dump:field2arg");
  memory->grow(field2index,n,"dump:field2index");
  memory->grow(argindex,n,"dump:argindex");
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpSurf::pack_compute(int n)
{
  double *vector = compute[field2index[n]]->vector_surf;
  double **array = compute[field2index[n]]->array_surf;
  int index = argindex[n];

  int *loc2glob;
  int nlocal = compute[field2index[n]]->surfinfo(loc2glob);
  
  if (index == 0)
    surf->collate_vec(nlocal,loc2glob,vector,1,&buf[n],size_one,0);
  else {
    int istride = compute[field2index[n]]->size_per_surf_cols;
    if (array)
      surf->collate_vec(nlocal,loc2glob,&array[0][index-1],istride,
                        &buf[n],size_one,0);
    else
      surf->collate_vec(nlocal,loc2glob,NULL,istride,
                        &buf[n],size_one,0);
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_fix(int n)
{
  double *vector = fix[field2index[n]]->vector_surf;
  double **array = fix[field2index[n]]->array_surf;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < nslocal; i++) {
      buf[n] = vector[i];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nslocal; i++) {
      buf[n] = array[i][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_variable(int n)
{
  double *vector = vbuf[field2index[n]];

  for (int i = 0; i < nslocal; i++) {
    buf[n] = vector[i];
    n += size_one;
  }
}

/* ----------------------------------------------------------------------
   one method for every attribute dump surf can output
   the surf property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

void DumpSurf::pack_id(int n)
{
  int *mysurfs = surf->mysurfs;

  for (int i = 0; i < nslocal; i++) {
    buf[n] = mysurfs[i] + 1;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v1x(int n)
{
  Surf::Point *pts = surf->pts;
  int *mysurfs = surf->mysurfs;

  if (dimension == 2) {
    Surf::Line *lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      buf[n] = pts[lines[mysurfs[i]].p1].x[0];
      n += size_one;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      buf[n] = pts[tris[mysurfs[i]].p1].x[0];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v1y(int n)
{
  Surf::Point *pts = surf->pts;
  int *mysurfs = surf->mysurfs;

  if (dimension == 2) {
    Surf::Line *lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      buf[n] = pts[lines[mysurfs[i]].p1].x[1];
      n += size_one;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      buf[n] = pts[tris[mysurfs[i]].p1].x[1];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v1z(int n)
{
  Surf::Point *pts = surf->pts;
  Surf::Tri *tris = surf->tris;
  int *mysurfs = surf->mysurfs;

  for (int i = 0; i < nslocal; i++) {
    buf[n] = pts[tris[mysurfs[i]].p1].x[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v2x(int n)
{
  Surf::Point *pts = surf->pts;
  int *mysurfs = surf->mysurfs;

  if (dimension == 2) {
    Surf::Line *lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      buf[n] = pts[lines[mysurfs[i]].p2].x[0];
      n += size_one;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      buf[n] = pts[tris[mysurfs[i]].p2].x[0];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v2y(int n)
{
  Surf::Point *pts = surf->pts;
  int *mysurfs = surf->mysurfs;

  if (dimension == 2) {
    Surf::Line *lines = surf->lines;
    for (int i = 0; i < nslocal; i++) {
      buf[n] = pts[lines[mysurfs[i]].p2].x[1];
      n += size_one;
    }
  } else if (dimension == 3) {
    Surf::Tri *tris = surf->tris;
    for (int i = 0; i < nslocal; i++) {
      buf[n] = pts[tris[mysurfs[i]].p2].x[1];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v2z(int n)
{
  Surf::Point *pts = surf->pts;
  Surf::Tri *tris = surf->tris;
  int *mysurfs = surf->mysurfs;

  for (int i = 0; i < nslocal; i++) {
    buf[n] = pts[tris[mysurfs[i]].p2].x[2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v3x(int n)
{
  Surf::Point *pts = surf->pts;
  Surf::Tri *tris = surf->tris;
  int *mysurfs = surf->mysurfs;

  for (int i = 0; i < nslocal; i++) {
    buf[n] = pts[tris[mysurfs[i]].p3].x[0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v3y(int n)
{
  Surf::Point *pts = surf->pts;
  Surf::Tri *tris = surf->tris;
  int *mysurfs = surf->mysurfs;

  for (int i = 0; i < nslocal; i++) {
    buf[n] = pts[tris[mysurfs[i]].p3].x[1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpSurf::pack_v3z(int n)
{
  Surf::Point *pts = surf->pts;
  Surf::Tri *tris = surf->tris;
  int *mysurfs = surf->mysurfs;

  for (int i = 0; i < nslocal; i++) {
    buf[n] = pts[tris[mysurfs[i]].p3].x[2];
    n += size_one;
  }
}

