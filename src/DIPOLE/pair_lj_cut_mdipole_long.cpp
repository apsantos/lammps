/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <cstring>
#include "pair_lj_cut_mdipole_long.h"
#include "atom.h"
#include "force.h"
#include "kspace.h"
#include "neighbor.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutMDipoleLong::PairLJCutMDipoleLong(LAMMPS *lmp) : 
  PairLJCutDipoleLong(lmp) {}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutMDipoleLong::init_style()
{
  if (!atom->mu_flag || !atom->torque_flag)
    error->all(FLERR,"Pair mdipole/cut requires atom attributes mu, torque");

  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");

  g_ewald = force->kspace->g_ewald;

  cut_coulsq = cut_coul * cut_coul;

  neighbor->request(this,instance_me);

  if (!atom->dipole_magnetic)
    error->all(FLERR,"Using magnetic dipole pair style with electric dipoles");

  // NOTE: warn if charges are actually set ??

  qqrd2e = 0.0;
  qmurd2e = 0.0;
  mumurd2e = force->bmubmur2e;
}
