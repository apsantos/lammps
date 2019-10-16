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
#include "pair_lj_long_mdipole_long.h"
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

PairLJLongMDipoleLong::PairLJLongMDipoleLong(LAMMPS *lmp) : 
  PairLJLongDipoleLong(lmp) {}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJLongMDipoleLong::init_style()
{
  if (!atom->mu_flag || !atom->torque_flag)
    error->all(FLERR,"Pair mdipole/long requires atom attributes mu, torque");

  // NOTE: warn if charges are actually set ??

  const char *style3[] = {"ewald/disp", NULL};
  const char *style6[] = {"ewald/disp", NULL};
  int i;

  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  // require an atom style with charge defined

  if (!atom->q_flag && (ewald_order&(1<<1)))
    error->all(FLERR,
        "Invoking coulombic in pair style lj/long/mdipole/long requires atom attribute q");
  if (!atom->mu && (ewald_order&(1<<3)))
    error->all(FLERR,"Pair lj/long/mdipole/long requires atom attributes mu, torque");
  if (!atom->torque && (ewald_order&(1<<3)))
    error->all(FLERR,"Pair lj/long/mdipole/long requires atom attributes mu, torque");

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;

  // ensure use of KSpace long-range solver, set g_ewald

  if (ewald_order&(1<<3)) {                             // r^-1 kspace
    if (force->kspace == NULL)
      error->all(FLERR,"Pair style requires a KSpace style");
    for (i=0; style3[i]&&strcmp(force->kspace_style, style3[i]); ++i);
    if (!style3[i])
      error->all(FLERR,"Pair style requires use of kspace_style ewald/disp");
  }
  if (ewald_order&(1<<6)) {                             // r^-6 kspace
    if (force->kspace == NULL)
      error->all(FLERR,"Pair style requires a KSpace style");
    for (i=0; style6[i]&&strcmp(force->kspace_style, style6[i]); ++i);
    if (!style6[i])
      error->all(FLERR,"Pair style requires use of kspace_style ewald/disp");
  }
  if (force->kspace) g_ewald = force->kspace->g_ewald;

  if (!atom->dipole_magnetic)
    error->all(FLERR,"Using magnetic dipole pair style with electric dipoles");

  qqrd2e = 0.0;
  qmurd2e = 0.0;
  mumurd2e = force->bmubmur2e;
}
