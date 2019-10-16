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
#include <stdlib.h>
#include "pair_lj_cut_mdipole_cut.h"
#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutMDipoleCut::PairLJCutMDipoleCut(LAMMPS *lmp) : 
  PairLJCutDipoleCut(lmp) {}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutMDipoleCut::init_style()
{
  if (!atom->mu_flag || !atom->torque_flag)
    error->all(FLERR,"Pair mdipole/cut requires atom attributes mu, torque");

  if (!atom->dipole_magnetic)
    error->all(FLERR,"Using magnetic dipole pair style with electric dipoles");

  // NOTE: warn if charges are actually set ??

  qqrd2e = 0.0;
  qmurd2e = 0.0;
  mumurd2e = force->bmubmur2e;

  neighbor->request(this,instance_me);
}
