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

#include <cmath>
#include <cstring>
#include <cstdlib>
#include "fix_bfield.h"
#include "atom.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixBfield::FixBfield(LAMMPS *lmp, int narg, char **arg) :
  FixEfield(lmp, narg, arg)
{
  // initialize prefactor constants for magnetic dipoles

  init_prefactors();
}

/* ---------------------------------------------------------------------- */

void FixBfield::init_prefactors()
{
  if (!atom->dipole_magnetic)
    error->all(FLERR,"Using fix bfield with electric dipoles");

  qe2f = 0.0;
  mue2f = force->bmub2f;
}

/* ---------------------------------------------------------------------- */

void FixBfield::init()
{
  FixEfield::init();
  qflag = 0;             // no magnetic monopoles, pending future Nobel prize

  // NOTE: warn if charges are actually set ??
}
