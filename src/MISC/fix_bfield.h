/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(bfield,FixBfield)

#else

#ifndef LMP_FIX_BFIELD_H
#define LMP_FIX_BFIELD_H

#include "fix_efield.h"

namespace LAMMPS_NS {

class FixBfield : public FixEfield {
 public:
  FixBfield(class LAMMPS *, int, char **);
  ~FixBfield() {}
  void init();

 private:
  void init_prefactors();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
