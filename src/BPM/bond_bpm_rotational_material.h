/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(bpm/rotational/material,BondBPMRotationalMaterial);
// clang-format on
#else

#ifndef LMP_BOND_BPM_ROTATIONAL_MATERIAL_H
#define LMP_BOND_BPM_ROTATIONAL_MATERIAL_H

#include "bond_bpm.h"

namespace LAMMPS_NS {

class BondBPMRotationalMaterial : public BondBPM {
 public:
  BondBPMRotationalMaterial(class LAMMPS *);
  ~BondBPMRotationalMaterial() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void settings(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *Kr, *Ks, *Kt, *Kb, *gnorm, *gslide, *groll, *gtwist;
  double *Fcr, *Fcs, *Tct, *Tcb;
  double *Kh, *Fch;
  int smooth_flag;
  int heat_flag;

  double elastic_forces(int, int, int, double, double, double, double *, double *, double *,
                        double *, double *, double *);
  void damping_forces(int, int, int, double *, double *, double *, double *, double *);
  double calculate_heat(int, double, double, double, double);
  double tempbreak(int, double, double);

  void allocate();
  void store_data();
  double store_bond(int, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
