//
// Created by shivam on 12/27/20.
//

#ifndef ENERGY_H
#define ENERGY_H
#include "system.h"

extern double RunningVdwEnergyff;
extern double RunningVdwEnergysf;
extern double RunningFourierEnergyff;
extern double RunningRealEwaldEnergyff;


void EnergyParticlePair(VECTOR *PosPtr, int i, int jb, double *EVdw, double *EShort, double *Vir);
void EnergyExternal(VECTOR *PosPtr, double *En);
void EnergySystem(void);
double hypergeometric(double a, double b, double c, double x);
double SteelePotential(double z, int aa);
double SlitPotential(VECTOR *PosPtr);
double CylindricalPotential(VECTOR *PosPtr);
double SphericalPotential(VECTOR *PosPtr);


#endif // ENERGY_H
