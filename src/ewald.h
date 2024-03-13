#ifndef EWALD_H
#define EWALD_H


#include <math.h>
#include "system.h"

extern double COULOMB_CONVERSION_FACTOR;
// Ewald Summation parameters

extern double Alpha;   // Gaussian width
extern int kmax;   // Max number of wave vectors in Fourier space
extern double ReciprocalCutOffSquared;  // the squared cutoff in fourier space
extern double EwaldSelfInteractionEnergy;  // Self Interaction Energy of one Adsorbate Molecule

void InitializeEwald(double precision);

double DotProduct(VECTOR a, VECTOR b);

double ErrorFunction(double x);

double ErrorFunctionComplement(double x);

double EwaldFourierEnergyAdsorbate(VECTOR *PosPtr);

double FourierEnergyDifferenceAdsorbate(int New, int Old, int i);

void SelfInteractionEwaldEnergy(void);

#endif
