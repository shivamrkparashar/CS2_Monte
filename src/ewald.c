#include <math.h>
#include <stdio.h>
#include "system.h"
#include "ewald.h"
#include "mc_move.h"

// The ratio of time to evaluate real and reciprocal space summation
// http://li.mit.edu/Archive/CourseWork/Ju_Li/Simulation/Ewald/ewald.c
# define TRTF 5.5
#define NINT(x) ((int)((x)>=0.0?((x)+0.5):((x)-0.5)) )
double COULOMB_CONVERSION_FACTOR;
double Alpha;   // Inversely proportional Gaussian width
int kmax;   // Max number of wave vectors in Fourier space
double ReciprocalCutOffSquared;
double EwaldSelfInteractionEnergy;  // Self Interaction Energy of one Adsorbate Molecule

/*
 * Initializes the alpha parameter and number of Wave vectors
 */
void InitializeEwald(double precision) {
    double tol; // initial guess value of Alpha * ChargeCutOff
    double tol1;
    double rc; // Charge Charge cutoff in Angstroms

    // http://li.mit.edu/Archive/CourseWork/Ju_Li/Simulation/Ewald/ewald.c
    // Alpha = sqrt(M_PI)*pow(NumberOfParticles*TRTF, 1./6.)/Box;

    // http://users.clas.ufl.edu/hping/teaching/Fall2012/USRMAN2.20.pdf
    // Equation 3.1
    //Alpha = sqrt(-log(precision*ChargeCutOff))/ChargeCutOff;

    // Taken from RASPA ewald.c
    rc = ChargeCutOff * Sigma[primaryatom];  // rc is cutoff in Angstrom
    tol = sqrt(fabs(log(precision * rc)));
    Alpha = sqrt(fabs(log(precision * rc * tol))) / rc;  // in Angstrom -1

    tol1 = sqrt(fabs(log(precision * rc * SQR(2.0 * tol * Alpha))));
    kmax = NINT(0.25 + Box * Sigma[primaryatom] * Alpha * tol1 / M_PI);
    ReciprocalCutOffSquared = SQR(1.05 * kmax);

    Alpha *= Sigma[primaryatom]; // in reduced units
    //Alpha = 0.20;

    fprintf(FileOutput, "Ewald Summation Parameters     \n");
    fprintf(FileOutput, "-------------------------------------------------------------------\n");
    fprintf(FileOutput, "Relative precession           :%lf\n", precision);
    fprintf(FileOutput, "Alpha (reduced units)         :%lf\n", Alpha);
    fprintf(FileOutput, "kmax                          :%d\n", kmax);


}

double DotProduct(VECTOR a, VECTOR b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

// COMPLEMENTARY ERROR FUNCTION; THIS IS A REALLY GOOD ONE !!!
double ErrorFunctionComplement(double x) {
    double t, u, y;
    const double PA = 3.97886080735226000, P0 = 2.75374741597376782e-1, P1 = 4.90165080585318424e-1;
    const double P2 = 7.74368199119538609e-1, P3 = 1.07925515155856677, P4 = 1.31314653831023098;
    const double P5 = 1.37040217682338167, P6 = 1.18902982909273333, P7 = 8.05276408752910567e-1;
    const double P8 = 3.57524274449531043e-1, P9 = 1.66207924969367356e-2, P10 = -1.19463959964325415e-1;
    const double P11 = -8.38864557023001992e-2, P12 = 2.49367200053503304e-3, P13 = 3.90976845588484035e-2;
    const double P14 = 1.61315329733252248e-2, P15 = -1.33823644533460069e-2, P16 = -1.27223813782122755e-2;
    const double P17 = 3.83335126264887303e-3, P18 = 7.73672528313526668e-3, P19 = -8.70779635317295828e-4;
    const double P20 = -3.96385097360513500e-3, P21 = 1.19314022838340944e-4, P22 = 1.27109764952614092e-3;

    t = PA / (PA + fabs(x));
    u = t - 0.5;
    y = (((((((((P22 * u + P21) * u + P20) * u + P19) * u + P18) * u + P17) * u + P16) * u + P15) * u + P14) * u +
         P13) * u + P12;
    y = ((((((((((((y * u + P11) * u + P10) * u + P9) * u + P8) * u + P7) * u + P6) * u + P5) * u + P4) * u + P3) * u +
           P2) * u + P1) * u + P0) * t * exp(-x * x);
    if (x < 0) y = 2.0 - y;
    return y;
}

double ErrorFunction(double x) {
    return 1.0 - ErrorFunctionComplement(x);
}

double EwaldFourierEnergyAdsorbate(VECTOR *PosPtr) {

    // 1. Calculates Fourier energy of one molecule at position PosPtr

    VECTOR Rk;
    double Rksq, exp_term;
    double energy_sum = 0;
    double EnergyFourier;
    double EnergyEwald;
    double ksq;

    for (int kx = -kmax; kx <= kmax; kx++) {
        Rk.x = 2 * M_PI * kx / Box;
        for (int ky = -kmax; ky <= kmax; ky++) {
            Rk.y = 2 * M_PI * ky / Box;
            for (int kz = -kmax; kz <= kmax; kz++) {
                Rk.z = 2 * M_PI * kz / Box;
                Rksq = SQR(Rk.x) + SQR(Rk.y) + SQR(Rk.z);
                ksq = SQR(kx) + SQR(ky) + SQR(kz);

                //  Ignoring when kx = ky = kz =0 as the exp_term will blow up
                //  And ignoring when Rksq is greater than the cutoff
                if ((ksq == 0) || (ksq > ReciprocalCutOffSquared)) continue;

                float SumCosine = 0;
                float SumSine = 0;

                // sum over all Coulombic atoms in old molecules
                for (int aa = 0; aa < NumberOfAdsorbateAtom; aa++) {
                    SumCosine += Charge[aa] * cos(DotProduct(Rk, PosPtr[aa]));
                    SumSine += Charge[aa] * sin(DotProduct(Rk, PosPtr[aa]));
                }

                exp_term = exp(-0.25 * Rksq / SQR(Alpha)) / Rksq;
                energy_sum += (SQR(SumCosine) + SQR(SumSine)) * exp_term;

            }
        }

    }
    // 1. Fourier energy
    EnergyFourier = COULOMB_CONVERSION_FACTOR * (2 * M_PI / CUBE(Box)) * energy_sum;

    // 2. Self-interaction energy of each gaussian cloud with itself

    // Total Fourier energy after discounting for self-interactions for one adsorbate molecule

    EnergyEwald = EnergyFourier - EwaldSelfInteractionEnergy;

    return EnergyEwald;

}

double FourierEnergyDifferenceAdsorbate(int New, int Old, int i) {
    // Used in all moves
    // New and Old are boolean = 0 or 1
    // i = index of randomly picked molecule for translation, rotation, and deletion.
    // No role of i in addition move

    if (!Charge[0]) {
        return 0.0;
    }
    // 1. Calculates Fourier energy of a molecule at position PosPtr


    VECTOR Rk;  // reciprocal space basis vector
    // Any vector in the reciprocal space is = l*Rk.x + m*Rk.y + n*Rk.z where l, m, n are integers
    double Rksq, exp_term;
    double DeltaUFourier = 0;
    double ksq;

    for (int kx = -kmax; kx <= kmax; kx++) {
        Rk.x = 2 * M_PI * kx / Box;
        for (int ky = -kmax; ky <= kmax; ky++) {
            Rk.y = 2 * M_PI * ky / Box;
            for (int kz = -kmax; kz <= kmax; kz++) {
                Rk.z = 2 * M_PI * kz / Box;
                Rksq = SQR(Rk.x) + SQR(Rk.y) + SQR(Rk.z);
                ksq = SQR(kx) + SQR(ky) + SQR(kz);

                //  Ignoring when kx = ky = kz =0 as the exp_term will blow up
                //  And ignoring when ksq is greater than the cutoff
                if ((ksq == 0) || (ksq > ReciprocalCutOffSquared)) continue;

                float SumCosine_old = 0;
                float SumSine_old = 0;

                // sum over all Coulombic atoms in old molecules
                if (Old) {
                    for (int aa = 0; aa < NumberOfAdsorbateAtom; aa++) {
                        SumCosine_old += Charge[aa] * cos(DotProduct(Rk, Positions[i][aa]));
                        SumSine_old += Charge[aa] * sin(DotProduct(Rk, Positions[i][aa]));
                    }
                }

                float SumCosine_new = 0;
                float SumSine_new = 0;


                // sum over all Coulombic atoms in new molecule
                if (New) {
                    for (int aa = 0; aa < NumberOfAdsorbateAtom; aa++) {
                        SumCosine_new += Charge[aa] * cos(DotProduct(Rk, NewPosition[aa]));
                        SumSine_new += Charge[aa] * sin(DotProduct(Rk, NewPosition[aa]));
                    }
                }

                exp_term = exp(-0.25 * Rksq / SQR(Alpha)) / Rksq;
                DeltaUFourier +=
                        exp_term * (SQR(SumCosine_new) + SQR(SumSine_new) - SQR(SumCosine_old) - SQR(SumSine_old));

            }
        }

    }
    DeltaUFourier *= COULOMB_CONVERSION_FACTOR * (2 * M_PI / CUBE(Box));

    // Taking in the self interaction term into account

    if (Old && New)    // Translation and Rotation move
        DeltaUFourier += 0;
    if (!Old && New)   // Addition move (Self interaction should be subtracted for N+1th molecule)
        DeltaUFourier -= EwaldSelfInteractionEnergy;
    if (Old && !New)   // Deletion move
        DeltaUFourier += EwaldSelfInteractionEnergy;

    return DeltaUFourier;
}

void SelfInteractionEwaldEnergy(void) {
    EwaldSelfInteractionEnergy = 0.0;
    double distance;

    for (int aa = 0; aa < NumberOfAdsorbateAtom; aa++) {
        for (int aa2 = 0; aa2 < NumberOfAdsorbateAtom; aa2++) {
            if (aa == aa2)  // self atomic energy
                EwaldSelfInteractionEnergy += COULOMB_CONVERSION_FACTOR * Alpha / sqrt(M_PI) * SQR(Charge[aa]);
            else   //  Energy of atoms (Gaussian clouds) within the same molecule
            {
                // Interatomic distance, convert into reduced units
                distance = fabs(AdsorbateAtomicDistance[aa] - AdsorbateAtomicDistance[aa2]) / Sigma[primaryatom];
                EwaldSelfInteractionEnergy +=
                        0.5 * COULOMB_CONVERSION_FACTOR * Charge[aa] * Charge[aa2] * ErrorFunction(Alpha * distance) /
                        distance;
            }
        }
    }

    printf("EwaldSelfInteractionEnergy = %lf \n", EwaldSelfInteractionEnergy);
}