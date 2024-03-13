
#include "system.h"
#include <math.h>
#include "ewald.h"
#include "mc_move.h"
#include "energy.h"


double RunningVdwEnergyff;
double RunningVdwEnergysf;
double RunningFourierEnergyff;
double RunningRealEwaldEnergyff;

// gives En- External potential energy of one particle due to solid-fluid interaction
void EnergyExternal(VECTOR *PosPtr, double *En) {

    double Enij;

    // For particles in a box
    if (PoreGeometry == 0) {
        Enij = 0.0;
    }

        // Cylindrical Pore
    else if (PoreGeometry == 1) {
        Enij = CylindricalPotential(PosPtr);
    }

        // Slit Pore
    else if (PoreGeometry == 2) {
        Enij = SlitPotential(PosPtr);
    } else if (PoreGeometry == 3) {
        Enij = SphericalPotential(PosPtr);
    }

    //passing the value to the input argument En
    *En = Enij;
}


// calculates the energy of particle i with particles j=(jb to NumberOfParticles except itself)
// This is not particle energy but the particle pair energy (That's why no factors of 0.5 in energies)
//  Fluid-fluid energy
void EnergyParticlePair(VECTOR *PosPtr, int i, int jb, double *EVdw, double *EShort, double *Vir) {
    //pos is the VECTOR of particle
    // PosPtr is a pointer to a Vector (to first atom of the molecule)
    // i is used just to avoid calculating energy when j=i

    double r2, Virij;
    double r2i, r6i;
    VECTOR dr;
    int j;
    double sigma_factor, epsilon_factor, sigma_factor_6;
    double UvdWff, URealEwaldff;

    UvdWff = URealEwaldff = 0;
    Virij = 0.0;

    // fluid-fluid Potential
    for (j = jb; j < NumberOfParticles; j++) {
        if (j != i) // Exclude intramolecular interactions
        {
            for (int aa = 0; aa < NumberOfAdsorbateAtom; aa++) {
                for (int aa2 = 0; aa2 < NumberOfAdsorbateAtom; aa2++) {
                    dr.x = PosPtr[aa].x - Positions[j][aa2].x;
                    dr.y = PosPtr[aa].y - Positions[j][aa2].y;
                    dr.z = PosPtr[aa].z - Positions[j][aa2].z;

                    // Why apply boundary conditions on Delta X?
                    // A- If dr > 0.5: dr should be = 1 - dr
                    //    If dr < 0.5: dr should be = dr
                    //    because of interactions with the periodic cells,
                    //    (I have assumed box length to be 1 in above case
                    // apply boundary conditions
                    dr.x -= Box * rint(dr.x / Box);
                    dr.y -= Box * rint(dr.y / Box);
                    dr.z -= Box * rint(dr.z / Box);

                    r2 = SQR(dr.x) + SQR(dr.y) + SQR(dr.z);

                    // calculate the LJ energy
                    if (r2 < SQR(VdWCutOff)) {
                        r2i = 1.0 / r2;
                        r6i = CUBE(r2i);

                        //  Precompute and tabulate these values
                        sigma_factor = 0.5 * (Sigma[aa] + Sigma[aa2]) / Sigma[primaryatom];
                        sigma_factor_6 = pow(sigma_factor, 6);
                        epsilon_factor = pow(Epsilon[aa] * Epsilon[aa2], 0.5) / Epsilon[primaryatom];

                        //  Precompute and tabulate these values
                        double Rc2 = SQR(1 / VdWCutOff);
                        double Rc6 = CUBE(Rc2);
                        PotentialatCutoff = 4.0 * epsilon_factor * (SQR(sigma_factor_6 * Rc6) - sigma_factor_6 * Rc6);

                        // Shifted potential
                        UvdWff += 4.0 * epsilon_factor * (SQR(sigma_factor_6 * r6i) - sigma_factor_6 * r6i) -
                                  PotentialatCutoff;
                        //Enij += 4.0 * epsilon_factor* (SQR(sigma_factor*r6i) - sigma_factor*r6i) - PotentialatCutoff;


                        Virij += 24.0 * (2 * SQR(r6i) - r6i);

                        //printf("VdW: %1.2lf \n", UvdWff);
                    }

                    // Calculate the Short-range Ewald energy (Real space)
                    if (r2 < SQR(ChargeCutOff)) {
                        double q1 = Charge[aa];
                        double q2 = Charge[aa2];
                        double rij = sqrt(r2);

                        URealEwaldff +=
                                COULOMB_CONVERSION_FACTOR * q1 * q2 * ErrorFunctionComplement(Alpha * rij) / rij;
                    }

                }
            }
        }
    }

    // Tail correction
    //Enij += (8.0/3.0)*M_PI*(NumberOfParticles/CUBE(Box))*
    //      ((1.0/3.0)*pow(1/VdWCutOff,9)-pow(1/VdWCutOff,3));

    *EShort = URealEwaldff;
    *EVdw = UvdWff;
    *Vir = Virij;
}

double SlitPotential(VECTOR *PosPtr) {
    /*
     * Takes the pointer to a vector which points to first atom of adsorbate molecules
     * Compute the slit pore potential for all atoms of that molecules
     * Slit pore is along z direction
     */
    double MoleculePotential = 0;
    for (int aa = 0; aa < NumberOfAdsorbateAtom; aa++) // loop over atoms within adsorbate
        MoleculePotential += SteelePotential(PosPtr[aa].z + Diameter / 2 - Box / 2, aa) + \
            SteelePotential(-PosPtr[aa].z + Diameter / 2 + Box / 2, aa);
    return MoleculePotential;

}

double SteelePotential(double z, int aa) {
    /*
     *  Returns one wall potential (Steele 10-4-3) gives the distance z and adsorbate atom index
     *  Wall is along z direction
     *  z is the distance between atom and wall (reduced units by primaryatom)
     *  rhos is the volumetric reduced density
     */

    if (z < 0) return 10e20; // potential on the other side of the wall

    double alpha = 0.61;  // alpha*Delta is the distance between the exposed layer and continuum slab
    double ss = Sigmasf[aa] / Sigma[primaryatom];
    double V = 2 * M_PI * rhos * Delta * SQR(ss) * Epsilonsf[aa] / Epsilon[primaryatom];
    V *= 2. / 5 * pow(ss / z, 10) - pow(ss / z, 4) - pow(ss, 4) / (3. * Delta * pow((z + alpha * Delta), 3));

    return V;

}


double CylindricalPotential(VECTOR *PosPtr) {
    /* Equation taken from- J. Chem. Phys. 135, 084703 (2011) Eq 6
     *
     * Returns Cylindrical potential energy for one adsorbate molecule
     * Cylinder is oriented along z direction
     * rhos is volumetric number density of solid atoms in reduced units
     */
    double MoleculePotential = 0.0;
    double R = Diameter / 2.0;
    double si6_term, si3_term, phi3_term;
    double alpha = 0.61;

    for (int aa = 0; aa < NumberOfAdsorbateAtom; aa++) {
        double r = sqrt(SQR(PosPtr[aa].x - Box / 2) + SQR(PosPtr[aa].y - Box / 2));
        if (r >= R) return 10e20; // potential outside the cylinder

        double ss = Sigmasf[aa] / Sigma[primaryatom];
        double V = 2 * M_PI * rhos * Delta * Epsilonsf[aa] / Epsilon[primaryatom] * SQR(ss);

        si6_term = M_PI * 63.0 / 64.0 * pow(ss / R / (1 - SQR(r / R)), 10) * hypergeometric(-4.5, -4.5, 1, SQR(r / R));
        si3_term = M_PI * 3.0 / 2.0 * pow(ss / R / (1 - SQR(r / R)), 4) * hypergeometric(-1.5, -1.5, 1, SQR(r / R));

        phi3_term = M_PI / 2.0 * pow(ss / (R + alpha * Delta) / (1 - SQR(r / (R + alpha * Delta))), 3) *
                    hypergeometric(-1.5, -0.5, 1, SQR(r / (R + alpha * Delta)));

        V *= si6_term - si3_term - ss / Delta * phi3_term;

        if (r < R) MoleculePotential += V; // Potential inside the cylinder
    }
    return MoleculePotential;
}

double SphericalPotential(VECTOR *PosPtr) {

    // Equation taken from- J. Chem. Phys. 135, 084703 (2011)
    // rhos is the Volumetric density number density of solid atoms in reduced units

    double MolecularPotential = 0.0;
    double alpha = 0.61;
    double R = Diameter / 2.0;


    for (int aa = 0; aa < NumberOfAdsorbateAtom; aa++) {
        double r = sqrt(SQR(PosPtr[aa].x - Box / 2) + SQR(PosPtr[aa].y - Box / 2) + SQR((PosPtr[aa].z - Box / 2)));
        if (r >= R) return 10e20; // potential outside the cylinder

        double ss = Sigmasf[aa] / Sigma[primaryatom];
        double AtomicPotential = 2 * M_PI * rhos * Delta * Epsilonsf[aa] / Epsilon[primaryatom] * SQR(ss);
        double sum1, sum2, sum4, sum5;
        sum1 = sum2 = sum4 = sum5 = 0.0;
        double term3;
        //double x = R-r; // Distance from the wall

        for (int i = 0; i <= 9; i++) {
            sum1 += 1.0 / (pow(R, i) * pow(R - r, 10 - i)) + 1.0 / (pow(R, i) * pow(R + r, 10 - i));
            if (i <= 3) {
                sum2 += 1.0 / (pow(R, i) * pow(R - r, 4 - i)) + 1.0 / (pow(R, i) * pow(R + r, 4 - i));
            }

            if (i == 1 || i == 2) {
                sum4 += 1.0 / (pow(R + alpha * Delta, i) * pow(R - r + alpha * Delta, 3 - i));
                sum5 += 1.0 / (pow(R + alpha * Delta, i) * pow(R + r + alpha * Delta, 3 - i));
            }
        }

        term3 = (CUBE(1.0 / (R - r + alpha * Delta)) + CUBE(1.0 / (R + r + alpha * Delta)));


        AtomicPotential *= 2.0 / 5.0 * pow(ss, 10) * sum1 - pow(ss, 4) * sum2;
        -pow(ss, 4) / 3.0 / Delta * (term3 + 1.5 * (sum4 + sum5));

        if (r < R) MolecularPotential += AtomicPotential; // Potential inside the sphere
    }
    return MolecularPotential;


}


// calculates total system energy
void EnergySystem(void) {
    double EniVdw, Viri;
    double EniShort;
    int i;
    double EnergyFourierSystem = 0;
    double Energyff = 0, Energysf = 0;

    TotalEnergy = 0.0;
    TotalVirial = 0.0;

    // add fluid-fluid energy
    for (i = 0; i < NumberOfParticles - 1; i++) {
        EnergyParticlePair(Positions[i], i, i + 1, &EniVdw, &EniShort, &Viri);
        Energyff += (EniVdw + EniShort);
        TotalEnergy += (EniVdw + EniShort);
        TotalVirial += Viri;
    }

    // add Adsorbate-Adsorbate Fourier energy
    for (i = 0; i < NumberOfParticles; i++) {
        EnergyFourierSystem += EwaldFourierEnergyAdsorbate(Positions[i]);
    }

    Energyff += EnergyFourierSystem;
    TotalEnergy += EnergyFourierSystem;

    fprintf(FileOutput, "\nNumber of Molecules         : %d \n", NumberOfParticles);
    fprintf(FileOutput, "Fourier Fluid-Fluid Energy  : %lf \n", EnergyFourierSystem);
    fprintf(FileOutput, "Total Fluid-Fluid Energy    : %lf \n", Energyff);


    // add solid-fluid energy
    for (i = 0; i < NumberOfParticles; i++) {
        EnergyExternal(Positions[i], &EniVdw);
        Energysf += EniVdw;
        TotalEnergy += EniVdw;

    }


    fprintf(FileOutput, "Total Solid-Fluid Energy    : %lf \n", Energysf);

    // add tail-correction
    //TotalEnergy+= NumberOfParticles*(8.0/3.0)*M_PI*(NumberOfParticles/CUBE(Box))*
    //      ((1.0/3.0)*pow(1/VdWCutOff,9)-pow(1/VdWCutOff,3));
}


double hypergeometric(double a, double b, double c, double x) {
    const double TOLERANCE = 1.0e-10;
    double term = a * b * x / c;
    double value = 1.0 + term;
    int n = 1;

    while (fabs(term) > TOLERANCE) {
        a++, b++, c++, n++;
        term *= a * b * x / c / n;
        value += term;
    }

    return value;
}
