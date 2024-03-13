#include <math.h>
#include <stdlib.h>
#include "ran_uniform.h"
#include "system.h"


// place NumberOfParticles particles on a simple cubic
// lattice with density 'Rho'
void Lattice(void) {
    int i, j, k, Count, N, NBoxMax;
    double Dx, Dy, Dz, spacing;

    // Initialization
    if (PoreGeometry != 0) // For cylindrical and slit pore
    {
        NumberOfParticles = 0;
        return;
    }

    // Maximum number of particles that can be fitted inside box of given dimension
    NBoxMax = CUBE(Box);
    // Start with empty pore
    /* Initialize NumberOfParticles*/
    //NumberOfParticles = RandomNumber()*MIN(NBoxMax,Npmax);
    NumberOfParticles = 0;
    return;

    // N is the Number of particles in 1 d
    N = (pow(NumberOfParticles, 1.0 / 3.0)) + 1;
    if (N == 0) N = 1;
    spacing = Box / N;
    Count = 0;

    Dx = 0.0;
    for (i = 0; i < N; i++) {
        Dy = 0.0;
        for (j = 0; j < N; j++) {
            Dz = 0.0;
            for (k = 0; k < N; k++) {
                if (Count < NumberOfParticles) {
                    Positions[Count][0].x = Dx;
                    Positions[Count][0].y = Dy;
                    Positions[Count][0].z = Dz;
                    Count++;

                    for (int aa = 1; aa < NumberOfAdsorbateAtom; aa++) {
                        Positions[Count][aa].x = Positions[Count][aa - 1].x + BondDistance;
                        Positions[Count][aa].y = Dy;
                        Positions[Count][aa].z = Dz;
                    }
                }
                Dz += spacing;
            }
            Dy += spacing;
        }
        Dx += spacing;
    }
    //fprintf(FileOutput, "Initialisation on a lattice: %d particles placed\n",Count);
    printf("Initialisation on a lattice: %d particles placed %d\n", Count, NumberOfParticles);
}

