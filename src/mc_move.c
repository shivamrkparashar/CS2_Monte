#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ewald.h"
#include "system.h"
#include "ran_uniform.h"
#include "mc_move.h"
#include "energy.h"

VECTOR NewPosition[MaxAdsorbateAtom];

// attempts to displace a randomly selected particle
void TranslationMoveAdsorbate(void)
{

    if(NumberOfParticles == 0) return;

  double EVdwNewff,VirialNew,EVdwOldff,VirialOld;
  double ERealEwaldNewff, ERealEwaldOldff;
  double EnergyNewsf, EnergyOldsf;
  VECTOR Center;
  double rndmx, rndmy, rndmz;
  int i;
  double DeltaUVdwff, DeltaUVdwsf, DeltaURealEwaldff, DeltaUFourierff, DeltaU;
 
  NumberOfDisplacementAttempts++;

  // choose a random particle
  i=NumberOfParticles*RandomNumber();

  // calculate old energy
    // passing the pointer struct of 1st atom of ith molecule
    EnergyParticlePair(Positions[i], i, 0, &EVdwOldff, &ERealEwaldOldff, &VirialOld);
  EnergyExternal(Positions[i], &EnergyOldsf);

  // give a random displacement
    rndmx = RandomNumber();
    rndmy = RandomNumber();
    rndmz = RandomNumber();
    for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
    {
        NewPosition[aa].x = Positions[i][aa].x + (rndmx - 0.5) * MaximumDisplacement;
        NewPosition[aa].y = Positions[i][aa].y + (rndmy - 0.5) * MaximumDisplacement;
        NewPosition[aa].z = Positions[i][aa].z + (rndmz - 0.5) * MaximumDisplacement;
    }

  // calculate new energy between NewParticle position and all j=0 to NumberOfParticles except j=i
    EnergyParticlePair(NewPosition, i, 0, &EVdwNewff, &ERealEwaldNewff, &VirialNew);
  EnergyExternal(NewPosition, &EnergyNewsf);

    DeltaUFourierff = FourierEnergyDifferenceAdsorbate(1, 1, i);
    DeltaURealEwaldff = ERealEwaldNewff - ERealEwaldOldff ;

    DeltaUVdwff = EVdwNewff  - EVdwOldff;
    DeltaUVdwsf = EnergyNewsf- EnergyOldsf;

  DeltaU = DeltaUVdwff + DeltaURealEwaldff + DeltaUFourierff + DeltaUVdwsf;

  if(RandomNumber()<exp(-Beta*(DeltaU)))
  {
      // accept
      NumberOfAcceptedDisplacementMoves++;
      RunningVdwEnergyff += DeltaUVdwff;
      RunningVdwEnergysf += DeltaUVdwsf;
      RunningRealEwaldEnergyff += DeltaURealEwaldff;
      RunningFourierEnergyff += DeltaUFourierff;
      RunningEnergy+=DeltaU;
      RunningVirial+=(VirialNew-VirialOld);

      Center = CenterOfMass(NewPosition);

      // put particle in simulation box
      if (Center.x < 0.0)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].x += Box;
      else if (Center.x >= Box)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].x -= Box;

      if (Center.y < 0.0)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].y += Box;
      else if (Center.y >= Box)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].y -= Box;

      if (Center.z < 0.0)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].z += Box;
      else if (Center.z >= Box)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].z -= Box;

      // update new position

      for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
      {
          Positions[i][aa].x = NewPosition[aa].x;
          Positions[i][aa].y = NewPosition[aa].y;
          Positions[i][aa].z = NewPosition[aa].z;
      }
  }
}
void ExchangeMoveAdsorbate(void)

{
    // attempts to add/remove a randomly selected particle

    double DeltaUVdwff, VirialChange;
    double DeltaURealEwaldff;
    double DeltaUVdwsf;
    int i;
    double arg;
    double DeltaUFourierff, DeltaU;

    NumberOfExchangeAttempts++;

    if (RandomNumber() < 0.5)
    {
        // Particle Removal

        if (NumberOfParticles == 0) return; // do nothing

        // choose a random particle
        i = NumberOfParticles*RandomNumber(); // i belongs integral values in [0 to N-1]



        // This DeltaUVdwff includes DeltaUVdw
        // DeltaUVdwff = U(N) - U(N-1)
        EnergyParticlePair(Positions[i], i, 0, &DeltaUVdwff, &DeltaURealEwaldff, &VirialChange);
        // Solid fluid energy for particle i
        EnergyExternal(Positions[i],&DeltaUVdwsf);
        // Fourier space difference in energy = U(New) - U(Old) = U(N-1) - U(N)
        DeltaUFourierff = FourierEnergyDifferenceAdsorbate(0, 1, i);

        // DeltaU =  U(New) - U(Old) = U(N-1) - U(N)
        DeltaU = -DeltaUVdwff - DeltaURealEwaldff - DeltaUVdwsf + DeltaUFourierff;
        arg = NumberOfParticles * Lambda3_Sigma3 / VolumeOfPore * exp(-Beta * DeltaU) * exp(-Beta * Chempot);

        if(RandomNumber() < arg)
        {
            //accept and Remove the Particle
            NumberOfAcceptedDeletionMoves++;
            RunningVdwEnergyff -= DeltaUVdwff;
            RunningVdwEnergysf -= DeltaUVdwsf;
            RunningRealEwaldEnergyff -= DeltaURealEwaldff;
            RunningFourierEnergyff += DeltaUFourierff;
            RunningEnergy += DeltaU;
            RunningVirial -= VirialChange;

            NumberOfParticles -= 1;
            // Move the last particle to the removed particle's position
            for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
                Positions[i][aa] = Positions[NumberOfParticles][aa];

        }

    }
    else
    {

        // Particle addition
        //RandomPosition(&NewPosition);
        RandomMoleculePosition(NewPosition);

        EnergyParticlePair(NewPosition, NumberOfParticles, 0, &DeltaUVdwff, &DeltaURealEwaldff, &VirialChange);
        EnergyExternal(NewPosition, &DeltaUVdwsf);
        // Difference in Fourier space energy
        DeltaUFourierff =  FourierEnergyDifferenceAdsorbate(1, 0, 0);

        DeltaU = DeltaUVdwff + DeltaURealEwaldff + DeltaUVdwsf + DeltaUFourierff;
        arg = VolumeOfPore/ Lambda3_Sigma3 / (NumberOfParticles + 1) * exp(-Beta * DeltaU) * exp(Beta * Chempot);

        if(RandomNumber()<arg)
        {
            //accept
            NumberOfAcceptedAdditionMoves++;
            RunningVdwEnergyff += DeltaUVdwff;
            RunningVdwEnergysf += DeltaUVdwsf;
            RunningRealEwaldEnergyff += DeltaURealEwaldff;
            RunningFourierEnergyff += DeltaUFourierff;
            RunningEnergy += DeltaU;
            RunningVirial += VirialChange;

            for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
                Positions[NumberOfParticles][aa] = NewPosition[aa];
            NumberOfParticles += 1;

        }

    }

}

void ExchangeMoveAdsorbateGaugeCell(void)

{
    // attempts to add/remove a randomly selected particle

    double DeltaUVdwff, VirialChange;
    double DeltaURealEwaldff;
    double DeltaUVdwsf;
    int i;
    double arg;
    double DeltaUFourierff, DeltaU;

    NumberOfExchangeAttempts++;

    if (RandomNumber() < 0.5)
    {
        // Particle Removal

        if (NumberOfParticles == 0) return; // do nothing

        // choose a random particle
        i = NumberOfParticles*RandomNumber(); // i belongs integral values in [0 to N-1]



        // This DeltaUVdwff includes DeltaUVdw
        // DeltaUVdwff = U(N) - U(N-1)
        EnergyParticlePair(Positions[i], i, 0, &DeltaUVdwff, &DeltaURealEwaldff, &VirialChange);
        // Solid fluid energy for particle i
        EnergyExternal(Positions[i],&DeltaUVdwsf);
        // Fourier space difference in energy = U(New) - U(Old) = U(N-1) - U(N)
        DeltaUFourierff = FourierEnergyDifferenceAdsorbate(0, 1, i);

        // DeltaU =  U(New) - U(Old) = U(N-1) - U(N)
        DeltaU = -DeltaUVdwff - DeltaURealEwaldff - DeltaUVdwsf + DeltaUFourierff;
        GaugeCellNumberOfParticles = Ntotal - NumberOfParticles;
        arg =  GaugeCellVolume*NumberOfParticles/(GaugeCellNumberOfParticles+1)/VolumeOfPore*exp(-Beta * DeltaU);

        if(RandomNumber() < arg)
        {
            //accept and Remove the Particle
            NumberOfAcceptedDeletionMoves++;
            RunningVdwEnergyff -= DeltaUVdwff;
            RunningVdwEnergysf -= DeltaUVdwsf;
            RunningRealEwaldEnergyff -= DeltaURealEwaldff;
            RunningFourierEnergyff += DeltaUFourierff;
            RunningEnergy += DeltaU;
            RunningVirial -= VirialChange;

            NumberOfParticles -= 1;
            // Move the last particle to the removed particle's position
            for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
                Positions[i][aa] = Positions[NumberOfParticles][aa];

        }

    }
    else
    {

        // Particle addition
        //RandomPosition(&NewPosition);
        RandomMoleculePosition(NewPosition);

        EnergyParticlePair(NewPosition, NumberOfParticles, 0, &DeltaUVdwff, &DeltaURealEwaldff, &VirialChange);
        EnergyExternal(NewPosition, &DeltaUVdwsf);
        // Difference in Fourier space energy
        DeltaUFourierff =  FourierEnergyDifferenceAdsorbate(1, 0, 0);

        DeltaU = DeltaUVdwff + DeltaURealEwaldff + DeltaUVdwsf + DeltaUFourierff;

        GaugeCellNumberOfParticles = Ntotal - NumberOfParticles;
        arg = GaugeCellNumberOfParticles * VolumeOfPore/(NumberOfParticles+1)/ GaugeCellVolume*exp(-Beta * DeltaU);

        if(RandomNumber()<arg)
        {
            //accept
            NumberOfAcceptedAdditionMoves++;
            RunningVdwEnergyff += DeltaUVdwff;
            RunningVdwEnergysf += DeltaUVdwsf;
            RunningRealEwaldEnergyff += DeltaURealEwaldff;
            RunningFourierEnergyff += DeltaUFourierff;
            RunningEnergy += DeltaU;
            RunningVirial += VirialChange;

            for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
                Positions[NumberOfParticles][aa] = NewPosition[aa];
            NumberOfParticles += 1;

        }

    }

}

void RandomPosition(VECTOR *NewPos)
{
    // Returns one random position vector within a pore

    if(PoreGeometry == 0)
    {
        // For Particles in a box
        NewPos->x = RandomNumber()*Box;
        NewPos->y = RandomNumber()*Box;
        NewPos->z = RandomNumber()*Box;

    }

    if(PoreGeometry == 1)
    {
        // For Cylindrical pore
        double angle = 2*M_PI* RandomNumber();
        // This radius is inner radius of cylinder
        // Bonddistance is also subtracted from radius so that O in CO2 do not overlap with solid.
        double radius = (Diameter - Sigmass/Sigma[primaryatom] - 2*BondDistance)/2*sqrt(RandomNumber());
        NewPos->x = Box/2 + radius*cos(angle);
        NewPos->y = Box/2 + radius*sin(angle);
        NewPos->z = RandomNumber()*Box;

    }

    if(PoreGeometry == 2)
    {
        // For Slit pore
        NewPos->x = RandomNumber()*Box;
        NewPos->y = RandomNumber()*Box;
        // Bonddistance is also subtracted from radius so that O in CO2 do not overlap with solid.
        NewPos->z = (Box - (Diameter - BondDistance) )/2 + RandomNumber()* (Diameter - BondDistance);

    }

    if(PoreGeometry == 3)
    {
        // Spherical pore

        double Radius = (Diameter- Sigmass/Sigma[primaryatom])/2.;
        double theta, phi, r, rnd;
        theta = RandomNumber()*2*M_PI;
        rnd = RandomNumber();
        phi = acos(2*rnd - 1);
        r = Radius*pow(RandomNumber(), 1/3.);

        NewPos->x = Box/2 + r*sin(phi)*cos(theta);
        NewPos->y = Box/2 + r*sin(phi)*sin(theta);
        NewPos->z = Box/2 + r*cos(phi);

    }

}
void RandomMoleculePosition(VECTOR *PosPtr)
{
    /*
     * Returns position of all the atoms of adsorbate molecules inserted at a random
     * location and random orientation
     * PosPtr is the Position Ptr of one Adsorbate molecule
     */

    // Generate a random position vector for center of the molecule
    VECTOR Center;
    RandomPosition(&Center);
    RandomMoleculeOrientation(Center, PosPtr);

}

void RandomMoleculeOrientation(VECTOR Center, VECTOR *PosPtr)
{
    // Returns PosPtr array of vector with Center vector

    VECTOR Orientation;
    double theta = RandomNumber()*2*M_PI;  // 0 < theta < 2pi
    double u = 2*RandomNumber() - 1;  // -1 < u < 1

    // Central atom
    PosPtr[1] = Center;

    // Random point on a unit sphere
    Orientation.x = sqrt(1-SQR(u))*cos(theta);
    Orientation.y = sqrt(1-SQR(u))*sin(theta);
    Orientation.z = u;

    //  atom left to the central atom coordinates (Oxygen left)
    PosPtr[0].x = PosPtr[1].x - Orientation.x * BondDistance;
    PosPtr[0].y = PosPtr[1].y - Orientation.y * BondDistance;
    PosPtr[0].z = PosPtr[1].z - Orientation.z * BondDistance;

    //  rightmost atom coordinates (Oxygen right)
    PosPtr[2].x = PosPtr[1].x + Orientation.x * BondDistance;
    PosPtr[2].y = PosPtr[1].y + Orientation.y * BondDistance;
    PosPtr[2].z = PosPtr[1].z + Orientation.z * BondDistance;

}

void RotationMoveAdsorbate(void)
{
    /*
     * The rotation move randomly rotates the adsorbate molecule along a randomly chosen vector.
     *
     */
    if (NumberOfParticles == 0 ) return;
    if (NumberOfAdsorbateAtom == 1) return;  // Rotation not required for monoatomic molecule

    double EVdwNewff, VirialNew, EVdwOldff, VirialOld;
    double EnergyNewsf, EnergyOldsf;
    double ERealEwaldNewff, ERealEwaldOldff;
    VECTOR Center;
    double DeltaUVdwff, DeltaUVdwsf, DeltaURealEwaldff, DeltaUFourierff, DeltaU;
    int i;

    NumberOfRotationAttempts++;
    // Choose a random Particle
    i = NumberOfParticles*RandomNumber();

    // Calculate Old Energy: ff and sf components
    EnergyParticlePair(Positions[i], i, 0, &EVdwOldff, &ERealEwaldOldff, &VirialOld);
    EnergyExternal(Positions[i], & EnergyOldsf);

    // Generate new position by randomly rotating the molecule about its center
    Center = Positions[i][1];
    RandomMoleculeOrientation(Center, &NewPosition[0]);

    // Calculate New Energy
    EnergyParticlePair(NewPosition, i, 0, &EVdwNewff, &ERealEwaldNewff, &VirialNew);
    EnergyExternal(NewPosition, &EnergyNewsf);

    DeltaUFourierff = FourierEnergyDifferenceAdsorbate(1, 1, i);

    DeltaUVdwff = EVdwNewff - EVdwOldff;
    DeltaUVdwsf = EnergyNewsf - EnergyOldsf;
    DeltaURealEwaldff = ERealEwaldNewff - ERealEwaldOldff;
    DeltaU = DeltaUVdwff + DeltaURealEwaldff + DeltaUFourierff + DeltaUVdwsf;

    if (RandomNumber() < exp(-Beta*DeltaU))
    {
        // Accept
        NumberOfAcceptedRotationMoves++;
        RunningVdwEnergyff += DeltaUVdwff;
        RunningVdwEnergysf += DeltaUVdwsf;
        RunningRealEwaldEnergyff += DeltaURealEwaldff;
        RunningFourierEnergyff += DeltaUFourierff;
        RunningEnergy+=DeltaU;
        RunningVirial+=(VirialNew-VirialOld);

        // put particle in simulation box
        if (Center.x < 0.0)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].x += Box;
        else if (Center.x >= Box)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].x -= Box;

        if (Center.y < 0.0)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].y += Box;
        else if (Center.y >= Box)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].y -= Box;

        if (Center.z < 0.0)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].z += Box;
        else if (Center.z >= Box)
          for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
              NewPosition[aa].z -= Box;

        // update new position
        for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
        {
          Positions[i][aa].x = NewPosition[aa].x;
          Positions[i][aa].y = NewPosition[aa].y;
          Positions[i][aa].z = NewPosition[aa].z;
        }
    }
}



VECTOR CenterOfMass(VECTOR *MoleculePos)
{
    VECTOR Center;
    Center.x = Center.y = Center.z = 0.0;

    for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
    {
        Center.x += MoleculePos[aa].x;
        Center.y += MoleculePos[aa].y;
        Center.z += MoleculePos[aa].z;
    }
    Center.x /= NumberOfAdsorbateAtom;
    Center.y /= NumberOfAdsorbateAtom;
    Center.z /= NumberOfAdsorbateAtom;
    return Center;
}