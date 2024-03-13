#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>  // for mkdir()
#include "system.h"
#include "ran_uniform.h"
#include "mc_move.h"
#include "energy.h"


int main(int argc, char *argv[]) {

    clock_t StartTime, EndTime;
    double CPUTimeUsed;
    StartTime = clock();
    int i, j, k;
    int NumberOfCycles, NumberOfSteps;
    int NparticleSumBlock = 0; // within a block
    double PressureSum, PressureCount, Pressure;
    double NParticleSum, NParticleCount, AverageNumberOfParticle;
    double EnergySquaredSum, EnergySum, EnergyCount;
    double rand;
    //VECTOR pos;

    // initialize system
    ReadInputFile(argc, argv);
    SetSimulationUnits();
    // Creates FileOutput, FileDens, FileResult, FileEnergy
    CreateOutputFiles();
    // Prints input parameters to output file
    PrintInputData();
    // Initialize Ewald summation parameters
    InitializeEwald(1e-6);
    // Calculate Self Interaction of one adsorbate molecule into EwaldSelfInteractionEnergy variable
    SelfInteractionEwaldEnergy();

    PressureSum = 0.0;
    PressureCount = 0.0;
    EnergySquaredSum = 0.0;
    EnergySum = 0.0;
    EnergyCount = 0;

    // total energy of the system

    fprintf(FileOutput, "\n\n =================================================================\n");
    fprintf(FileOutput, "\nBefore Starting the Simulation\n");
    fprintf(FileOutput, "-------------------------------------------------------------------\n");
    EnergySystem();
    fprintf(FileOutput, " Total Energy Initial Configuration: %lf\n", TotalEnergy);
    fprintf(FileOutput, " Total Virial Initial Configuration: %lf\n", TotalVirial);
    fprintf(FileOutput, " Total Initial Number of Particles %d\n", NumberOfParticles);

    RunningEnergy = TotalEnergy;
    RunningVirial = TotalVirial;
    RunningFourierEnergyff = 0;
    RunningRealEwaldEnergyff = 0;
    RunningVdwEnergyff = 0;

    // start MC-cycles

    for (i = 1; i <= 2; i++) {
        // i=1 equilibration
        // i=2 production

        if (i == EQUILIBRATION) {
            NumberOfCycles = NumberOfEquilibrationCycles;
            if (NumberOfCycles != 0) {
                fprintf(FileOutput, "\n\n =================================================================\n");
                fprintf(FileOutput, "\nStart Equilibration\n");
                fprintf(FileOutput, "-------------------------------------------------------------------\n");
            }
        } else {
            //NumberOfCycles=NumberOfProductionCycles;
            NumberOfCycles = 10000;
            if (NumberOfCycles != 0) {

                fprintf(FileOutput, "\n\n =================================================================\n");
                fprintf(FileOutput, "\nStart Production\n");
                fprintf(FileOutput, "-------------------------------------------------------------------\n");
            }
        }
        printf("Number of production cycles  = %d \n", NumberOfProductionCycles);
        printf("i = %d \n", i);

        NumberOfDisplacementAttempts = 0;
        NumberOfRotationAttempts = 0;
        NumberOfAcceptedDisplacementMoves = 0;
        NumberOfAcceptedAdditionMoves = 0;
        NumberOfAcceptedDeletionMoves = 0;
        NumberOfAcceptedRotationMoves = 0;

        // initialize the subroutine that adjust the maximum displacement
        Adjust();

        for (j = 1; j < NumberOfCycles + 1; j++) {
            NumberOfSteps = MAX(MinimumInnerSteps, NumberOfParticles);

            for (k = 0; k < NumberOfSteps; k++) {
                rand = RandomNumber();
                // Monte Carlo moves
                if (rand <= ProbabilityOfExchangeMove) {
                    if (strcasecmp(SimulationMethod, "GCMC") == 0)
                        ExchangeMoveAdsorbate();
                    else if (strcasecmp(SimulationMethod, "IGGC") == 0)
                        ExchangeMoveAdsorbateGaugeCell();
                } else if (rand > ProbabilityOfExchangeMove &&
                           rand <= ProbabilityOfExchangeMove + ProbabilityOfDisplacementMove)
                    TranslationMoveAdsorbate();
                else if (rand > ProbabilityOfDisplacementMove + ProbabilityOfExchangeMove)
                    RotationMoveAdsorbate();
            }

            if (i == PRODUCTION) {
                // sample averages

                if ((j % SamplingFrequency) == 0) {
                    //Sample(j, RunningEnergy, RunningVirial, &Pressure, FileEnergy);
                    if (strcasecmp(SimulationMethod, "IGGC") == 0) {
                        int NumberOfParticleGaugeCell = Ntotal - NumberOfParticles;
                        GaugeCellParticleHistogram[NumberOfParticleGaugeCell]++;
                    }

                    fprintf(FileDens, "%d %lf \n", j, NumberOfParticles / VolumeOfPore);
                    PressureSum += Pressure;
                    PressureCount += 1.0;
                    NParticleSum += NumberOfParticles;
                    NparticleSumBlock += NumberOfParticles;
                    NParticleCount += 1.0; // Number of times particles have been counted

                    EnergySum += RunningEnergy;
                    EnergySquaredSum += SQR(RunningEnergy);
                    EnergyCount += 1;
                }

                if ((j % MovieFrequency) == 0) WritePdb(FileMovie);

                if ((j % (NumberOfCycles / 5)) == 0) // Block averages
                {
                    fprintf(FileOutput, "---------->> Done %d out of %d Cycles\n", j, NumberOfCycles);
                    fprintf(FileOutput, "Block Average:\n");
                    fprintf(FileOutput, "    NumberOfParticles = %lf \n", NparticleSumBlock / (NumberOfCycles / 5.));

                    // adjust maximum displacements
                    Adjust();
                    NparticleSumBlock = 0; // Restart with 0 for next block
                }
            }
        }

        // For each of Equilibration and Production:
        if (NumberOfCycles != 0) {
            if (NumberOfDisplacementAttempts != 0) {
                fprintf(FileOutput, "\n\n =================================================================\n");
                fprintf(FileOutput, "\nSummary\n");
                fprintf(FileOutput, "-------------------------------------------------------------------\n");
                fprintf(FileOutput, "Number Of Displacement move attempts  : %d\n", NumberOfDisplacementAttempts);
                fprintf(FileOutput, "Success: %d (%.2lf %%)\n", NumberOfAcceptedDisplacementMoves,
                        100.0 * NumberOfAcceptedDisplacementMoves / NumberOfDisplacementAttempts);
            }

            if (NumberOfRotationAttempts != 0) {
                fprintf(FileOutput, "Number Of Rotation move attempts  : %d\n", NumberOfRotationAttempts);
                fprintf(FileOutput, "Success: %d (%.2lf %%)\n", NumberOfAcceptedRotationMoves,
                        100.0 * NumberOfAcceptedRotationMoves / NumberOfRotationAttempts);
            }
            if (NumberOfExchangeAttempts != 0) {
                fprintf(FileOutput, "Number Of Exchange move attempts  : %d \n", NumberOfExchangeAttempts);
                fprintf(FileOutput, "Addition Success : %d (%.2lf %%) \n", NumberOfAcceptedAdditionMoves,
                        100.0 * NumberOfAcceptedAdditionMoves / NumberOfExchangeAttempts);
                fprintf(FileOutput, "Deletion Success : %d (%.2lf %%) \n", NumberOfAcceptedDeletionMoves,
                        100.0 * NumberOfAcceptedDeletionMoves / NumberOfExchangeAttempts);
            }

            // test total energy
            EnergySystem();

            if (fabs(TotalEnergy - RunningEnergy) > 1.0e-6)
                fprintf(FileOutput, "--------- Problems Energy ---------\n");
            if (fabs(TotalVirial - RunningVirial) > 1.0e-6)
                fprintf(FileOutput, "--------- Problems Virial ---------\n");

            fprintf(FileOutput, "\nTotal Energy End Of Simulation          : %lf\n", TotalEnergy);
            fprintf(FileOutput, "    Running ff VdW Energy               : %lf\n", RunningVdwEnergyff);
            fprintf(FileOutput, "    Running ff Ewald Short-rangeEnergy  : %lf\n", RunningRealEwaldEnergyff);
            fprintf(FileOutput, "    Running ff Fourier                  : %lf\n", RunningFourierEnergyff);
            fprintf(FileOutput, "    Running sf VdW Energy               : %lf\n", RunningVdwEnergysf);
            fprintf(FileOutput, "    Running Energy Total                 : %lf\n", RunningEnergy);
            fprintf(FileOutput, "    Difference                           : %lf\n", TotalEnergy - RunningEnergy);

            // print chemical potential and pressure
            if (i == PRODUCTION) {
                AverageNumberOfParticle = (double) NParticleSum / NParticleCount;
                fprintf(FileOutput, "\nIdeal gas pressure (Pa)       : %lf\n", IdealGasPressureInSI);
                fprintf(FileOutput, "Average Number Of Particles       : %lf\n", AverageNumberOfParticle);
                fprintf(FileOutput, "Loading (mmol/cm3)                : %lf\n",
                        1e3 * AverageNumberOfParticle / (AVOGADRO_CONSTANT * VolumeOfPoreInSI * 1e6));
            }
        }
    }
    Store(); // generates restart file that contains x,y,z coordinates of last frame

    fclose(FileMovie);
    fclose(FileEnergy);

    EndTime = clock();
    CPUTimeUsed = (EndTime - StartTime) / CLOCKS_PER_SEC / 60.; // in minutes

    WriteResult(AverageNumberOfParticle, CPUTimeUsed);
    return 0;
}
