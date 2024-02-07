#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "system.h"


void WritePdb(FILE *FilePtr)
{
  int i;
  static int Countmodel=0,CountMolecule=0;

  Countmodel++;

  fprintf(FilePtr,"MODEL %d\n", Countmodel);
  fprintf(FilePtr,"CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n", Sigma[primaryatom]*Box, Sigma[primaryatom]*Box, Sigma[primaryatom]*Box, 90.0, 90.0, 90.0);

  for(i=0;i<NumberOfParticles;i++)
  {
    CountMolecule++;
    for (int aa=0; aa<NumberOfAdsorbateAtom;aa++)
        fprintf(FilePtr,"ATOM%7d  %2s%12d   %8.3lf%8.3lf%8.3lf\n",
            (i+1),Symbol[aa],aa, Sigma[primaryatom]*Positions[i][aa].x,Sigma[primaryatom]*Positions[i][aa].y,Sigma[primaryatom]*Positions[i][aa].z);
  }
  fprintf(FilePtr,"%s\n","ENDMDL");
}

void Store(void)
{
  int i;
  FILE *FileRestart;
  char buffer[256];

    if (strcasecmp(SimulationMethod, "GCMC") ==0) {
        sprintf(buffer, "P_%.3f/restart.dat", IdealGasPressureInSI);
        FileRestart = fopen(buffer, "w");
    }
    else if (strcasecmp(SimulationMethod, "IGGC")==0){
        sprintf(buffer, "NTotal_%d/restart.dat", Ntotal);
        FileRestart = fopen(buffer, "w");
    }
    fprintf(FileRestart,"%lf\n",Box);
    fprintf(FileRestart,"%d\n",NumberOfParticles);
    fprintf(FileRestart,"%lf\n",MaximumDisplacement);

    for(i=0;i<NumberOfParticles;i++)
    {

        for (int aa=0; aa<NumberOfAdsorbateAtom;aa++)
            fprintf(FileRestart, "%s %lf %lf %lf\n",
               Symbol[aa], Positions[i][aa].x, Positions[i][aa].y, Positions[i][aa].z);
    }
    fclose(FileRestart);
}

void WriteResult(double AverageNumberOfParticle, double CPUTimeUsed)
{
    if (strcasecmp(SimulationMethod,"GCMC")==0) {
        double Loading_g_cm3;
        double MoleculeMass = 0;
        for (int aa = 0; aa < NumberOfAdsorbateAtom; aa++)
            MoleculeMass += Mass[aa];

        Loading_g_cm3 = (AverageNumberOfParticle * MoleculeMass) / (AVOGADRO_CONSTANT * VolumeOfPoreInSI * 1e6);

        fprintf(FileResult, "%20s %20s %20s %20s %20s %20s %20s \n", "Chempot", "P(kPa)", "<N>", "Loading(g/cm3)",
                "Loading(mumol/m2)", "Loading(mol/m3)", "Time(min)");
        fprintf(FileResult, "%20.3lf %20.3lf %20.3lf %20.3lf %20.3lf %20.3lf %20.3lf\n",
                Chempot, IdealGasPressureInSI / 1000, AverageNumberOfParticle, Loading_g_cm3,
                1e6 * AverageNumberOfParticle / (AVOGADRO_CONSTANT * SQR(Box * SigmaInSI)),
                AverageNumberOfParticle / (AVOGADRO_CONSTANT * VolumeOfPoreInSI), CPUTimeUsed);

        printf("%20s %20s %20s %20s %20s %20s %20s \n", "Chempot", "P(kPa)", "<N>", "Loading(g/cm3)",
               "Loading(mumol/m2)", "Loading(mol/m3)", "Time(min)");
        printf("%20.3lf %20.3lf %20.3lf %20.3lf %20.3lf %20.3lf %20.3lf\n",
               Chempot, IdealGasPressureInSI / 1000, AverageNumberOfParticle, Loading_g_cm3,
               1e6 * AverageNumberOfParticle / (AVOGADRO_CONSTANT * SQR(Box * SigmaInSI)),
               AverageNumberOfParticle / (AVOGADRO_CONSTANT * VolumeOfPoreInSI), CPUTimeUsed);
    }

  if (strcasecmp(SimulationMethod, "IGGC") ==0) {
      double ChemicalPotentialAvg;
      // Print Gauge Cell Number of Particles histogram to a file
      fprintf(FileHistogram, "%20s  %20s\n", "NumberOfParticlesGaugeCell", "Frequency");
      for (int i = 0; i < Ntotal + 1; i++) {
          fprintf(FileHistogram, "%d  %d\n", i, GaugeCellParticleHistogram[i]);
      }

      fclose(FileHistogram);

      // Print results on screen
      int Ngauge, Nsim;
      fprintf(FileIGGCresults, "%20s %20s %20s %20s %20s \n", "MaxNparticles", "Npore", "Ngauge", "Chempot", "Time(min)");
      printf("%20s %20s %20s %20s %20s \n", "MaxNparticles", "Npore", "Ngauge", "Chempot", "Time(min)");

      for (Nsim = 0; Nsim < Ntotal; Nsim++) {
          Ngauge = Ntotal - Nsim;
          if (GaugeCellParticleHistogram[Ngauge] !=0 && GaugeCellParticleHistogram[Ngauge-1] !=0 && Ngauge !=0)
          {
              // Number of particles in gauge cell is non-zero

              ChemicalPotentialCanonical[Nsim] = -Temp * log(GaugeCellVolume / Lambda3_Sigma3 / Ngauge) +
                                                 Temp * log((double) GaugeCellParticleHistogram[Ngauge] /
                                                            GaugeCellParticleHistogram[Ngauge - 1]);
              fprintf(FileIGGCresults, "%20d %20d %20d %20.3lf %20.3f\n", Ntotal, Nsim,
                      Ntotal - Nsim, ChemicalPotentialCanonical[Nsim], CPUTimeUsed);
              printf("%20d %20d %20d %20.3lf %20.3f\n", Ntotal, Nsim,
                     Ntotal - Nsim, ChemicalPotentialCanonical[Nsim], CPUTimeUsed);
          }
      }

      // print results.dat
      fprintf(FileResult, "%20s %20s %20s %20s %20s \n", "Ntotal", "Npore_avg", "Ngauge_avg", "Chempot", "Time(min)");
      ChemicalPotentialAvg = -Temp * log(GaugeCellVolume / Lambda3_Sigma3 / (Ntotal-AverageNumberOfParticle + 1) );
      fprintf(FileResult, "%20d %20.3lf %20.3lf %20.3lf %20.3f\n", Ntotal, AverageNumberOfParticle,
              Ntotal-AverageNumberOfParticle, ChemicalPotentialAvg, CPUTimeUsed);
  }

    fclose(FileResult);
}