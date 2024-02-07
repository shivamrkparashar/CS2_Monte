#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "system.h"
#include <time.h>
#include <string.h>
#include "ran_uniform.h"
#include "ewald.h"

void ReadInputFile(int argc, char *argv[])
{
    /*
     * Easier way of reading input script compared to RASPA's method
     */

    // By default, "input" is the InputFileName
   if (argc == 1) {
        InputFileName = "input";
        printf("FileInputPtr Filename : %s \n", InputFileName);
    }

    if (argc == 2) {
        InputFileName = argv[1];
        printf("FileInputPtr Filename : %s \n", InputFileName);
    }

    if (argc > 2) {
        printf("Number of arguments cannot exceed 2 \n");
        exit(0);

    }

    FILE *FileInputPtr;
    FileInputPtr = fopen(InputFileName, "r");
    if (FileInputPtr == NULL) {
        printf("Could not open file: %s \n", InputFileName);
        return;
    }

    char line[MAXCHAR];
    char AdsorbateFileName[1024];
    while (fgets(line, MAXCHAR, FileInputPtr) != NULL)
    {
        char keyword[1024] = "", arguments[1024] = "";
        // Skip comments
        if (line[0] == '#') continue;

        //sscanf(line,"%s%[^\n]",keyword,arguments);
        sscanf(line,"%s %s \n",keyword,arguments);
        //printf("Keyword: %-30s\t Argument: %-20s\n", keyword, arguments);

        if(strcasecmp("SimulationMethod",keyword)==0) sscanf(arguments, "%s", SimulationMethod);
        if(strcasecmp("PoreGeometry",keyword)==0) sscanf(arguments, "%d", &PoreGeometry);
        if(strcasecmp("PoreDiameter",keyword)==0) sscanf(arguments, "%lf", &Diameter);
        if(strcasecmp("BoxLength",keyword)==0) sscanf(arguments, "%lf", &Box);
        if(strcasecmp("GaugeCellBoxLength",keyword)==0) sscanf(arguments, "%lf", &GaugeCellBoxLength);
        if(strcasecmp("TotalNumberOfParticles",keyword)==0) sscanf(arguments, "%d", &Ntotal);
        if(strcasecmp("Epsilonss",keyword)==0) sscanf(arguments, "%lf", &Epsilonss);
        if(strcasecmp("Sigmass",keyword)==0) sscanf(arguments, "%lf", &Sigmass);
        if(strcasecmp("Rhos",keyword)==0) sscanf(arguments, "%lf", &rhos_Ang);
        if(strcasecmp("Delta",keyword)==0) sscanf(arguments, "%lf", &Delta_Ang);
        if(strcasecmp("NumberOfEquilibrationCycles",keyword)==0) sscanf(arguments, "%d", &NumberOfEquilibrationCycles);
        if(strcasecmp("NumberOfProductionCycles",keyword)==0) sscanf(arguments, "%d", &NumberOfProductionCycles);
        //printf("Number of production cycles  = %d \n", NumberOfProductionCycles);
        if(strcasecmp("SampleEvery",keyword)==0) sscanf(arguments, "%d", &SamplingFrequency);
        if(strcasecmp("VdWCutoff",keyword)==0) sscanf(arguments, "%lf", &VdWCutOff);
        if(strcasecmp("ChargeCutOff",keyword)==0) sscanf(arguments, "%lf", &ChargeCutOff);
        if(strcasecmp("ExternalTemperature",keyword)==0) sscanf(arguments, "%lf", &Temp);
        if(strcasecmp("ExternalChemicalPotential",keyword)==0) sscanf(arguments, "%lf", &Chempot);
        if(strcasecmp("WriteMoviesEvery",keyword)==0) sscanf(arguments, "%d", &MovieFrequency);
        if(strcasecmp("TranslationProbability",keyword)==0) sscanf(arguments, "%lf", &ProbabilityOfDisplacementMove);
        if(strcasecmp("SwapProbability",keyword)==0) sscanf(arguments, "%lf", &ProbabilityOfExchangeMove);
        if(strcasecmp("RotationProbability",keyword)==0) sscanf(arguments, "%lf", &ProbabilityOfRotationMove);
        if(strcasecmp("Restart",keyword)==0) sscanf(arguments, "%d", &Restart);
        if(strcasecmp("AdsorbateDefinition",keyword)==0) sscanf(arguments, "%s", AdsorbateFileName);


    }
    fclose(FileInputPtr);
    MaximumDisplacement = 0.09;

    if (strcasecmp(SimulationMethod,"GCMC")==0)
        printf("Running GCMC simulation \n");
    else if (strcasecmp(SimulationMethod, "IGGC")==0)
        printf("Running Ideal Gas Gauge Cell \n");
    else {
        printf("Error: SimulationMethod can be 'GCMC' or 'IGGC' ");
        exit(0);
    }

    printf("AdsorbateFileName: %s\n", AdsorbateFileName);
    // Reading Adsorbate Definition file

    FILE *AdsorbateInput;
    AdsorbateInput =  fopen(AdsorbateFileName,"r");
    if (AdsorbateInput == NULL)
    {
        printf("Could not open file: %s \n", AdsorbateFileName);
        exit(0);
    }

    if (ProbabilityOfExchangeMove + ProbabilityOfDisplacementMove + ProbabilityOfRotationMove != 1)
    {
        printf("Sum of probabilities of moves is not equal to 1 \n");
        exit(0);
    }

    int aa = 0;
    //  Read adsorbate.def file
    while (fgets(line, MAXCHAR, AdsorbateInput) != NULL)
    {
        char keyword[1024], arguments[1024];
        // Skip comments
        if (line[0] == '#') continue;

        sscanf(line,"%s %[^\n]",keyword,arguments);

        if(strcasecmp("NumberOfAdsorbateAtom",keyword)==0) sscanf(arguments, "%d", &NumberOfAdsorbateAtom);
        else if(strcasecmp("BondDistance",keyword)==0) sscanf(arguments, "%lf", &BondDistance);
        else {
            //char sym[10];
            //double a, b, c;

            //printf("Keyword: %-30s\t Argument: %-20s\n", keyword, arguments);
            sscanf(arguments,"%s %lf %lf %lf %lf %lf", Symbol[aa], &Epsilon[aa], &Sigma[aa], &Mass[aa],&Charge[aa], &AdsorbateAtomicDistance[aa]);
            //sscanf(arguments,"%s %lf %lf %lf", sym, &a, &b, &c);
            aa++;

            //printf("Symbol:%s Epsilon:%lf Sigma:%lf Mass:%lf \n", Symbol[aa], Epsilon[aa], Sigma[aa], Mass[aa]);

            }

    }

    fclose(AdsorbateInput);

    InitializeRandomNumberGenerator(time(0l));
  if(ProbabilityOfExchangeMove + ProbabilityOfDisplacementMove + ProbabilityOfRotationMove != 1.0)
  {
      fprintf(FileOutput,"ProbabilityOfExchangeMove = %lf \n ",ProbabilityOfDisplacementMove);
      fprintf(FileOutput,"Error: Probabilities addition = %lf \n", ProbabilityOfExchangeMove + ProbabilityOfDisplacementMove);
      exit(1);

  }
    // calculate parameters:
    Beta=1.0/Temp;
    //Box = 2 * VdWCutOff + Diameter;

  if(Restart==0)
  {
      // generate configuration from lattice
      Lattice();
  }
  else if (Restart == 1)
  {
      FILE *FileRestartPtr;
      double Boxf, Rhof, Rho;

    //fprintf(FileOutput," Read confs from disk\n");
    printf(" Read confs from disk\n");
    FileRestartPtr=fopen("restart.dat","r");
    fscanf(FileRestartPtr,"%lf\n",&Boxf);
    fscanf(FileRestartPtr,"%d\n",&NumberOfParticles);
    fscanf(FileRestartPtr,"%lf\n",&MaximumDisplacement);
    Rhof=NumberOfParticles/CUBE(Boxf);
    //if(fabs(Boxf-Box)>1.0e-6)
    //fprintf(FileOutput,"%lf %lf\n",Rho,Rhof);
    for(int i=0;i<NumberOfParticles;i++)
    {
        for (int aa=0; aa<NumberOfAdsorbateAtom;aa++)
        {
            fscanf(FileRestartPtr, "%lf %lf %lf\n", &Positions[i][aa].x, &Positions[i][aa].y, &Positions[i][aa].z);
            Positions[i][aa].x *= Box / Boxf;
            Positions[i][aa].y *= Box / Boxf;
            Positions[i][aa].z *= Box / Boxf;
        }
    }
    fclose(FileRestartPtr);
  }

}



void PrintInputData (void)
{
    if (FileOutput == NULL)
    {
        printf("Could not open file: Output.txt \n");
        return;
    }

    // write input data
    fprintf(FileOutput, "Input file data \n ");
    fprintf(FileOutput, "------------------------------------------------------------------\n \n");

    fprintf(FileOutput,"  Simulation Method                          : %s\n",SimulationMethod);
    fprintf(FileOutput, "  Type of Simulation                         : %d\n", PoreGeometry);
    fprintf(FileOutput,"  Number Of Equilibration Cycles             : %d\n",NumberOfEquilibrationCycles);
    fprintf(FileOutput,"  Number Of Production Cycles                : %d\n",NumberOfProductionCycles);
    fprintf(FileOutput,"  Sample Frequency                           : %d\n",SamplingFrequency);
    fprintf(FileOutput,"  Probability To Attempt Displacement Move    : %lf\n", ProbabilityOfDisplacementMove);
    fprintf(FileOutput,"  Probability To Attempt Rotation Move       : %lf\n", ProbabilityOfRotationMove);
    fprintf(FileOutput,"  Probability To Attempt Exchange Move       : %lf\n", ProbabilityOfExchangeMove);
    fprintf(FileOutput,"  Maximum Displacement                       : %lf\n",MaximumDisplacement);
    fprintf(FileOutput,"  Chemical Potential (reduced units)         : %lf\n",Chempot);
    fprintf(FileOutput,"  Temperature (reduced units)                : %lf\n",Temp);
    fprintf(FileOutput,"  Box Length (in Sigma)                      : %lf\n",Box);
    fprintf(FileOutput,"  Gauge Cell Box Length (in Sigma)           : %lf\n",GaugeCellBoxLength);
    fprintf(FileOutput, "  Total Number of Particles (Gauge + Pore)   : %d\n", Ntotal);
    fprintf(FileOutput, "\n===================================================================\n");
    fprintf(FileOutput,"  \n Fluid-Fluid Parameters\n");
    fprintf(FileOutput, "-------------------------------------------------------------------\n");
    fprintf(FileOutput,"  Number of Adsorbate Atoms: %d\n", NumberOfAdsorbateAtom);
    for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
    {
        fprintf(FileOutput, "  Symbol: %s\n", Symbol[aa]);
        fprintf(FileOutput, "    Sigma (in Angstrom)                      : %lf\n", Sigma[aa]);
        fprintf(FileOutput, "    Epsilon (in Kelvin)                      : %lf\n", Epsilon[aa]);
        fprintf(FileOutput, "    Mass (in g/mol)                          : %lf\n", Mass[aa]);
        fprintf(FileOutput, "    Charge (electronic)                      : %lf\n", Charge[aa]);
    }
    fprintf(FileOutput, "    VdWCutOff (in Sigma units)                  : %lf\n", VdWCutOff);
    fprintf(FileOutput, "    ChargeCutOff (in Sigma units)               : %lf\n", ChargeCutOff);
    fprintf(FileOutput, "\n===================================================================\n");
    fprintf(FileOutput,"  \n Solid-Solid Parameters\n");
    fprintf(FileOutput, "-------------------------------------------------------------------\n");
    fprintf(FileOutput,"    Sigma ss (in Angstrom)                      : %lf\n",Sigmass);
    fprintf(FileOutput,"    Epsilon ss (in Kelvin)                      : %lf\n",Epsilonss);
    fprintf(FileOutput,"    Surface density (Ang^3)                     : %lf\n", rhos_Ang);
    fprintf(FileOutput,"    Delta (Ang)                                 : %lf\n", Delta_Ang);
    fprintf(FileOutput,"    Outer Diameter/Width  (reduced units)       : %lf\n",Diameter);

    fprintf(FileOutput, "\n\n ===================================================================\n \n");

}

void SetSimulationUnits(void)
{
    primaryatom = 0;

    // Converts Coulomb potential energy into reduced units
    COULOMB_CONVERSION_FACTOR = SQR(ELECTRONIC_CHARGE) / (4 * M_PI * ELECTRIC_CONSTANT);
    COULOMB_CONVERSION_FACTOR /= (BOLTZMANN_CONSTANT *Epsilon[primaryatom] * Sigma[primaryatom] * ANGSTROM);
    //printf("COULOMB_CONVERSION_FACTOR = %lf \n", COULOMB_CONVERSION_FACTOR);


    // Convert BondDistance into reduced units
    BondDistance /= Sigma[primaryatom];

    EpsilonInSI = Epsilon[primaryatom]*BOLTZMANN_CONSTANT;
    SigmaInSI = Sigma[primaryatom]*ANGSTROM;
    double MoleculeMass = 0;
      for (int aa=0; aa<NumberOfAdsorbateAtom; aa++)
          MoleculeMass += Mass[aa];

    printf("Molecular Mass = %lf AU \n", MoleculeMass);
    MassInSI = MoleculeMass*ATOMIC_MASS_UNIT;
    TempInSI = Epsilon[primaryatom]*Temp;
    Lambda3_Sigma3  = CUBE(sqrt(SQR(PLANCK_CONSTANT) / (2 * M_PI * MassInSI * BOLTZMANN_CONSTANT * TempInSI * SQR(SigmaInSI))));

    IdealGasPressure = Temp*exp(Chempot/Temp)/Lambda3_Sigma3;
    IdealGasPressureInSI = EpsilonInSI/CUBE(SigmaInSI)*IdealGasPressure;

    rhos = rhos_Ang * CUBE(Sigma[primaryatom]); // Volumetric number density
    Delta = Delta_Ang /Sigma[primaryatom];

    // Lorentz Berthelot mixing rule (For fluid and solid)
    for (int aa=0; aa<NumberOfAdsorbateAtom;aa++)
    {
        Epsilonsf[aa] = pow(Epsilon[aa] * Epsilonss, 0.5);
        Sigmasf[aa] = 0.5*(Sigma[aa] + Sigmass);
    }
    
    if(PoreGeometry == 0)
    {
        // GCMC in a box
        printf("Running GCMC in a box \n");
        VolumeOfPore = CUBE(Box);
        VolumeOfPoreInSI = VolumeOfPore*CUBE(SigmaInSI);
        printf("Volume of Simulation box (in A3)    : %lf\n",VolumeOfPoreInSI*pow(10,30));
    }
    else if(PoreGeometry == 1)
    {
        // GCMC in cylindrical pore
        //VolumeOfPore = M_PI*SQR((Diameter-Sigmass/Sigma[primaryatom])/2)*Box;
        printf("Running GCMC in a cylindrical pore \n");
        VolumeOfPore = M_PI*SQR(Diameter/2)*Box;
        VolumeOfPoreInSI = VolumeOfPore*CUBE(SigmaInSI);
        printf("  VolumeOfCylindricalPore (in A3)    : %lf\n",VolumeOfPoreInSI*pow(10,30));
    }
    else if(PoreGeometry == 2)
    {
        // GCMC in slit pore
        //VolumeOfPore = (Diameter-Sigmass/Sigma[primaryatom])*SQR(Box);
        printf("Running GCMC in a Slit pore\n");
        VolumeOfPore = Diameter*SQR(Box);
        VolumeOfPoreInSI = VolumeOfPore*CUBE(SigmaInSI);
        printf("  VolumeOfSlitPore (in A3)    : %lf\n",VolumeOfPoreInSI*pow(10,30));
    }

    else if(PoreGeometry == 3)
    {
        // GCMC in a spherical pore
        printf("Running GCMC in a spherical pore \n");
        VolumeOfPore = 4/3.0 *M_PI *CUBE(Diameter/2.0);
        VolumeOfPoreInSI = VolumeOfPore*CUBE(SigmaInSI);
        printf("  VolumeOfSphericalPore (in A3)    : %lf\n",VolumeOfPoreInSI*pow(10,30));

    }

    if(strcmp(SimulationMethod, "IGGC") ==0)
    {
        GaugeCellVolume = CUBE(GaugeCellBoxLength);
    }




    printf("Temp(K): %lf \n" ,TempInSI);
    //printf("Lambda3_Sigma3: %lf \n" , Lambda3_Sigma3);
    printf("Ideal Gas set pressure (Pa) : %lf \n", IdealGasPressureInSI);
}

void CreateOutputFiles(void)
{
    char buffer[256];


    if (strcasecmp(SimulationMethod,"GCMC")==0) {
        // 1. Create a folder for that Pressure
        sprintf(buffer, "P_%.3f", IdealGasPressureInSI);
        mkdir(buffer, S_IRWXU);
        // 2. Energy.dat file in FileEnergy Variable
        sprintf(buffer, "P_%.3f/Energy.dat", IdealGasPressureInSI);
        FileEnergy = fopen(buffer, "w");
        fprintf(FileEnergy, "Energy_Per_Particle   Pressure_(reduced) \n");
        // 3. Movie.pdb file in FileMovie Variable
        sprintf(buffer, "P_%.3f/Movie.pdb", IdealGasPressureInSI);
        FileMovie = fopen(buffer, "w");
        // 4. Density.dat file in FileDens Variable
        sprintf(buffer, "P_%.3f/Density.dat", IdealGasPressureInSI);
        FileDens = fopen(buffer, "w");
        fprintf(FileDens, "NCycles  Density\n");
        // 5. results.dat file keeps track of all pressure results; located outside any pressure folder
        FileResult = fopen("results.dat", "a");
        // 6. Output.txt file
        sprintf(buffer, "P_%.3f/Output.txt", IdealGasPressureInSI);
        FileOutput = fopen(buffer, "w");

        fprintf(FileOutput, "===================================================================\n \n");
    }
    if(strcasecmp(SimulationMethod, "IGGC")==0){
        // 1. Create a folder for that total number of particles
        sprintf(buffer, "NTotal_%d", Ntotal);
        mkdir(buffer, S_IRWXU);
        // 2. Energy.dat file in FileEnergy Variable
        sprintf(buffer, "NTotal_%d/Energy.dat", Ntotal);
        FileEnergy = fopen(buffer, "w");
        fprintf(FileEnergy, "Energy_Per_Particle   Pressure_(reduced) \n");
        // 3. Movie.pdb file in FileMovie Variable
        sprintf(buffer, "NTotal_%d/Movie.pdb", Ntotal);
        FileMovie = fopen(buffer, "w");
        // 4. Density.dat file in FileDens Variable
        sprintf(buffer, "NTotal_%d/Density.dat", Ntotal);
        FileDens = fopen(buffer, "w");
        fprintf(FileDens, "NCycles  Density\n");
        // 5. results.dat file keeps track of all pressure results; located outside any pressure folder
        FileResult = fopen("results.dat", "a");
        // 6. Output.txt file
        sprintf(buffer, "NTotal_%d/Output.txt", Ntotal);
        FileOutput = fopen(buffer, "w");
        // 7. Histogram Gauge Cell file
        sprintf(buffer, "NTotal_%d/HistogramGaugeCell.dat", Ntotal);
        FileHistogram = fopen(buffer, "w");
        // 8. IGGC result file
        sprintf(buffer, "NTotal_%d/IGGC_results.dat", Ntotal);
        FileIGGCresults = fopen(buffer, "w");

        fprintf(FileOutput, "===================================================================\n \n");
        }
}
