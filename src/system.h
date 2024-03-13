#ifndef SYSTEM_H
#define SYSTEM_H

#include <stdio.h>

#define Npmax 10000
#define MaxAdsorbateAtom 3
#define MinimumInnerSteps 50
#define MaxCharSymbol 7
#define MAXCHAR 1000

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MIN(x, y) ((x)<(y)?(x):(y))
#define MAX(x, y) ((x)>(y)?(x):(y))
#define ANGSTROM 1e-10
#define ATOMIC_MASS_UNIT 1.6605402e-27   //kg
#define PLANCK_CONSTANT 6.6260687652e-34 // J.s
#define BOLTZMANN_CONSTANT 1.380650324e-23    // J K^-1
#define AVOGADRO_CONSTANT 6.0221419947e23     // mol^-1
#define ELECTRONIC_CHARGE 1.60217662e-19 // Coulombs
#define ELECTRIC_CONSTANT  8.8541878176e-12     // C^2/(N.m^2)  // Epsilon naught

enum {
    EQUILIBRATION = 1, PRODUCTION = 2
};

typedef struct {
    double x;
    double y;
    double z;
} VECTOR;

extern int PoreGeometry; // 0: Box, 1: Cylindrical pore, 2: Slit pore, 3: Spherical pore
extern char SimulationMethod[10];   // 'GCMC' or 'IGGC'
extern int Restart;

extern char *InputFileName;

extern FILE *FileOutput; // Main output file
extern FILE *FileMovie, *FileEnergy, *FileDens, *FileResult, *FileHistogram, *FileIGGCresults;

extern double Box;    // Simulation Box Length
extern double Temp;   // Temperature
extern double Beta;   // 1/Temp
extern double Chempot; // Set chemical potential

extern int NumberOfParticles;
extern int GaugeCellNumberOfParticles;
extern int Ntotal; // Ngauge + Npore
extern double GaugeCellBoxLength;
extern double GaugeCellVolume;
extern int GaugeCellParticleHistogram[Npmax];
extern double ChemicalPotentialCanonical[Npmax];
//extern int Occupancy[Npmax]; // is 0 if particle is not present

extern int NumberOfAdsorbateAtom;
extern double Epsilon[MaxAdsorbateAtom];   // Epsilon ff
extern double Sigma[MaxAdsorbateAtom];     // Sigma ff
extern char Symbol[MaxAdsorbateAtom][MaxCharSymbol];
extern int primaryatom; // for sigma, epsilon reference
extern double Mass[MaxAdsorbateAtom];      // Mass Of The atoms
extern double Charge[MaxAdsorbateAtom];    // Charge on the atoms
extern double AdsorbateAtomicDistance[MaxAdsorbateAtom];  //  Atomic distance matrix
extern double Epsilonss; // Epsilon Solid-solid
extern double Sigmass;   // Sigma Solid-solid
extern double Epsilonsf[MaxAdsorbateAtom]; // Epsilon Solid-fluid
extern double Sigmasf[MaxAdsorbateAtom];   // Sigma Solid-fluid
extern double Diameter;  // Diameter/ Width of the pore
extern double VolumeOfPore; // Inner Volume of pore
extern double rhos_Ang;      // Surface number density for cylindrical pore and volume density for slit pore
extern double rhos;       // reduced units
extern double Delta_Ang;
extern double Delta;

extern VECTOR Positions[Npmax][MaxAdsorbateAtom];
extern double BondDistance; // Adsorbate bond distance
extern double VdWCutOff;    // Cut-Off Radius for the VdW Potential
extern double ChargeCutOff;   // Cut-Off radius for Coulombic potential

extern double EpsilonInSI;
extern double SigmaInSI;
extern double MassInSI;
extern double TempInSI;
extern double VolumeOfPoreInSI;
extern double Lambda3_Sigma3;   // de broglie wavelength

extern double TotalEnergy;
extern double TotalVirial;

extern double RunningEnergy;
extern double RunningVirial;

extern int NumberOfDisplacementAttempts;       // number of attemps that have been performed to displace a particle
extern int NumberOfAcceptedDisplacementMoves;  // number of successful attemps to displace a particle
extern int NumberOfRotationAttempts;
extern int NumberOfAcceptedRotationMoves;
extern double MaximumDisplacement; // maximum displacement

extern int NumberOfExchangeAttempts; // number of attempts that have been performed to exchange a particle
extern int NumberOfAcceptedDeletionMoves;
extern int NumberOfAcceptedAdditionMoves;
extern double ProbabilityOfExchangeMove;
extern double ProbabilityOfRotationMove;
extern double ProbabilityOfDisplacementMove;


extern int NumberOfEquilibrationCycles;
extern int NumberOfProductionCycles;
extern int SamplingFrequency;
extern int MovieFrequency;

extern double IdealGasPressure, IdealGasPressureInSI;

extern double PotentialatCutoff;

void Adjust(void);

void ReadInputFile(int argc, char *argv[]);

void PrintInputData(void);


void Sample(int i, double En, double Vir, double *Pressure, FILE *FilePtr);

void WritePdb(FILE *FilePtr);

void Store(void);

void Lattice(void);

void SetSimulationUnits(void);

void CreateOutputFiles(void);

void WriteResult(double AverageNumberOfParticle, double CPUTimeUsed);

#endif