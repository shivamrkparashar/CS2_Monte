#include <stdio.h>
#include <math.h>
#include "system.h"

int PoreGeometry;
char SimulationMethod[10];
int Restart; // 0 if not restarting, 1 if restarting

char *InputFileName;
FILE *FileOutput; // Main output file
FILE *FileMovie, *FileEnergy, *FileDens, *FileResult, *FileHistogram, *FileIGGCresults;

double Box;    // Simulation Box Length
double Temp;   // Temperature
double Beta;   // 1/Temp
double Chempot; // Chemical Potential


int NumberOfParticles;
int NumberOfAdsorbateAtom;
VECTOR Positions[Npmax][MaxAdsorbateAtom];
//MoleculeProperty Adsorbate[MaxAdsorbateAtom];

int GaugeCellNumberOfParticles;
int Ntotal; // Ngauge + Npore
double GaugeCellBoxLength;
double GaugeCellVolume;
int GaugeCellParticleHistogram[Npmax];
double ChemicalPotentialCanonical[Npmax];  // Chempot obtained by Gauge cell method

double TotalEnergy;
double TotalVirial;

double RunningEnergy;
double RunningVirial;

double Epsilon[MaxAdsorbateAtom];   // Epsilon  (in Kelvin)
double Sigma[MaxAdsorbateAtom];     // Sigma    (in Angstrom)
double Mass[MaxAdsorbateAtom];      // Mass Of The Molecules (in Atomic Mass Unit)
double Charge[MaxAdsorbateAtom];
double AdsorbateAtomicDistance[MaxAdsorbateAtom];
char Symbol[MaxAdsorbateAtom][MaxCharSymbol];
int primaryatom;// for sigma, epsilon reference
double Epsilonsf[MaxAdsorbateAtom];
double Sigmasf[MaxAdsorbateAtom];
double Sigmass;
double Epsilonss;
double Diameter; // Outer Diameter
double rhos_Ang;
double rhos;   // reduced units
double Delta_Ang;
double Delta;
double BondDistance; // Adsorbate Bond distance
double VdWCutOff;    // Cut-Off Radius Of The Potential
double ChargeCutOff;


double EpsilonInSI;
double SigmaInSI;
double MassInSI;
double TempInSI;
double VolumeOfPore;
double VolumeOfPoreInSI;
double Lambda3_Sigma3;   // de Broglie Wavelength cubed

int NumberOfDisplacementAttempts;
int NumberOfAcceptedDisplacementMoves;
int NumberOfRotationAttempts;
int NumberOfAcceptedRotationMoves;
double MaximumDisplacement;

int NumberOfExchangeAttempts;
int NumberOfAcceptedDeletionMoves;
int NumberOfAcceptedAdditionMoves;
double ProbabilityOfExchangeMove;
double ProbabilityOfRotationMove;
double ProbabilityOfDisplacementMove;

int NumberOfEquilibrationCycles;
int NumberOfProductionCycles;
int SamplingFrequency;
int MovieFrequency;

double IdealGasPressure, IdealGasPressureInSI;
double PotentialatCutoff;

