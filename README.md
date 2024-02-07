# CS2_Monte
## A Monte Carlo Code for Gas Adsorption Simulations in Cylindrical, Slit and Spherical (CS2) pore geometry.

CS2_Monte allows monte carlo simulations of gas adsorption in two ensembles: grand canonical and mesocanonical (or gauge cell).

This code can model adsorbate of upto 3 nodes. For examples, single site CH4, Ar, TraPPE N2 and CO2 can be modelled. 


## Input File Format
```
SimulationMethod             IGGC       # IGGC or GCMC
SimulationType               2          # 0: Box, 1: Cylindrical, 2: Slit, 3: Spherical
BoxLength                    14.41      # Box length/sigma_ff
GaugeCellBoxLength           150        # Gauge cell length/sigma_ff
TotalNumberOfParticles       10         # Ntotal (Ngauge + Npore)
NumberOfEquilibrationCycles  5000       # Number of equilibration cycle. 1 Cycle = max(50, Number of particles in pore cell)
NumberOfProductionCycles     10000      # Number of production cycle
SampleEvery                  1          # Frequency of writing output 
WriteMoviesEvery             5000       # Frequency of printing coordinates of adsorbate
Restart                      0          # 0- start from scrath, 1- start from a restart file. The restart file should be named 'restart.dat'.
VdWCutoff                    5.00       # Cutoff distance for van der walls interactions in sigma_ff 
ChargeCutOff                 5.00       # Cutoff distance for real space ewald interactsion in sigma_ff
ExternalTemperature          0.73       # External temperature/epsilon_ff
AdsorbateDefinition          Ar.def     # Definition of adsorbate file
TranslationProbability       0.5        # Target proability of translation move
SwapProbability              0.5        # Target probability of exchange move
RotationProbability          0          # Target probability of rotation move
Sigmass                      3.405      # LJ parameter sigma for solid-solid in Angstroms
Epsilonss                    28.047     # LJ parameter epsilon for solid-solid in Kelvin
PoreDiameter                 4.41       # Pore diameter in sigma_ff
Rhos                         0.114      # Surface number density for cylindrical pore and volume density for slit pore (A^2 or A^3)
Delta                        3.354      # Distance between layers of solid atoms in stelle potential (Ang)
```


## How to build the c code
Make sure that file `CMakeLists.txt` is inside the `src` directory
```bash
cd src
mkdir build
cd build
cmake ../
cmake --build .
```
This will generate an executable `cs2_monte` in the `build` directory

