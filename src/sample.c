#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

// Write Quantities (Pressure And Energy) To File
// Ener (Input) : Total Energy
// Vir  (Input) : Total Virial

void Sample(int i,double En,double Vir,double *Pressure,FILE *FilePtr)
{
  double Enp,Vol;
 
  if(NumberOfParticles>0)
  {
    Enp=En/NumberOfParticles;
    Vol=CUBE(Box);
    *Pressure=(NumberOfParticles/Vol)/Beta+Vir/(3.0*Vol);
    *Pressure+=(16.0/3.0)*M_PI*SQR(NumberOfParticles/Vol)*
        ((2.0/3.0)*pow(1 / VdWCutOff, 9) - pow(1 / VdWCutOff, 3));
  }
  else
  {
    Enp=0.0;
    *Pressure=0.0;
  }
  fprintf(FilePtr,"%lf %lf\n",Enp,*Pressure);
}
