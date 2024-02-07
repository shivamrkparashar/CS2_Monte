#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

// adjusts maximum displacement such that 50% of the
// moves will be accepted
void Adjust(void)
{
  static int Attempp=0,Naccp=0;
  double Dro,Frac;

  //printf("Adjusting: %d %d\n", NumberOfDisplacementAttempts, NumberOfAcceptedDisplacementMoves);
  if(NumberOfDisplacementAttempts == 0 || Attempp >= NumberOfDisplacementAttempts)
  {
    Naccp=NumberOfAcceptedDisplacementMoves;
    Attempp=NumberOfDisplacementAttempts;
  }
  else
  {
    Frac= (double)(NumberOfAcceptedDisplacementMoves - Naccp) / ((double)(NumberOfDisplacementAttempts - Attempp));
    Dro=MaximumDisplacement;
    MaximumDisplacement*=fabs(Frac/0.5);

    // limit the change
    if(MaximumDisplacement/Dro>1.5) MaximumDisplacement=Dro*1.5;
    if(MaximumDisplacement/Dro<0.5) MaximumDisplacement=Dro*0.5;
    if(MaximumDisplacement>Box/4.0) MaximumDisplacement=Box/4.0;
    //printf(" Max. Displ. Set To : %lf\n",MaximumDisplacement);
    //printf(" (Old : %lf)\n",Dro);
    fprintf(FileOutput," Frac. Acc.: %lf\n",Frac);
    fprintf(FileOutput, " Attempts: %d\n", NumberOfDisplacementAttempts - Attempp);
    fprintf(FileOutput," Success: %d\n\n", NumberOfAcceptedDisplacementMoves - Naccp);

    //printf("Number Of Particles:  %d \n", NumberOfParticles);
  }
}
