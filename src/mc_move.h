//
// Created by shivam on 12/25/20.
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ewald.h"
#include "ran_uniform.h"
#ifndef SOURCE_POLYATOMIC_MC_MOVE_H
#define SOURCE_POLYATOMIC_MC_MOVE_H

#endif //SOURCE_POLYATOMIC_MC_MOVE_H

extern VECTOR NewPosition[MaxAdsorbateAtom];


void RandomMoleculeOrientation(VECTOR Center, VECTOR *PosPtr);
void TranslationMoveAdsorbate(void);
void ExchangeMoveAdsorbate(void);
void ExchangeMoveAdsorbateGaugeCell(void);
void RotationMoveAdsorbate(void);
void RandomPosition(VECTOR *NewPos);
void RandomMoleculePosition(VECTOR *PosPtr);
VECTOR CenterOfMass(VECTOR *MoleculePos);