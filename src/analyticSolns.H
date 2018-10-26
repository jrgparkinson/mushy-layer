
#ifndef _ANALYTICSOLNS_H_
#define _ANALYTICSOLNS_H_

#include "IntVect.H"
#include "RealVect.H"
#include "MushyLayerParams.h"
#include "LevelData.H"
#include "FArrayBox.H"
#include "BoxIterator.H"
#include "MushyLayerUtils.H"
#include "BoxIterator.H"
//#include "spline.h"
#include "ChomboSpline.h"
#include "Box.H"

#include "NamespaceHeader.H"


/// Produce a simple channel which gives some sense of a mushy layer solution
void channelFieldMushyLayer(LevelData<FArrayBox>& enthalpyAnalytic,
                                      LevelData<FArrayBox>& bulkConcentrationAnalytic,
                                      Real a_domainHeight, Real a_domainWidth, Real a_dx,
                                      MushyLayerParams a_parameters);


void analyticSolnSolidificationNoFlow(LevelData<FArrayBox>& enthalpyAnalytic,
		LevelData<FArrayBox>& bulkConcentrationAnalytic,
		Real a_domainLength, Real a_dx,
		MushyLayerParams a_params);


#endif
