
#ifndef _ANALYTICSOLNS_H_
#define _ANALYTICSOLNS_H_

#include "BoxIterator.H"
#include "FArrayBox.H"
#include "IntVect.H"
#include "LevelData.H"
#include "MushyLayerParams.h"
#include "MushyLayerUtils.H"
#include "RealVect.H"
#include "Box.H"
#include "ChomboSpline.h"

#include "NamespaceHeader.H"

/// Produce a simple channel which gives some sense of a mushy layer solution
void channelFieldMushyLayer(LevelData<FArrayBox> &enthalpyAnalytic,
                            LevelData<FArrayBox> &bulkConcentrationAnalytic,
                            Real a_domainHeight, Real a_domainWidth, Real a_dx,
                            MushyLayerParams a_parameters);

void analyticSolnSolidificationNoFlow(
    LevelData<FArrayBox> &enthalpyAnalytic,
    LevelData<FArrayBox> &bulkConcentrationAnalytic, Real a_domainLength,
    Real a_dx, MushyLayerParams a_params);

#endif
