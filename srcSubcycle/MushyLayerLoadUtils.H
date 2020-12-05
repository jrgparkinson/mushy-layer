
#ifndef _MUSHYLAYERLOADUTILS_H_
#define _MUSHYLAYERLOADUTILS_H_

#include <cmath>

#include "AMR.H"
#include "AMRLevelMushyLayer.H"
#include "ParmParse.H"
#include "memtrack.H"
#include "parstream.H"

#include "NamespaceHeader.H"

void getAMRHierarchy(string inFile, Vector<AMRLevelMushyLayer *> &amrLevels,
                     int &finest_level, HDF5HeaderData &header);

void addExtraParams(string runInputs, ParmParse &pp);

#include "NamespaceFooter.H"

#endif
