#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CELLTOEDGE2_H_
#define _CELLTOEDGE2_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

enum CellToEdgeAveragingMethod
{
    arithmeticAveraging,
    geometricAveraging
};

///
void CellToEdge2(const LevelData<FArrayBox> &a_cellData,
                 LevelData<FluxBox> &a_edgeData,
                 CellToEdgeAveragingMethod method = arithmeticAveraging);

///
void CellToEdge2(const FArrayBox &a_cellData, FluxBox &a_edgeData,
                 CellToEdgeAveragingMethod method = arithmeticAveraging);

/// this one does a single direction
void CellToEdge2(const FArrayBox &a_cellData, FArrayBox &a_edgeData,
                 const int a_dir,
                 CellToEdgeAveragingMethod method = arithmeticAveraging);

/// this one does a single direction and component
void CellToEdge2(const FArrayBox &a_cellData, const int a_cellComp,
                 FArrayBox &a_edgeData, const int a_edgeComp, const int a_dir,
                 CellToEdgeAveragingMethod method = arithmeticAveraging);

#include "NamespaceFooter.H"

#endif
