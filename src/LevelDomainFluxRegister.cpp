#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelDomainFluxRegister.h"
#include "computeSum.H"

LevelDomainFluxRegister::LevelDomainFluxRegister ()
: m_fineFR(nullptr), m_coarseFR(nullptr), m_defined(false),  m_refRat(-1), m_dx(-1), m_numComp(-1)
{ }

LevelDomainFluxRegister::~LevelDomainFluxRegister ()
{
}

void LevelDomainFluxRegister::define(const ProblemDomain a_domain,
                                     const DisjointBoxLayout a_grids,
                                     const int a_refRat,
                                     const Real a_dx,
                                     LevelDomainFluxRegister* a_fineFR,
                                     LevelDomainFluxRegister* a_coarseFR,
                                     int a_numComp)
{
  m_probDomain = a_domain;
  m_grids = a_grids;
  m_defined = true;

  m_refRat = a_refRat;
  m_dx = a_dx;

  m_numComp = a_numComp;

  m_fineFR = a_fineFR;
  m_coarseFR = a_coarseFR;

  m_fluxHi.resize(m_numComp);
  m_fluxLo.resize(m_numComp);

  for (int comp=0; comp<m_numComp; comp++)
  {
    m_fluxHi[comp] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grids, SpaceDim, IntVect::Zero));
    m_fluxLo[comp] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grids, SpaceDim, IntVect::Zero));
  }

  setToZero();

}

void LevelDomainFluxRegister::incrFlux(const LevelData<FluxBox>& a_flux,
                                       const Real a_scale,
                                       const int a_fluxCompStart,
                                       const int a_localCompStart,
                                       const int a_numComps)

{
  for (int comp = 0; comp < a_numComps; comp++)
  {

    for (int fluxDir = 0; fluxDir < SpaceDim; fluxDir++)
    {
      Box probDom = m_probDomain.domainBox();
      Box nodeBox = probDom.surroundingNodes(fluxDir);

      Box bottomFluxBox = ::bdryLo(nodeBox, fluxDir, 1);
      Box topFluxBox = ::bdryHi(nodeBox, fluxDir, 1);

      LevelData<FArrayBox>& localFluxHi =  (*m_fluxHi[a_localCompStart + comp]);
      LevelData<FArrayBox>& localFluxLo =  (*m_fluxLo[a_localCompStart + comp]);

      // Set flux to 0 everywhere except top/bottom edges
      for (DataIterator mlDit = a_flux.dataIterator(); mlDit.ok(); ++mlDit)
      {
        const FluxBox& thisFlux = a_flux[mlDit];
        const FArrayBox& thisFluxDir = thisFlux[fluxDir];

        Box fluxBox = thisFluxDir.box();

        Box localTopBox = fluxBox;
        localTopBox &= topFluxBox;

        Box localBottomBox = fluxBox;
        localBottomBox &= bottomFluxBox;

        for (BoxIterator bit(localTopBox); bit.ok(); ++ bit)
        {
          IntVect iv = bit();
          (localFluxHi)[mlDit](iv - BASISV(fluxDir), fluxDir) += a_scale*thisFluxDir(iv, a_fluxCompStart + comp);
        }

        for (BoxIterator bit(localBottomBox); bit.ok(); ++ bit)
        {
          IntVect iv = bit();
          (localFluxLo)[mlDit](iv, fluxDir) += a_scale*thisFluxDir(iv, a_fluxCompStart + comp);
        }

      } // end loop over boxes

    } // end loop over flux directions

  } // end loop over local components

}



Real LevelDomainFluxRegister::getFluxHierarchy(const int a_dir,
                                               const Side::LoHiSide a_side,
                                               const Real a_scale,
                                               const int a_localComp)
{
  Real flux;

  // get coarsest FR
  LevelDomainFluxRegister* fr = this;
  while (fr->m_coarseFR)
  {
    fr = fr->m_coarseFR;
  }

  Vector<int> a_refRat;
  Vector<LevelData<FArrayBox> * > a_fluxes;

  while (fr)
  {
    a_refRat.push_back(fr->m_refRat);

    if (a_side == Side::Lo)
    {
      a_fluxes.push_back(&(*fr->m_fluxLo[a_localComp]));
    }
    else
    {
      a_fluxes.push_back(&(*fr->m_fluxHi[a_localComp]));
    }

    fr = fr->m_fineFR;
  }
  Real vol;
  flux = a_scale * computeSum(vol, a_fluxes, a_refRat, m_dx, Interval(a_dir, a_dir), 0);

//  flux = flux/vol;

  return flux;

}

//Real LevelDomainFluxRegister::getFluxLevel(const int a_dir,
//                                           const Side::LoHiSide a_side,
//                                           const Real a_scale)
//{
//
//  Real flux;
//
//  if (a_side == Side::Lo)
//  {
//    flux= a_scale * computeSum(m_fluxLo, nullptr, -1, m_dx, Interval(a_dir, a_dir), 0);
//  }
//  else
//  {
//    flux= a_scale * computeSum(m_fluxHi, nullptr, -1,  m_dx, Interval(a_dir, a_dir), 0);
//  }
//
//  return flux;
//
//}

void LevelDomainFluxRegister::setToZero()
{
  for (int comp=0; comp < m_numComp; comp++)
  {
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      (*m_fluxHi[comp])[dit].setVal(0.0);
      (*m_fluxLo[comp])[dit].setVal(0.0);
    }
  }
}
