#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Read in a chk file and write another, with some slight modifications

#include <iostream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "BoxLayout.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"
#include <iostream>
#include <dirent.h>
#include <stdio.h>
#include "AMRLevelMushyLayer.H"
#include "Diagnostics.h"
#include "MushyLayerLoadUtils.H"

// One more function for MPI
void dumpmemoryatexit();

void getBoxes(Vector<Box>& outBoxes, Box newDomain,
              Real block_factor, Real box_size)
{
  //Copied from AMR

  if (box_size < block_factor)
  {
    MayDay::Abort("Base grid size must be greater than blocking factor");
  }



  Tuple<Vector<int>,SpaceDim> box_sizes;
  IntVect num_grids;
  IntVect base_size;

  for (int d = 0; d < SpaceDim; ++d)
  {
    int num_div = 1;
    int domain_size = newDomain.size(d);
    while (num_div * box_size < domain_size)
    {
      ++num_div;
    }

    // int(x/y) +(x%y)?1:0 is integer division with rounding upwards, for x,y>0
    base_size[d] = int(domain_size/num_div) +((domain_size%num_div) ? 1 : 0);
    box_sizes[d].resize(num_div, base_size[d]);
    box_sizes[d][num_div-1] = domain_size -(num_div - 1) * base_size[d];
    num_grids[d] = num_div;
  }

  Box b(IntVect::Zero,num_grids - IntVect::Unit);
  const IntVect& domain_hi = newDomain.bigEnd();
  const IntVect& domain_lo = newDomain.smallEnd();

  BoxIterator bit(b);
  for (bit.begin(); bit.ok(); ++bit)
  {
    const IntVect& iv = bit();
    IntVect lo = domain_lo + iv * base_size;
    IntVect hi = min(lo + base_size - IntVect::Unit, domain_hi);
    Box grid(lo,hi);

    pout() << "New box (" << lo[0] << "," << lo[1] << ") - (" << hi[0] << "," << hi[1] << ")" << endl;


    // I have no idea why we were previously refining these boxes
//    Box refinedBox = refine(grid, block_factor);
//    lo = refinedBox.smallEnd();
//    hi = refinedBox.bigEnd();
//    pout() << "Refined box (" << lo[0] << "," << lo[1] << ") - (" << hi[0] << "," << hi[1] << ")" << endl;

//    outBoxes.push_back(refinedBox );

    outBoxes.push_back(grid );
  }

}

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  // setChomboMPIErrorHandler();
  MPI_Barrier(Chombo_MPI::comm);  // Barrier #1
#endif

  // ------------------------------------------------
  // parse command line and input file
  // ------------------------------------------------
  // infile must be first
  if (argc < 2)
  {
    cerr << "  need inputs file" << endl;
    abort();
  }

  char* in_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, in_file);



#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  pout() << endl;
  pout() << "Initializing..." << endl;

  // declare variable to store hierarchy
  string    inFile;
  string    outFile;
  string    runInputs ;
  int changeLeft = 0;
  int changeRight = 0;
  int changeTop = 0;
  int changeBottom = 0;
  int xShift = 0;
  Real box_size = -1;
  int refinement = 1;
  Real dtReductionFactor = 10;
  Real smoothing = 0.0;
  //int Ny;
  //  Real H;
  pp.get("run_inputs",runInputs);
  pp.get("inFile" ,inFile);
  pp.get("outFile" ,outFile);
  pp.query("deltaNxLeft" ,changeLeft);
  pp.query("deltaNxRight" ,changeRight);
  pp.query("deltaNzBottom" ,changeBottom);
  pp.query("deltaNzTop" ,changeTop);
  pp.query("smoothing", smoothing);

  pp.query("xShift" ,xShift);
  pp.query("refinement", refinement);
  pp.query("box_size" ,box_size);
  pp.query("dt_reduction_factor" ,dtReductionFactor);
  //  pp.get("Ny" ,Ny);
  //  pp.get("H" ,H);

  Real newLev0Dx = -1;
  pp.query("newDx", newLev0Dx);

  addExtraParams(runInputs, pp);


  // Get the AMR hierarchy
  Vector<AMRLevelMushyLayer*> amrlevels;
  int finest_level;
  HDF5HeaderData header;
  getAMRHierarchy(inFile, amrlevels, finest_level, header);


  // Get parameters
  int block_factor;
  ParmParse ppMain("main");
  ppMain.get("block_factor", block_factor);


  // Make changes

  // E.g. change domain size
  const ProblemDomain oldLev0Domain = amrlevels[0]->problemDomain();

  int dir = 0;

  // Define new boxes if we want to change the domain width or
  // enforce a certain box_size
  if ((  changeLeft != 0 || changeRight != 0
      || changeTop != 0 || changeBottom != 0)
      || (box_size > 0 && refinement == 1)     )
  {

    if ((changeLeft != 0 || changeRight != 0 || changeTop != 0 || changeBottom != 0) != 0
        && refinement != 1)
    {
      pout() << "Warning - not doing refinement" << endl;
    }

    // Initially level 0 domain, then get's refined
    Box newDomainBox(oldLev0Domain.domainBox());



    ProblemDomain newDomain(oldLev0Domain);




    // Grow/shrink domain
    //    newDomain.grow(dir, deltaNx);
    newDomain.growLo(0, changeLeft);
    newDomain.growHi(0, changeRight);

    newDomain.growLo(1, changeBottom);
    newDomain.growHi(1, changeTop);




    for (int level = 0; level <= finest_level; level++)
    {
      AMRLevelMushyLayer* ml = amrlevels[level];

      Box oldDomain = amrlevels[level]->problemDomain().domainBox();

      const DisjointBoxLayout curDBL = ml->grids();

      pout() << "Processing level " << level << "..." << endl;

      Vector<Box> outBoxes;
      Vector<int> outProcs;

      pout() << "  Creating data holders..." << endl;



      // Define new grids with specified box_size
      if (box_size > 0)
      {
        getBoxes(outBoxes, newDomain.domainBox(), block_factor, box_size);

        //        outBoxes = DisjointBoxLayout(outBoxes, outProcs);
        //makeGrids(outBoxes, newDomain);
      }
      else
      {
        // adjust current grids to fit new domain;

        MayDay::Error("Need to know box size!");

      }

      LoadBalance(outProcs,outBoxes);

      DisjointBoxLayout outGrids(outBoxes,outProcs,newDomain);

      // Need to copy all data on this level to the new grids
      // scalar vars, vector vars and m_advVel are written out

      // could  maybe just used regrid here?
      pout() << "Reshaping data for new grids" << endl;
      ml->reshapeData(outGrids, newDomain);
      pout() << "Finished reshaping data for new grids" << endl;

      // Ready for next level
      newDomain.refine(ml->refRatio());
    }

  } // end if changing domain
  else if (refinement != 1)
  {
    pout() << "Refining data with refinement factor: " << refinement << endl;
    // Refine current data/domain
    for (int level = 0; level <= finest_level; level++)
    {
      Vector<Box> outBoxes;
      Vector<int> outProcs;

      AMRLevelMushyLayer* ml = amrlevels[level];

      ProblemDomain domain = ml->problemDomain();

      domain.refine(refinement);

      pout() << "Generating boxes on level: " << level << endl;
      getBoxes(outBoxes, domain.domainBox(), block_factor, box_size);

      LoadBalance(outProcs, outBoxes);

      DisjointBoxLayout outGrids(outBoxes, outProcs, domain);

      // Need to interpolate data onto newgrids
      ml->refine(refinement, outGrids, domain);
      //      ml->regrid(outBoxes);
      //      ml->reshapeData(outGrids, domain);

      // Update dx
      ml->dx(ml->dx()/refinement);

    }
  }

  // Shift data if asked
  if (xShift != 0)
  {   
    for (int level = 0; level <= finest_level; level++)
    {
      pout() << "Shifting data on level: " << level << endl;
      AMRLevelMushyLayer* ml = amrlevels[level];

      ml->shiftData(dir, xShift);
    }
  }

  // Do smoothing if required
  if (smoothing > 0)
  {
    bool smoothVel = true;
    bool smoothScalar = true;
    pp.query("smoothVel", smoothVel);
    pp.query("smoothScalar", smoothScalar);


    pout() << "Doing smoothing with coefficient " << smoothing << endl;
    // Only have to call this from the coarsest level
      AMRLevelMushyLayer* mlCoarsest = amrlevels[0];
      mlCoarsest->setSmoothingCoeff(smoothing);
      mlCoarsest->doPostRegridSmoothing(smoothVel, smoothScalar);

      // Now ensure U=0 in the right places
      for (int level = 0; level <= finest_level; level++)
      {
        Real porosity_cap = 0.01;
        ppMain.query("ccvel_porosity_cap", porosity_cap);

        AMRLevelMushyLayer* ml = amrlevels[level];
        ml->updateEnthalpyVariables();
        ml->setCCVelZero(porosity_cap);
      }

  }

  // Write new file

  HDF5Handle handleOut;
  {
    CH_TIME("AMR::writeCheckpointFile.openFile");
    handleOut.open(outFile.c_str(), HDF5Handle::CREATE);
  }

  // Reset iteration and time

  // write amr data
  //  HDF5HeaderData header;
  //  header.m_int ["max_level"]  = max_level;
  //  header.m_int ["num_levels"] = finest_level + 1;
  header.m_int ["iteration"]  = 0;
  header.m_real["time"]       = 0.0;

  // Set periodicity info
  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  //bool isPeriodic[SpaceDim];
  Vector<int> isPeriodic(SpaceDim, -1);

  pp.queryarr("periodic", isPeriodic, 0, SpaceDim);

  D_TERM(
       if (isPeriodic[0] >= 0)
         header.m_int ["is_periodic_0"] = isPeriodic[0];
        ,

         if (isPeriodic[1] >= 0)
           header.m_int ["is_periodic_1"] = isPeriodic[1];
         ,

           if (isPeriodic[2] >= 0)
             header.m_int ["is_periodic_2"] = isPeriodic[2];
           );


  //
  //  for (int level = 0; level < regrid_intervals.size(); ++level)
  //  {
  //    char headername[100];
  //    sprintf(headername, "regrid_interval_%d", level);
  //    header.m_int[headername] = regrid_intervals[level];
  //  }
  //
  //  for (int level = 0; level < max_level; ++level)
  //  {
  //    char headername[100];
  //    sprintf(headername, "steps_since_regrid_%d", level);
  //    header.m_int[headername] = steps_since_regrid[level];
  //  }

  // should steps since regrid be in the checkpoint file?
  {
    CH_TIME("writeHeader");
    header.writeToFile(handleOut);
  }

  // write physics class data
  Real lev0dt = amrlevels[0]->dt();
  amrlevels[0]->dt(lev0dt/dtReductionFactor); // need this

  // might have changed dx
  if (newLev0Dx > 0)
  {
    amrlevels[0]->dx(newLev0Dx);
    Real prevDx = newLev0Dx;
    for (int level = 1; level <= finest_level; ++level)
    {
      Real levDx = prevDx/amrlevels[level]->refRatio();
      amrlevels[level]->dx(levDx);

      prevDx = levDx;
    }
  }

  pout() << " Write checkpoint header" << endl;
  amrlevels[0]->writeCheckpointHeader(handleOut);


  for (int level = 0; level <= finest_level; ++level)
  {
    pout() << " Write checkpoint data for level: " << level << endl;
    amrlevels[level]->writeCheckpointLevel(handleOut);
  }

  handleOut.close();



#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main



