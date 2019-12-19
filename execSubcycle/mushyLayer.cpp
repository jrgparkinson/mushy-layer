
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"


#include "AMR.H"
#include "AMRLevelMushyLayerFactory.H"
#include "DebugDump.H"
#include "memtrack.H"
#include "CH_Attach.H"
#include "FABView.H"
#include "MushyLayerSubcycleUtils.H"
#include "MushyLayerUtils.H"
#include "parstream.H"

#include <string>
#include <fstream>
#include <streambuf>

#include "UsingNamespace.H"

/***************/
void mushyLayer(const Real& a_stopTime,
                const int&  a_nstop,
                const Vector<int>& a_refRat)
{
  CH_TIME("mushyLayer");

//  printRepoVersion();

  // Define the problem domain (size, resolution, and periodicity)
  ProblemDomain prob_domain;
  getProblemDomain(prob_domain);

  // Define the object which
  RefCountedPtr<AMRLevelMushyLayerFactory>  amrg_fact;
  getAMRFactory(amrg_fact);

  AMR amr;
  defineAMR(amr, amrg_fact, prob_domain, a_refRat);

  //
  setupAMRForAMRRun(amr, prob_domain);

  // run the simulation
  amr.run(a_stopTime,a_nstop);

  // Output last pltfile and statistics.
  // Also cleanup everything.
  amr.conclude();

}
/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif
  { //scoping trick



    // Check for an input file
    char* inFile = NULL;

    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
    // Parse the command line and the input file (if any)
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

    ParmParse ppMain("main");


    bool stopTimeOrStep = false;

    Real stopTime = 0.0;
    if (ppMain.contains("max_time"))
    {
      ppMain.get("max_time",stopTime);
      stopTimeOrStep = true;
    }
    else
    {
      stopTime = 100000.0;
      pout() << "No max_time given, using max_time = " << stopTime << endl;

    }

    Real nstop = 0;
    if (ppMain.contains("max_step"))
    {
      ppMain.get("max_step",nstop);
      stopTimeOrStep = true;
    }
    else
    {
      nstop = 9999999;
      pout() << "No max_step given, using max_step = " << nstop << endl;
    }
    int nstop_int = round(nstop);

    if (!stopTimeOrStep)
    {
      pout() << "Neither max_step or max_time given, please define one of these to decide when to stop the simulation." << endl;
      MayDay::Error("Quitting as no max_step or max_time");
    }

    int max_level = 0;
    if (ppMain.contains("max_level"))
    {
      ppMain.get("max_level",max_level);
    }
    else
    {
      pout() << "No max_level give, using max_level = " << max_level << endl;
    }

    int num_read_levels = Max(max_level,1);
    Vector<int> ref_ratios; // (num_read_levels,1);

    // Only require ref_ratio to be defined for AMR simulations
    if (max_level > 0)
    {
      ppMain.getarr("ref_ratio",ref_ratios,0,num_read_levels);
    }
    // AMR expects max_level + 1 values of ref_ratio
    // so add a dummy value
    ref_ratios.push_back(1);



    mushyLayer(stopTime, nstop_int, ref_ratios);
  }
#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  return  0;
}
