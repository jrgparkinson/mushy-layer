//===========================================================================
// driver.cpp
//
//===========================================================================
#include <iostream>

#include "../srcNonSubcycle/amrMushyLayer.H"
#include "ParmParse.H"
#include "AMRIO.H"
#include "SPMD.H"


//===========================================================================
// Simulation of solidification of a binary alloy
//
//===========================================================================
int main(int argc, char* argv[]) {

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif
    
    pout() << "mushyLayer starting" << endl;

    if(argc < 2) 
      { std::cerr << " usage: " << argv[0] << " <input_file>\n"; exit(0); }
    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file);
    
    amrMushyLayer amrObject;

    // set up initial grids, initialize data, etc.
    amrObject.initialize();

    ParmParse pp2("main");
    int maxStep;
    Real maxTime;
    //Real startTime;
    pp2.get("max_time", maxTime);
    pp2.get("max_step", maxStep);
    
    amrObject.run(maxTime, maxStep);

  }  // end nested scope
  
  CH_TIMER_REPORT();


#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
