#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// This code reads a plotfile and computes diagnostics from it

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

  string inFolder; // full path to folder containing inputs and checkpoint file(s)
  string outFileName = "cpp_diagnostics.out";
  string prefix;

  pp.get("prefix" ,prefix);
  pp.get("in_folder" ,inFolder);
  pp.query("out_file",outFileName);

  string runInputs  = "inputs"; // default inputs filename

  pp.query("run_inputs",runInputs);
  string inputsFullPath = inFolder;
  inputsFullPath += runInputs;
  addExtraParams(inputsFullPath, pp);

  // Overwrite original params with any of those in our file
  addExtraParams(in_file, pp);

  // If we've specified a particular file, use it.
  // Else, search for files in the directory
  Vector<string> fileNames;

  pout() << "Finding files..." << endl;

  // Get all the files in the directory with the correct prefix
  DIR           *dirp;
  struct dirent *directory;

  dirp = opendir(inFolder.c_str());
  if (dirp)
  {
    while ((directory = readdir(dirp)) != NULL)
    {
      string thisFilename = directory->d_name;
      if (thisFilename.find(prefix.c_str()) != std::string::npos)
      {
        fileNames.push_back(thisFilename);
        //        printf("%s\n", thisFilename.c_str());
      }

    }

    closedir(dirp);
  }

  // Try and sort files

  fileNames.sort();
  //std::sort (fileNames.begin(), fileNames.end(), fileNames);

  // All the things we're going to measure
  int numFiles = fileNames.size();
  Vector<Real> Fs_top(numFiles);
  Vector<Real> Fs_bottom(numFiles);
  Vector<Real> Fs_av(numFiles);
  Vector<Real> time(numFiles);
  Vector<Real> FsVertDiffusionL2(numFiles), FsVertFluidL2(numFiles), FsVertFrameL2(numFiles);
  Vector<Real> FsVertDiffusionL1(numFiles), FsVertFluidL1(numFiles), FsVertFrameL1(numFiles);
  Vector<Real> FsVertDiffusionL0(numFiles), FsVertFluidL0(numFiles), FsVertFrameL0(numFiles);
  Vector<Real> channelWidth(numFiles);
  Vector<Real> channelSpacing(numFiles);

  // Write info to file
  ofstream myfile;
  myfile.open (outFileName.c_str());
  myfile << "Time,Fs_average,Fs_top,Fs_bottom,channel_width,channel_spacing," <<
      "L2FsVertDiffusion,L2FsVertFluid,L2FsVertFrame," <<
      "L1FsVertDiffusion,L1FsVertFluid,L1FsVertFrame," <<
      "L0FsVertDiffusion,L0FsVertFluid,L0FsVertFrame" << endl;

  for (int file_i = 0; file_i < numFiles; file_i++)
  {

    string chkFile = inFolder;
    chkFile += fileNames[file_i];

    // Load this file


    // read physics class header data
    Vector<AMRLevelMushyLayer*> amrlevels;
    int finest_level;
    HDF5HeaderData header;
    getAMRHierarchy(chkFile, amrlevels, finest_level, header);

    time[file_i] = amrlevels[0]->time();

    const ProblemDomain oldLev0Domain = amrlevels[0]->problemDomain();

    // If we haven't loaded the advection velocity, calculate it
    if (!amrlevels[0]->loadedAdvVel())
    {

      // Make sure we're using the eutectic point for nondimensionalisation when we compute velocities (else the calculation fails)
      for (int lev = finest_level; lev >=0; lev--)
      {
        amrlevels[lev]->setDimensionlessReferenceEutectic();
      }

      // Init pressure
      for (int lev = finest_level; lev >=0; lev--)
      {
        amrlevels[lev]->postInitialize();
      }

      // Initialise velocities
      for (int lev = 0; lev<=finest_level; lev++)
      {
      amrlevels[lev]->computeAllVelocities(false);
      }

    }

    // Now shift to the nondimensionalisation relative to the initial point for computing diagnostics
    for (int lev = finest_level; lev >=0; lev--)
    {
      amrlevels[lev]->setDimensionlessReferenceInitial();
    }

    // Compute diagnostics
    for (int lev = finest_level; lev >=0 ; lev--)
    {
      // Make sure we're set up to compute diagnostics
      amrlevels[lev]->set_compute_diagnostics(true);
      amrlevels[lev]->computeDiagnostics();
    }

    Diagnostics& diag = amrlevels[0]->m_diagnostics;

    Fs_top[file_i] = diag.getDiagnostic(DiagnosticNames::diag_soluteFluxTop, time[file_i]);
    Fs_bottom[file_i] = diag.getDiagnostic(DiagnosticNames::diag_soluteFluxBottom, time[file_i]);
    Fs_av[file_i] = diag.getDiagnostic(DiagnosticNames::diag_averageVerticalSaltFlux, time[file_i]);

    FsVertDiffusionL2[file_i] = diag.getDiagnostic(DiagnosticNames::diag_L2FsVertDiffusion, time[file_i]);
    FsVertFluidL2[file_i] = diag.getDiagnostic(DiagnosticNames::diag_L2FsVertFluid, time[file_i]);
    FsVertFrameL2[file_i] = diag.getDiagnostic(DiagnosticNames::diag_L2FsVertFrame, time[file_i]);

    FsVertDiffusionL1[file_i] = diag.getDiagnostic(DiagnosticNames::diag_L1FsVertDiffusion, time[file_i]);
        FsVertFluidL1[file_i] = diag.getDiagnostic(DiagnosticNames::diag_L1FsVertFluid, time[file_i]);
        FsVertFrameL1[file_i] = diag.getDiagnostic(DiagnosticNames::diag_L1FsVertFrame, time[file_i]);

    FsVertDiffusionL0[file_i] = diag.getDiagnostic(DiagnosticNames::diag_L0FsVertDiffusion, time[file_i]);
       FsVertFluidL0[file_i] = diag.getDiagnostic(DiagnosticNames::diag_L0FsVertFluid, time[file_i]);
       FsVertFrameL0[file_i] = diag.getDiagnostic(DiagnosticNames::diag_L0FsVertFrame, time[file_i]);

    //    amrlevels[0]->computeChimneyDiagnostics();

    //    channelWidth[file_i] = diag.getDiagnostic(DiagnosticNames::diag_chimneyWidth, time[file_i]);
    //    channelSpacing[file_i] = diag.getDiagnostic(DiagnosticNames::diag_chimneySpacing, time[file_i]);

    // clean up memory

    myfile << setprecision(10) << time[file_i] << "," << Fs_av[file_i] << "," << Fs_top[file_i] << "," <<
        Fs_bottom[file_i] << "," << channelWidth[file_i] << ","
        << FsVertDiffusionL2[file_i] << ","<< FsVertFluidL2[file_i] << ","<< FsVertFrameL2[file_i] << ","
        << FsVertDiffusionL1[file_i] << ","<< FsVertFluidL1[file_i] << ","<< FsVertFrameL1[file_i] << ","
        << FsVertDiffusionL0[file_i] << ","<< FsVertFluidL0[file_i] << ","<< FsVertFrameL0[file_i] << ","
        << channelSpacing[file_i] << endl;


    // Clean up by deleting mushy layer levels
    for (int lev = finest_level; lev >=0 ; lev--)
        {
          delete amrlevels[lev];
          amrlevels[lev] = NULL;
        }
  }


  myfile.close();


  pout() << "Done..." << endl;
  pout() << endl;

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main

