#include "AMRLevelMushyLayer.H"
#include "CH_HDF5.H"

void AMRLevelMushyLayer::writeAMRHierarchy(string filename)
{
  CH_assert(m_level==0);

  // For now, we're just writing out one vector variable
//  int numVars = 1;
//  int numComps = SpaceDim*m_numVectorVars + m_numScalarVars;
  int numComps = SpaceDim*m_outputVectorVars.size() + m_outputScalarVars.size();

  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  int nlevels = hierarchy.size();

  Vector<AMRLevelMushyLayer*> a_hierarchy; // don't actually need this, except to fill the arguments in getHierarchyAndGrids(...)
  Vector<DisjointBoxLayout> a_vectGrids;
  Vector<LevelData<FArrayBox>* > a_vectData;
  Vector<string> a_vectNames(numComps);
  ProblemDomain a_baseDomain;
  Real a_baseDx;
  Real a_dt;
  Real a_time;
  Vector<int> a_vectRatio;
  Interval a_levels;

  // Get data
  a_vectData.resize(nlevels);

  AMRLevelMushyLayer* ml = this;
  for (int lev=0; lev < nlevels; lev++)
  {
    ml->getExtraPlotFields();

    a_vectData[lev] = new LevelData<FArrayBox>(ml->m_grids, numComps);

    int startComp = 0;

    for (int vari = 0; vari < m_outputScalarVars.size(); vari++)
    {
      int thisVar = m_outputScalarVars[vari];
      (*ml->m_scalarNew[thisVar]).copyTo(Interval(0, 0), *a_vectData[lev], Interval(startComp, startComp));
      startComp = startComp+1;
    }

    for (int vari = 0; vari < m_outputVectorVars.size(); vari++)
    {
      int thisVar = m_outputVectorVars[vari];
      (*ml->m_vectorNew[thisVar]).copyTo(Interval(0, SpaceDim-1), *a_vectData[lev], Interval(startComp, startComp+SpaceDim-1));
      startComp = startComp+SpaceDim;
    }

    ml = ml->getFinerLevel();
  }

  for (int vari = 0; vari < m_outputScalarVars.size(); vari++)
  {
    a_vectNames[vari] = m_scalarVarNames[m_outputScalarVars[vari]];
  }

  for (int vari = 0; vari < m_outputVectorVars.size(); vari++)
  {
    for (int i=0; i < SpaceDim; i++)
    {
      char vectName[100];
      string component;
      if (i==0)
      {
        component = "x";
      }
      else if(i==1)
      {
        component = "y";
      }
      else
      {
        component = "z";
      }

      sprintf(vectName, "%s%s", component.c_str(), m_vectorVarNames[m_outputVectorVars[vari]].c_str());
      a_vectNames[m_outputScalarVars.size() + vari*SpaceDim + i] = string(vectName);
    }
  }

  getHierarchyAndGrids(a_hierarchy,
                             a_vectGrids, a_vectRatio,
                             a_baseDomain, a_baseDx);



  a_dt = m_dt;
  a_time = m_time;
  a_levels = Interval(0, nlevels-1);

  string output_dir = "";
  ParmParse ppMain("main");
  ppMain.query("output_folder", output_dir);
//  string fullFilename = sprintf("%s/%s", output_dir, filename);

  char fullFilename[300];
  sprintf(fullFilename, "%s/%s", output_dir.c_str(), filename.c_str());



  WritePartialAMRHierarchyHDF5(/// file to send output to
                               fullFilename,
                               /// grids at each level
                               a_vectGrids,
                               /// data indexed by level
                                a_vectData,
                               /// names of variables in <i>a_vectData</i>
                               a_vectNames,
                               /// domain at base level given by <i>a_levels</i>.begin()
                               a_baseDomain.domainBox(),
                               /// grid spacing at base level
                               a_baseDx,
                               /// time step at base level
                               a_dt,
                               /// time
                               a_time,
                               /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                               a_vectRatio,
                               /// indices of levels to write out
                               a_levels);
}


#ifdef CH_USE_HDF5

/*******/

void
AMRLevelMushyLayer::
writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::writeCheckpointHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;

  int numChkScalarComps = m_chkScalarVars.size();
  int numChkVectorComps = m_chkVectorVars.size()*SpaceDim;
  int numChkComps = numChkScalarComps + numChkVectorComps;

  header.m_int["num_components"] = numChkComps;

  // Setup the component names
  char compStr[30];

  // First, get all the variable names (even the ones we're not writing out)
  Vector<string> scalarVarNames(m_numScalarVars);

  getScalarVarNames(scalarVarNames, m_numScalarVars, m_scalarVarNames);

  Vector<string> chkVectorVarNames(m_numVectorVars);
  getChkVectorVarNames(chkVectorVarNames, m_chkVectorVars, m_vectorVarNames);

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMushyLayer::writeCheckpointHeader - get scalar var names" << endl;
    }


  // Now write out the ones we want
  // First do scalar comps
  for (int comp = 0; comp < numChkScalarComps; ++comp)
  {
    sprintf(compStr,"component_%d",comp);

    // Work out which variable this is that we want to write out
    int scalarVarComp = m_chkScalarVars[comp];

    header.m_string[compStr] = scalarVarNames[scalarVarComp];
  }

  if (s_verbosity >= 3)
     {
       pout() << "AMRLevelMushyLayer::writeCheckpointHeader - get vector var names" << endl;
     }

  // Now do vector comps. Slightly more complicated because
  // each variable has x, y, (maybe z) components
  for (int i = 0; i < numChkVectorComps; ++i)
  {

    sprintf(compStr,"component_%d",i+numChkScalarComps);
    header.m_string[compStr] = chkVectorVarNames[i];
  }


  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
     {
       pout() << "AMRLevelMushyLayer::writeCheckpointHeader - finished" << endl;
     }

}

/*******/
void
AMRLevelMushyLayer::
writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::writeCheckpointLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"] = m_ref_ratio;
  header.m_real["dx"] = m_dx;
  header.m_real["dt"] = m_dt;
  header.m_real["time"] = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();
  header.m_int ["tag_buffer_size"] = m_opt.tagBufferSize;

  // Setup the periodicity info
  D_TERM(
      if (m_problem_domain.isPeriodic(0))
        header.m_int ["is_periodic_0"] = 1;
      else
        header.m_int ["is_periodic_0"] = 0; ,

        if (m_problem_domain.isPeriodic(1))
          header.m_int ["is_periodic_1"] = 1;
        else
          header.m_int ["is_periodic_1"] = 0; ,

          if (m_problem_domain.isPeriodic(2))
            header.m_int ["is_periodic_2"] = 1;
          else
            header.m_int ["is_periodic_2"] = 0; );

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMushyLayer::writeCheckpointLevel - write header" << endl;
    }

  // Write the header for this level
  header.writeToFile(a_handle);

  // Write the data for this level
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMushyLayer::writeCheckpointLevel - write dbl" << endl;
    }
  write(a_handle,m_scalarNew[0]->boxLayout());

  //  Vector<string> scalarVarNames, vectorVarNames;

  int numChkScalarComps = m_chkScalarVars.size();
  int numChkVectorComps = m_chkVectorVars.size();

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMushyLayer::writeCheckpointLevel - write scalar fields" << endl;
    }

  for (int i=0; i<numChkScalarComps; i++)
  {
    int var = m_chkScalarVars[i];

    write(a_handle,*m_scalarNew[var],m_scalarVarNames[var]);
  }


  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMushyLayer::writeCheckpointLevel - write vector fields" << endl;
    }

  for (int i=0; i<numChkVectorComps; i++)
  {
    int var = m_chkVectorVars[i];
    write(a_handle,*m_vectorNew[var],m_vectorVarNames[var]);
  }

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMushyLayer::writeCheckpointLevel - write adv vel" << endl;
    }

  write(a_handle, m_advVel, "advVel");

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMushyLayer::writeCheckpointLevel - finished" << endl;
    }
}

/*******/
void
AMRLevelMushyLayer::
readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::readCheckpointHeader" << endl;
  }

  // Reader the header
  HDF5HeaderData header;
  header.readFromFile(a_handle);

}

/*******/
void
AMRLevelMushyLayer::
readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::readCheckpointLevel (level " << m_level << ")" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
  {
    MayDay::Error("AMRLevelMushyLayer::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelMushyLayer::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_opt.tagBufferSize= header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_opt.tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelMushyLayer::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRLevelMushyLayer::readCheckpointLevel: file does not contain dt");
  }

  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
  {
    pout() << "read dt = " << m_dt << endl;
  }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
  {
    MayDay::Error("AMRLevelMushyLayer::readCheckpointLevel: file does not contain time");
  }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
  {
    pout() << "read time = " << m_time << endl;
  }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRLevelMushyLayer::readCheckpointLevel: file does not contain prob_domain");
  }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
    isPeriodic[0] = (header.m_int["is_periodic_0"] == 1);
  else
    isPeriodic[0] = false; ,

    if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
      isPeriodic[1] = (header.m_int["is_periodic_1"] == 1);
    else
      isPeriodic[1] = false; ,

      if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
        isPeriodic[2] = (header.m_int["is_periodic_2"] == 1);
      else
        isPeriodic[2] = false;);

  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle,grids);

  if (gridStatus != 0)
  {
    MayDay::Error("AMRLevelMushyLayer::readCheckpointLevel: file does not contain a Vector<Box>");
  }

  // Create level domain
  Vector<int> procs;
  LoadBalance(procs, grids);
  m_grids = DisjointBoxLayout(grids,procs, m_problem_domain);

  // Indicate/guarantee that the indexing below is only for reading
  // otherwise an error/assertion failure occurs
  const DisjointBoxLayout& constGrids = m_grids;

  LayoutIterator lit = constGrids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
  {
    const Box& b = constGrids[lit()];
    m_level_grids.push_back(b);
  }

  if (s_verbosity >= 4)
  {
    pout() << "read level domain: " << endl;
    LayoutIterator lit = m_grids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = m_grids[lit()];
      pout() << lit().intCode() << ": " << b << endl;
    }
    pout() << endl;
  }

  int dataStatus = -1;

  // Resize arrays of data structures
  initDataStructures();
  createDataStructures();

  LevelData<FArrayBox> tempScalar(m_grids, 1);
  LevelData<FArrayBox> tempVector(m_grids, SpaceDim);

  if (s_verbosity >= 10)
  {
    pout() << "AMRLevelMushyLayer::readCheckpointLevel() - read in data" << endl;
  }

  for (int var=0; var<m_numScalarVars; var++)
  {
    if (s_verbosity >= 10)
    {
      pout() << "     reading in variable: " << m_scalarVarNames[var] << endl;
    }

    // Check that this is a variable we have been writing out

    std::vector<int> v = m_chkScalarVars.stdVector();

    if(std::find(v.begin(), v.end(), var) != v.end())
    {

      dataStatus = read<FArrayBox>(a_handle,tempScalar,m_scalarVarNames[var], m_grids);
      tempScalar.copyTo(Interval(0,0), *m_scalarNew[var], Interval(0,0));
      m_scalarNew[var]->exchange();

      if (dataStatus != 0)
      {
        if (s_verbosity >= 10)
        {
          pout() << "       ... variable not found in checkpoint file, continuing anyway " << endl;
        }
        // This is no longer an error - just read in what we can and deal with it
        //      MayDay::Error("AMRLevelMushyLayer::readCheckpointLevel: file does not contain state data");
      }
      else if (s_verbosity >= 10)
      {
        pout() << "       ... done " << endl;
      }

    }
    else
    {
      if (s_verbosity >= 10)
      {
        pout() << "       ... skipping " << endl;
      }
    }
  }


  for (int var=0; var<m_numVectorVars; var++)
  {
    if (s_verbosity >= 10)
    {
      pout() << "     reading in variable: " << m_vectorVarNames[var] << endl;
    }

    std::vector<int> v = m_chkVectorVars.stdVector();

    if(std::find(v.begin(), v.end(), var) != v.end())
    {
      dataStatus =  read<FArrayBox>(a_handle, tempVector,m_vectorVarNames[var], m_grids);
      tempVector.copyTo(Interval(0, SpaceDim-1), *m_vectorNew[var], Interval(0, SpaceDim-1) );
      if (dataStatus != 0)
      {
        if (s_verbosity >= 10)
        {
          pout() << "       ... variable not found in checkpoint file, continuing anyway " << endl;
        }
        // This is no longer an error - just read in what we can and deal with it
        //      MayDay::Error("AMRLevelMushyLayer::readCheckpointLevel: file does not contain state data");
      }
      else if (s_verbosity >= 10)
      {
        pout() << "       ... done " << endl;
      }

    }
    else
    {
      if (s_verbosity >= 10)
      {
        pout() << "       ... skipping " << endl;
      }
    }
  }

//  LevelData<FluxBox> tempFluxBox(m_grids, 1);
//
//  if (s_verbosity >= 10)
//  {
//    pout() << "     reading in advVel" << endl;
//  }
//  dataStatus = read<FluxBox>(a_handle, tempFluxBox, "advVel", m_grids);
//  tempFluxBox.copyTo(Interval(0,0), m_advVel, Interval(0,0));
//  if (s_verbosity >= 10)
//  {
//    pout() << "       ... done " << endl;
//  }

  // Take care of changing dimensionless units if necessary
  ParmParse pp("parameters");
  ParmParse ppMain("main");

  // Default reference temperature/salinity
  Real refTemp = -23;
  Real refSalinity = 230;

  // Check to see if we've specified other reference values
  pp.query("referenceTemperature", refTemp);
  pp.query("referenceSalinity", refSalinity);

  // Default reference values
  Real prevRefTemp = refTemp;
  Real prevRefSalinity = refSalinity;

  // In case we've specified other reference values
  pp.query("prevReferenceTemperature", prevRefTemp);
  pp.query("prevReferenceSalinity", prevRefSalinity);

  if (refTemp != prevRefTemp && refSalinity != prevRefSalinity)
  {
    // Transform H and S fields
    // For now do this very crudely

    // The only case we'll currently handle - transform
    // from Katz to Wells units
    if (refTemp == -3.5 && prevRefTemp == -23
        && refSalinity == 35 && prevRefSalinity == 230)
    {
      // H = H - 1
      // S = S + 1
      pout() << "Converting Katz to Wells units (H=H-1; S=S+1)" << endl;

      for (DataIterator dit = m_scalarNew[0]->dataIterator(); dit.ok(); ++dit)
      {
        (*m_scalarNew[ScalarVars::m_enthalpy])[dit].plus(-1.0);
        (*m_scalarOld[ScalarVars::m_enthalpy])[dit].plus(-1.0);
        (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit].plus(1.0);
        (*m_scalarOld[ScalarVars::m_bulkConcentration])[dit].plus(1.0);
      }
    }
    else
    {
      pout() << "Don't know how to transform these units" << endl;
    }
  }

  Real dtReductionFactor = 1.0;
  ppMain.query("restart_dtReduction", dtReductionFactor);
  m_dt = m_dt/dtReductionFactor;
  pout() << "dt reduced by factor " << dtReductionFactor << " to " << m_dt << endl;

  // Calculate diagnostic variables from (H, S)
  updateEnthalpyVariables();

  // Fill old time states too
  copyNewToOldStates();

  // Set up data structures
  levelSetup();
}

/*******/
void
AMRLevelMushyLayer::
writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::writePlotHeader" << endl;
  }

  ParmParse pp("main");
  bool writeChkFile = false;
    pp.query("writePltAsChkFile", writeChkFile);

    if (writeChkFile)
    {
      writeCheckpointHeader(a_handle);
      return;
    }

  // Setup the number of components -- include space for error
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numOutputComps;

  // Setup the component names
  char compStr[30];
  Vector<string> varNames(m_numOutputComps);
  getVarNames(varNames, m_outputScalarVars, m_outputVectorVars,
              m_scalarVarNames, m_vectorVarNames);

  int comp = 0;
  for (comp = 0; comp < m_numOutputComps; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = varNames[comp];
  }

  // Write out the key parameters
  m_parameters.writeToHDF5(header);



  // Write the header
  header.writeToFile(a_handle);




}



/*******/
void
AMRLevelMushyLayer::
writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::writePlotLevel (level " << m_level << ")" << endl;
  }

  // Add some debugging options
  bool writeBoxes = true;
  bool writeScalars = true;
  bool writeVectors = true;

  ParmParse pp("main");
  pp.query("writeBoxes", writeBoxes);
  pp.query("writeScalars", writeScalars);
  pp.query("writeVectors", writeVectors);

  bool writeSeparateData = false;
  pp.query("writePltFileDataSeperately", writeSeparateData);

  bool writeChkFile = false;
  pp.query("writePltAsChkFile", writeChkFile);

  if (writeChkFile)
  {
    writeCheckpointLevel(a_handle);
    return;
  }


  if (writeSeparateData && s_verbosity >= 5)
    {
    pout() << "Writing data separately " << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"] = m_ref_ratio;
  header.m_real["dx"] = m_dx;
  header.m_real["dt"] = m_dt;
  header.m_real["time"] = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();
  header.m_int ["timestepFailed"] = m_timestepFailed;

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::writePlotLevel - done setup level header information" << endl;
  }

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::writePlotLevel - done write to handle" << endl;
  }

  const DisjointBoxLayout& levelGrids = m_scalarNew[0]->getBoxes();
    if (s_verbosity >= 5) { pout() << "AMRLevelMushyLayer::writePlotLevel - got levelGrids" << endl; }

    if (writeBoxes)
    {
    write(a_handle,levelGrids);
    }

  // Write the data for this level
  LevelData<FArrayBox> plotData(levelGrids, m_numOutputComps);

  // first copy data to plot data holder
  if (writeScalars)
  {

    for (int i=0; i < m_outputScalarVars.size(); i++)
    {

      if (writeSeparateData)
      {
        int var = m_outputScalarVars[i];
        write(a_handle,*m_scalarNew[var],m_scalarVarNames[var]);
      }
      else
      {
        m_scalarNew[m_outputScalarVars[i]]->copyTo(m_scalarNew[i]->interval(), plotData,
                                                   Interval(i, i));

        // Subtract 1 from lambda before writing out
        if (m_outputScalarVars[i] == m_lambda)
        {
          for (DataIterator dit = plotData.dataIterator(); dit.ok(); ++dit)
          {
            plotData[dit].plus(-1, i);
          }
        }
      }

    }

  }

  if (writeVectors)
  {

    Interval vecSrcComps(0, SpaceDim-1);
    for (int i = 0; i<m_outputVectorVars.size(); i++)
    {
      int startComp = m_outputScalarVars.size() + (i*SpaceDim);
      Interval destComps(startComp, startComp + SpaceDim-1);

      if (writeSeparateData)
      {
          int var = m_outputVectorVars[i];
          write(a_handle,*m_vectorNew[var],m_vectorVarNames[var]);
      }
      else
      {


      m_vectorNew[m_outputVectorVars[i]]->copyTo(vecSrcComps, plotData, destComps);
      }
    }
  }

//  plotData.exchange(plotData.interval());

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::writePlotLevel - made plotData" << endl;
  }




  if (s_verbosity >= 5) {  pout() << "AMRLevelMushyLayer::writePlotLevel - wrote levelGrids" << endl;}

  if(!writeSeparateData)
  {
  write(a_handle,plotData,"data");
  }

  if (s_verbosity >= 5) { pout() << "AMRLevelMushyLayer::writePlotLevel - data written" << endl; }

//  write(a_handle, m_advVel, "advVel");

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::writePlotLevel - finished all" << endl;
  }

}

void AMRLevelMushyLayer::writePlotFile(int iter)
{

  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::writePlotFile() " << endl;
  }
  string  m_plotfile_prefix;
  ParmParse pp("main");
  pp.query("plot_prefix", m_plotfile_prefix);

  int m_cur_step = m_time/m_dt;

  string iter_str = m_plotfile_prefix;

  char suffix[100];
  sprintf(suffix,"%06dlev%diter%d.%dd.hdf5",m_cur_step,m_level,iter,SpaceDim);

  iter_str += suffix;

  HDF5Handle handle(iter_str.c_str(), HDF5Handle::CREATE);

  // write amr data
  HDF5HeaderData header;
  header.m_int ["max_level"]  = 0;
  header.m_int ["num_levels"] = 1;
  header.m_int ["iteration"]  = m_cur_step;
  header.m_real["time"]       = m_time;

  // should steps since regrid be in the checkpoint file?
  header.writeToFile(handle);

  // write physics class header data
  writePlotHeader(handle);

  // write physics class per-level data
  writePlotLevel(handle);


  handle.close();
}


#endif

void AMRLevelMushyLayer::getExtraPlotFields()
{


  EdgeToCell(m_advVel, *m_vectorNew[VectorVars::m_advectionVel]);

  // Copy pi and phi from projection
  if (m_projection.isInitialized())
  {
    m_projection.Pi().copyTo(*m_scalarNew[ScalarVars::m_pi]);
    m_projection.phi().copyTo(*m_scalarNew[ScalarVars::m_phi]);

    m_projection.Pi().copyTo(*m_scalarNew[ScalarVars::m_MACBC]);

    m_projection.MACrhs().copyTo(*m_scalarNew[ScalarVars::m_MACrhs]);
    m_projection.CCrhs().copyTo(*m_scalarNew[ScalarVars::m_CCrhs]);

    m_projection.CCcorrection().copyTo(*m_vectorNew[VectorVars::CCcorrection]);
    EdgeToCell( m_projection.MACcorrection(), *m_vectorNew[VectorVars::m_MACcorrection]);

    EdgeToCell(m_projection.grad_eLambda(), *m_vectorNew[VectorVars::m_freestreamCorrection]);

    Divergence::levelDivergenceMAC(*m_scalarNew[ScalarVars::m_divUcorr], m_projection.grad_eLambda(), m_dx);

    bool scaleMACBCWithChi = false;
    Real BCscale = 1.0;
    ParmParse pp("projection");
    pp.query("scaleMACBCWithChi", scaleMACBCWithChi);

    pp.query("MACbcScale", BCscale);

//    Real phiScale = m_projection.getPhiScale(m_dt);

    AMRLevelMushyLayer* coarsestML = this->getCoarsestLevel();
    Real coarsestDt = coarsestML->dt();

    for (DataIterator dit = m_scalarNew[ScalarVars::m_MACBC]->dataIterator(); dit.ok(); ++dit)
    {
      (*m_scalarNew[ScalarVars::m_MACBC])[dit].mult(0.5*m_dt);

//      (*m_scalarNew[ScalarVars::m_CCrhs])[dit].mult(m_dt);

//      (*m_scalarNew[ScalarVars::m_MACrhs])[dit].mult(phiScale);


//      if (scaleMACBCWithChi)
//      {
//  //      (*m_scalarNew[ScalarVars::m_MACBC])[dit].mult((*m_scalarNew[ScalarVars::m_porosity])[dit]);
//        (*m_scalarNew[ScalarVars::m_MACBC])[dit].mult(BCscale);
//      }

      // Scale MAC correction
      (*m_vectorNew[VectorVars::m_MACcorrection])[dit].divide(m_dt/coarsestDt);
    }

  }

}

void AMRLevelMushyLayer::computeVorticityStreamfunction()
{

  CH_TIME("AMRLevelMushyLayer::computeVorticityStreamfunction");

  // First vorticity = curl (u)
  computeVorticity();

  // Now streamfunction psi - given by the solution to lap(psi) = - vorticity

  Vector<AMRLevelMushyLayer*> hierarchy;
  Vector<DisjointBoxLayout> grids;
  Vector<int> refRat;
  ProblemDomain lev0Dom;
  Real lev0Dx;
  getHierarchyAndGrids(hierarchy, grids, refRat, lev0Dom, lev0Dx);
  int nLevels = grids.size();

  // First let's get a poisson op
  AMRPoissonOpFactory op;
  Real alpha = 0.0;
  Real beta = 1.0;
  //  op.define(m_problem_domain, m_grids, m_dx, m_physBCPtr->streamFunctionBC(), maxDepth, alpha, beta);

  op.define(lev0Dom, grids, refRat, lev0Dx, m_physBCPtr->streamFunctionBC(), alpha, beta);



  // Need a bottom solver
//  RelaxSolver<LevelData<FArrayBox> > bottomSolver;
  BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_verbosity = m_opt.verbosity_multigrid;

  // And a multigrid solver
  AMRMultiGrid<LevelData<FArrayBox> > solver;
  AMRLevelOpFactory<LevelData<FArrayBox> >& castFact =
      (AMRLevelOpFactory<LevelData<FArrayBox> >&) op;

  solver.define(lev0Dom, castFact,
                &bottomSolver, nLevels);

  solver.m_verbosity = m_opt.verbosity_multigrid;

  // We were spending a lot of time doing this solve
  // Fix these to try and make the solve quicker
  solver.m_eps = 10e-10;
  solver.m_pre = 4;
  solver.m_post = 4;

  // Construct RHS, including coarse BCs if necessary
  Vector<LevelData<FArrayBox>* > amrPsi;
  Vector<LevelData<FArrayBox>* > amrVorticity;

  for (int lev = 0; lev < nLevels; lev++)
  {
    AMRLevelMushyLayer* ml = hierarchy[lev];

    amrPsi.push_back(&(*ml->m_scalarNew[ScalarVars::m_streamfunction]));
    amrVorticity.push_back(&(*ml->m_scalarNew[ScalarVars::m_vorticity]));
  }

  // now solve on this level
  solver.solve(amrPsi, amrVorticity, m_level,
               m_level, false); // don't initialize to zero


}

void AMRLevelMushyLayer::computeVorticity()
{
  CH_TIME("AMRLevelMushyLayer::computeVorticity");



  // this breaks the "const"-ness of the function,
  // but is necessary to ensure that boundary
  // conditions are properly set.
  LevelData<FArrayBox>* velPtr = &(*m_vectorNew[VectorVars::m_fluidVel]);
  LevelData<FArrayBox>& vel =
      *(const_cast<LevelData<FArrayBox>*>(&*velPtr));

  const DisjointBoxLayout& grids = vel.getBoxes();

  if (m_level > 0)
  {
    // do quadratic C/F BC's
    // for now, assume that BC's are with new velocity
    // (may need to be changed for fine-level regridding)
    LevelData<FArrayBox>& crseVel = (*getCoarserLevel()->m_vectorNew[VectorVars::m_fluidVel]);
    const DisjointBoxLayout& crseGrids = crseVel.getBoxes();
    int nRefCrse = getCoarserLevel()->refRatio();

    QuadCFInterp interpolator(grids, &crseGrids, m_dx, nRefCrse,
                              SpaceDim, m_problem_domain);

    interpolator.coarseFineInterp(vel, crseVel);
  }

  Interval velComps(0,SpaceDim-1);
  vel.exchange(velComps);

  DataIterator dit = vel.dataIterator();

  // at the moment, do extrap at physical boundaries
  // (same effect as one-sided differencing)
  BCHolder vortVelBC = m_physBCPtr->extrapFuncBC();  //m_physBCPtr->extrapFuncBC(1);

  for (dit.reset(); dit.ok(); ++dit)
  {
    // apply physical BC's
    vortVelBC(vel[dit], grids[dit],
              m_problem_domain, m_dx,
              true); // homogeneous
    if (SpaceDim == 3)
    {
      for (int dir=0; dir<SpaceDim; dir++)
      {
        //todo - 3D: vorticity needs more components in 3D, can't really store it in m_scalarNew.
        // will need to make a new data structure which has the write number of componenets for the dimensionality of the problem

        // and then compute vorticity
//        FORT_COMPUTEVORT(CHF_FRA1((*m_scalarNew[ScalarVars::m_vorticity])[dit],dir),
//                         CHF_CONST_FRA(vel[dit]),
//                         CHF_BOX(grids[dit]),
//                         CHF_CONST_REAL(m_dx),
//                         CHF_CONST_INT(dir));
      }
    }
    else if (SpaceDim == 2)
    {
      int dir = 2;
      FORT_COMPUTEVORT(CHF_FRA1((*m_scalarNew[ScalarVars::m_vorticity])[dit],0),
                       CHF_CONST_FRA(vel[dit]),
                       CHF_BOX(grids[dit]),
                       CHF_CONST_REAL(m_dx),
                       CHF_CONST_INT(dir));
    } // end dimensionality choice
  } // end loop over grids

}

