#include "MushyLayerLoadUtils.H"
#include "AMRLevelMushyLayerFactory.H"
#include "MushyLayerSubcycleUtils.H"

void getAMRHierarchy(string inFile, Vector<AMRLevelMushyLayer*>& amrlevels, int& finest_level, HDF5HeaderData& header)
{
  // Get variables for initializing amrlevelmushylayer objects

  Vector<int> steps_since_regrid; //regrid_intervals,
  std::vector<int> regrid_intervals; // (num_read_levels,1);

  std::vector<bool> periodic;

  ParmParse ppMain("main");

  int max_level = 0;
  ppMain.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  ppMain.getarr("regrid_interval",regrid_intervals,0,num_read_levels);

  ppMain.getarr("periodic_bc",periodic,0,SpaceDim-1);

  HDF5Handle handle(inFile,HDF5Handle::OPEN_RDONLY);
  // read from checkpoint file

  // Shamelessly copied from AMR::setupForRestart
  handle.setGroup("/");
  header.readFromFile(handle);

  // read finest level
  if (header.m_int.find("num_levels") == header.m_int.end())
  {
    MayDay::Error("AMR::restart: checkpoint file does not contain num_levels");
  }
  int num_levels = header.m_int ["num_levels"];
  finest_level = num_levels - 1;



  //int cur_step = header.m_int ["iteration"];

  steps_since_regrid.resize(max_level+1, 0);
  for (int level = 0; level < max_level; ++level)
  {
    char headername[100];
    sprintf(headername, "steps_since_regrid_%d", level);
    if (header.m_int.find(headername) != header.m_int.end())
    {
      steps_since_regrid[level] = header.m_int [headername];
    }
  }

  // Need level 0 header to get domain box

  // Setup the level string
    char levelStr[20];
    sprintf(levelStr,"%d",0);
    const std::string label = std::string("level_") + levelStr;

    // Read the header for this level
    handle.setGroup(label);

    HDF5HeaderData lev0header;
    lev0header.readFromFile(handle);

  Box domainBox = lev0header.m_box["prob_domain"];

    // Get the periodicity info -- this is more complicated than it really
    // needs to be in order to preserve backward compatibility
    bool isPeriodic[SpaceDim];
    D_TERM(if (!(lev0header.m_int.find("is_periodic_0") == lev0header.m_int.end()))
      isPeriodic[0] = (lev0header.m_int["is_periodic_0"] == 1);
    else
      isPeriodic[0] = false; ,

      if (!(lev0header.m_int.find("is_periodic_1") == lev0header.m_int.end()))
        isPeriodic[1] = (lev0header.m_int["is_periodic_1"] == 1);
      else
        isPeriodic[1] = false; ,

        if (!(lev0header.m_int.find("is_periodic_2") == lev0header.m_int.end()))
          isPeriodic[2] = (lev0header.m_int["is_periodic_2"] == 1);
        else
          isPeriodic[2] = false;);

  ProblemDomain problemDomain = ProblemDomain(domainBox,isPeriodic);

  // read physics class header data
  amrlevels.resize(num_levels);

  RefCountedPtr<AMRLevelMushyLayerFactory> mlFact;
  getAMRFactory(mlFact);

  for (int level = 0; level <= finest_level; ++level)
  {

    amrlevels[level] = (static_cast <AMRLevelMushyLayer*> (mlFact->new_amrlevel()) );

    AMRLevelMushyLayer* crsePtr = NULL;
    if (level > 0)
    {
      crsePtr =  amrlevels[level-1];
    }
    // Do these actually need to be correct?
    int refRatio = 1;

    amrlevels[level]->define(crsePtr, problemDomain, level, refRatio);

    // reset to root group to read header info
    handle.setGroup("/");
    amrlevels[level]->readCheckpointHeader(handle);

    amrlevels[level]->readCheckpointLevel(handle);

    amrlevels[level]->defineSolvers(amrlevels[level]->time());
  }

  handle.close();
}


void addExtraParams(string runInputs, ParmParse& pp)
{
  pout() << "addExtraParams() from " << runInputs << endl;
  const char* pFile = runInputs.c_str();
    std::ifstream t(pFile);
    std::string str((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());
    pp.addEntries(str);
}
