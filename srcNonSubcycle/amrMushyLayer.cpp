/* _______              __
  / ___/ /  ___  __ _  / /  ___
 / /__/ _ \/ _ \/  ' \/ _ \/ _ \
 \___/_//_/\___/_/_/_/_.__/\___/ 
 */
//
// This software is copyright (C) by the Lawrence Berkeley
// National Laboratory.  Permission is granted to reproduce
// this software for non-commercial purposes provided that
// this notice is left intact.
// 
// It is acknowledged that the U.S. Government has rights to
// this software under Contract DE-AC03-765F00098 between
// the U.S.  Department of Energy and the University of
// California.
//
// This software is provided as a professional and academic
// contribution for joint exchange. Thus it is experimental,
// is provided ``as is'', with no warranties of any kind
// whatsoever, no support, no promise of updates, or printed
// documentation. By using this software, you acknowledge
// that the Lawrence Berkeley National Laboratory and
// Regents of the University of California shall have no
// liability with respect to the infringement of other
// copyrights by any part of this software.
//
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>

using std::ifstream; 
using std::ios;

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include "amrMushyLayer.H"

#include "parstream.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "computeNorm.H"
#include "utils_F.H"

//#include "FourthOrderFillPatch.H"

// small parameter defining when times are equal
#define TIME_EPS 1.0e-12

int amrMushyLayer::s_verbosity = 1;

bool 
amrMushyLayer::isDefined() const
{
  return m_is_defined;
}


amrMushyLayer::amrMushyLayer()
{
  setDefaults();
}

amrMushyLayer::~amrMushyLayer()
{
  // clean up memory
  for (int lev=0; lev<m_fluidAdv.size(); lev++)
  {



    // Vector vars
    for (int a_var = 0; a_var<m_numVectorVars; a_var++)
    {

      if (m_vectorOld[a_var][lev] != NULL)
      {
        delete m_vectorOld[a_var][lev];
        m_vectorOld[a_var][lev] = NULL;
      }

      if (m_vectorNew[a_var][lev] != NULL)
      {
        delete m_vectorNew[a_var][lev];
        m_vectorNew[a_var][lev] = NULL;
      }

      if (m_dVector[a_var][lev] != NULL)
      {
        delete m_dVector[a_var][lev];
        m_dVector[a_var][lev] = NULL;
      }

    }

    // Advection
    if (m_fluidAdv[lev] != NULL)
    {
      delete m_fluidAdv[lev];
      m_fluidAdv[lev] = NULL;
    }
    if (m_frameAdv[lev] != NULL)
    {
      delete m_frameAdv[lev];
      m_frameAdv[lev] = NULL;
    }


  } // end loop over currently defined levels

  //	for (int lev=0; lev<=m_max_level; lev++)
  //	{
  //		if (m_projection[lev] != NULL)
  //		{
  //			delete m_projection[lev];
  //			m_projection[lev] = NULL;
  //		}
  //	}


}

void
amrMushyLayer::run(Real a_max_time, int a_max_step)
{
  CH_TIME("amrMushyLayer::run()");

  if (s_verbosity > 3)
  {
    pout() << "amrMushyLayer::run -- max_time= " << a_max_time
        << ", max_step = " << a_max_step << endl;
  }

  // only call computeInitialDt if we're not doing restart
  if (!m_do_restart)
  {
    computeInitialDt();
  }
  else
  {
    computeDt();
  }

  if (m_printAnalyticSoln)
  {
    //Enforce analytic solution and regrid to a single level
    enforceAnalyticSolution();
    m_max_level=0;
    regrid();

    //Make sure we skip all timestepping and go straight to writing out the solution
    a_max_step=0;
  }


  // advance solution until done
  while ((        (a_max_time > m_time)
      && (m_cur_step < a_max_step) //Stop if we reach max_step
      && (m_dt > TIME_EPS) // Stop if dt becomes stupidly small
      && (!m_steady_state) //Stop if we reach steady state
  )
      ||  (m_cur_step < m_minTimestep) // Don't stop until we've completed the min number of timesteps
  )
  {

    // dump plotfile before regridding
    if ((m_cur_step%m_plot_interval == 0) && m_plot_interval > 0)
    {
      writePlotFile();
    }

    /* Regrid if:
     * Not the first timestep
     * AND
     * Time step is integer multiple of the regrid interval
     * AND
     * Regrid interval is > 0
     * AND
     * Haven't just restarted
     */
    if ((m_cur_step != 0)
        && (m_cur_step%m_regrid_interval ==0)
        && m_regrid_interval > 0
        && m_restart_step != m_cur_step)
    {
      regrid();
    }

    // compute dt after regridding in case number of levels has changed
    computeDt();
    if ((a_max_time - m_time) < m_dt)
    {
      m_dt = a_max_time - m_time;
    }


    if ((m_cur_step%m_check_interval == 0) && (m_check_interval > 0)
        && (m_cur_step != m_restart_step))
    {
      writeCheckpointFile();
    }

    timeStep();

    //Check if we've reached steady state
    convergedToSteadyState();
  }

  //Put some logging in to determine why we finished
  std::ostringstream oss;
  oss << "amrMushyLayer::run() - Simulation finished. \n " <<
      "============================ \n" <<
      "a_max_time = " << a_max_time << ", m_time = " << m_time << "\n" <<
      "m_cur_step = " << m_cur_step << ", a_max_step = " << a_max_step << "\n" <<
      "m_dt = " << m_dt << ", TIME_EPS = " << TIME_EPS << "\n" <<
      "m_steady_state = " << m_steady_state << "\n" <<
      "total_num_iterations = " << m_total_num_iterations << "\n"
      "============================ \n";

  string msg = oss.str();
  logMessage(1, msg);

  // dump out final plotfile, if appropriate
  if (m_plot_interval >= 0)
  {
    writePlotFile();
    writeCheckpointFile();
    //		writeSolutionToTextFile(0);
  }

}

void amrMushyLayer::
averageCoarseToFineSolutions()
{
  CH_TIME("amrMushyLayer::averageCoarseToFineSolutions()");

  if (m_finest_level != 0)
  {
    for (int lev=m_finest_level; lev>0; lev--)
    {
      // probably will make sense to predefine these and
      // keep them around rather than re-allocating them
      // every timestep
      CoarseAverage avgDownScalar(m_amrGrids[lev], 1, m_refinement_ratios[lev-1]);
      CoarseAverage avgDownVector(m_amrGrids[lev], SpaceDim, m_refinement_ratios[lev-1]);

      for (int a_var=0; a_var<m_numVars; a_var++)
      {
        avgDownScalar.averageToCoarse(*m_scalarNew[a_var][lev-1], *m_scalarNew[a_var][lev]);
        avgDownScalar.averageToCoarse(*m_scalarOld[a_var][lev-1], *m_scalarOld[a_var][lev]);
      }

      for (int a_var=0; a_var<m_numVectorVars; a_var++)
      {
        avgDownVector.averageToCoarse(*m_vectorNew[a_var][lev-1], *m_vectorNew[a_var][lev]);
        avgDownVector.averageToCoarse(*m_vectorOld[a_var][lev-1], *m_vectorOld[a_var][lev]);
      }
    }
  }

}


void amrMushyLayer::
updateEnthalpyVariables()
{
  CH_TIME("amrMushyLayer::updateEnthalpyVariables");

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    //Get level data's from the vector of all scalar variables to de-clutter the remaining code.
    LevelData<FArrayBox>& HC = *m_HC[lev];

    LevelData<FArrayBox>& theta = *m_scalarNew[m_theta][lev];
    LevelData<FArrayBox>& compositionLiquid = *m_scalarNew[m_compositionLiquid][lev];
    LevelData<FArrayBox>& compositionSolid = *m_scalarNew[m_compositionSolid][lev];
    LevelData<FArrayBox>& porosity = *m_scalarNew[ScalarVars::m_porosity][lev];

    LevelData<FArrayBox>& solidus = *m_scalarNew[ScalarVars::m_enthalpySolidus][lev];
    LevelData<FArrayBox>& liquidus = *m_scalarNew[ScalarVars::m_enthalpyLiquidus][lev];
    LevelData<FArrayBox>& eutectic = *m_scalarNew[ScalarVars::m_enthalpyEutectic][lev];



    ::updateEnthalpyVariables(HC, theta, compositionLiquid, compositionSolid, porosity,
                              solidus, liquidus, eutectic,
                              m_parameters);

//    LevelData<FArrayBox>& solidFraction = *m_scalarNew[m_solidFraction][lev];
//    for (DataIterator dit = theta.dataIterator(); dit.ok(); ++dit)
//    {
//      solidFraction[dit].setVal(1.0);
//      solidFraction[dit].minus(porosity[dit]);
//    }
  }
}

/*
 * Requires both 'from' and 'to' to be defined across m_amrGrids
 */
void
amrMushyLayer::activeLevelCopy(Vector<RefCountedPtr<LevelData<FArrayBox> > > from,
                               Vector<RefCountedPtr<LevelData<FArrayBox> > > to)
{
  CH_TIME("amrMushyLayer::activeLevelCopy()");

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    //This doesn't copy ghostcells
    //from[lev]->copyTo(*to[lev]);

    //This method should include ghost cells
    for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      (*to[lev])[dit].copy((*from[lev])[dit]);
    }
  }
}


// do regridding
void
amrMushyLayer::regrid()
{

  if (s_verbosity > 3)
  {
    pout() << "amrMushyLayer::regrid" << endl;
  }

  // only do any of this if the max level > 0
  if (m_max_level > 0)
  {

    // in this code, lbase is always 0
    int lbase =0;

    // first generate tags
    Vector<IntVectSet> tagVect(m_max_level);
    tagCells(tagVect);

    // now generate new boxes
    int top_level = min(m_finest_level, m_max_level-1);
    Vector<Vector<Box> > old_grids(m_finest_level+1);
    Vector<Vector<Box> > new_grids;

    // this is clunky, but i don't know of a better way to turn
    // a DisjointBoxLayout into a Vector<Box>
    for (int lev=0; lev<= m_finest_level; lev++)
    {
      const DisjointBoxLayout& levelDBL = m_amrGrids[lev];
      old_grids[lev].resize(levelDBL.size());
      LayoutIterator lit = levelDBL.layoutIterator();
      int boxIndex = 0;
      for (lit.begin(); lit.ok(); ++lit, ++boxIndex)
      {
        old_grids[lev][boxIndex] = levelDBL[lit()];
      }
    }


    int new_finest_level;

    //if (procID() == uniqueProc(SerialTask::compute))
    //{

    BRMeshRefine meshrefine(m_amrDomains[0], m_refinement_ratios,
                            m_fill_ratio, m_block_factor,
                            m_nesting_radius, m_max_box_size);



    new_finest_level = meshrefine.regrid(new_grids, tagVect,
                                         lbase, top_level,
                                         old_grids);

    //}
    //broadcast(new_finest_level, uniqueProc(SerialTask::compute));
    //broadcast(new_grids, uniqueProc(SerialTask::compute));

    //int numLevels = Min(new_finest_level, m_max_level)+1;

    // now loop through levels and redefine if necessary
    logMessage(8,  "    amrMushyLayer::regrid - redefine levels");


    for (int lev=lbase+1; lev<= new_finest_level; ++lev)
    {
      int numGridsNew = new_grids[lev].size();
      Vector<int> procIDs(numGridsNew);
      LoadBalance(procIDs, new_grids[lev]);
      const DisjointBoxLayout newDBL(new_grids[lev], procIDs,
                                     m_amrDomains[lev]);
      m_amrGrids[lev] = newDBL;

      // may eventually want to do post-regrid smoothing on this!
      FineInterp interpolatorScalar(newDBL, 1,
                                    m_refinement_ratios[lev-1],
                                    m_amrDomains[lev]);

      //Do this for each scalar variable
      logMessage(20,  "    amrMushyLayer::regrid - build new storage");
      for (int a_var=0; a_var < m_numVars; a_var++)
      {
        // get old storage
        RefCountedPtr<LevelData<FArrayBox> > old_oldDataPtr = m_scalarOld[a_var][lev];
        RefCountedPtr<LevelData<FArrayBox> > old_newDataPtr = m_scalarNew[a_var][lev];

        // build new storage
        RefCountedPtr<LevelData<FArrayBox> > new_oldDataPtr
        (new LevelData<FArrayBox>(newDBL, 1,m_ghostVect));
        RefCountedPtr<LevelData<FArrayBox> > new_newDataPtr
        (new LevelData<FArrayBox>(newDBL, 1,m_ghostVect));


        //We don't have to do this, but it makes debugging easier
        for (DataIterator dit=new_newDataPtr->dataIterator(); dit.ok(); ++dit)
        {
          (*new_newDataPtr)[dit].setVal(0);
          (*new_oldDataPtr)[dit].setVal(0);
        }

        // For m_dScalar, we just need new grids (don't worry about them being populated)
        RefCountedPtr<LevelData<FArrayBox> > new_dDataPtr
        (new LevelData<FArrayBox>(newDBL, 1,m_ghostVect));

        // fill with interpolated data from coarser level
        interpolatorScalar.interpToFine(*new_newDataPtr, *m_scalarNew[a_var][lev-1]);
        interpolatorScalar.interpToFine(*new_oldDataPtr, *m_scalarOld[a_var][lev-1]);

        // now copy old-grid data into new holder
        if (old_newDataPtr != NULL)
        {
          Interval dataComps = new_newDataPtr->interval();
          old_newDataPtr->copyTo(dataComps, *new_newDataPtr, dataComps);
          old_oldDataPtr->copyTo(dataComps, *new_oldDataPtr, dataComps);

          // can now delete old data
          //This actually seems to be a bad idea - don't call delete on a ref counted ptr
          //					delete old_oldDataPtr;
          //					delete old_newDataPtr;
        }


        // now copy new holders into multilevel arrays
        m_scalarOld[a_var][lev] = new_oldDataPtr;
        m_scalarNew[a_var][lev] = new_newDataPtr;
        m_dScalar[a_var][lev] = new_dDataPtr;

      }


      //Do the same for vector variables

      // may eventually want to do post-regrid smoothing on this!
      FineInterp interpolatorVector(newDBL, SpaceDim,
                                    m_refinement_ratios[lev-1],
                                    m_amrDomains[lev]);

      logMessage(20,  "    amrMushyLayer::regrid - build new storage (vector vars)");
      for (int a_var=0; a_var < m_numVectorVars; a_var++)
      {
        // build new storage
        LevelData<FArrayBox>* old_oldDataPtr = m_vectorOld[a_var][lev];
        LevelData<FArrayBox>* old_newDataPtr = m_vectorNew[a_var][lev];
        LevelData<FArrayBox>* old_dDataPtr = m_dVector[a_var][lev];

        LevelData<FArrayBox>* new_oldDataPtr = new LevelData<FArrayBox>(newDBL, SpaceDim,m_ghostVect);
        LevelData<FArrayBox>* new_newDataPtr = new LevelData<FArrayBox>(newDBL, SpaceDim,m_ghostVect);

        // For m_dScalar, we just need new grids (don't worry about them being populated)
        LevelData<FArrayBox>* new_dDataPtr =	new LevelData<FArrayBox>(newDBL, SpaceDim,m_ghostVect);


        interpolatorVector.interpToFine(*new_newDataPtr, *m_vectorNew[a_var][lev-1]);
        interpolatorVector.interpToFine(*new_oldDataPtr, *m_vectorOld[a_var][lev-1]);

        // now copy old-grid data into new holder
        if (old_newDataPtr != NULL)
        {
          Interval dataComps = new_newDataPtr->interval();
          old_newDataPtr->copyTo(dataComps, *new_newDataPtr, dataComps);
          old_oldDataPtr->copyTo(dataComps, *new_oldDataPtr, dataComps);

          // can now delete old data
          delete old_oldDataPtr;
          delete old_newDataPtr;
          delete old_dDataPtr;
        }


        // now copy new holders into multilevel arrays
        m_vectorOld[a_var][lev] = new_oldDataPtr;
        m_vectorNew[a_var][lev] = new_newDataPtr;
        m_dVector[a_var][lev] = new_dDataPtr;
      }




      //Now fill ghost cells from coarser level.
      //			PiecewiseLinearFillPatch filler(m_amrGrids[lev], //this level
      //					m_amrGrids[lev-1], //coarser level
      //					1, //1 component
      //					m_amrDomains[lev-1], //domain on the coarser level
      //					m_refinement_ratios[lev-1], //refinement ratio of the coarser level
      //					m_num_ghost);

      QuadCFInterp quadFillerScalar(m_amrGrids[lev], //this level
                                    &(m_amrGrids[lev-1]), //coarser level
                                    m_amrDx[lev], //fine dx
                                    m_refinement_ratios[lev-1], //refinement ratio of the coarser level
                                    1,  // 1 component
                                    m_amrDomains[lev]); //problem domain on the fine level



      QuadCFInterp quadFillerVector(m_amrGrids[lev], //this level
                                    &(m_amrGrids[lev-1]), //coarser level
                                    m_amrDx[lev], //fine dx
                                    m_refinement_ratios[lev-1], //refinement ratio of the coarser level
                                    SpaceDim,  // multiple components
                                    m_amrDomains[lev]); //problem domain on the fine level




      logMessage(20,  "    amrMushyLayer::regrid - fill ghost cells");
      for (int a_var=0; a_var < m_numVars; a_var++)
      {



        //This fills exterior ghost cells
        applyBCs(a_var, lev);

        quadFillerScalar.coarseFineInterp(*m_scalarNew[a_var][lev], *m_scalarNew[a_var][lev-1]);
        quadFillerScalar.coarseFineInterp(*m_scalarOld[a_var][lev], *m_scalarOld[a_var][lev-1]);
        //				fillerScalar.fillInterp(*m_scalarNew[a_var][lev], *m_scalarNew[a_var][lev-1],0,0,1);
        //				fillerScalar.fillInterp(*m_scalarOld[a_var][lev], *m_scalarOld[a_var][lev-1],0,0,1);


      }

      for (int a_var =0; a_var<m_numVectorVars; a_var++)
      {
        quadFillerVector.coarseFineInterp(*m_vectorNew[a_var][lev], *m_vectorNew[a_var][lev-1]);
        quadFillerVector.coarseFineInterp(*m_vectorOld[a_var][lev], *m_vectorOld[a_var][lev-1]);

        //Fill physical boundary ghost cells
        applyVectorBCs(a_var, lev);
      }

    } // end loop over currently defined levels

    // Re-populate advection fluxboxes (these doesn't change, but the grids that they're defined on do)
    // We do this AFTER we've filled ghost cells
    for (int lev=0; lev <= new_finest_level; lev++)
    {
      getFrameAdvection(lev);
      getFluidAdvection(lev);
    }
    logMessage(20,  "    amrMushyLayer::regrid - got fluid/frame advection");

    for (int lev=m_scalarOld[0].size()-1; lev > new_finest_level; lev--)
    {
      //We should actually delete these! (Maybe not)
      DisjointBoxLayout emptyDBL;
      m_amrGrids[lev] = emptyDBL;
      //			m_amrGrids.pop_back();
    }

    m_finest_level = new_finest_level;


    // finally, set up covered_level flags
    if (s_verbosity > 4)
    {
      pout() << "    amrMushyLayer::regrid - new finest level = "  << m_finest_level << endl;
    }
    m_covered_level.resize(m_max_level+1, 0);

    // note that finest level can't be covered.
    for (int lev=m_finest_level-1; lev>=0; lev--)
    {

      // if the next finer level is covered, then this one is too.
      if (m_covered_level[lev+1] == 1)
      {
        m_covered_level[lev] = 1;
      }
      else
      {
        // see if the grids finer than this level completely cover it
        IntVectSet fineUncovered(m_amrDomains[lev+1].domainBox());
        const DisjointBoxLayout& fineGrids = m_amrGrids[lev+1];

        LayoutIterator lit = fineGrids.layoutIterator();
        for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& thisBox = fineGrids.get(lit());
          fineUncovered.minus_box(thisBox);
        }

        if (fineUncovered.isEmpty())
        {
          m_covered_level[lev] = 1;
        }
      }
    } // end loop over levels to determine covered levels
  } // end if max level > 0 in the first place

  //Regrid flux registers
  logMessage(8, "    amrMushyLayer::regrid - redefine flux registers");


  for (int lev=0; lev<m_finest_level; lev++)
  {
    for (int a_var=0; a_var<m_numVars; a_var++)
    {
      m_fluxRegister[a_var][lev] = RefCountedPtr<LevelFluxRegister>
      (new LevelFluxRegister(m_amrGrids[lev+1], m_amrGrids[lev], m_amrDomains[lev+1], m_refinement_ratios[lev], 1));
      m_fluxRegister[a_var][lev]->setToZero();
    }
    for (int a_var=0; a_var<m_numVectorVars; a_var++)
    {
      m_vectorFluxRegister[a_var][lev] = RefCountedPtr<LevelFluxRegister>
      (new LevelFluxRegister(m_amrGrids[lev+1], m_amrGrids[lev], m_amrDomains[lev+1], m_refinement_ratios[lev], SpaceDim));
      m_vectorFluxRegister[a_var][lev]->setToZero();
    }
  }

  //Regrid projection operators
  //	setupProjectionOperators();

  logMessage(8,  "    amrMushyLayer::regrid - finished regridding");

} 



void 
amrMushyLayer::tagCells(Vector<IntVectSet>& a_tags)
{

  if (s_verbosity > 3)
  {
    pout() << "amrMushyLayer::tagCells" << endl;
  }


  int top_level = a_tags.size();
  top_level = min(top_level-1, m_finest_level);
  // loop over levels
  for (int lev=0; lev<=top_level; lev++)
  {
    IntVectSet& levelTags = a_tags[lev];
    tagCellsLevel(levelTags, lev);
  }

  if (s_verbosity > 5)
  {
    pout() << "amrMushyLayer::tagCells finished" << endl;
  }

}

void
amrMushyLayer::tagCellsLevel(IntVectSet& a_tags, int a_level)
{

  if (s_verbosity > 4)
  {
    pout() << "amrMushyLayer::tagCellsLevel " << a_level << endl;
  }

  {
    IntVectSet local_tags;

    if (m_tagging_vector_var > -1 and m_tagging_scalar_var < 0)
    {
      DataIterator dit = m_vectorNew[m_tagging_vector_var][a_level]->dataIterator();
      LevelData<FArrayBox>& levelPhi = *m_vectorNew[m_tagging_vector_var][a_level];
      const DisjointBoxLayout& levelGrids = m_amrGrids[a_level];

      // need to ensure that ghost cells are set properly
      levelPhi.exchange(levelPhi.interval());

      for (dit.begin(); dit.ok(); ++dit)
      {

        BoxIterator bit(levelGrids[dit()]);
        for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          for (int idir=0; idir<SpaceDim; idir++)
          {
            // Tag if this component of the vector field is greater than the threshold value * (level+1)
            if (abs(levelPhi[dit](iv, idir)) > m_tagging_val*(a_level+1))
            {
              local_tags |= iv;
            }
          }
        } // end loop over cells
      }

    }
    else if (m_tagging_scalar_var > -1)
    {
      DataIterator dit = m_scalarNew[m_tagging_scalar_var][a_level]->dataIterator();
      LevelData<FArrayBox>& levelPhi = *m_scalarNew[m_tagging_scalar_var][a_level];
      const DisjointBoxLayout& levelGrids = m_amrGrids[a_level];

      // need to ensure that ghost cells are set properly
      levelPhi.exchange(levelPhi.interval());

      for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox gradPhi(levelGrids[dit()], 1);

        if (s_verbosity > 5)
        {
          pout() << "Tagging cells on box " << levelGrids[dit()].smallEnd() << "->" << levelGrids[dit()].bigEnd() << endl;
        }

        for (int dir=0; dir<SpaceDim; dir++)
        {
          const Box b = levelGrids[dit()];
          const Box bcenter = b & grow ( m_amrDomains[a_level], -BASISV(dir) );
          const Box blo = b & adjCellLo( bcenter, dir );
          const Box bhi = b & adjCellHi( bcenter, dir );
          const int haslo = ! blo.isEmpty();
          const int hashi = ! bhi.isEmpty();
          FORT_UNDIVIDEDGRAD ( CHF_FRA(gradPhi),
                               CHF_CONST_FRA(levelPhi[dit()]),
                               CHF_BOX(bcenter),
                               CHF_BOX(blo),
                               CHF_BOX(bhi),
                               CHF_CONST_INT(dir),
                               CHF_CONST_INT(haslo),
                               CHF_CONST_INT(hashi));


          // now tag cells based on values
          BoxIterator bit(levelGrids[dit()]);
          for (bit.begin(); bit.ok(); ++bit)
          {
            const IntVect& iv = bit();
            if (abs(gradPhi(iv)) > m_tagging_val)
              local_tags |= iv;
          } // end loop over cells
        } // end loop over directions

      } // end loop over grids
    }
    else
    {
      // No refinement
    }

    // now buffer tags
    local_tags.grow(m_tags_grow);
    local_tags &= m_amrDomains[a_level];

    // gather all tags into one place

    Vector<IntVectSet> all_tags;

    const int dest_proc = uniqueProc (SerialTask::compute);

    gather(all_tags, local_tags, dest_proc);
    if (procID()==uniqueProc (SerialTask::compute)) {
      for (int i=0; i!=all_tags.size(); ++i)
        a_tags |= all_tags[i];
    }

  }
}

void
amrMushyLayer::tagCellsInit(Vector<IntVectSet>& a_tags)
{

  if (s_verbosity > 3)
  {
    pout() << "amrMushyLayer::tagCellsInit" << endl;
  }


  // for initial time, refine entire domain for all levels

  if (m_refine_initial_domain)
  {
    for (int lev=0; lev<m_max_level; lev++)
    {
      a_tags[lev].define(m_amrDomains[lev].domainBox());
    }
  }
  else
  {
    tagCells(a_tags);
  }

}


// compute timestep
void
amrMushyLayer::computeDt()
{
  logMessage(10, "amrMushyLayer::computeDt");

  // Ignoring advection and flow for now
//  setupAdvectionSolvers(); // We need these to get dt
//  Real maxFluidVel = m_godunovEnthalpy[0]->getMaxWaveSpeed(*m_scalarNew[ScalarVars::m_enthalpy][0], *m_fluidAdv[0]);
  Real maxFluidVel  = 0.0;

  //Get the dt for the finest level based on frame advection
  Real finest_frameAdv_dt = m_amrDx[m_finest_level] /  m_parameters.nonDimVel;

  //Just do this on level 0 (base level)

  Real finest_fluidAdv_dt = m_amrDx[m_finest_level] / maxFluidVel;

  Real finest_dt = min(finest_frameAdv_dt, finest_fluidAdv_dt);

  // Multiply by initial CFL and grow for each timestep we've done
  Real cfl = min(m_initial_cfl * pow(m_max_dt_grow, m_cur_step), m_cfl);
  Real grown_dt = finest_dt * cfl;

  //Check we haven't got a dt greater than allowed by our max CFL number


  //Finally check we haven't exceeded our dt cap - this only really matters when fluid velocities are tiny or non-existant
  m_dt = min(grown_dt, m_max_dt*cfl);


  // This is so we can fix the same dt for runs on different grids, for
  // comparison of the time dependent solution
  ParmParse pp("main");
  Real fixed_dt = -1;
  pp.query("fixed_dt", fixed_dt);
  if (fixed_dt > 0)
  {
    m_dt = fixed_dt;
  }




}

void
amrMushyLayer::computeInitialDt()
{

  if (s_verbosity > 3)
  {
    pout() << "amrMushyLayer::computeInitialDt" << endl;
  }

  // for now, just call computeDt;
  computeDt();

}

void
amrMushyLayer::convergedToSteadyState()
{
  ParmParse pp("main");

  //following Katz & Worster:
  // if max(|theta^(n+1) - theta^n|)/(delta t) < 1e-3 then stop
  Real tolerance = 1e-3;

  pp.query("convergence_criteria", tolerance);

  Vector<LevelData<FArrayBox>* > dthetadt, dThetadt;
  dthetadt.resize(m_finest_level+1, NULL);
  dThetadt.resize(m_finest_level+1, NULL);

  for (int lev=0; lev <=m_finest_level; lev++)
  {
    dthetadt[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_ghost*IntVect::Zero);
    dThetadt[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_ghost*IntVect::Zero);

    for (DataIterator dit = m_scalarNew[m_theta][lev]->dataIterator(); dit.ok(); ++dit)
    {
      (*dthetadt[lev])[dit()].setVal(0);
      (*dthetadt[lev])[dit()] += (*m_scalarNew[m_theta][lev])[dit()];
      (*dthetadt[lev])[dit()] -= (*m_scalarOld[m_theta][lev])[dit()];
      (*dthetadt[lev])[dit()] /= m_dt;

      (*dThetadt[lev])[dit()].setVal(0);
      (*dThetadt[lev])[dit()] += (*m_scalarNew[ScalarVars::m_bulkConcentration][lev])[dit()];
      (*dThetadt[lev])[dit()] -= (*m_scalarOld[ScalarVars::m_bulkConcentration][lev])[dit()];
      (*dThetadt[lev])[dit()] /= m_dt;
    }
  }

  Interval comps(0,0);
  Real thetaMaxDiff = computeNorm(dthetadt, m_refinement_ratios, m_amrDx[0], comps, 0);
  Real ThetaMaxDiff = computeNorm(dThetadt, m_refinement_ratios, m_amrDx[0], comps, 0);

  pout() << "(delta T)/(delta t) = " << thetaMaxDiff << ", (delta Theta)/(delta t) = " << ThetaMaxDiff <<endl;
  pout() << "Needs to be < " << tolerance << " for convergence" << endl;

  if (thetaMaxDiff < tolerance and
      ThetaMaxDiff < tolerance)
  {
    pout() << "Reached steady state!" << endl;
    m_steady_state = true;
  }
  else
  {
    // This is important - we may enforce a minimum number of timesteps, and if a 'steady state' was reached during them but
    // this is no longer the case at some later time, we want to keep running the model
    m_steady_state = false;
  }

  //Clean up memory
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    if (dThetadt[lev] != NULL)
    {
      delete dThetadt[lev];
      dThetadt[lev] = NULL;
    }

    if (dthetadt[lev] != NULL)
    {
      delete dthetadt[lev];
      dthetadt[lev] = NULL;
    }
  }

}

void amrMushyLayer::
ptrToRefCountedPtr(Vector<LevelData<FArrayBox>* >& ptr,
                   Vector<RefCountedPtr<LevelData<FArrayBox> > >& refPtr)
{
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    refPtr[lev] = RefCountedPtr<LevelData<FArrayBox> >(ptr[lev]);
  }
}

void amrMushyLayer::
refCountedPtrToPtr(Vector<RefCountedPtr<LevelData<FArrayBox> > >& refPtr,
                   Vector<LevelData<FArrayBox>* >& ptr)
{
  for (int lev=0; lev<=m_finest_level; lev++)
  {

    //		LevelData<FArrayBox>& levDat = *refPtr[lev];
    ptr[lev] = &(*refPtr[lev]);

  }
}

Vector<LevelData<FArrayBox>*> amrMushyLayer::
timeCenteredScalar(const int a_var)
{
  Vector<LevelData<FArrayBox>*> av(m_finest_level+1,NULL);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    const DisjointBoxLayout& grids = m_amrGrids[lev];
    av[lev] = new LevelData<FArrayBox>(grids, 1,m_ghostVect);

    DataIterator dit = av[lev]->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      (*av[lev])[dit].setVal(0);
      (*av[lev])[dit] += (*m_scalarNew[a_var][lev])[dit];
      (*av[lev])[dit] += (*m_scalarOld[a_var][lev])[dit];
      (*av[lev])[dit].divide(2);
    }
  }

  return av;
}


int amrMushyLayer::getRefinementRatio(int coarseLevel)
{
  int refRatio;

  if(coarseLevel >= 0 and coarseLevel < m_max_level)
  {
    refRatio = m_refinement_ratios[coarseLevel];
  }
  else
  {
    //Default value
    refRatio = 2;
  }

  return refRatio;
}

void amrMushyLayer::interpToFineSolutions(int lev,
                                          const std::vector<int> scalarVars,
                                          const std::vector<int> vectorVars)
{
  bool averageFromDest = true;
  //For consistency, get level > 0 vars by interpolation
  FineInterp interpScalar(m_amrGrids[lev], 1, m_refinement_ratios[lev-1], m_amrDomains[lev]);
  FineInterp interpVector(m_amrGrids[lev], SpaceDim, m_refinement_ratios[lev-1], m_amrDomains[lev]);

  if (scalarVars.size() > 0)
  {
    // If we've specified the vars to interpolate, do them

    for (int var_i = 0; var_i < scalarVars.size(); var_i++)
    {
      int a_var = scalarVars[var_i];
      interpScalar.interpToFine(*m_scalarNew[a_var][lev],*m_scalarNew[a_var][lev-1], averageFromDest);
      interpScalar.interpToFine(*m_scalarOld[a_var][lev],*m_scalarOld[a_var][lev-1], averageFromDest);
    }
  }
  else
  {
    //Otherwise, interpolate all variables

    for (int a_var = 0; a_var < m_numVars; a_var++)
    {
      interpScalar.interpToFine(*m_scalarNew[a_var][lev],*m_scalarNew[a_var][lev-1], averageFromDest);
      interpScalar.interpToFine(*m_scalarOld[a_var][lev],*m_scalarOld[a_var][lev-1], averageFromDest);
    }
  }

  if (vectorVars.size() > 0)
  {
    for (int var_i = 0; var_i < vectorVars.size(); var_i++)
    {
      int a_var = vectorVars[var_i];
      interpVector.interpToFine(*m_vectorNew[a_var][lev],*m_vectorNew[a_var][lev-1], averageFromDest);
      interpVector.interpToFine(*m_vectorOld[a_var][lev],*m_vectorOld[a_var][lev-1], averageFromDest);

    }
  }
  else
  {
    for (int a_var = 0; a_var < m_numVectorVars; a_var++)
    {
      interpVector.interpToFine(*m_vectorNew[a_var][lev],*m_vectorNew[a_var][lev-1], averageFromDest);
      interpVector.interpToFine(*m_vectorOld[a_var][lev],*m_vectorOld[a_var][lev-1], averageFromDest);
    }
  }

}
