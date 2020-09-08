#include "AMRLevelMushyLayer.H"

void AMRLevelMushyLayer::setSmoothingCoeff(Real a_coeff)
{
  // Need this so the setup_new_run program can set this coefficient

  s_regrid_smoothing_coeff = a_coeff;
}

void AMRLevelMushyLayer::setSmoothingDone(bool a_smoothingDone)
{
  m_regrid_smoothing_done = a_smoothingDone;
}


void AMRLevelMushyLayer::doPostRegridSmoothing(bool a_smoothVel, bool a_smoothScalar)
{
  // first compute Laplacians of diffused quantities
  // for later use in smoothing ops.
  // make use of the fact that all levels should be
  // at the same time, so that we can store this info
  // in the old-data registers for later use

  // We only hit this function on level 0, so can just assume that smoothing hasn't yet been done
  m_regrid_smoothing_done = false;

  if (!m_regrid_smoothing_done)
  {
    // since we need a composite laplacian here _and_
    // we need to do averaging down, do this all at
    // once at beginning of regrid cycle for velocity
    // and for all diffused scalars (anything that's hit
    // with a diffusion operator).
    // for now, do this the easy way using AMRSolvers (ouch!)
    // and doing composite ops.  In all cases, we store
    // the Laplacian in the old-data holder (since we
    // shouldn't need it until the next timestep -- all
    // levels should be at the same time)

    AMRLevelMushyLayer* thisMLPtr = this;
    while (!thisMLPtr->finestLevel())
    {
      thisMLPtr = thisMLPtr->getFinerLevel();
    }
    CH_assert (thisMLPtr->finestLevel());
    int finest_level = thisMLPtr->m_level;

    // this looks a lot like what we do in post-timestep
    thisMLPtr = this;
    // if this level > 0, we will actually need up to two
    // more coarser levels for C/F BCs (recall that we need to
    // apply the operator to the level from which we will be
    // interpolating data, so we need one level coarser than
    // that one for a C/F bc.
    int startLev = m_level;
    // this moves startLev to the level on which we need to apply
    // the operator
    if (startLev > 0)
    {
      startLev = startLev -1;
      thisMLPtr = thisMLPtr->getCoarserLevel();
    }
    // this moves startLev to the level from which we get C/F BCs for
    // that level
    if (startLev > 0)
    {
      startLev = startLev - 1;
      thisMLPtr = thisMLPtr->getCoarserLevel();
    }

    // now set up multilevel stuff
    Vector<LevelData<FArrayBox>* > amrS(finest_level+1,nullptr);
    Vector<LevelData<FArrayBox>* > amrLapS(finest_level+1);
    Vector<DisjointBoxLayout> amrGrids(finest_level+1);
    Vector<int> amrRefRatios(finest_level+1,0);
    Vector<Real> amrDx(finest_level+1,0);
    Vector<ProblemDomain> amrDomains(finest_level+1);
    // will use this to do averaging down
    Vector<CoarseAverage*> amrAvgDown(finest_level+1, nullptr);

    // loop over levels, allocate temp storage for velocities,
    // set up for amrsolves

    for (int lev=startLev; lev<=finest_level; lev++)
    {
      const DisjointBoxLayout& levelGrids = thisMLPtr->m_grids;
      // since AMRSolver can only handle one component at a
      // time, need to allocate temp space to copy stuff
      // into, and then back out of to compute laplacian
      IntVect ghostVect(D_DECL(1,1,1));
      amrS[lev] = new LevelData<FArrayBox>(levelGrids, 1,
                                           ghostVect);
      amrLapS[lev] = new LevelData<FArrayBox>(levelGrids,1,
                                              ghostVect);

      amrGrids[lev] = levelGrids;
      amrRefRatios[lev] = thisMLPtr->refRatio();
      amrDx[lev] = thisMLPtr->m_dx;
      amrDomains[lev] = thisMLPtr->problemDomain();
      thisMLPtr = thisMLPtr->getFinerLevel();
      if (lev>= m_level && lev > 0)
      {
        amrAvgDown[lev] = new CoarseAverage(levelGrids,1,
                                            amrRefRatios[lev-1]);
      }
    } // end loop over levels

    AMRPoissonOpFactory localPoissonOpFactory;
    // allocate solver -- do this in a separate function in
    // order to make sure the two steps which involve
    // elliptic operators are consistent.

    // reset startLevel if necessary -- at this point, startLev
    // will be the coarsest level on which we apply the operator.
    if (m_level > 0) startLev=m_level-1;

    // note that we use startLev instead of m_level
    // as lBase for the solver because we will need
    // to compute operator on startLev
    defineRegridAMROp(localPoissonOpFactory, amrGrids,
                      amrDomains, amrDx,
                      amrRefRatios, startLev);

    int numLevs = finest_level - startLev + 1;
    Vector<AMRPoissonOp*> localPoissonOpPtrVec(numLevs);
    for (int lev=startLev; lev<=finest_level; lev++)
    {
      // virtual  AMRLevelOp<LevelData< FArrayBox> >*
      localPoissonOpPtrVec[lev - startLev] =
          (AMRPoissonOp*)
          localPoissonOpFactory.AMRnewOp(amrDomains[lev]);
    }


    if (a_smoothVel)
    {

      // now loop over velocity components
      for (int dir=0; dir<SpaceDim; ++dir)
      {
        Interval velComps(dir,dir);
        Interval tempComps(0,0);
        // for each velocity component, copy all levels
        // of data into temp storage.  then, loop over
        // all levels and compute laplacian of this comonent.
        // finally, copy lap(vel) into old_vel storage space.
        thisMLPtr = this;
        if (startLev < m_level) thisMLPtr = thisMLPtr->getCoarserLevel();

        for (int lev=startLev; lev<=finest_level; lev++)
        {
          thisMLPtr->m_vectorNew[VectorVars::m_fluidVel]->copyTo(velComps, *amrS[lev],
                                                     tempComps);
          thisMLPtr = thisMLPtr->getFinerLevel();
        }

        // now loop over levels and apply multilevel operator.
        // also do averaging down here if necessary

        for (int lev=startLev; lev<=finest_level; lev++)
        {
          LevelData<FArrayBox>& levelS = *amrS[lev];
          LevelData<FArrayBox>& levelLapS = *amrLapS[lev];
          int indVec = lev - startLev;
          if (lev == 0)
          { // no coarser level
            if (lev == finest_level)
            { // no finer level:  this is all there is
              localPoissonOpPtrVec[indVec]->
              applyOp(levelLapS, levelS);
            }
            else
            { // finer level exists
              localPoissonOpPtrVec[indVec]->
              AMROperatorNC(levelLapS,
                            *amrS[lev + 1],
                            levelS, false,
                            localPoissonOpPtrVec[indVec+1]);
            }
          }
          else if (lev < finest_level - 1)
          { // three-level operator
            localPoissonOpPtrVec[indVec]->
            AMROperator(levelLapS,
                        *amrS[lev + 1],
                        levelS,
                        *amrS[lev - 1], false,
                        localPoissonOpPtrVec[indVec+1]);
          }
          else
          { // no finer level
            localPoissonOpPtrVec[indVec]->
            AMROperatorNF(levelLapS,
                          levelS,
                          *amrS[lev - 1], false);
          }
        }

        // need to do average down from finest level on down
        // to speed this up eventually, may want to do this
        // _after_ we copy into multicomponent velocity LDF
        for (int lev= finest_level; lev>startLev; lev--)
        {
          amrAvgDown[lev]->averageToCoarse(*amrLapS[lev-1],
                                           *amrLapS[lev]);
        }

        // finally, loop over levels and copy to old-vel storage
        thisMLPtr = this;
        if (m_level > 0) thisMLPtr = thisMLPtr->getCoarserLevel();

        for (int lev=startLev; lev<=finest_level; lev++)
        {
          amrLapS[lev]->copyTo(tempComps,*(thisMLPtr->m_vectorNew[VectorVars::m_fluidVel]),
                               velComps);
          thisMLPtr = thisMLPtr->getFinerLevel();
        }

      } // end loop over velocity components

    } // end if smoothing velocity

    // since we're done with velocity, can clear up storage space
    for (int lev=0; lev<=finest_level; lev++)
    {
      if (amrS[lev] != nullptr)
      {
        delete amrS[lev];
        amrS[lev] = nullptr;
      }
      if (amrLapS[lev] != nullptr)
      {
        delete amrLapS[lev];
        amrLapS[lev] = nullptr;
      }
    }

    // now do scalars -- since (at the moment),
    // we only have single-component scalars, can
    // avoid copies and work directly with the scalars
    // themselves
    if (a_smoothScalar)
    {


      for (int scalComp=0; scalComp < m_numScalarVars; scalComp++)
      {
        // only do all this if scalar is diffused
        if (m_scalarDiffusionCoeffs[scalComp] > 0)
        {
          thisMLPtr = this;
          if (startLev < m_level)
          {
            thisMLPtr = thisMLPtr->getCoarserLevel();
            // may need even coarser level for C/F BCs
            if (startLev > 0)
            {
              startLev = startLev -1;
              thisMLPtr = thisMLPtr->getCoarserLevel();
            }
          }

          for (int lev=startLev; lev<=finest_level; lev++)
          {
            amrS[lev] = thisMLPtr->m_scalarNew[scalComp];
            //                                 amrLapS[lev] = thisMLPtr->m_scalarOld[scalComp];
            amrLapS[lev] = new LevelData<FArrayBox> (thisMLPtr->m_grids, 1, IntVect::Unit);
            thisMLPtr = thisMLPtr->getFinerLevel();
          }

          if (m_level > 0) startLev = m_level-1;

          // now compute laplacian and average down if necessary
          for (int lev=startLev; lev<=finest_level; lev++)
          {
            LevelData<FArrayBox>& levelS = *amrS[lev];
            LevelData<FArrayBox>& levelLapS = *amrLapS[lev];
            int indVec = lev - startLev;
            if (lev == 0)
            { // no coarser level
              if (lev == finest_level)
              { // no finer level:  this is all there is
                localPoissonOpPtrVec[indVec]->
                applyOp(levelLapS, levelS);
              }
              else
              { // finer level exists
                localPoissonOpPtrVec[indVec]->
                AMROperatorNC(levelLapS,
                              *amrS[lev + 1],
                              levelS, false,
                              localPoissonOpPtrVec[indVec+1]);
              }


            }
            else if (lev < finest_level - 1)
            { // three-level operator
              localPoissonOpPtrVec[indVec]->
              AMROperator(levelLapS,
                          *amrS[lev + 1],
                          levelS,
                          *amrS[lev - 1], false,
                          localPoissonOpPtrVec[indVec+1]);
            }
            else
            { // no finer level
              localPoissonOpPtrVec[indVec]->
              AMROperatorNF(levelLapS,
                            levelS,
                            *amrS[lev - 1], false);
            }
            if (lev > m_level)
            {
              amrAvgDown[lev]->averageToCoarse(*amrLapS[lev-1],
                                               *amrLapS[lev]);
            }
          } // end loop over levels


          for (int lev=startLev; lev <= finest_level; ++lev)
          {
            LevelData<FArrayBox>& levelS = *amrS[lev]; //NB: amrS points to scalarNew
            LevelData<FArrayBox>& levelLapS = *amrLapS[lev]; // NB: amrLapS points to scalarOld
            DataIterator levelDit = levelS.dataIterator();
            for (levelDit.reset(); levelDit.ok(); ++levelDit)
            {
              // levelLapS contains (1 - delta t * smooth factor * Lap) S^{n+1}
              // subtract off solution
              levelLapS[levelDit] -= levelS[levelDit] ;



              // Scale by diffusion/viscosity (was previously scaled by viscosity coefficient)
              levelLapS[levelDit] *= m_scalarDiffusionCoeffs[scalComp];
              if (m_parameters.m_viscosityCoeff > 0.0)
              {
                levelLapS[levelDit] *= 1.0/m_parameters.m_viscosityCoeff;
              }

              // Add to solution
              levelS[levelDit] += levelLapS[levelDit];


              // Basic smoothing
//              Box b = levelS[levelDit].box();
//              b.grow(-1);
//              for (BoxIterator bit = BoxIterator(b); bit.ok(); bit.begin())
//              {
//                IntVect iv = bit();
//                levelS[levelDit](iv) = levelS[levelDit](iv) + this->s_regrid_smoothing_coeff* (levelS[levelDit](iv + BASISV(0))   + levelS[levelDit](iv - BASISV(0)));
//              }

            }

//                      int temp=0;


          } // end loop over levels
        } // end if scalar is diffused
      } // end loop over scalar components

    } // end if smoothing scalars

    // finally, clean up storage for coarseAverages
    for (int lev=0; lev<=finest_level; ++lev)
    {
      if (amrAvgDown[lev] != nullptr)
      {
        delete amrAvgDown[lev];
        amrAvgDown[lev] = nullptr;
      }


    }
    for (int lev=startLev; lev<=finest_level; lev++)
    {
      delete localPoissonOpPtrVec[lev - startLev];
    }
    thisMLPtr = this;
    if (startLev < m_level) thisMLPtr = thisMLPtr->getCoarserLevel();
    for (int lev=startLev; lev <= finest_level; ++lev)
    {
      thisMLPtr->m_regrid_smoothing_done = true;
      thisMLPtr = thisMLPtr->getFinerLevel();
    }
  } // end if post-regridding smoothing hasn't been done yet


  CH_assert(m_level==0);

  //  this->smoothVelocityField(0);
}

void
AMRLevelMushyLayer::defineRegridAMROp(AMRPoissonOpFactory& a_factory,
                                      const Vector<DisjointBoxLayout>& a_grids,
                                      const Vector<ProblemDomain>& a_domains,
                                      const Vector<Real>& a_amrDx,
                                      const Vector<int>& a_refRatios,
                                      const int& a_lBase)
{
  // want to use dt from lbase
  Real dtLBase = m_dt;
  if (a_lBase > m_level)
  {
    AMRLevelMushyLayer* thisMLPtr = getFinerLevel();
    while (a_lBase > thisMLPtr->m_level)
    {
      thisMLPtr = thisMLPtr->getFinerLevel();
    }
    dtLBase = thisMLPtr->m_dt;
  }
  else if (a_lBase < m_level)
  {
    AMRLevelMushyLayer* thisMLPtr = getCoarserLevel();
    while (a_lBase < thisMLPtr->m_level)
    {
      thisMLPtr = thisMLPtr->getCoarserLevel();
    }
    dtLBase = thisMLPtr->m_dt;
  }

  // define coefficient
  Real nu = m_parameters.m_viscosityCoeff;  //m_parameters.prandtl;
  Real mu = -m_opt.regrid_smoothing_coeff*dtLBase;
  if (nu > 0.0)
  {
    mu = mu*nu;
  }


  // Would like to use extrap BC's, since they're probably the safest

  a_factory.define(a_domains[0],
                   a_grids,
                   a_refRatios,
                   a_amrDx[0],
                   m_physBCPtr->BasicGradPressureFuncBC(), // this is just extrap BCs
                   1., // alpha
                   mu); // beta
}


void AMRLevelMushyLayer::tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::tagCells " << m_level << endl;
  }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.

  // Create tags based on criteria
  IntVectSet localTags;

  if (m_opt.min_regrid_time >= 0 && m_time < m_opt.min_regrid_time)
  {
    a_tags = localTags;
    pout() << "Not tagging any cells as current time (" << m_time << ") is less than " << m_opt.min_regrid_time << endl;
    return;
  }

  AMRLevelMushyLayer* ml = this;
  while(ml->getFinerLevel())
  {
    ml = ml->getFinerLevel();
  }
  int finestLevel = ml->m_level;


  // Compute liquid, mushy and inbetween cells - these will be useful
  IntVectSet mushyCells;
  IntVectSet shrunkMushyCells;

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    for (BoxIterator bit(m_grids[dit]); bit.ok(); ++bit)
    {
      Real porosity = (*m_scalarNew[ScalarVars::m_porosity])[dit](bit());
      if (porosity < m_opt.taggingMarginalPorosityLimit
          && porosity > 0.0)
      {
        mushyCells |= bit();
      }
    }
  }

  shrunkMushyCells = mushyCells;
  IntVectSet shiftedTags = shrunkMushyCells;
  IntVect shiftVect = IntVect::Zero;
  shiftVect[SpaceDim-1] = m_opt.porousCellsShrink; // distance to shift
  shiftedTags.shift(shiftVect);
  shrunkMushyCells &= shiftedTags;

  if (!shrunkMushyCells.isEmpty())
  {
    // Recover original cells
    shrunkMushyCells.shift(-shiftVect);
    for (int i=0; i < m_opt.porousCellsShrink-1; i++)
    {
      shrunkMushyCells.growHi(SpaceDim-1);

      // Grow horizontally to get around some edge effects
      for (int horizontal_dir = 0; horizontal_dir < SpaceDim-1; horizontal_dir++)
      {
        shrunkMushyCells.growHi(horizontal_dir);
      }
    }
  }

  if (m_opt.refinementMethod == RefinementMethod::tagSpeed) // m_opt.tag_velocity
  {

    if (s_verbosity >= 2)
    {
      pout() << "AMRLevelMushyLayer::tagCells - tag fluid speed >  " << m_opt.vel_thresh << endl;
    }

    tagCellsVar(localTags, m_opt.vel_thresh, -1, m_fluidVel, TaggingMethod::Magnitude);
  }
  else if (m_opt.refinementMethod == RefinementMethod::tagMushChannelsCompositeCriteria)
  {
    if (s_verbosity >= 2)
       {
         pout() << "AMRLevelMushyLayer::tagCells - using composite criteria to find mush and channels >  " << endl;
       }

    DataIterator dit = m_grids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& b = m_grids[dit()];

          FArrayBox& porosity = (*m_scalarNew[m_porosity])[dit];
          FArrayBox& velocity = (*m_vectorNew[m_fluidVel])[dit];
          FArrayBox& bulkC = (*m_scalarNew[m_bulkConcentration])[dit];

          FArrayBox taggingMetricFab(b, 1);
          FArrayBox gradFab(b, SpaceDim);

          // empirical parameter weighting for channels
          // want to avoid
          Real alpha=0.0;

          // Calculated undivided gradient
          for (int dir = 0; dir < SpaceDim; ++dir)
          {
            const Box bCenter = b & grow(m_problem_domain, -BASISV(dir));
            const Box bLo = b & adjCellLo(bCenter, dir);
            const int hasLo = !bLo.isEmpty();
            const Box bHi = b & adjCellHi(bCenter, dir);
            const int hasHi = !bHi.isEmpty();
            FORT_GETGRADF(CHF_FRA1(gradFab, dir), CHF_CONST_FRA1(porosity, 0),
                          CHF_CONST_INT(dir), CHF_BOX(bLo), CHF_CONST_INT(hasLo),
                          CHF_BOX(bHi), CHF_CONST_INT(hasHi), CHF_BOX(bCenter));
          }

          FORT_MAGNITUDEF(CHF_FRA1(taggingMetricFab, 0), CHF_CONST_FRA(gradFab),
                          CHF_BOX(b));


          for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
          {
            IntVect iv = bit();
            RealVect loc;
            ::getLocation(iv, loc, m_dx);

            if (taggingMetricFab(iv)*(1 - alpha*min(velocity(iv, SpaceDim-1) * (m_parameters.compositionRatio-bulkC(iv)), 0.0) ) > m_opt.refineThresh)
            {
              localTags |= iv;
            }
          }
        }

  }
  else if (m_opt.refinementMethod == RefinementMethod::tagChannels) // m_opt.tag_channels
  {
    Real porosity_limit = 0.5;
    Real vel_limit = 0.5;
    Real porosity_criteria = 0.9; // how low the horizontally averaged porosity must be for us to decide we're in a mushy layer
    Real y_limit =  m_domainHeight - computeMushDepth(porosity_criteria);

    DataIterator dit = m_grids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& b = m_grids[dit()];

      FArrayBox& porosity = (*m_scalarNew[m_porosity])[dit];
      FArrayBox& velocity = (*m_vectorNew[m_fluidVel])[dit];

      for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        ::getLocation(iv, loc, m_dx);

        if (loc[1] > y_limit
            && (porosity(iv) > porosity_limit
                || abs(velocity(iv, 1)) > vel_limit ))
        {
          localTags |= iv;
        }
      }
    }

    // Shrink the regrow tags invertical direction to get rid of individual odd cells
    //    int grow_dist = 5;
    //    localTags.grow(1, -grow_dist);
    //    localTags.grow(1, grow_dist);

    int vertical_dir = SpaceDim-1;

    // Only keep tags which overlap with the same tags shifted up by shift_dist cells
    int shift_dist = 5;
    IntVectSet shiftedTags = localTags;

    // Shift in vertical direction (which has an index given by number of spatial dimensions - 1) by shift_dist
    IntVect shiftDir = IntVect::Zero;
    shiftDir[vertical_dir] = shift_dist;
    shiftedTags.shift(shiftDir);

    localTags &= shiftedTags;

    localTags.grow(vertical_dir, shift_dist);

    // Grow further downwards
    // Can only do this by shifting down then growing in both directions
    int growLo = 4;

    // Create new shift vector
    shiftDir = IntVect::Zero;
    shiftDir[vertical_dir] = -growLo;
    localTags.shift(shiftDir);
    localTags.grow(vertical_dir, growLo);

  }
  else if (m_opt.refinementMethod == RefinementMethod::tagPlumeMush
      || m_opt.refinementMethod == RefinementMethod::tagMushChannels)  // m_opt.tag_plume_mush
  {
    // Trying to refine plumes

    if (s_verbosity >= 2)
    {
      if (m_opt.refinementMethod == RefinementMethod::tagPlumeMush)
      {
        pout() << "AMRLevelMushyLayer::tagCells - refine plume mush - " << m_level << endl;
      }
      else
      {
        pout() << "AMRLevelMushyLayer::tagCells - refine mush channels - " << m_level << endl;
      }
    }


    // Place finest resolution around channels
    if (m_level == finestLevel - 1)
    {

      if (m_opt.refinementMethod == RefinementMethod::tagMushChannels)
      {

        // Create new field, bulk concentration * vertical velocity
        LevelData<FArrayBox> concentrationVelocity(m_grids, 1, IntVect::Zero);
        m_scalarNew[m_bulkConcentration]->copyTo(concentrationVelocity);

        for (DataIterator dit= m_grids.dataIterator(); dit.ok(); ++dit)
        {
          concentrationVelocity[dit].plus(1.0);
          concentrationVelocity[dit].mult((*m_vectorNew[m_fluidVel])[dit], SpaceDim-1, 0, 1);

          Box b = concentrationVelocity[dit].box();

          BoxIterator bit(b);
          for (bit.begin(); bit.ok(); ++bit)
          {
            const IntVect& iv = bit();

            if (concentrationVelocity[dit](iv) < -m_opt.refineThresh)
            {
              localTags |= iv;
            }
          }
        }

      }
      else
      {


        // Refine cells that are already refined + have large salinity gradients

        IntVectSet downflowCells;


        IntVectSet porosityGradientCells;
        tagCellsVar(porosityGradientCells, m_opt.refineThresh, m_porosity, -1, m_opt.taggingMethod);

        // Grow these cells to go further into the channel,
        // where the salinity is high and fluid velocities are large and negative
        porosityGradientCells.grow(2*m_opt.tagBufferSize);

        // Tag regions of downflow
        tagCellsVar(downflowCells,
                    m_opt.plumeVelThreshold, // refine thresh
                    -1, // don't consider a scalar field
                    m_fluidVel, TaggingMethod::CompareLargerThanNegative, 1); // y component of velocity

        IntVectSet saltyCells;
        tagCellsVar(saltyCells, m_opt.plumeSalinityThreshold,
                    ScalarVars::m_bulkConcentration,
                    -1, TaggingMethod::CompareLargerThan);

        // Only take the cells which have the porosity gradient
        // and either downflow or very salty

        // Combine downflow and salty
        downflowCells |= saltyCells;

        // Only take downflow cells which are also in the set of cells satisfying the porosity gradient condition
        downflowCells &= porosityGradientCells;

        if (s_verbosity >= 5)
        {
          pout() << "Downflow & porous gradient cells: " << downflowCells << endl;
        }

        // Add to the local tag set
        localTags |= downflowCells;
      }


    }
    else
    {
      // Make sure we tag mushy regions on level 0
      if (m_level == 0)
      {
        if (s_verbosity >= 5)
        {
          pout() << "Refining on level " << m_level << " where porosity < 1 (all mushy cells)" << endl;
        }

        localTags |= mushyCells;
      }
      else
      {


        if (s_verbosity >= 5)
        {
          pout() << "Refining on level " << m_level << " where porosity gradients > refine thresh (" << m_opt.refineThresh << ")" << endl;
        }

        // Refine mushy regions which may be about generate channels
        IntVectSet mushyCells;
        tagCellsVar(mushyCells,
                    m_opt.refineThresh, // refine thresh for gradients
                    m_porosity, -1,
                    TaggingMethod::UndividedGradient);

        localTags = mushyCells;


      }


    }

  }
  else if (m_opt.refinementMethod == RefinementMethod::tagScalar)
  {
    if (s_verbosity >= 2)
    {

      pout() << "AMRLevelMushyLayer::tagCells - refining on variable - " << m_scalarVarNames[m_opt.taggingVar] << endl;

    }

    tagCellsVar(localTags, m_opt.refineThresh, m_opt.taggingVar, -1, m_opt.taggingMethod);
  }
  else if (m_opt.refinementMethod == RefinementMethod::tagVector)
  {
    if (s_verbosity >= 2)
    {

      pout() << "AMRLevelMushyLayer::tagCells - refining on variable - " << m_vectorVarNames[m_opt.taggingVectorVar] << endl;

    }

    tagCellsVar(localTags, m_opt.refineThresh, -1, m_opt.taggingVectorVar, m_opt.taggingMethod);
  }
  else
  {
    if (s_verbosity >= 2)
    {
      pout() << "AMRLevelMushyLayer::tagCells - custom tagging criteria, tagging mushy cells " << endl;
    }
    localTags |= mushyCells;
  }

  if (m_opt.onlyTagPorousCells)
  {
    localTags &= shrunkMushyCells;
  }


  if (m_opt.tagMLboundary)
  {
    if (s_verbosity >= 2)
    {
      pout() << "AMRLevelMushyLayer::tagCells - also tag Mush-Liquid boundary " << endl;
    }



    // Tag liquid-mush boundary and add to existing tags
    tagMushLiquidBoundary(localTags);

    // Tags cells with porosity < 1
    //    tagMushyCells(localTags);
  }

  if (m_opt.tagDomainBoundary)
  {
    if (s_verbosity >= 2)
    {
      pout() << "AMRLevelMushyLayer::tagCells - also tag domain boundary " << endl;
    }

    tagBoundaryLayerCells(localTags);
  }

  if (m_opt.tagCenterBoxSize > 0)
  {
    if (s_verbosity >= 2)
    {
      pout() << "AMRLevelMushyLayer::tagCells - only tag the middle of the domain " << endl;
    }

    localTags = IntVectSet();

    if (m_time >= m_opt.tagCenterBoxRegridTime)
    {

      tagCenterCells(localTags, m_opt.tagCenterBoxSize);
    }
  }

  localTags.grow(m_opt.tagBufferSize);

  // Another option - when regridding, move to specified gridfile
  if (m_opt.fixed_grid_time >= 0)
  {
    // Define new intvect set
    localTags = IntVectSet();
    if (m_time > m_opt.fixed_grid_time)
    {
      if (s_verbosity >= 2)
      {
        pout() << "AMRLevelMushyLayer::tagCells - tag all cells covered by specified grid file " <<  endl;
      }

      // Get whole hierarchy of grids
      Vector<Vector<Box> > amrGrids;
      getFixedGrids(amrGrids,  m_problem_domain, "regrid_gridfile");

      // Now get the grids for the level finer than this one, to
      // decide what cells to tag on this level
      int numLevels = amrGrids.size();
      int levelForGrids = m_level+1;
      if (levelForGrids < numLevels)
      {


        Vector<Box> levelGrids = amrGrids[levelForGrids];

        for (int i=0; i < levelGrids.size(); i++)
        {
          localTags |= levelGrids[i];
        }

        // Now need to coarsen tags to current level
        localTags.coarsen(m_ref_ratio);
      }
    }
  }

  // Ensure tags are contained within problem domain
  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;


  // This is some code which will force any refined levels to coarsen again,
  // setting up a cycle of refine->coarsen->refine to artifically test how
  // the code handles such procedures (e.g. testing for memory leaks)
  if (m_opt.testRegridCoarsening)
  {
    if (s_verbosity >= 2)
    {
      pout() << "AMRLevelMushyLayer::tagCells - testing regrid coarsening" << endl;
    }

    bool hasValidFinerLevel = false;
    if (m_hasFiner)
    {
      AMRLevelMushyLayer* fine =getFinerLevel();
      if (fine  != nullptr)
      {
        Vector<Box> boxes = fine->boxes();
        if (boxes.size() > 0)
        {
          hasValidFinerLevel = true;
        }
      }
    }

    if (m_level>0 || hasValidFinerLevel)
    {
      IntVectSet blankTags;
      a_tags = blankTags;
      return;
    }
  }

  if (s_verbosity >= 3)
  {
    pout() << "Final local tags: " << localTags << endl;
    pout() << "  Num points: " << localTags.numPts() << endl;
    pout() << "  Min enclosing box: " << localTags.minBox() << endl;
  }

  a_tags = localTags;
}

void AMRLevelMushyLayer::tagMushyCells(IntVectSet& localTags)
{
  IntVectSet mushyCells;

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    Box b = m_grids[dit];

    FArrayBox& porosityFAB = (*m_scalarNew[ScalarVars::m_porosity])[dit];
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();


      if (porosityFAB(iv) < 1)
      {
        mushyCells |=iv;
      }
    }
  }

  localTags |= mushyCells;
}

void AMRLevelMushyLayer::tagBoundaryLayerCells(IntVectSet& localTags)
{
  IntVectSet boundaryCells;

  Box domBox = m_problem_domain.domainBox();
  Box interiorBox(domBox);
  interiorBox.grow(-1); // contains all interior cells (ones we don't want to refine on)

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    Box b = m_grids[dit];

    IntVectSet ivs(b);

    // remove all ivs in the interior of the domain, leaving just those which border the domain boundary
    ivs -= interiorBox;

    //b -= interiorBox;

    boundaryCells |= ivs;

  } // end loop over grids

  localTags |= boundaryCells;

}

void AMRLevelMushyLayer::tagMushLiquidBoundary(IntVectSet& localTags)
{
  IntVectSet boundaryCells;

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    Box b = m_grids[dit];
    b.grow(-1); // shrink by 1 in all directions so we can search adjacent cells

    FArrayBox& porosityFAB = (*m_scalarNew[ScalarVars::m_porosity])[dit];
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      bool nearbyCellIsLiquid = false;
      // Check adjacent cells to see if one of them is liquid
      for (int dir=0; dir<SpaceDim; dir++)
      {
        IntVect ivUp =   iv + BASISV(dir);
        IntVect ivDown = iv - BASISV(dir);

        if (porosityFAB(ivUp) == 1.0 || porosityFAB(ivDown) == 1.0)
        {
          nearbyCellIsLiquid = true;
        }
      }

      Real thisPorosity = porosityFAB(iv);
      if (thisPorosity < 1.0 && nearbyCellIsLiquid)
      {
        boundaryCells |= iv;
      }
    }
  }

  localTags |= boundaryCells;

}

void AMRLevelMushyLayer::tagCenterCells(IntVectSet& localTags, Real radius)
{
  IntVectSet centerCells;

  int shrinkRadius = round(m_numCells.min()*radius);

  Box centerBox = m_problem_domain.domainBox();
  centerBox.grow(-shrinkRadius);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    Box b = m_grids[dit];
    IntVectSet boxCells(b);

    b &= centerBox;

    centerCells |= b;
  }

  localTags |= centerCells;

}


void AMRLevelMushyLayer::tagCellsVar(IntVectSet& localTags, Real refineThresh,
                                     int taggingVar,
                                     int taggingVectorVar,
                                     TaggingMethod taggingMethod,int comp)
{
  const DisjointBoxLayout& levelDomain = m_scalarNew[0]->disjointBoxLayout();

  // If there is a coarser level, interpolate to fill undefined ghost cells
  if (m_hasCoarser)
  {
    const AMRLevelMushyLayer* amrGodCoarserPtr = getCoarserLevel();

    if (taggingVar >= 0)
    {

      PiecewiseLinearFillPatch pwl(levelDomain,
                                   amrGodCoarserPtr->m_scalarNew[0]->disjointBoxLayout(),
                                   1, amrGodCoarserPtr->m_problem_domain,
                                   amrGodCoarserPtr->m_ref_ratio, 1);

      pwl.fillInterp(*m_scalarNew[taggingVar],
                     *amrGodCoarserPtr->m_scalarNew[taggingVar],
                     *amrGodCoarserPtr->m_scalarNew[taggingVar], 1.0, 0, 0, 1);
    }
    if (taggingVectorVar >=0)
    {
      PiecewiseLinearFillPatch pwl(levelDomain,
                                   amrGodCoarserPtr->m_vectorNew[0]->disjointBoxLayout(),
                                   SpaceDim, amrGodCoarserPtr->m_problem_domain,
                                   amrGodCoarserPtr->m_ref_ratio, 1);

      pwl.fillInterp(*m_vectorNew[taggingVectorVar],
                     *amrGodCoarserPtr->m_vectorNew[taggingVectorVar],
                     *amrGodCoarserPtr->m_vectorNew[taggingVectorVar], 1.0, 0, 0, SpaceDim);
    }

    if (taggingVectorVar < 0 && taggingVar < 0)
    {
      MayDay::Error("No valid tagging variable");
    }
  }

  if (taggingVar >= 0)
  {
    m_scalarNew[taggingVar]->exchange(Interval(0, 0));
  }

  if (taggingVectorVar > 0)
  {
    m_vectorNew[taggingVectorVar]->exchange(Interval(0, SpaceDim - 1));
  }

  // Compute undivided gradient
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& b = levelDomain[dit()];
    FArrayBox gradFab(b, SpaceDim);
    FArrayBox UFab((*m_scalarNew[0])[dit()].box(), 1);
    UFab.setVal(0.0);

    if (taggingVar >= 0)
    {
      UFab.plus((*m_scalarNew[taggingVar])[dit()]);
    }
    else
    {
      if (comp == -1)
      {
        // Get scalar magnitude of vector field
        FORT_MAGNITUDEF(CHF_FRA1(UFab, 0), CHF_CONST_FRA((*m_vectorNew[taggingVectorVar])[dit()]),
                        CHF_BOX(UFab.box()));
      }
      else
      {
        UFab.plus((*m_vectorNew[taggingVectorVar])[dit()], comp, 0);
      }
    }

    FArrayBox taggingMetricFab(b, 1);

    if (taggingMethod == TaggingMethod::UndividedGradient)
    {
      // Calculated undivided gradient
      for (int dir = 0; dir < SpaceDim; ++dir)
      {
        const Box bCenter = b & grow(m_problem_domain, -BASISV(dir));
        const Box bLo = b & adjCellLo(bCenter, dir);
        const int hasLo = !bLo.isEmpty();
        const Box bHi = b & adjCellHi(bCenter, dir);
        const int hasHi = !bHi.isEmpty();
        FORT_GETGRADF(CHF_FRA1(gradFab, dir), CHF_CONST_FRA1(UFab, 0),
                      CHF_CONST_INT(dir), CHF_BOX(bLo), CHF_CONST_INT(hasLo),
                      CHF_BOX(bHi), CHF_CONST_INT(hasHi), CHF_BOX(bCenter));
      }

      FORT_MAGNITUDEF(CHF_FRA1(taggingMetricFab, 0), CHF_CONST_FRA(gradFab),
                      CHF_BOX(b));
    }
    else if (taggingMethod == TaggingMethod::Magnitude)
    {
      // Calculate magnitude
      FORT_MAGNITUDESIGN(CHF_FRA1(taggingMetricFab, 0), CHF_CONST_FRA(UFab),
                         CHF_BOX(b));
      // Scale it?
      // Let's stop doing this - so regions of with high values will also be at max refinement
      //      taggingMetricFab.divide(m_level + 1);

    }
    else if (taggingMethod == TaggingMethod::CompareLargerThan)
    {
      // Just compare the actual field
      taggingMetricFab.copy(UFab);

    }
    else if (taggingMethod == TaggingMethod::CompareLargerThanNegative)
    {
      // Just compare the actual field, but make negative (so we refine where the field is negative).
      taggingMetricFab.copy(UFab);
      taggingMetricFab.mult(-1);
    }

    // Tag where gradient exceeds threshold
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      if (taggingMetricFab(iv) >= refineThresh)
      {
        localTags |= iv;
      }
    }
  }
}

/*******/
void AMRLevelMushyLayer::tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

/*******/
void AMRLevelMushyLayer::regrid(const Vector<Box>& a_newGrids)
{
  // Always print when we regrid
  if (s_verbosity >= 0)
  {
    pout() << "AMRLevelMushyLayer::regrid (level " << m_level << ")" << endl;
  }

  if (s_verbosity >= 4)
  {
    for (int i=0; i < a_newGrids.size(); i++)
    {
      const Box& b = a_newGrids[i];
      pout() << "  Box spanning (" << b.smallEnd() << ") -  (" << b.bigEnd() << ")" << endl;
    }

  }

  // Are new grids any different?
  m_newGrids_different = true;

  Vector<int> new_procs;
  LoadBalance(new_procs, a_newGrids);
  DisjointBoxLayout new_grids = DisjointBoxLayout(a_newGrids, new_procs, m_problem_domain);
  if (new_grids.sameBoxes(m_grids))
  {
    m_newGrids_different = false;
  }

  pout() << "  New grids are different (time = " << m_time << ")? " << m_newGrids_different << endl;

  if (m_newGrids_different)
  {
    m_level_grids = a_newGrids;
  }
  else
  {
    return;
  }

  // Save original grids and load balance
  Vector<int> procs;
  LoadBalance(procs, m_level_grids);
  m_grids = DisjointBoxLayout(m_level_grids, procs, m_problem_domain);

  // Save data for later
  Vector<RefCountedPtr<LevelData<FArrayBox> > > scalarOld_OldGrids, scalarNew_OldGrids,
  vectorOld_OldGrids, vectorNew_OldGrids;

  IntVect iv_ghost = m_scalarNew[0]->ghostVect();
  DisjointBoxLayout old_grids = m_scalarNew[0]->disjointBoxLayout();

  DataIterator ditOld = old_grids.dataIterator();

  scalarNew_OldGrids.resize(m_numScalarVars);
  vectorNew_OldGrids.resize(m_numVectorVars);
  scalarOld_OldGrids.resize(m_numScalarVars);
  vectorOld_OldGrids.resize(m_numVectorVars);

  for (int scalarVar = 0; scalarVar < m_numScalarVars; scalarVar++)
  {
    scalarNew_OldGrids[scalarVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(old_grids, 1, iv_ghost));

    scalarOld_OldGrids[scalarVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(old_grids, 1, iv_ghost));

    for (ditOld.reset(); ditOld.ok(); ++ditOld)
    {
      (*scalarNew_OldGrids[scalarVar])[ditOld].copy((*m_scalarNew[scalarVar])[ditOld]);
      (*scalarOld_OldGrids[scalarVar])[ditOld].copy((*m_scalarOld[scalarVar])[ditOld]);
    }

    //    m_scalarNew[scalarVar]->copyTo(m_scalarNew[scalarVar]->interval(),
    //                                      *scalarNew_OldGrids[scalarVar],
    //                                      scalarNew_OldGrids[scalarVar]->interval());
    //    m_scalarOld[scalarVar]->copyTo(m_scalarOld[scalarVar]->interval(),
    //                                   *scalarOld_OldGrids[scalarVar],
    //                                   scalarOld_OldGrids[scalarVar]->interval());
  }

  for (int vectorVar = 0; vectorVar < m_numVectorVars; vectorVar++)
  {
    vectorNew_OldGrids[vectorVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(old_grids, SpaceDim, iv_ghost));

    vectorOld_OldGrids[vectorVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(old_grids, SpaceDim, iv_ghost));

    for (ditOld.reset(); ditOld.ok(); ++ditOld)
    {
      (*vectorNew_OldGrids[vectorVar])[ditOld].copy((*m_vectorNew[vectorVar])[ditOld]);
      (*vectorOld_OldGrids[vectorVar])[ditOld].copy((*m_vectorOld[vectorVar])[ditOld]);
    }

    //    m_vectorNew[vectorVar]->copyTo(m_vectorNew[vectorVar]->interval(),
    //                                       *vectorNew_OldGrids[vectorVar],
    //                                       vectorNew_OldGrids[vectorVar]->interval());
    //    m_vectorOld[vectorVar]->copyTo(m_vectorOld[vectorVar]->interval(),
    //                                   *vectorOld_OldGrids[vectorVar],
    //                                   vectorOld_OldGrids[vectorVar]->interval());
  }


  // Create new grids
  createDataStructures();

  // Set up data structures
  // We do this in postregrid
  levelSetup(); // still need to do this to have m_fineInterp defined
  //  defineSolvers();

  // Interpolate from coarser level
  if (m_hasCoarser)
  {
    AMRLevelMushyLayer* amrMushyLayerCoarserPtr = getCoarserLevel();

    // Crucial to use 4th order scheme here
    // Lower order doesn't produce a sufficiently smooth fine level solution

    FourthOrderFineInterp scalarInterp4, vectorInterp4;
    FineInterp scalarInterp, vectorInterp;

    int crseRefRat = amrMushyLayerCoarserPtr->m_ref_ratio;
    scalarInterp4.define(m_grids, 1, crseRefRat, m_problem_domain);
    vectorInterp4.define(m_grids, SpaceDim, crseRefRat, m_problem_domain);

    scalarInterp.define(m_grids, 1, crseRefRat, m_problem_domain);
    vectorInterp.define(m_grids, SpaceDim, crseRefRat, m_problem_domain);


    for (int scalarVar = 0; scalarVar < m_numScalarVars;
        scalarVar++)
    {
      if (m_opt.scalarHOinterp)
      {
        scalarInterp4.interpToFine(*m_scalarNew[scalarVar],
                                   *(amrMushyLayerCoarserPtr->m_scalarNew[scalarVar]));
        scalarInterp4.interpToFine(*m_scalarOld[scalarVar],
                                   *(amrMushyLayerCoarserPtr->m_scalarOld[scalarVar]));
      }
      else
      {
        scalarInterp.interpToFine(*m_scalarNew[scalarVar],
                                  *(amrMushyLayerCoarserPtr->m_scalarNew[scalarVar]));
        scalarInterp.interpToFine(*m_scalarOld[scalarVar],
                                  *(amrMushyLayerCoarserPtr->m_scalarOld[scalarVar]));
      }

      ParmParse pp("regrid");
      bool initLambdaOne = false;
      pp.query("init_lambda_one", initLambdaOne);
      if (scalarVar == ScalarVars::m_lambda && initLambdaOne)
      {
        setValLevel(*m_scalarNew[scalarVar], 1.0);
      }

    }

    for (int vectorVar = 0; vectorVar < m_numVectorVars;
        vectorVar++)
    {

      if (m_opt.vectorHOinterp)
      {
        vectorInterp4.interpToFine(*m_vectorNew[vectorVar],
                                   *(amrMushyLayerCoarserPtr->m_vectorNew[vectorVar]));
        vectorInterp4.interpToFine(*m_vectorOld[vectorVar],
                                   *(amrMushyLayerCoarserPtr->m_vectorOld[vectorVar]));
      }
      else
      {
        vectorInterp.interpToFine(*m_vectorNew[vectorVar],
                                  *(amrMushyLayerCoarserPtr->m_vectorNew[vectorVar]));
        vectorInterp.interpToFine(*m_vectorOld[vectorVar],
                                  *(amrMushyLayerCoarserPtr->m_vectorOld[vectorVar]));
      }

    }

  }

  // Copy from old grids into new grids

  for (int scalarVar = 0; scalarVar < m_numScalarVars; scalarVar++)
  {

    scalarNew_OldGrids[scalarVar]->copyTo(
        scalarNew_OldGrids[scalarVar]->interval(),
        *m_scalarNew[scalarVar],
        m_scalarNew[scalarVar]->interval());

    scalarOld_OldGrids[scalarVar]->copyTo(
        scalarOld_OldGrids[scalarVar]->interval(),
        *m_scalarOld[scalarVar],
        m_scalarOld[scalarVar]->interval());
  }

  for (int vectorVar = 0; vectorVar < m_numVectorVars; vectorVar++)
  {
    vectorNew_OldGrids[vectorVar]->copyTo(
        vectorNew_OldGrids[vectorVar]->interval(),
        *m_vectorNew[vectorVar],
        m_vectorNew[vectorVar]->interval());

    vectorOld_OldGrids[vectorVar]->copyTo(
        vectorOld_OldGrids[vectorVar]->interval(),
        *m_vectorOld[vectorVar],
        m_vectorOld[vectorVar]->interval());
  }

  // Re calculate analytic solns on new grid
  // May want to turn this off if it starts taking a lot of time
  calculateAnalyticSolns(false);

  fillFrameVelocity();

  // Ensure the new theta, porosity etc. are consistent with the regridded H and S
  // This was introducing errors during regridding for some reason - turn it off
  //  updateEnthalpyVariables();

  //    int temp=0;
}


// Note - for levels > 0, we need the lower level to be defined consistently for doing interpolation
void AMRLevelMushyLayer::refine(Real ref_ratio, DisjointBoxLayout a_grids, ProblemDomain a_domain)
{
  // Need data structures to backup old data
  RefCountedPtr<LevelData<FArrayBox> > previousScal, previousVect;
  LevelData<FluxBox> prevAdvVel;

  IntVect ivGhost = 4*IntVect::Unit;

  Interval scalInterval(0,0);
  Interval vectInterval(0, SpaceDim-1);

  FineInterp scalarInterp, vectorInterp;
  scalarInterp.define(a_grids, 1, ref_ratio, a_domain);
  vectorInterp.define(a_grids, SpaceDim, ref_ratio, a_domain);

//  prevAdvVel.define(m_grids, 1, advectionGhost);
  previousScal = RefCountedPtr<LevelData<FArrayBox> >(
      new LevelData<FArrayBox>(m_grids, 1, ivGhost));
  previousVect = RefCountedPtr<LevelData<FArrayBox> >(
      new LevelData<FArrayBox>(m_grids, SpaceDim, ivGhost));

  // FineInterp doesn't handle fluxboxes, so can't refine advection velocity for now
//    m_advVel.copyTo(scalInterval, prevAdvVel, scalInterval);
//    m_advVel.define(newGrids, 1, advectionGhost); // reshape
  //  prevAdvVel.copyTo(scalInterval,m_advVel, scalInterval); // copy back

  for (int scalarVar = 0; scalarVar < m_numScalarVars; scalarVar++)
  {
//    m_scalarNew[scalarVar]->copyTo(scalInterval, *previousScal, scalInterval);
    fillScalars(*previousScal, m_time, scalarVar, true, true);
    m_scalarNew[scalarVar]->define(a_grids, 1, ivGhost); //reshape
    if (m_opt.regrid_linear_interp)
    {
      scalarInterp.interpToFine(*m_scalarNew[scalarVar], *previousScal);
    }
    else
    {
      scalarInterp.pwcinterpToFine(*m_scalarNew[scalarVar], *previousScal);
    }
  }

  for (int vectorVar = 0; vectorVar < m_numVectorVars; vectorVar++)
  {
    m_vectorNew[vectorVar]->copyTo(vectInterval, *previousVect, vectInterval);
    m_vectorNew[vectorVar]->define(a_grids, SpaceDim, ivGhost); //reshape
    if (m_opt.regrid_linear_interp)
    {
      vectorInterp.interpToFine(*m_vectorNew[vectorVar], *previousVect);
    }
    else
    {
      vectorInterp.pwcinterpToFine(*m_vectorNew[vectorVar], *previousVect);
    }
  }

  //  ProblemDomain oldDomain = m_problem_domain;
  m_problem_domain = a_domain;
  m_grids = a_grids;
}

void AMRLevelMushyLayer::postRegrid(int a_base_level)
{
  int methodVersion = 1;
  ParmParse pp("regrid");
  pp.query("post_method", methodVersion);

  // Create data structures for new grids
  levelSetup();

  if (methodVersion == 0)
  {
    postRegridOld(a_base_level);
  }
  else
  {
    postRegridNew(a_base_level);
  }
}

void AMRLevelMushyLayer::postRegridNew(int a_base_level)
{

  // Don't run method if this isn't the finest level, or this is the base level
  if (!finestLevel() || m_level == a_base_level)
  {
    return;
  }

  pout() << "postRegridNew between base level " << a_base_level << " and current level " << m_level << endl;
  pout() << "     level " << m_level << ": time=" << m_time << ", dt=" << m_dt << endl;
  pout() << "     level " << getCoarserLevel()->m_level << ": time=" << getCoarserLevel()->m_time << ", dt=" << getCoarserLevel()->m_dt << endl;

  // Get the mushy layer object for the base level
  AMRLevelMushyLayer* ml_base = this;
  while (ml_base->level() > a_base_level)
  {
    ml_base = ml_base->getCoarserLevel();
  }

  CH_assert(ml_base->level() == a_base_level);

  // Fill AMR data
  int numLevels = m_level + 1;

  Vector<LevelData<FArrayBox>*> amrVel(numLevels);
  Vector<LevelData<FArrayBox>*> amrLambda(numLevels);
  Vector<RefCountedPtr<LevelData<FluxBox> > > amrPorosityFace(numLevels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > amrPorosity(numLevels);

  fillAMRVelPorosity(amrVel, amrPorosityFace, amrPorosity);
  fillAMRLambda(amrLambda);

  Vector<int> nRef;
  AMRLevelMushyLayer* mlTemp = getCoarsestLevel();
  while (mlTemp)
  {
    nRef.push_back(mlTemp->m_ref_ratio);
    mlTemp = mlTemp->getFinerLevel();
  }

  Projector& base_projection = ml_base->m_projection;
  Real base_dt = ml_base->dt();

  // Fill BCs
  AMRLevelMushyLayer* thisLevelData = ml_base;

  VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_opt.viscousBCs));

  for (int lev = a_base_level; lev <= m_level; lev++)
  {
    const ProblemDomain& levelDomain = thisLevelData->problemDomain();
    Real levelDx = thisLevelData->m_dx;
    LevelData<FArrayBox>& thisAmrVel = *amrVel[lev];
    const DisjointBoxLayout& thisLevelGrids = thisAmrVel.getBoxes();
    velBC.applyBCs(thisAmrVel, thisLevelGrids, levelDomain, levelDx,
                   false); // inhomogeneous

    // This deals with periodic BCs
    thisAmrVel.exchange();

    if (thisLevelData->m_finer_level_ptr != nullptr)
    {
      thisLevelData =
          dynamic_cast<AMRLevelMushyLayer*>(thisLevelData->m_finer_level_ptr);
    }
  }

  bool homogeneousCFBC = true;
  base_projection.initialVelocityProject(amrVel, amrPorosityFace, amrPorosity, homogeneousCFBC);

  Real maxLambda = ::computeNorm(amrLambda, nRef, getCoarsestLevel()->dx(), Interval(0,0), 0, getCoarsestLevel()->m_level) - 1.0;
  m_diagnostics.addDiagnostic(DiagnosticNames::diag_postRegridLambda, m_time+m_dt, maxLambda);
  Real crseNewTime =getCoarsestLevel()->m_time + getCoarsestLevel()->m_dt;
  getCoarsestLevel()->m_diagnostics.addDiagnostic(DiagnosticNames::diag_postRegridLambda, crseNewTime, maxLambda);
  pout() << "PostRegrid(level " << m_level <<", time " << m_time <<") - max|Lambda| = " << maxLambda << endl;

  base_projection.doPostRegridOps(amrLambda,amrPorosityFace,
                                  base_dt,m_time, 1.0);

  // Reset BCs
  thisLevelData = ml_base;
  for (int lev = a_base_level; lev <= m_level; lev++)
    {
      const ProblemDomain& levelDomain = thisLevelData->problemDomain();
      Real levelDx = thisLevelData->m_dx;
      LevelData<FArrayBox>& thisAmrVel = *amrVel[lev];
      const DisjointBoxLayout& thisLevelGrids = thisAmrVel.getBoxes();
      velBC.applyBCs(thisAmrVel, thisLevelGrids, levelDomain, levelDx,
                     false); // inhomogeneous

      // This deals with periodic BCs
      thisAmrVel.exchange();

      if (thisLevelData->m_finer_level_ptr != nullptr)
      {
        thisLevelData =
            dynamic_cast<AMRLevelMushyLayer*>(thisLevelData->m_finer_level_ptr);
      }
    }

  if (m_opt.initialize_pressures)
  {
    ml_base->initializeGlobalPressureNew();
  }

  // Sanity checking
  pout() << "postRegridNew (level " << m_level << ") finished. time=" << m_time << ", dt=" << m_dt << endl;


}

void AMRLevelMushyLayer::postRegridOld(int a_base_level)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::postRegrid (level " << m_level << ")" << endl;
  }

  // If the new grids are no different to old grids, don't do anything
  // Need to check across all levels

  AMRLevelMushyLayer* mlPtr = this;
  while(mlPtr->getCoarserLevel())
  {
    mlPtr =  mlPtr->getCoarserLevel();
  }

  bool allGridsTheSame = !mlPtr->m_newGrids_different;

  while(mlPtr)
  {
    allGridsTheSame = allGridsTheSame && !mlPtr->m_newGrids_different;
    mlPtr = mlPtr->getFinerLevel();
  }

  if (allGridsTheSame)
  {
    return;
  }

  // Only do this on level 0
  if (m_level == 0 && m_opt.do_postRegrid_smoothing)
  {
    doPostRegridSmoothing();
  }

  //Need to ensure this is done on the base level
  levelSetup();
  defineSolvers(m_time); // not entirely sure time this should be

  if (m_level == 0 && m_opt.doProjection)
  {
    AMRLevelMushyLayer* thisLevelData = this;

    thisLevelData = this;
    while (thisLevelData->hasFinerLevel())
    {
      thisLevelData = thisLevelData->getFinerLevel();
    }
    //          CH_assert(thisLevelData->finestLevel());
    int numLevels = thisLevelData->m_level + 1;

    Vector<LevelData<FArrayBox>*> amrVel(numLevels);
    Vector<LevelData<FArrayBox>*> amrLambda(numLevels);
    Vector<RefCountedPtr<LevelData<FluxBox> > > amrPorosityFace(numLevels);
    Vector<RefCountedPtr<LevelData<FArrayBox> > > amrPorosity(numLevels);

    fillAMRVelPorosity(amrVel, amrPorosityFace, amrPorosity);
    fillAMRLambda(amrLambda);
    Projector& level0Proj = m_projection;

    // set physical boundary conditions on velocity

    thisLevelData = this;

    VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_opt.viscousBCs));

    for (int lev = 0; lev < numLevels; lev++)
    {
      const ProblemDomain& levelDomain = thisLevelData->problemDomain();
      Real levelDx = thisLevelData->m_dx;
      LevelData<FArrayBox>& thisAmrVel = *amrVel[lev];
      const DisjointBoxLayout& thisLevelGrids = thisAmrVel.getBoxes();
      velBC.applyBCs(thisAmrVel, thisLevelGrids, levelDomain, levelDx,
                     false); // inhomogeneous

      // This deals with periodic BCs
      thisAmrVel.exchange();

      //      amrLambda[lev] = &(*thisLevelData->m_scalarNew[ScalarVars::m_lambda]);

      if (thisLevelData->m_finer_level_ptr != nullptr)
      {
        thisLevelData =
            dynamic_cast<AMRLevelMushyLayer*>(thisLevelData->m_finer_level_ptr);
      }
    }

    if (m_opt.project_initial_vel)
    {
      // NB this doesn't calculate pressure, just corrects velocity
      if (m_opt.makeRegridPlots)
      {
        writeAMRHierarchy("regrid1.hdf5");
      }

      level0Proj.initialVelocityProject(amrVel, amrPorosityFace, amrPorosity,
                                        false); // inhomogenous BC
      if (m_opt.makeRegridPlots)
      {
        writeAMRHierarchy("regrid2.hdf5");
      }
    }

    setVelBCs(numLevels, amrVel, velBC);

    updateEnthalpyVariables();
    updateEnthalpyVariablesOld();

    int finest_level = numLevels-1;
    Real dtInit = computeDtInit(finest_level);

    if (m_opt.regrid_init_pressure)
    {
      if (m_opt.makeRegridPlots)
      {
        writeAMRHierarchy("regrid3.hdf5");
      }

      dtInit *= m_opt.regrid_dt_scale;

      if (solvingFullDarcyBrinkman())
      {
        initializeGlobalPressure(dtInit, false);
      }

      if (!(solvingFullDarcyBrinkman() && doVelocityAdvection()))
      {

        // Should do some initialisation here?
        thisLevelData = this;

        for (int lev = 0; lev < numLevels; lev++)
        {
          int num_passes = 1;
          initTimeIndependentPressure(thisLevelData, num_passes);
          thisLevelData = thisLevelData->getFinerLevel();
        }

      }

    }

    if (m_opt.makeRegridPlots)
    {
      writeAMRHierarchy("regrid4.hdf5");
    }

    if(m_opt.initLambdaPostRegrid)
    {
      thisLevelData = this;

      for (int lev = 0; lev < numLevels; lev++)
      {
        thisLevelData->setFluxRegistersZero();
        thisLevelData = thisLevelData->getFinerLevel();
      }

      // Finally, let's try and compute the freestream preservation carefully (i.e. use subcycling)

      // Save dt for each level so it can be reset later
      Vector<Real> dtSave(numLevels, 1.0e8);
      thisLevelData = this;
      for (int lev = m_level; lev < numLevels; lev++)
      {
        dtSave[lev] = thisLevelData->m_dt;
        thisLevelData->backupTimestep();

        thisLevelData = thisLevelData->getFinerLevel();
      }


      if (m_opt.regrid_advect_before_freestream)
      {
        if (s_verbosity > 2)
        {
          pout() << "AMRLevelMushyLayer::postRegrid - advecting lambda on all levels to compute correction" << endl;
        }

        // This is a fake subcycled advance
        // Only currently works for 2 levels?
        thisLevelData = this;

        Real dtLev = dtInit;
        for (int lev = 0; lev < numLevels; lev++)
        {
          Real tlev = 0;
          while(tlev < dtInit)
          {

            thisLevelData->dt(dtLev);
            thisLevelData->computeAllVelocities(true);
            thisLevelData->advectLambda(true);

            tlev = tlev + dtLev;

          }

          if (m_opt.regrid_freestream_subcycle)
          {
            dtLev = dtLev/thisLevelData->m_ref_ratio;
          }

          thisLevelData = thisLevelData->getFinerLevel();

        }

      } // end if doing advection of lambda

      if (m_opt.makeRegridPlots)
      {
        writeAMRHierarchy("regrid5.hdf5");
      }

      // Reflux Lambda to compute VD correction
      if (m_opt.regrid_reflux_lambda)
      {
        AMRRefluxLambda();
      }

      if (m_opt.makeRegridPlots)
      {
        writeAMRHierarchy("regrid6.hdf5");
      }

      fillAMRVelPorosity(amrVel, amrPorosityFace, amrPorosity);
      fillAMRLambda(amrLambda);

      Vector<int> nRef;
      AMRLevelMushyLayer* mlTemp = getCoarsestLevel();
      while (mlTemp)
      {
        nRef.push_back(mlTemp->m_ref_ratio);
        mlTemp = mlTemp->getFinerLevel();
      }

      Real maxLambda = ::computeNorm(amrLambda, nRef, getCoarsestLevel()->dx(), Interval(0,0), 0, getCoarsestLevel()->m_level)-1.0;
      m_diagnostics.addDiagnostic(DiagnosticNames::diag_postRegridLambda, m_time, maxLambda);
      Real crseNewTime =getCoarsestLevel()->m_time + getCoarsestLevel()->m_dt;
      getCoarsestLevel()->m_diagnostics.addDiagnostic(DiagnosticNames::diag_postRegridLambda, crseNewTime, maxLambda);

      level0Proj.doPostRegridOps(amrLambda,amrPorosityFace,dtInit,m_time, m_opt.regrid_eta_scale);

      if (m_opt.makeRegridPlots)
      {
        writeAMRHierarchy("regrid7.hdf5");
      }

      // Reset dt
      thisLevelData = this;
      for (int lev = m_level;  lev < numLevels; lev++)
      {
        thisLevelData->dt(dtSave[lev]);
        thisLevelData->restartTimestepFromBackup(true); // true - don't replace pressure
        thisLevelData = thisLevelData->getFinerLevel();
      }

      if (m_opt.makeRegridPlots)
      {
        writeAMRHierarchy("regrid8.hdf5");
      }

    } // end if initialising lambda correction

  } // end if level 0 and doing projection

}

void AMRLevelMushyLayer::fillAMRLambda(Vector<LevelData<FArrayBox>*>& amrLambda)
{
  //  CH_assert(m_level==0);
  AMRLevelMushyLayer* thisLevelData = getCoarsestLevel();

  int numLevels = amrLambda.size();

  thisLevelData = getCoarsestLevel();

  for (int lev = 0; lev < numLevels; lev++)
  {

    amrLambda[lev] = &(*thisLevelData->m_scalarNew[ScalarVars::m_lambda]);

    if (thisLevelData->m_finer_level_ptr != nullptr)
    {
      thisLevelData =
          dynamic_cast<AMRLevelMushyLayer*>(thisLevelData->m_finer_level_ptr);
    }
  }
}
