/*
 * AMRLevelMushyLayerBiogeochemistry.cpp
 *
 *  Created on: 20 Aug 2019
 *      Author: parkinsonjl
 *
 *      source file to contain methods for doing biogeochemistry
 */
#include "AMRLevelMushyLayer.H"
#include "VCAMRPoissonOp2.H"


void AMRLevelMushyLayer::computeRadianceIntensity()
{
  // Compute radiation intensity by Beer-Lambert law (or something better)

  // First need to compute transmittance at each point in space
  LevelData<FArrayBox> attenuation(m_grids, 1, IntVect::Unit);
  LevelData<FArrayBox> attenuationSum(m_grids, 1, IntVect::Unit);

 Real liquid_attenuation = 0.01;
  Real solid_attenuation = 1.0; // nondimensionalised

  // T= (1-chi)*Ts + chi*Tl
  // T = Ts + chi*(Tl-Ts)

  for (DataIterator dit = attenuation.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox& porosity = (*m_scalarNew[m_porosity])[dit];
    attenuation[dit].setVal(liquid_attenuation-solid_attenuation);
    attenuation[dit].mult(porosity, 0, 0, 1);
    attenuation[dit].plus(solid_attenuation);
  }

  attenuation.exchange();

  // Integrate in columns
  // intensity = 10^{ - int_{top}^{z} T dz}
  attenuation.copyTo(attenuationSum); // fill top row
  Box domBox = m_problem_domain.domainBox();
  domBox.growHi(SpaceDim-1, -1);

  int lo_i = domBox.smallEnd()[0];
  int hi_i = domBox.bigEnd()[0];

  for (int z_i=domBox.bigEnd()[SpaceDim-1]; z_i >= domBox.smallEnd()[SpaceDim-1]; z_i--)
  {
    // need to find which box contains this z_i

    for (DataIterator dit = attenuation.dataIterator(); dit.ok(); ++dit)
       {
         Box b = m_grids[dit];


//         b.growLo(SpaceDim-1, -1);

         IntVect iv_lo = IntVect(lo_i, z_i);

         if (b.contains(iv_lo))
         {

           IntVect iv_hi = IntVect(hi_i, z_i);
           Box horizBox = Box(iv_lo, iv_hi);

         for (BoxIterator bit = BoxIterator(horizBox); bit.ok(); ++bit)
         {
           IntVect iv = bit();

           IntVect ivUp = iv + BASISV(SpaceDim-1);
           attenuationSum[dit](iv) = attenuationSum[dit](ivUp) + attenuation[dit](iv)*m_dx;

         }
         }
       }

    attenuationSum.exchange();

  }




  // Compute intensity
  for (DataIterator dit = attenuation.dataIterator(); dit.ok(); ++dit)
  {
    Box b = m_grids[dit];
    for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
    {
      IntVect iv= bit();
      (*m_scalarNew[ScalarVars::m_lightIntensity])[dit](iv) = pow(10, -attenuationSum[dit](iv));
    }
   }


}


void AMRLevelMushyLayer::advectTracer(int a_tracerVar, LevelData<FArrayBox>& a_src)
{
  LevelFluxRegister* coarserFRPtr = NULL;
  LevelFluxRegister* finerFRPtr = NULL;
  LevelData<FArrayBox>* coarserDataOldPtr = NULL;
  LevelData<FArrayBox>* coarserDataNewPtr = NULL;
  Real tCoarserOld, tCoarserNew;

  getCoarseScalarDataPointers(a_tracerVar,
                              &coarserDataOldPtr,  &coarserDataNewPtr,  // we don't use these two, but need them as dummy arguments
                              &coarserFRPtr, &finerFRPtr, // get the flux registers for the thing we're updating, a_scalarVar
                              tCoarserOld, tCoarserNew); // don't need these either, they're just dummy arguments

//  if (!doFRupdates)
//  {
//    coarserFRPtr = NULL;
//    finerFRPtr = NULL;
//  }

  DataIterator dit(m_grids);

//  LevelData<FArrayBox> liquid_tracer_conc(m_grids, 1, this->m_numGhostAdvection*IntVect::Unit);
//  for (DataIterator dit = liquid_tracer_conc.dataIterator(); dit.ok(); ++dit)
//  {
//    liquid_tracer_conc[dit].copy((*m_scalarNew[a_tracerVar])[dit]);
//    liquid_tracer_conc[dit].divide((*m_scalarNew[m_porosity])[dit]);
//  }
//
//  liquid_tracer_conc.exchange();
  LevelData<FArrayBox> liquid_tracer_conc;
  computeScalarConcInLiquid(liquid_tracer_conc, a_tracerVar);


  // should have one less ghost cell than src
  LevelData<FluxBox> flux(m_grids, 1, a_src.ghostVect()-IntVect::Unit);

  // Get the flux of a_advectionVar, i.e. u*a_advectionVar
  //    computeScalarAdvectiveFlux(flux, a_advectionVar, -1, a_advVel, m_time-m_dt, m_dt); // -1 means no diffusive src
  computeScalarAdvectiveFlux(flux, liquid_tracer_conc, a_src, m_advVel, a_tracerVar, m_time-m_dt, m_dt); // -1 means no diffusive src

  // Need to get tracer var/chi to advect

  // Make the source term, div(u*a_advectionVar)
  LevelData<FArrayBox> update(m_grids, 1);
  Divergence::levelDivergenceMAC(update, flux, m_dx);

  // Add the source term to the old time solution
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    update[dit].mult(m_dt);
    (*m_scalarNew[a_tracerVar])[dit] -= update[dit];

  }

  // Flux register updates

    Real scale = m_dt;
    updateScalarFluxRegister(a_tracerVar, flux, scale);


}


void AMRLevelMushyLayer::computeScalarConcInLiquid(LevelData<FArrayBox>& liquid_tracer_conc, int a_tracerVar)
{
  liquid_tracer_conc.define(m_grids, 1, this->m_numGhostAdvection*IntVect::Unit);
  for (DataIterator dit = liquid_tracer_conc.dataIterator(); dit.ok(); ++dit)
  {
    liquid_tracer_conc[dit].copy((*m_scalarNew[a_tracerVar])[dit]);
    liquid_tracer_conc[dit].divide((*m_scalarNew[m_porosity])[dit]);
  }

  liquid_tracer_conc.exchange();
}

void AMRLevelMushyLayer::advectPassiveTracer()
{

  LevelData<FArrayBox> a_src(m_grids, 1, 1*IntVect::Unit);

  computeScalarDiffusiveSrc(ScalarVars::m_passiveScalar, a_src);

//  setValLevel(a_src, 0);
  advectTracer(ScalarVars::m_passiveScalar, a_src);

}

void AMRLevelMushyLayer::computeScalarDiffusiveSrc(int a_scalarBulkConc, LevelData<FArrayBox>& a_src)
{
  BCHolder bc;
  getScalarBCs(bc, a_scalarBulkConc, false);

  // Op should actually be a variable coefficient op with porosity as the variable coefficient

  AMRPoissonOpFactory* op = new AMRPoissonOpFactory();
  op->define(m_problem_domain, m_grids, m_dx, bc);
  RefCountedPtr<AMRPoissonOpFactory> OpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(op);
  RefCountedPtr<AMRPoissonOp> amrpop =  RefCountedPtr<AMRPoissonOp>(
      (AMRNonLinearMultiCompOp*) OpFact->AMRnewOp(m_problem_domain) );

  // todo - finish this
  VCAMRPoissonOp2* vcop = new VCAMRPoissonOp2();
  //op->define(m_problem_domain, m_grids, m_dx, bc, -1, 0.0, 1);
//  vcop->define(m_grids, a_gridsFiner, a_gridsCoarser, a_dxLevel, a_refRatio, a_refRatioFiner, a_domain, a_bc, a_exchange, a_cfregion, a_nComp)
//

  LevelData<FArrayBox> *crseVar = NULL;

//  AMRLevelMushyLayer* crseML = getCoarserLevel();
//  if (crseML)
//  {
//    crseVar = &(*crseML->m_scalarNew[a_var]);
//  }

  amrpop->setAlphaAndBeta(0, m_scalarDiffusionCoeffs[a_scalarBulkConc]);

  LevelData<FArrayBox> liquidConc;
  computeScalarConcInLiquid(liquidConc, a_scalarBulkConc);

  // This just calls applyOpI if crseHC = NULL, else does CF interpolation
  amrpop->applyOpMg(a_src, liquidConc, crseVar, false);

}


void AMRLevelMushyLayer::advectActiveTracer()
{
  // src needs one more int vect than flux
  LevelData<FArrayBox> a_src(m_grids, 1, IntVect::Unit);

  computeActiveTracerSourceSink(a_src);
  a_src.exchange();
  advectTracer(ScalarVars::m_activeScalar, a_src);

}

void AMRLevelMushyLayer::computeActiveTracerSourceSink(LevelData<FArrayBox>& a_srcSink)
{
  // Choose this function yourself

//  setValLevel(a_srcSink, 0.0);
  computeScalarDiffusiveSrc(ScalarVars::m_activeScalar, a_srcSink);

  // Example: grow tracer wherever we're a) in the ice (porosity < 1) and b) have enough light
  // growth rate proportional to light intensity
  // Also add some saturation so we don't grow any more if concentration > saturation_value
  // This is very simple, but demonstrates the sort of model we can consider

  // You may also want to consider a source which varies with other fields e.g. temperature

  Real critical_intensity = 0.1;
  Real saturation_value = 1.0;
  for (DataIterator dit = a_srcSink.dataIterator(); dit.ok(); ++dit)
  {
    Box b = a_srcSink[dit].box();

    FArrayBox& porosity = (*m_scalarNew[m_porosity])[dit];
    FArrayBox& intensity = (*m_scalarNew[ScalarVars::m_lightIntensity])[dit];
    FArrayBox& tracer_conc = (*m_scalarNew[m_activeScalar])[dit];

    for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
    {
      IntVect iv= bit();

      if (porosity(iv) < 1 && intensity(iv) > critical_intensity)
      {
        a_srcSink[dit](iv) =  a_srcSink[dit](iv) + (intensity(iv)-critical_intensity)*max(0.0, saturation_value-tracer_conc(iv));
      }

    }
  }

}



