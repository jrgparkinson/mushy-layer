#include "analyticSolns.H"



void channelFieldMushyLayer(LevelData<FArrayBox>& enthalpyAnalytic,
                                      LevelData<FArrayBox>& bulkConcentrationAnalytic,
                                      Real a_domainHeight, Real a_domainWidth, Real a_dx,
                                      MushyLayerParams a_parameters)
{

  // Enthalpy:
  for (DataIterator dit = enthalpyAnalytic.dataIterator(); dit.ok(); ++dit)
  {
    Box b = enthalpyAnalytic[dit].box();
    FArrayBox& H = enthalpyAnalytic[dit];
    FArrayBox& S= bulkConcentrationAnalytic[dit];
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      ::getLocation(iv, loc, a_dx);

      Real x = loc[0]; Real z = loc[1];


      Real HBottom = a_parameters.bcValEnthalpyLo[1];

      Real HTop = a_parameters.bcValEnthalpyHi[1] + 0.5*a_parameters.stefan;

      // Vertical salinity structure
      Real S_vert = exp(-200*pow(z - 0.5*a_domainHeight, 2));

      H(iv) = HBottom + (HTop - HBottom)*pow(z/a_domainHeight, 4);
      S(iv) = a_parameters.bcValBulkConcentrationLo[1] + (0.2) * exp(-pow(100*(x-(a_domainWidth/2)), 2))*S_vert;

    }
  }

}


void analyticSolnSolidificationNoFlow(LevelData<FArrayBox>& enthalpyAnalytic,
                                      LevelData<FArrayBox>& bulkConcentrationAnalytic,
                                      Real a_domainLength, Real a_dx,
                                      MushyLayerParams a_parameters)
{

  // First calculate z = z(theta) for finely spaced theta
  vector<Real> thetaGrid, thetaDefined;
  vector<Real> zCalc;
  int numPoints = (int) a_domainLength/a_dx;
  numPoints = numPoints * 4; //Use an even finer grid to make interpolation more accurate

  Real concRatio = a_parameters.compositionRatio;
  Real h;

  Real thetaBottom = a_parameters.bcValTemperatureLo[1];
  Real thetaTop = a_parameters.bcValTemperatureHi[1];

  Real porosityTop = a_parameters.bcValPorosityHi[1];

  Real vel = a_parameters.nonDimVel/a_parameters.m_heatDiffusionCoeff;

  Real zEutectic = a_domainLength; // Position of the eutectic. For now fix this at the top of the domain



  if (thetaTop < a_parameters.thetaEutectic)
  {
    Real HBottom = a_parameters.bcValEnthalpyLo[1];
    Real HTop = a_parameters.bcValEnthalpyHi[1];
    zEutectic = (HBottom-a_parameters.thetaEutectic)/(HBottom-HTop) * a_domainLength;
  }

  pout() << "zEutectic =  " << zEutectic << endl;

  if (abs(vel) > 0)
  {
    //Use a fixed point iteration to find thetaInf

    Real hotBoundary = max(thetaBottom, thetaTop);
    // Initial guesses for upper and lower bounds of theta infinity
    Real lower = hotBoundary; // theta inf can't be lower than this!
    Real upper = hotBoundary*100.0; // make this sufficiently large that we should definitely capture hot boundary

    pout() << setiosflags(ios::scientific) << setprecision(10);
    pout() << "Initial thetaInf guess:  " << lower << " < thetaInf < " << upper << endl;

    // Initial guess
    a_parameters.thetaInterface = a_parameters.thetaEutectic + 1; // Strictly: ( a_parameters.lewis * a_parameters.ThetaInf - a_parameters.thetaInf) / (a_parameters.lewis -1 );

    Real tolerance = 1e-8; // don't make this any smaller of we'll be inconsistent with double vs single precision
    Real hGuess, hotBoundaryGuess; //halfway, thetaInfGuess,
    hotBoundaryGuess = 0;

    int counter = 0;

    while(sqrt(abs(pow(thetaBottom-hotBoundaryGuess, 2))) > tolerance && counter < 100)
    {

      // Guess that theta infinity is halfway between our two bounds
      a_parameters.thetaInf =  (upper+lower)/2;

      // Find the predicted mush-liquid interface for this theta_infinity
      hGuess = a_parameters.directionalSolidificationMushyZ(a_parameters.thetaInterface, zEutectic);

      // Find the predicted value at the hot boundary for this theta_infinity
      hotBoundaryGuess = a_parameters.thetaInf + 	(1-a_parameters.thetaInf)*exp(vel*(0-hGuess));

      // If we guessed too high a temperature, reduce thetaInf
      if (hotBoundaryGuess > thetaBottom)
      {
        upper = a_parameters.thetaInf;
      }
      else
      {
        // We guessed too low, increase thetaInf
        lower = a_parameters.thetaInf;
      }

      pout() << setiosflags(ios::scientific) << setprecision(10);
      pout() << lower << " < thetaInf < " << upper << endl;

      counter++;
      if (counter % 10 == 0)
      {
        pout() << setiosflags(ios::scientific) << setprecision(10);
        pout() << counter << " iterations, thetaInf = " << a_parameters.thetaInf << endl;
      }
    }

    a_parameters.thetaInterface = ( - a_parameters.lewis * a_parameters.ThetaInf - a_parameters.thetaInf) / (a_parameters.lewis -1 );
    //	Real thetaInterface = a_parameters.thetaInterface;
    Real deltaTheta = (a_parameters.thetaInf -  a_parameters.thetaEutectic);

    //Constants needed for finding z(theta)

    h = a_parameters.directionalSolidificationMushyZ(a_parameters.thetaInterface, zEutectic); // position of the mush-liquid interface

    pout() << "Calculating analytic solution for solidification without flow: " << endl;
    pout() << "zEutectic was set to " << zEutectic << endl;
    pout() << "h = " << h << ", thetaInf = " << a_parameters.thetaInf << ", thetaInterface = " << a_parameters.thetaInterface << endl;

    // Need to start with low temperatures i.e. cold boundary
    for (int i = numPoints; i >= 0; i--)
    {
      Real theta;

      theta = a_parameters.thetaEutectic + deltaTheta*((float)i / (float)numPoints);

      thetaGrid.push_back(theta);

      Real z = a_parameters.directionalSolidificationMushyZ(theta, zEutectic);

      if (z > h)
      {
        zCalc.push_back(z);
        thetaDefined.push_back(theta);
      }
    }
  }
  else
  {
    // No frame advection
    a_parameters.thetaInf = thetaBottom;

    a_parameters.thetaInterface = ( a_parameters.lewis * a_parameters.ThetaInf - a_parameters.thetaInf) / (a_parameters.lewis -1 );

    h = (a_parameters.thetaInterface - thetaBottom)/(thetaTop - thetaBottom);
  }

  //Now interpolate to find theta(z)

  for (DataIterator dit = enthalpyAnalytic.dataIterator(); dit.ok(); ++dit)
  {
    Box a_box = enthalpyAnalytic[dit()].box();
    int nComp = enthalpyAnalytic.nComp();

    BoxIterator bit(a_box);

    FArrayBox thetaAnalytic(a_box, nComp);
    FArrayBox analyticSolidFrac(a_box, nComp);

    Real zero = 0.0;

    for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      Real z = (iv[SpaceDim-1]+0.5)*a_dx;

      if (abs(vel) == 0)
      {

        Real theta = thetaBottom + (thetaTop-thetaBottom ) * z;
        Real porosity = 1;
        if ( z >= h)
        {
          // assume linear chi profile (roughly correct)

          porosity = max(zero, 1 + (porosityTop - 1)*(z-h)/(1-h));
        }
        enthalpyAnalytic[dit](iv) = porosity*a_parameters.stefan + theta ;
        bulkConcentrationAnalytic[dit](iv) = 1;

      }
      else
      {
        Real theta=0.0, solidFraction=0.0;

        // If stefan=0 we're always in liquid
        if (z < h or a_parameters.stefan == 0 or zCalc.size() == 0)
        {
          //We're in the liquid
          theta = a_parameters.thetaInf + (a_parameters.thetaInterface - a_parameters.thetaInf)* exp(vel * (z-h));

          solidFraction = 0;
        }
        else if (z >= zEutectic)
        {
          // Liquid phase
          solidFraction = 1;
          Real linearTgradient = (thetaTop - a_parameters.thetaEutectic )/(a_domainLength-zEutectic);
          theta = a_parameters.thetaEutectic + linearTgradient*(z-zEutectic);
        }
        else if (abs(vel) > 0)
        {
          //Spline interpolation to get mush

           if (zCalc.size() > 2)
           {
             spline spl;
             spl.set_points(zCalc, thetaDefined);
             theta = spl(z);
           }
           else
           {
             theta = thetaDefined[0];
           }

           solidFraction = (1-theta)/(concRatio - theta);

        }

        else
        {
          // At the liquidus
          theta = 1;
          solidFraction = 0;
        }

        enthalpyAnalytic[dit](iv) = (1-solidFraction)*a_parameters.stefan + theta;
        bulkConcentrationAnalytic[dit](iv) = a_parameters.bcValBulkConcentrationLo[1];

      }

    } //end loop over intvects in box

  } // end loop over boxes on level

}


#include "NamespaceFooter.H"
