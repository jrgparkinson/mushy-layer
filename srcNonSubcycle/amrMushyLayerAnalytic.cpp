#include "amrMushyLayer.H"

#include "spline.h"



/*
 * Calculate the analytic solution for this problem (if possible)
 */
void
amrMushyLayer::
calculateAnalyticSolution()
{
	CH_TIME("amrMushyLayer::calculateAnalyticSolution()");

	if (s_verbosity > 3)
	{
		pout() << "amrMushyLayer::calculateAnalyticSolution" << endl;
	}


	if (m_parameters.physicalProblem == m_bmDiffusiveSolidification)
	{
		calculateSolnBm1();
	}
	else if(m_parameters.physicalProblem == m_bmHRL)
	{
		calculateSolnBm2();
	}
	else if(m_parameters.physicalProblem == m_bmDiffusion)
	{
		calculateSolnDiffusion();
	}
	else if(m_parameters.physicalProblem == m_bmAdvection)
	{
		//Do nothing
	}


}

//Just set theta = z
void amrMushyLayer::
calculateSolnDiffusion()
{
	for (int lev=0; lev<=m_finest_level; lev++)
	{
		for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
		{
			Box a_box = m_amrGrids[lev][dit()];
			BoxIterator bit(a_box);

			FArrayBox& analytictheta = (*m_scalarNew[m_thetaTrue][lev])[dit()];

			for (bit.begin(); bit.ok(); ++bit)
			{
				IntVect iv = bit();
				Real z = (iv[1]+0.5)*m_amrDx[lev];
				Real theta = m_parameters.thetaBottom +
						(m_parameters.thetaTop - m_parameters.thetaBottom) * z/m_domainSize[1];
				analytictheta.set(iv, 0, theta);
			}

			//Porosity = 1, H = S + theta
			(*m_scalarNew[m_enthalpyAnalytic][lev])[dit].setVal(m_parameters.stefan);
			(*m_scalarNew[m_enthalpyAnalytic][lev])[dit] += (*m_scalarNew[m_theta][lev])[dit];

			//This doesn't really matter
			(*m_scalarNew[m_ThetaAnalytic][lev])[dit].setVal(0);
		}
	}

}

void amrMushyLayer::
calculateSolnBm2()
{
	for (int lev=0; lev<=m_finest_level; lev++)
	{
		for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
		{
			Box a_box = (*m_scalarNew[m_thetaTrue][lev])[dit()].box();
			BoxIterator bit(a_box);

			FArrayBox& analytictheta = (*m_scalarNew[m_thetaTrue][lev])[dit()];

			for (bit.begin(); bit.ok(); ++bit)
			{
				IntVect iv = bit();
				Real z = (iv[1]+0.5)*m_amrDx[lev];
				Real theta = m_parameters.thetaBottom +
						(m_parameters.thetaTop - m_parameters.thetaBottom) * z/m_domainSize[1];
				analytictheta.set(iv, 0, theta);
			}
		}
	}

}


void amrMushyLayer::
calculateSolnBm1()
{

	// First calculate z = z(theta) for finely spaced theta
	vector<Real> thetaGrid, thetaDefined;
	vector<Real> zCalc;
	int numPoints = (int) m_domainSize[1]/m_amrDx[m_finest_level];
	numPoints = numPoints * 4; //Use an even finer grid to make interpolation more accurate

	Real thetaEutectic = m_parameters.thetaEutectic;
	Real thetaInf = m_parameters.thetaInf;
	Real thetaInterface = m_parameters.thetaInterface;
	Real deltaTheta = (thetaInf - thetaEutectic);


	//Constants needed for finding z(theta)
	Real concRatio = m_parameters.compositionRatio;

	Real h = m_parameters.directionalSolidificationMushyZ(thetaInterface);

	for (int i = 0; i <= numPoints; i++)
	{
		Real theta = thetaEutectic + deltaTheta*((float)i / (float)numPoints);
		thetaGrid.push_back(theta);

		//Calculate z(theta)
		//			Real z = (1/m_parameters.nonDimVel) * (((alpha - concRatio) / (alpha-beta)) * log((alpha)/(alpha-theta)) +
		//													 ((concRatio-beta)/(alpha-beta)) * 		log((beta)/(beta-theta)));
		Real z = m_parameters.directionalSolidificationMushyZ(theta);

		if (z < h)
		{
			zCalc.push_back(z);
			thetaDefined.push_back(theta);
		}

		//pout() << "theta = " << theta << ", z = " << zVal << "\n";
	}


	ofstream myfile;
	myfile.open("analyticSoln.data");

	if (!myfile.is_open())
	{
		MayDay::Error("Couldn't open file for writing out analytic solution");
	}


	//Now interpolate to find theta(z)
	for (int lev=0; lev<=m_finest_level; lev++)
	{
		for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
		{
//			Box a_box = m_amrGrids[lev][dit()];
			Box a_box = (*m_scalarNew[m_Tanalytic][lev])[dit()].box();
			BoxIterator bit(a_box);

			FArrayBox& analyticTemp = (*m_scalarNew[m_Tanalytic][lev])[dit()];
			FArrayBox& analyticSolidFrac = (*m_scalarNew[m_solidFractionTrue][lev])[dit()];

			for (bit.begin(); bit.ok(); ++bit)
			{
				IntVect iv = bit();
				Real z = (iv[1]+0.5)*m_amrDx[lev];

				Real theta, solidFraction, ThetaL;


				// If stefan=0 we're always in liquid
				if (z > h or m_parameters.stefan == 0)
				{
					//We're in the liquid
					theta = thetaInf + (thetaInterface - thetaInf)* exp(- m_parameters.nonDimVel * (z-h));
					ThetaL = 1;

					solidFraction = 0;
				}
				else
				{
					//First find the two z values either side of this one
					int i = 1;
					while (zCalc[i] < z)
					{
						i++;
					}

					//Spline interpolation

					tk::spline thetaSpline;
					thetaSpline.set_points(zCalc, thetaDefined);
					theta = thetaSpline(z);

					ThetaL = theta;

					solidFraction = (1-theta)/(concRatio - theta);

				}


				if (iv[0]==1)
				{
					myfile << z << ", " << theta << ", " << solidFraction << "\n";
				}

				//				thetaKatz = thetaMatlabKatz[iv[1]];
				//				solidFractionKatz = solidFracMatlab[iv[1]];
				//				ThetaLKatz = ThetaLMatlabKatz[iv[1]];


				//Convert theta to T
				Real T = m_parameters.bottomTemp +
						m_parameters.deltaTemp * theta;

				analyticTemp.set(iv, 0, T);


				analyticSolidFrac.set(iv, 0, solidFraction);
				(*m_scalarNew[m_thetaTrue][lev])[dit()].set(iv, 0, theta);
				(*m_scalarNew[m_ThetaLAnalytic][lev])[dit()].set(iv, 0, ThetaL);


			} //end loop over intvects in box

			//H = (1-phi)*S + theta
			(*m_scalarNew[m_enthalpyAnalytic][lev])[dit].setVal(1);
			(*m_scalarNew[m_enthalpyAnalytic][lev])[dit] -= (*m_scalarNew[m_solidFractionTrue][lev])[dit];
			(*m_scalarNew[m_enthalpyAnalytic][lev])[dit].mult(m_parameters.stefan);
			(*m_scalarNew[m_enthalpyAnalytic][lev])[dit] += (*m_scalarNew[m_thetaTrue][lev])[dit];

			FArrayBox& thisH = (*m_scalarNew[m_enthalpyAnalytic][lev])[dit];

			//Theta = 1
			(*m_scalarNew[m_ThetaAnalytic][lev])[dit].setVal(1);

		} // end loop over boxes on level

	}// end loop over levels

	myfile.close();
	//matlabAnalyticSoln.close();

	//	// Ensure these are set at all timesteps
	//	activeLevelCopy(m_scalarNew[m_thetaTrue], m_scalarOld[m_thetaTrue]);
	//	activeLevelCopy(m_scalarNew[m_solidFractionTrue], m_scalarOld[m_solidFractionTrue]);
	//	activeLevelCopy(m_scalarNew[m_ThetaLAnalytic], m_scalarOld[m_ThetaLAnalytic]);


}





void amrMushyLayer::enforceAnalyticSolution()
{
	if (s_verbosity > 3)
	{
		pout() << "amrMushyLayer::enforceAnalyticSolution()" << endl;
	}

	if (m_parameters.physicalProblem == m_bmDiffusiveSolidification)
	{
		//Specific to directional solidification benchmark
		activeLevelCopy(m_scalarNew[m_thetaTrue], m_scalarNew[m_theta]);
		activeLevelCopy(m_scalarNew[m_solidFractionTrue], m_scalarNew[m_solidFraction]);
		activeLevelCopy(m_scalarNew[m_ThetaLAnalytic], m_scalarNew[m_compositionLiquid]);
		activeLevelCopy(m_scalarNew[m_ThetaAnalytic], m_scalarNew[m_composition]);

		int ilev = 0;
		for (DataIterator dit = m_amrGrids[ilev].dataIterator(); dit.ok(); ++dit)
		{
			(*m_scalarNew[m_porosity][ilev])[dit()].setVal(1);
			(*m_scalarNew[m_porosity][ilev])[dit()] -= (*m_scalarNew[m_solidFraction][ilev])[dit()];
		}

		calculateEnthalpy();

		vector<int> ignore;
//		ignore.push_back(m_theta);
//		ignore.push_back(m_porosity);
//		ignore.push_back(m_compositionLiquid);
		updateEnthalpyVariables(ignore);
	}

	else if (m_parameters.physicalProblem == m_bmDiffusion)
	{
		activeLevelCopy(m_scalarNew[m_thetaTrue], m_scalarNew[m_theta]);
	}
	else if (m_parameters.physicalProblem == m_bmHRL)
	{
		activeLevelCopy(m_scalarNew[m_thetaTrue], m_scalarNew[m_theta]);

		calculateEnthalpy();

		vector<int> ignore;
		updateEnthalpyVariables(ignore);

	}

	else
	{
		//Standard

		activeLevelCopy(m_scalarNew[m_ThetaAnalytic], m_scalarNew[m_composition]);
		activeLevelCopy(m_scalarNew[m_enthalpyAnalytic], m_scalarNew[m_enthalpy]);

		// Get porosity, theta, Theta_l, Theta_s
		vector<int> ignore;
		updateEnthalpyVariables(ignore);
	}




	//This is done at the start of timestep() anyway
	//Copy new to old
	//	activeLevelCopy(m_scalarNew[m_enthalpy], m_scalarOld[m_enthalpy]);
	//	activeLevelCopy(m_scalarNew[m_theta],  m_scalarOld[m_theta]);
	//	activeLevelCopy(m_scalarNew[m_solidFraction], m_scalarOld[m_solidFraction]);
	//	activeLevelCopy(m_scalarNew[m_compositionLiquid], m_scalarOld[m_compositionLiquid]);


}
