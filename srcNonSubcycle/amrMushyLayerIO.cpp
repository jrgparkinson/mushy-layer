#include "amrMushyLayer.H"

#include "AMRIO.H"
#include "CH_HDF5.H"
#include "ParmParse.H"

#include<sys/stat.h>
#include<sys/types.h>

void amrMushyLayer::
writeSolutionToTextFile(int lev)
{
	if (s_verbosity >= 3)
	{
		pout() << "amrMushyLayer:writeSolutionToTextFile()" << endl;
	}

	ofstream output;
	char (outputFile[100]);
	sprintf(outputFile, "finalSoln%d.data", m_ncells[1]);
	output.open(outputFile);

	for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
	{
		Box a_box = m_amrGrids[lev][dit()];
		BoxIterator bit(a_box);

		FArrayBox& theta = (*m_scalarNew[m_theta][lev])[dit()];
		FArrayBox& solidFrac = (*m_scalarNew[m_solidFraction][lev])[dit()];

		for (bit.begin(); bit.ok(); ++bit)
		{
			IntVect iv = bit();
			Real z = (iv[1]+0.5)*m_amrDx[lev];
			if (iv[0] == m_ncells[0]/2)
			{
				output << z <<"," << theta.get(iv, 0) << "," << solidFrac.get(iv, 0) << endl;
			}
		}
	}

	output.close();
}

void amrMushyLayer::
logMessage(int priority, string message)
{
	CH_TIME("amrMushyLayer::logMessage()");

	if (s_verbosity >= priority)
	{
		pout() << message << endl;
	}
}

#ifdef CH_USE_HDF5


void
amrMushyLayer::
writePlotHeader(HDF5Handle& a_handle) const
{
	if (s_verbosity >= 3)
	{
		pout() << "AMRNavierStokes::writePlotHeader " << endl;
	}

	HDF5HeaderData header;

	int numcomp = m_numVars + (m_numVectorVars*SpaceDim);

	header.m_int ["num_components"] = numcomp;
	char comp_str[30];
	int comp=0;

	//First do scalar vars
	for (int a_var=0; a_var < m_numVars; a_var++)
	{
		sprintf(comp_str, "component_%d", comp);
		header.m_string[comp_str] = m_varNames[a_var];
		comp++;
	}

	//Now vector vars
	Vector<string> dirString(SpaceDim);
	dirString[0] = string("x");
	if (SpaceDim == 2)
	{
		dirString[1] = string("z");
	}
	else if(SpaceDim == 3)
	{
		dirString[1] = string("y");
		dirString[2] = string("z");
	}

	for (int a_var=0; a_var < m_numVectorVars; a_var++)
	{
		for (int dir=0; dir<SpaceDim; dir++)
		{
			sprintf(comp_str, "component_%d", comp);
			header.m_string[comp_str] = dirString[dir]+m_varNames[a_var];
			comp++;
		}
	}


	header.writeToFile(a_handle);

	if (s_verbosity >= 3)
	{
		pout () << header << endl;
	}
}

void amrMushyLayer::getVarNames(Vector<string>& varNames)
{
	//First do scalar vars
	for (int a_var = 0; a_var<m_numVars; a_var++)
	{
		varNames[a_var] = m_varNames[a_var];
	}

	//Then vector vars
	Interval vecSrcComps(0, SpaceDim-1);
	for (int a_var = 0; a_var<m_numVectorVars; a_var++)
	{
		int startComp = m_numVars + (a_var*2);
		Vector<string> dirString(SpaceDim);
		dirString[0] = string("x");
		if (SpaceDim == 2)
		{
			dirString[1] = string("y");
		}
		else if(SpaceDim == 3)
		{
			dirString[1] = string("y");
			dirString[2] = string("z");
		}

		for (int dir = 0; dir<SpaceDim; dir++)
		{
			//varNames[startComp+dir] = m_vectorVarNames[a_var] + string(" ") + dirString[dir];
			varNames[startComp+dir] = dirString[dir] + m_vectorVarNames[a_var];
		}
	}
}


void
amrMushyLayer::writePlotFile()
{
	writePlotFile(-1);
}

/// write hdf5 plotfile to the standard location
void
amrMushyLayer::writePlotFile(int iteration)
{

	logMessage(3, "amrMushyLayer::writePlotFile" );

	int numPlotComps = m_numVars + (m_numVectorVars*SpaceDim);

	Box domain = m_amrDomains[0].domainBox();
	Real dt = 1.;
	int numLevels = m_finest_level +1;

	//variable names
	Vector<string> varNames(numPlotComps);
	getVarNames(varNames);

	// compute plot data
	Vector<LevelData<FArrayBox>* > plotData(m_scalarNew[0].size(), nullptr);
	for (int lev=0; lev<numLevels; lev++)
	{
		// first allocate storage
		plotData[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
				numPlotComps,
				m_ghostVect);

		//First do scalar vars
		for (int a_var = 0; a_var<m_numVars; a_var++)
		{
			Interval srcComps = m_scalarNew[a_var][lev]->interval();
			Interval destComps(a_var, a_var);
			m_scalarNew[a_var][lev]->copyTo(srcComps, *(plotData[lev]), destComps);
		}

		//Then vector vars
		Interval vecSrcComps(0, SpaceDim-1);
		for (int a_var = 0; a_var<m_numVectorVars; a_var++)
		{
			int startComp = m_numVars + (a_var*2);
			Interval destComps(startComp, startComp + SpaceDim-1);
			m_vectorNew[a_var][lev]->copyTo(vecSrcComps, *(plotData[lev]), destComps);
		}


	} // end loop over levels for computing plot data

	//Print out all var names and their index for matlab
//	if (m_time == 0)
//	{
//		for (int i=0; i<varNames.size(); i++)
//		{
//			pout() << varNames[i] << " = " << i+1 << ";\n";
//		}
//	}

	string iteration_folder = "";

	// usually, name the file by the timestep
	int file_index = m_cur_step;
	//If this is the plotfile during an iteration, put it in a separate folder
	if (iteration > -1 && iteration % m_iteration_plot_interval == 0)
	{
		char iterationFolder[100];
		sprintf(iterationFolder, "iter_%d/", m_cur_step);

		char newDir[100];
		sprintf(newDir, "%s/iter_%d", m_output_folder.c_str(), m_cur_step);
		mkdir(newDir, 0777);
		iteration_folder = string(iterationFolder);

		//If we're plotting at iteration, name the file on this
		file_index = iteration;
	}

	// generate plotfile name
	char iter_str[500];
	if (file_index < 10)
	{
		sprintf(iter_str, "%s/%s%s0000%d.%dd.hdf5", m_output_folder.c_str(), iteration_folder.c_str(), m_plot_prefix.c_str(),
				file_index, SpaceDim);
	}
	else if (file_index < 100)
	{
		sprintf(iter_str, "%s/%s%s000%d.%dd.hdf5", m_output_folder.c_str(), iteration_folder.c_str(), m_plot_prefix.c_str(),
				file_index, SpaceDim);
	}
	else if (file_index < 1000)
	{
		sprintf(iter_str, "%s/%s%s00%d.%dd.hdf5", m_output_folder.c_str(), iteration_folder.c_str(), m_plot_prefix.c_str(),
				file_index, SpaceDim);
	}
	else if (file_index < 10000)
	{
		sprintf(iter_str, "%s/%s%s0%d.%dd.hdf5", m_output_folder.c_str(), iteration_folder.c_str(), m_plot_prefix.c_str(),
				file_index, SpaceDim);
	}
	else
	{
		sprintf(iter_str, "%s/%s%s%d.%dd.hdf5", m_output_folder.c_str(), iteration_folder.c_str(), m_plot_prefix.c_str(),
				file_index, SpaceDim);
	}

	if (m_printAnalyticSoln)
	{
		sprintf(iter_str, "%s/%s-analytic.%dd.hdf5", m_output_folder.c_str(), m_plot_prefix.c_str(),
				SpaceDim);
	}


	string filename(iter_str);
	//	logMessage(1, iter_str);
	WriteAMRHierarchyHDF5(filename, m_amrGrids, plotData, varNames,
			domain, m_amrDx[0], dt, m_time, m_refinement_ratios,
			numLevels);

	// need to delete plotData
	for (int lev=0; lev<numLevels; lev++)
	{
		if (plotData[lev] != nullptr)
		{
			delete plotData[lev];
			plotData[lev] = nullptr;
		}
	}
}

/// write checkpoint file out for later restarting
void
amrMushyLayer::writeCheckpointFile() const
{

	if (s_verbosity > 3)
	{
		pout() << "amrMushyLayer::writeCheckpointfile" << endl;
	}

#ifdef CH_USE_HDF5

	// generate checkpointfile name
	char (iter_str[200]);

	if (m_cur_step < 10)
	{
		sprintf(iter_str, "%s/%s000%d.%dd.hdf5", m_output_folder.c_str(), m_check_prefix.c_str(),
				m_cur_step, SpaceDim);
	}
	else if (m_cur_step < 100)
	{
		sprintf(iter_str, "%s/%s00%d.%dd.hdf5", m_output_folder.c_str(), m_check_prefix.c_str(),
				m_cur_step, SpaceDim);
	}
	else if (m_cur_step < 1000)
	{
		sprintf(iter_str, "%s/%s0%d.%dd.hdf5", m_output_folder.c_str(), m_check_prefix.c_str(),
				m_cur_step, SpaceDim);
	}
	else
	{
		sprintf(iter_str, "%s/%s%d.%dd.hdf5", m_output_folder.c_str(), m_check_prefix.c_str(),
				m_cur_step, SpaceDim);
	}

	// overwrite the same checkpoint file, rather than re-writing them
	//	sprintf(iter_str, "%s.%d.hdf5", m_check_prefix.c_str(), SpaceDim);

	if (s_verbosity > 3)
	{
		pout() << "checkpoint file name = " << iter_str << endl;
	}

	HDF5Handle handle(iter_str, HDF5Handle::CREATE);

	// write amr data -- only dump out things which are essential
	// to restarting the computation (i.e. max_level, finest_level,
	// time, refinement ratios, etc.).  Other paramters (regrid
	// intervals, block-factor, etc can be changed by the inputs
	// file of the new run.
	// At the moment, the maximum level is not allowed to change,
	// although in principle, there is no real reason why it couldn't
	//
	HDF5HeaderData header;
	header.m_int["max_level"] = m_max_level;
	header.m_int["finest_level"] = m_finest_level;
	header.m_int["current_step"] = m_cur_step;
	header.m_real["time"] = m_time;
	header.m_real["dt"] = m_dt;
	header.m_int["m_numVars"] = m_numVars;
	header.m_int["m_numVectorVars"] = m_numVectorVars;
	// at the moment, save cfl, but it can be changed by the inputs
	// file if desired.
	header.m_real["cfl"] = m_cfl;

	// periodicity info
	D_TERM(
			if (m_amrDomains[0].isPeriodic(0))
				header.m_int["is_periodic_0"] = 1;
			else
				header.m_int["is_periodic_0"] = 0; ,

				if (m_amrDomains[0].isPeriodic(1))
					header.m_int["is_periodic_1"] = 1;
				else
					header.m_int["is_periodic_1"] = 0; ,

					if (m_amrDomains[0].isPeriodic(2))
						header.m_int["is_periodic_2"] = 1;
					else
						header.m_int["is_periodic_2"] = 0;
	);


	// set up component names
	char compStr[30];
	string compName;
	for (int a_var=0; a_var < m_numVars; a_var++)
	{
		sprintf(compStr, "var_%d", a_var);
		header.m_string[compStr] = m_varNames[a_var];
	}

	header.writeToFile(handle);

	// now loop over levels and write out each level's data
	// note that we loop over all allowed levels, even if they
	// are not defined at the moment.
	for (int lev=0; lev<= m_max_level; lev++)
	{
		// set up the level string
		char levelStr[20];
		sprintf(levelStr, "%d", lev);
		const std::string label = std::string("level_") + levelStr;

		handle.setGroup(label);

		// set up the header info
		HDF5HeaderData levelHeader;
		if (lev < m_max_level)
		{
			levelHeader.m_int["ref_ratio"] = m_refinement_ratios[lev];
		}
		else
		{
			levelHeader.m_int["ref_ratio"] = 2;
		}
		levelHeader.m_real["dx"] = m_amrDx[lev];
		levelHeader.m_box["prob_domain"] = m_amrDomains[lev].domainBox();


		levelHeader.writeToFile(handle);

		// now write the data for this level
		// only try to write data if level is defined.
		if (lev <= m_finest_level)
		{
			write(handle, m_amrGrids[lev]);

			for (int a_var=0; a_var<m_numVars; a_var++)
			{
				write(handle, *(m_scalarNew[a_var][lev]), m_varNames[a_var]);
			}

			for (int a_var=0; a_var<m_numVectorVars; a_var++)
			{
				write(handle, *(m_vectorNew[a_var][lev]), m_vectorVarNames[a_var]);
			}

			write(handle, *(m_fluidAdv[lev]), "FluidAdvection");
			write(handle, *(m_frameAdv[lev]), "FrameAdvection");
		}
	}// end loop over levels

	handle.close();
#endif
}


/// read checkpoint file for restart
void
amrMushyLayer::readCheckpointFile(HDF5Handle& a_handle)
{

	if (s_verbosity > 3)
	{
		pout() << "amrMushyLayer::readCheckpointFile" << endl;
	}

#ifndef CH_USE_HDF5
	MayDay::Error("code must be compiled with HDF5 to read checkpoint files");
#endif

#ifdef CH_USE_HDF5
	HDF5HeaderData header;
	header.readFromFile(a_handle);

	if (s_verbosity >= 3)
	{
		pout() << "hdf5 header data: " << endl;
		pout() << header << endl;
	}

	// read max level
	if (header.m_int.find("max_level") == header.m_int.end())
	{
		MayDay::Error("checkpoint file does not contain max_level");
	}
	// we can change max level upon restart
	int max_level_check = header.m_int["max_level"];
	if (max_level_check != m_max_level)
	{
		if (s_verbosity > 0)
		{
			pout() << "Restart file has a different max level than inputs file"
					<< endl;
			pout() << "     max level from inputs file = "
					<< m_max_level << endl;
			pout() << "     max level in checkpoint file = "
					<< max_level_check << endl;
			pout() << "Using max level from inputs file" << endl;
		}
	}
	// read finest level
	if (header.m_int.find("finest_level") == header.m_int.end())
	{
		MayDay::Error("checkpoint file does not contain finest_level");
	}

	m_finest_level = header.m_int["finest_level"];

	if (m_finest_level > m_max_level)
	{
		MayDay::Error("finest level in restart file > max allowable level!");
	}

	// read current step
	if (header.m_int.find("current_step") == header.m_int.end())
	{
		MayDay::Error("checkpoint file does not contain current_step");
	}

	m_cur_step = header.m_int["current_step"];
	m_restart_step = m_cur_step;

	// read time
	if (header.m_real.find("time") == header.m_real.end())
	{
		MayDay::Error("checkpoint file does not contain time");
	}

	m_time = header.m_real["time"];

	// read timestep
	if (header.m_real.find("dt") == header.m_real.end())
	{
		MayDay::Error("checkpoint file does not contain dt");
	}

	m_dt = header.m_real["dt"];

	// read num scalar vars
	if (header.m_int.find("m_numVars") == header.m_int.end())
	{
		MayDay::Error("checkpoint file does not contain m_numVars");
	}
//	int nVars  = header.m_int["m_numVars"];

	// read num vector vars
	if (header.m_int.find("m_numVectorVars") == header.m_int.end())
	{
		MayDay::Error("checkpoint file does not contain m_numVectorVars");
	}
	m_numVectorVars = header.m_int["m_numVectorVars"];


	// read cfl
	if (header.m_real.find("cfl") == header.m_real.end())
	{
		MayDay::Error("checkpoint file does not contain cfl");
	}

	Real check_cfl = header.m_real["cfl"];
	ParmParse ppCheck("main");

	if (ppCheck.contains("cfl"))
	{
		// check for consistency and warn if different
		if (check_cfl != m_cfl)
		{
			if (s_verbosity > 0)
			{
				pout() << "CFL in checkpoint file different from inputs file"
						<< endl;
				pout() << "     cfl in inputs file = " << m_cfl << endl;
				pout() << "     cfl in checkpoint file = " << check_cfl
						<< endl;
				pout() << "Using cfl from inputs file" << endl;
			}
		}  // end if cfl numbers differ
	} // end if cfl present in inputs file
	else
	{
		m_cfl = check_cfl;
	}

	// read periodicity info
	// Get the periodicity info -- this is more complicated than it really
	// needs to be in order to preserve backward compatibility
	bool isPeriodic[SpaceDim];
	D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
		isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
	else
		isPeriodic[0] = false; ,

		if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
			isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
		else
			isPeriodic[1] = false; ,

			if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
				isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
			else
				isPeriodic[2] = false;);

	// now resize stuff
	m_amrDomains.resize(m_max_level+1);
	m_amrGrids.resize(m_max_level+1);
	m_amrDx.resize(m_max_level+1);

	m_frameAdv.resize(m_max_level+1, nullptr);
	m_fluidAdv.resize(m_max_level+1, nullptr);

	m_fluxRegister.resize(m_numVars);


	m_scalarOld.resize(m_numVars);
	m_scalarNew.resize(m_numVars);
	m_dScalar.resize(m_numVars);

	m_vectorOld.resize(m_numVectorVars);
	m_vectorNew.resize(m_numVectorVars);
	m_dVector.resize(m_numVectorVars);


	for (int a_var=0; a_var<m_numVars; a_var++)
	{
		m_scalarOld[a_var].resize(m_max_level+1);
		m_scalarNew[a_var].resize(m_max_level+1);
		m_dScalar[a_var].resize(m_max_level+1);
		m_fluxRegister[a_var].resize(m_max_level+1);
	}

	for (int a_var=0; a_var<m_numVectorVars; a_var++)
	{
		m_vectorOld[a_var].resize(m_max_level+1);
		m_vectorNew[a_var].resize(m_max_level+1);
		m_dVector[a_var].resize(m_max_level+1);
	}

	// now read in level-by-level data
	for (int lev=0; lev<= m_max_level; lev++)
	{
		// set up the level string
		char levelStr[20];
		sprintf(levelStr, "%d", lev);
		const std::string label = std::string("level_") + levelStr;

		a_handle.setGroup(label);

		// read header info
		HDF5HeaderData header;
		header.readFromFile(a_handle);

		if (s_verbosity >= 3)
		{
			pout() << "level " << lev << " header data" << endl;
			pout() << header << endl;
		}

		// Get the refinement ratio
		if (lev < m_max_level)
		{
			int checkRefRatio;
			if (header.m_int.find("ref_ratio") == header.m_int.end())
			{
				MayDay::Error("checkpoint file does not contain ref_ratio");
			}
			checkRefRatio = header.m_int["ref_ratio"];

			// check for consistency
			if (checkRefRatio != m_refinement_ratios[lev])
			{
				MayDay::Error("inputs file and checkpoint file ref ratios inconsistent");
			}
		}

		// read dx
		if (header.m_real.find("dx") == header.m_real.end())
		{
			MayDay::Error("checkpoint file does not contain dx");
		}

		m_amrDx[lev] = header.m_real["dx"];

		// read problem domain box
		if (header.m_box.find("prob_domain") == header.m_box.end())
		{
			MayDay::Error("checkpoint file does not contain prob_domain");
		}
		Box domainBox = header.m_box["prob_domain"];

		m_amrDomains[lev] = ProblemDomain(domainBox, isPeriodic);


		// the rest is only applicable if this level is defined
		if (lev <= m_finest_level)
		{
			// read grids
			Vector<Box> grids;
			const int grid_status = read(a_handle, grids);
			if (grid_status != 0)
			{
				MayDay::Error("checkpoint file does not contain a Vector<Box>");
			}
			// do load balancing
			int numGrids = grids.size();
			Vector<int> procIDs(numGrids);
			LoadBalance(procIDs, grids);
			DisjointBoxLayout levelDBL(grids, procIDs, m_amrDomains[lev]);
			m_amrGrids[lev] = levelDBL;

			getFrameAdvection(lev);

			IntVect ghostVect(m_num_ghost*IntVect::Unit);


			//read this level's scalar data
			for (int a_var=0; a_var<m_numVars; a_var++)
			{
				// allocate storage
				m_scalarOld[a_var][lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(levelDBL, 1,
						ghostVect));
				m_scalarNew[a_var][lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(levelDBL, 1,
						ghostVect));
				m_dScalar[a_var][lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(levelDBL, 1,
						ghostVect));


				LevelData<FArrayBox>& new_data = *m_scalarNew[a_var][lev];
				const int dataStatus = read<FArrayBox>(a_handle,
						new_data,
						m_varNames[a_var],
						levelDBL);

				if (dataStatus != 0)
				{
					char err[80];
					strcpy(err, "checkpoint file does not contain data for var ");
					strcat(err, m_varNames[a_var].c_str());
					MayDay::Error(err);
				}


			} //finished reading scalar data

			//read this level's vector data
			for (int a_var=0; a_var<m_numVectorVars; a_var++)
			{

				//allocate storage
				m_vectorOld[a_var][lev] = new LevelData<FArrayBox>(levelDBL, SpaceDim, ghostVect);
				m_vectorNew[a_var][lev] = new LevelData<FArrayBox>(levelDBL, SpaceDim, ghostVect);
				m_dVector[a_var][lev] = new LevelData<FArrayBox>(levelDBL, SpaceDim, ghostVect);

				LevelData<FArrayBox>& new_vector_data = *m_vectorNew[a_var][lev];
				const int dataStatus = read<FArrayBox>(a_handle,
						new_vector_data,
						m_vectorVarNames[a_var],
						levelDBL);

				if (dataStatus != 0)
				{
					char err[80];
					strcpy(err, "checkpoint file does not contain data for var ");
					strcat(err, m_vectorVarNames[a_var].c_str());
					MayDay::Error(err);
				}
			}

			//Read advection fluxboxes
			m_frameAdv[lev] = new LevelData<FluxBox>(levelDBL, SpaceDim, m_ghostVect);
			m_fluidAdv[lev] = new LevelData<FluxBox>(levelDBL, SpaceDim, m_ghostVect);

			LevelData<FluxBox>& new_fluidAdv_data = *m_fluidAdv[lev];
			const int dataStatus = read<FluxBox>(a_handle,
					new_fluidAdv_data,
					"FluidAdvection",
					levelDBL);

			if (dataStatus != 0)
			{
				MayDay::Error("Couldn't read in fluid advection data");
			}


			const int dataStatusFrame = read<FluxBox>(a_handle,
					*m_frameAdv[lev],
					"FrameAdvection",
					levelDBL);

			if (dataStatusFrame != 0)
			{
				MayDay::Error("Couldn't read in frame advection data");
			}

		} // end if this level is defined
	} // end loop over levels

	a_handle.close();
#endif

}

/// set up for restart
void
amrMushyLayer::restart(string& a_restart_file)
{
	if (s_verbosity > 3)
	{
		pout() << "amrMushyLayer::restart" << endl;
	}

	HDF5Handle handle(a_restart_file, HDF5Handle::OPEN_RDONLY);
	// first read in data from checkpoint file
	readCheckpointFile(handle);

	// don't think I need to do anything else, do I?


}

#endif

