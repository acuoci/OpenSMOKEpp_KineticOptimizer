/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2018  Alberto Cuoci                                      |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_BatchReactorExperiment_H
#define OpenSMOKE_BatchReactorExperiment_H

// Utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/cema/OnTheFlyCEMA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// Standard batch reactors
#include "idealreactors/batch/BatchReactor"
#include "Grammar_BatchReactorExperiment.h"

// Optimization rules
#include "OptimizationRules_BatchReactorExperiment.h"

namespace OpenSMOKE
{
	
	class BatchReactorExperiment
	{
	public:

		void Setup(const std::string input_file_name,
			OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
			OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML);

		void Solve(const bool verbose = false);
		double Solve(const double x);

		double norm2_abs_error() const { return norm2_abs_error_; }
		double norm2_rel_error() const { return norm2_rel_error_; }

		const OpenSMOKE::OptimizationRules_BatchReactorExperiment* optimization() const { return optimization_; }

	private:

		OpenSMOKE::BatchReactor_Type type_;
		OpenSMOKE::BatchReactor_NonIsothermal_ConstantPressure* batch_nonisothermal_constantp_;
		OpenSMOKE::BatchReactor_Isothermal_ConstantPressure* batch_isothermal_constantp_;
		OpenSMOKE::BatchReactor_NonIsothermal_ConstantVolume* batch_nonisothermal_constantv_;
		OpenSMOKE::BatchReactor_Isothermal_ConstantVolume* batch_isothermal_constantv_;
		OpenSMOKE::BatchReactor_NonIsothermal_UserDefinedVolume* batch_nonisothermal_userdefinedv_;


		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML_;

		OpenSMOKE::OptimizationRules_BatchReactorExperiment*	optimization_;
		OpenSMOKE::PolimiSoot_Analyzer*							polimi_soot_;
		OpenSMOKE::OnTheFlyPostProcessing*						on_the_fly_post_processing_;
		OpenSMOKE::OnTheFlyROPA*								onTheFlyROPA_;
		OpenSMOKE::BatchReactor_Options*						batch_options_;
		OpenSMOKE::ODE_Parameters*								ode_parameters_;
		OpenSMOKE::SensitivityAnalysis_Options*					sensitivity_options_;
		OpenSMOKE::IgnitionDelayTimes_Analyzer*					idt;
		OpenSMOKE::OnTheFlyCEMA*								onTheFlyCEMA;

		double tStart_;
		double tEnd_;

		double norm2_abs_error_;
		double norm2_rel_error_;
	};
}

#include "BatchReactorExperiment.hpp"

#endif // OpenSMOKE_BatchReactorExperiment_H