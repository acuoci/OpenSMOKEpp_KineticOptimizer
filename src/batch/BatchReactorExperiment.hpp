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

namespace OpenSMOKE
{
	void BatchReactorExperiment::Setup(const std::string input_file_name,
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML)
	{
		// Pointers
		thermodynamicsMapXML_ = thermodynamicsMapXML;
		kineticsMapXML_ = kineticsMapXML;

		// Defines the grammar rules
		OpenSMOKE::Grammar_BatchReactorExperiment grammar_BatchReactor;

		// Define the dictionaries
		std::string main_dictionary_name_ = "BatchReactor";
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(main_dictionary_name_).SetGrammar(grammar_BatchReactor);



		// Read thermodynamics and kinetics maps
		double T, P_Pa;
		OpenSMOKE::OpenSMOKEVectorDouble omega;
		tEnd_ = 0.;
		tStart_ = 0.;                       // default 0
		double volume = 1.;					// default value [1 m3]

		// Read initial conditions
		{
			std::string name_of_gas_status_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@InitialStatus") == true)
				dictionaries(main_dictionary_name_).ReadDictionary("@InitialStatus", name_of_gas_status_subdictionary);

			GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary), *thermodynamicsMapXML, T, P_Pa, omega);
		}

		// Read end time
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@EndTime") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@EndTime", value, units);
				if (units == "s")		  tEnd_ = value;
				else if (units == "ms")   tEnd_ = value / 1000.;
				else if (units == "min")  tEnd_ = value * 60.;
				else if (units == "h")    tEnd_ = value * 3600.;
				else OpenSMOKE::FatalErrorMessage("Unknown time units");
			}
		}

		// Read start time
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@StartTime") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@StartTime", value, units);
				if (units == "s")         tStart_ = value;
				else if (units == "ms")   tStart_ = value / 1000.;
				else if (units == "min")  tStart_ = value * 60.;
				else if (units == "h")    tStart_ = value * 3600.;
				else OpenSMOKE::FatalErrorMessage("Unknown time units");
			}
		}

		// Read volume
		{
			if (dictionaries(main_dictionary_name_).CheckOption("@Volume") == true)
			{
				double value;
				std::string units;

				dictionaries(main_dictionary_name_).ReadMeasure("@Volume", value, units);
				if (units == "m3")		  volume = value;
				else if (units == "dm3")  volume = value / 1.e3;
				else if (units == "cm3")  volume = value / 1.e6;
				else if (units == "mm3")  volume = value / 1.e9;
				else if (units == "l")    volume = value / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown volume units");
			}
		}

		// Read exchange area
		double exchange_area = 0.;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@ExchangeArea") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@ExchangeArea", value, units);
				if (units == "m2")        exchange_area = value;
				else if (units == "dm2")  exchange_area = value / 1.e2;
				else if (units == "cm2")  exchange_area = value / 1.e4;
				else if (units == "mm2")  exchange_area = value / 1.e6;
				else OpenSMOKE::FatalErrorMessage("Unknown area units");
			}
		}

		// Read global thermal exchange coefficient
		double global_thermal_exchange_coefficient = 0.;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@GlobalThermalExchangeCoefficient") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@GlobalThermalExchangeCoefficient", value, units);
				if (units == "W/m2/K")			global_thermal_exchange_coefficient = value;
				else if (units == "W/m2/C")		global_thermal_exchange_coefficient = value;
				else if (units == "kcal/m2/K")		global_thermal_exchange_coefficient = value * 4186.8;
				else if (units == "kcal/m2/C")		global_thermal_exchange_coefficient = value * 4186.8;
				else OpenSMOKE::FatalErrorMessage("Unknown global thermal exchange coefficient units");
			}
		}

		// Environment temperature
		double T_environment = T;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@EnvironmentTemperature") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@EnvironmentTemperature", value, units);
				if (units == "K")		T_environment = value;
				else if (units == "C")	T_environment = value + 273.15;
				else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
			}
		}


		//Type
		{
			std::string value;
			if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
			{
				dictionaries(main_dictionary_name_).ReadString("@Type", value);
				if (value == "Isothermal-ConstantVolume")				type_ = OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV;
				else if (value == "Isothermal-ConstantPressure")		type_ = OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP;
				else if (value == "NonIsothermal-ConstantVolume")		type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV;
				else if (value == "NonIsothermal-ConstantPressure")		type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP;
				else if (value == "NonIsothermal-UserDefinedVolume")	type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME;
				else OpenSMOKE::FatalErrorMessage("Unknown batch reactor type: " + value);

			}
		}

		// Options
		OpenSMOKE::BatchReactor_Options* batch_options;
		{
			batch_options = new OpenSMOKE::BatchReactor_Options();
			if (dictionaries(main_dictionary_name_).CheckOption("@Options") == true)
			{
				std::string name_of_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@Options", name_of_options_subdictionary);
				batch_options->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			}
		}

		// ODE Parameters
		OpenSMOKE::ODE_Parameters*	ode_parameters;
		{
			ode_parameters = new OpenSMOKE::ODE_Parameters();
			if (dictionaries(main_dictionary_name_).CheckOption("@OdeParameters") == true)
			{
				std::string name_of_ode_parameters_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@OdeParameters", name_of_ode_parameters_subdictionary);
				ode_parameters->SetupFromDictionary(dictionaries(name_of_ode_parameters_subdictionary));
			}
		}

		// Sensitivity Options
		OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;
		{
			sensitivity_options = new OpenSMOKE::SensitivityAnalysis_Options();
			if (dictionaries(main_dictionary_name_).CheckOption("@SensitivityAnalysis") == true)
			{
				std::string name_of_sensitivity_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@SensitivityAnalysis", name_of_sensitivity_options_subdictionary);

				batch_options->SetSensitivityAnalysis(true);
				sensitivity_options->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
			}
		}

		// On the fly ROPA
		OpenSMOKE::OnTheFlyROPA* onTheFlyROPA = new OpenSMOKE::OnTheFlyROPA(*thermodynamicsMapXML, *kineticsMapXML);
		if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyROPA") == true)
		{
			std::string name_of_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyROPA", name_of_options_subdictionary);
			// No ROPA (disabled)
			//onTheFlyROPA->SetupFromDictionary(dictionaries(name_of_options_subdictionary), path_kinetics_output);
		}

		// On the fly CEMA
		onTheFlyCEMA = new OpenSMOKE::OnTheFlyCEMA(*thermodynamicsMapXML, *kineticsMapXML, batch_options->output_path());
		if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyCEMA") == true)
		{
			std::string name_of_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyCEMA", name_of_options_subdictionary);
			onTheFlyCEMA->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
		}

		// On the fly PostProcessing
		OpenSMOKE::OnTheFlyPostProcessing*		on_the_fly_post_processing;
		{
			on_the_fly_post_processing = new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML, *kineticsMapXML, batch_options->output_path());

			if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyPostProcessing") == true)
			{
				std::string name_of_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyPostProcessing", name_of_options_subdictionary);
				on_the_fly_post_processing->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			}
		}

		// Ignition Delay Times
		{
			idt = new OpenSMOKE::IgnitionDelayTimes_Analyzer();
			if (dictionaries(main_dictionary_name_).CheckOption("@IgnitionDelayTimes") == true)
			{
				if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV ||
					type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
					OpenSMOKE::FatalErrorMessage("The @IgnitionDelayTimes can be used only for NonIsothermal reactors");

				std::string name_of_idt_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@IgnitionDelayTimes", name_of_idt_subdictionary);
				idt->SetupFromDictionary(dictionaries(name_of_idt_subdictionary), *thermodynamicsMapXML);
			}
		}

		OpenSMOKE::BatchReactor_VolumeProfile* batchreactor_volumeprofile;
	{
			// Read pressure coefficient
			if (dictionaries(main_dictionary_name_).CheckOption("@PressureCoefficient") == true)
			{
				batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();

				if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile can be used only for NonIsothermal-UserDefinedVolume reactors");

				double value;
				std::string units;

				dictionaries(main_dictionary_name_).ReadMeasure("@PressureCoefficient", value, units);
				if (units == "Pa/s")		 batchreactor_volumeprofile->SetPressureCoefficient(value);
				else if (units == "bar/s")   batchreactor_volumeprofile->SetPressureCoefficient(value*1.e5);
				else if (units == "atm/s")   batchreactor_volumeprofile->SetPressureCoefficient(value*101325.);
				else if (units == "Pa/ms")	 batchreactor_volumeprofile->SetPressureCoefficient(value*1000.);
				else if (units == "bar/ms")  batchreactor_volumeprofile->SetPressureCoefficient(value*1.e5*1000.);
				else if (units == "atm/ms")  batchreactor_volumeprofile->SetPressureCoefficient(value*101325.*1000.);
				else OpenSMOKE::FatalErrorMessage("Unknown pressure coefficient units. Available: Pa/s || Pa/ms || bar/s || bar/ms || atm/s || atm/ms");
			}

			// Read volume profile
			if (dictionaries(main_dictionary_name_).CheckOption("@VolumeProfile") == true)
			{
				batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();

				std::string name_of_profile_subdictionary;

				if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile can be used only for NonIsothermal-UserDefinedVolume reactors");

				dictionaries(main_dictionary_name_).ReadDictionary("@VolumeProfile", name_of_profile_subdictionary);

				OpenSMOKE::OpenSMOKEVectorDouble x, y;
				std::string x_variable, y_variable;
				GetXYProfileFromDictionary(dictionaries(name_of_profile_subdictionary), x, y, x_variable, y_variable);

				if (x_variable != "time")
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile must be defined versus the time");
				if (y_variable != "volume")
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile must be defined versus the volume");

				if (y[1] != volume)
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile and the @Volume options must be consistent");

				batchreactor_volumeprofile->SetProfile(x, y);
			}

			// Read pressure profile
			if (dictionaries(main_dictionary_name_).CheckOption("@PressureProfile") == true)
			{
				batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();

				std::string name_of_profile_subdictionary;

				if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
					OpenSMOKE::FatalErrorMessage("The @PressureProfile can be used only for NonIsothermal-UserDefinedVolume reactors");

				dictionaries(main_dictionary_name_).ReadDictionary("@PressureProfile", name_of_profile_subdictionary);

				OpenSMOKE::OpenSMOKEVectorDouble x, y;
				std::string x_variable, y_variable;
				GetXYProfileFromDictionary(dictionaries(name_of_profile_subdictionary), x, y, x_variable, y_variable);

				if (x_variable != "time")
					OpenSMOKE::FatalErrorMessage("The @PressureProfile must be defined versus the time");
				if (y_variable != "pressure")
					OpenSMOKE::FatalErrorMessage("The @PressureProfile must be defined versus the volume");

				if (y[1] != P_Pa)
					OpenSMOKE::FatalErrorMessage("The @PressureProfile and the initial pressure of mixture must be consistent");

				batchreactor_volumeprofile->SetPressureProfile(x, y);
			}
		}

		// Polimi soot
		OpenSMOKE::PolimiSoot_Analyzer* polimi_soot = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML);
		{
			std::string name_of_polimisoot_analyzer_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@PolimiSoot") == true)
			{
				dictionaries(main_dictionary_name_).ReadDictionary("@PolimiSoot", name_of_polimisoot_analyzer_subdictionary);
				polimi_soot->SetupFromDictionary(dictionaries(name_of_polimisoot_analyzer_subdictionary));
			}
		}

		// Optimization
		{
			optimization_ = new OpenSMOKE::OptimizationRules_BatchReactorExperiment;

			std::string name_of_optimization_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@Optimization") == true)
			{
				dictionaries(main_dictionary_name_).ReadDictionary("@Optimization", name_of_optimization_subdictionary);
				optimization_->SetupFromDictionary(dictionaries(name_of_optimization_subdictionary), dictionaries);
			}
		}


		// Solve the ODE system: NonIsothermal, Constant Volume
		if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
		{
			batch_nonisothermal_constantv_ = new OpenSMOKE::BatchReactor_NonIsothermal_ConstantVolume(*thermodynamicsMapXML, *kineticsMapXML,
					*ode_parameters, *batch_options, *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot, volume, T, P_Pa, omega,
					global_thermal_exchange_coefficient, exchange_area, T_environment);
			batch_nonisothermal_constantv_->Solve(tStart_, tEnd_);
		}

		// Solve the ODE system: NonIsothermal, Volume assigned according to a specified law
		if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
		{
			batch_nonisothermal_userdefinedv_ = new OpenSMOKE::BatchReactor_NonIsothermal_UserDefinedVolume(*thermodynamicsMapXML, *kineticsMapXML,
					*ode_parameters, *batch_options, *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot, volume, T, P_Pa, omega, tStart_,
					global_thermal_exchange_coefficient, exchange_area, T_environment);

			batch_nonisothermal_userdefinedv_->SetVolumeProfile(*batchreactor_volumeprofile);

			batch_nonisothermal_userdefinedv_->Solve(tStart_, tEnd_);
		}

		// Solve the ODE system: Isothermal, Constant Volume
		if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV)
		{
			batch_isothermal_constantv_ = new OpenSMOKE::BatchReactor_Isothermal_ConstantVolume(*thermodynamicsMapXML, *kineticsMapXML,
					*ode_parameters, *batch_options, *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot, volume, T, P_Pa, omega);
			batch_isothermal_constantv_->Solve(tStart_, tEnd_);
		}

		// Solve the ODE system: NonIsothermal, Constant Pressure
		if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP)
		{
			batch_nonisothermal_constantp_ = new OpenSMOKE::BatchReactor_NonIsothermal_ConstantPressure(*thermodynamicsMapXML, *kineticsMapXML,
				*ode_parameters, *batch_options, *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot, volume, T, P_Pa, omega,
				global_thermal_exchange_coefficient, exchange_area, T_environment);
		}

		// Solve the ODE system: Isothermal, Constant Pressure
		if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
		{
			batch_isothermal_constantp_ = new OpenSMOKE::BatchReactor_Isothermal_ConstantPressure(*thermodynamicsMapXML, *kineticsMapXML,
				*ode_parameters, *batch_options, *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot, volume, T, P_Pa, omega);
		}
	}

	void BatchReactorExperiment::Solve(const bool verbose)
	{

		// Solve the ODE system: NonIsothermal, Constant Volume
		if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
		{
			batch_nonisothermal_constantv_->Solve(tStart_, tEnd_);
			
		}

		// Solve the ODE system: NonIsothermal, Volume assigned according to a specified law
		if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
		{
			batch_nonisothermal_userdefinedv_->Solve(tStart_, tEnd_);
		}

		// Solve the ODE system: Isothermal, Constant Volume
		if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV)
		{
			batch_isothermal_constantv_->Solve(tStart_, tEnd_);
		}

		// Solve the ODE system: NonIsothermal, Constant Pressure
		if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP)
		{
			batch_nonisothermal_constantp_->Solve(tStart_, tEnd_);
		}

		// Solve the ODE system: Isothermal, Constant Pressure
		if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
		{
			batch_isothermal_constantp_->Solve(tStart_, tEnd_);
		}

		double tau = idt->temperature_slope_tau();
		if (optimization_->slope_definition() == false)
			tau = idt->temperature_increase_tau();

		const double abs_error = (optimization_->tau() - tau);
		const double rel_error = (optimization_->tau() - tau) / optimization_->tau();
		norm2_abs_error_ = abs_error * abs_error;
		norm2_rel_error_ = rel_error * rel_error;

		if (verbose == true)
		{
			std::cout << "Norm2(abs): " << norm2_abs_error_ << std::endl;
			std::cout << "Norm2(rel): " << norm2_rel_error_ << std::endl;
		}

		idt->Reset();
	}

	double BatchReactorExperiment::Solve(const double x)
	{
		return 0.;
	}

}
