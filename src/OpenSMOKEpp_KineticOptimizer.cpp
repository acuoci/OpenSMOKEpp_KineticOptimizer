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

// Include OpenMP Header file
#if defined(_OPENMP)
#include <omp.h>
#endif

// Math library
#include <math.h>

// NLopt library
#if NLOPT_LIB == 1
#include <nlopt.h>
#endif

// Optim library
#if OPTIM_LIB == 1
#include "optim/optim.hpp"
#endif

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Ideal reactors/flames
#include "batch/BatchReactorExperiment.h"
#include "pfr/PlugFlowReactorExperiment.h"
#include "premixed1d/Premixed1DFlameExperiment.h"
#include "counterflow1d/CounterFlow1DFlameExperiment.h"

// Optimizer
#include "Grammar_KineticOptimizer.h"
#include "math/native-minimization-solvers/MinimizationSimplex.h"

const double UNFEASIBLE_BIG_NUMBER = 1.e16;

void CorrectMinMaxValues(Eigen::VectorXd& bMin, Eigen::VectorXd& bMax);
void FromMinimizationParametersToRealParameters(const Eigen::VectorXd& b, Eigen::VectorXd& parameters);
void FromRealParametersToMinimizationParameters(const Eigen::VectorXd& parameters, Eigen::VectorXd& b, Eigen::VectorXd& bMin, Eigen::VectorXd& bMax);
void WriteFinalResult(const Eigen::VectorXd& real_parameters_0, const Eigen::VectorXd& real_parameters);

double OptFunction(const Eigen::VectorXd &b);

#if NLOPT_LIB == 1
double NLOptFunction(unsigned n, const double *x, double *grad, void *my_func_data);
#endif

#if OPTIM_LIB == 1
double OPTIMFunction(const arma::vec& x, arma::vec* grad_out, void* opt_data);
#endif

OpenSMOKE::BatchReactorExperiment* batch_reactors;
OpenSMOKE::PlugFlowReactorExperiment* plugflow_reactors;
OpenSMOKE::Premixed1DFlameExperiment* premixed1d_flames;
OpenSMOKE::CounterFlow1DFlameExperiment* counterflow1d_flames;

unsigned int nExpBatch;
unsigned int nExpPlug;
unsigned int nExpPremixed1D;
unsigned int nExpCounterFlow1D;

bool obj_function_relative_errors;

int number_parameters;
bool central_difference;

double fobj_rel_best;
double fobj_abs_best;
std::ofstream fMonitoring;
std::ofstream fMonitoringBest;

std::vector<int> list_of_target_A;
std::vector<int> list_of_target_Beta;
std::vector<int> list_of_target_E_over_R;

std::vector<double> list_of_min_abs_A;
std::vector<double> list_of_max_abs_A;
std::vector<double> list_of_min_rel_A;
std::vector<double> list_of_max_rel_A;

std::vector<double> list_of_min_abs_Beta;
std::vector<double> list_of_max_abs_Beta;
std::vector<double> list_of_min_rel_Beta;
std::vector<double> list_of_max_rel_Beta;

std::vector<double> list_of_min_abs_E_over_R;
std::vector<double> list_of_max_abs_E_over_R;
std::vector<double> list_of_min_rel_E_over_R;
std::vector<double> list_of_max_rel_E_over_R;

bool violated_uncertainty;
std::vector<int> list_of_target_uncertainty_factors;
std::vector<double> list_of_uncertainty_factors;
std::vector<double> k0_500;
std::vector<double> k0_1000;
std::vector<double> k0_1500;
std::vector<double> k0_2000;

// Thermodynamics and kinetics maps
OpenSMOKE::ThermodynamicsMap_CHEMKIN*		thermodynamicsMapXML;
OpenSMOKE::KineticsMap_CHEMKIN*				kineticsMapXML;
OpenSMOKE::TransportPropertiesMap_CHEMKIN*	transportMapXML;

unsigned int numberOfGradientEvaluations;
unsigned int numberOfFunctionEvaluations;

int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_KineticOptimizer", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	std::string input_file_name_ = "input.dic";
	std::string main_dictionary_name_ = "Optimizer";

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the OpenSMOKEpp_VirtualChemistryOptimizer");
		description.add_options()
			("help", "print help messages")
			("input", po::value<std::string>(), "name of the file containing the main dictionary (default \"input.dic\")")
			("dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"PlugFlowReactor\")");

		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, description), vm); // can throw 

			if (vm.count("help"))
			{
				std::cout << "Basic Command Line Parameters" << std::endl;
				std::cout << description << std::endl;
				return OPENSMOKE_SUCCESSFULL_EXIT;
			}

			if (vm.count("input"))
				input_file_name_ = vm["input"].as<std::string>();

			if (vm.count("dictionary"))
				main_dictionary_name_ = vm["dictionary"].as<std::string>();

			po::notify(vm); // throws on error, so do after help in case  there are any problems 
		}
		catch (po::error& e)
		{
			std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
			std::cerr << description << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}
	}

	// Defines the grammar rules
	OpenSMOKE::Grammar_KineticOptimizer grammar_kinetic_optimizer;

	// Define the dictionaries
	OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
	dictionaries.ReadDictionariesFromFile(input_file_name_);
	dictionaries(main_dictionary_name_).SetGrammar(grammar_kinetic_optimizer);


	// Kinetic scheme
	boost::filesystem::path path_kinetics_output;
	if (dictionaries(main_dictionary_name_).CheckOption("@KineticsFolder") == true)
	{
		dictionaries(main_dictionary_name_).ReadPath("@KineticsFolder", path_kinetics_output);
		OpenSMOKE::CheckKineticsFolder(path_kinetics_output);
	}
	else
	{
		std::string name_of_rapid_kinetics_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@KineticsPreProcessor") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@KineticsPreProcessor", name_of_rapid_kinetics_subdictionary);

		OpenSMOKE::Grammar_RapidKineticMechanism grammar_rapid_kinetics;
		dictionaries(name_of_rapid_kinetics_subdictionary).SetGrammar(grammar_rapid_kinetics);


		boost::filesystem::path path_input_thermodynamics;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Thermodynamics") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Thermodynamics", path_input_thermodynamics);

		boost::filesystem::path path_input_kinetics;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Kinetics") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Kinetics", path_input_kinetics);

		boost::filesystem::path path_input_transport;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Transport") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Transport", path_input_transport);

		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Output") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Output", path_kinetics_output);

		OpenSMOKE::RapidKineticMechanismWithTransport(path_kinetics_output, path_input_transport.c_str(), path_input_thermodynamics.c_str(), path_input_kinetics.c_str());
	}

	// Import thermodynamic and kinetic maps
	{
		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_output / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, doc);
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);

		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
	}

	// List of experiments
	std::vector<std::string> list_of_batch_experiments;
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfBatchExperiments") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfBatchExperiments", list_of_batch_experiments);
	
	std::vector<std::string> list_of_plugflow_experiments;
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfPlugFlowExperiments") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfPlugFlowExperiments", list_of_plugflow_experiments);
	
	std::vector<std::string> list_of_premixed1d_experiments;
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfPremixed1DExperiments") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfPremixed1DExperiments", list_of_premixed1d_experiments);
	
	std::vector<std::string> list_of_counterflow1d_experiments;
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfCounterFlow1DExperiments") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfCounterFlow1DExperiments", list_of_counterflow1d_experiments);
	

	// List of optimization parameters
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfTarget_A") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfTarget_A", list_of_target_A);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfTarget_Beta") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfTarget_Beta", list_of_target_Beta);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfTarget_E_over_R") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfTarget_E_over_R", list_of_target_E_over_R);

	// List of absolute maximum parameters
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMaxAbs_A") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMaxAbs_A", list_of_max_abs_A);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMaxAbs_Beta") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMaxAbs_Beta", list_of_max_abs_Beta);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMaxAbs_E_over_R") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMaxAbs_E_over_R", list_of_max_abs_E_over_R);

	// List of absolute minimum parameters
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMinAbs_A") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMinAbs_A", list_of_min_abs_A);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMinAbs_Beta") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMinAbs_Beta", list_of_min_abs_Beta);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMinAbs_E_over_R") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMinAbs_E_over_R", list_of_min_abs_E_over_R);

	// List of relative maximum parameters
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMaxRel_A") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMaxRel_A", list_of_max_rel_A);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMaxRel_Beta") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMaxRel_Beta", list_of_max_rel_Beta);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMaxRel_E_over_R") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMaxRel_E_over_R", list_of_max_rel_E_over_R);

	// List of relative minimum parameters
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMinRel_A") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMinRel_A", list_of_min_rel_A);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMinRel_Beta") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMinRel_Beta", list_of_min_rel_Beta);
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfMinRel_E_over_R") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfMinRel_E_over_R", list_of_min_rel_E_over_R);

	// List of target reactions for uncertainty factors
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfTargetUncertaintyFactors") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfTargetUncertaintyFactors", list_of_target_uncertainty_factors);

	// List of uncertainty factors
	if (dictionaries(main_dictionary_name_).CheckOption("@ListOfUncertaintyFactors") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ListOfUncertaintyFactors", list_of_uncertainty_factors);

	// Objective function on relative errors
	obj_function_relative_errors = true;
	if (dictionaries(main_dictionary_name_).CheckOption("@RelativeErrors") == true)
		dictionaries(main_dictionary_name_).ReadBool("@RelativeErrors", obj_function_relative_errors);

	// Algorithm
	std::string algorithm = "OpenSMOKEpp-Simplex";
	if (dictionaries(main_dictionary_name_).CheckOption("@Algorithm") == true)
		dictionaries(main_dictionary_name_).ReadString("@Algorithm", algorithm);

	// Variant
	std::string variant = "";
	if (dictionaries(main_dictionary_name_).CheckOption("@Variant") == true)
		dictionaries(main_dictionary_name_).ReadString("@Variant", variant);

	// Central/forward gradient
	central_difference = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@CentralGradient") == true)
		dictionaries(main_dictionary_name_).ReadBool("@CentralGradient", central_difference);

	// Max number iterations
	int max_eval = 1000000;
	if (dictionaries(main_dictionary_name_).CheckOption("@MaxIterations") == true)
		dictionaries(main_dictionary_name_).ReadInt("@MaxIterations", max_eval);

	// Checks
	{
		if (list_of_max_abs_A.size() != list_of_target_A.size() && list_of_max_abs_A.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMaxAbs_A option with wrong number of arguments");
		if (list_of_min_abs_A.size() != list_of_target_A.size() && list_of_min_abs_A.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMinAbs_A option with wrong number of arguments");
		if (list_of_max_rel_A.size() != list_of_target_A.size() && list_of_max_rel_A.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMaxRel_A option with wrong number of arguments");
		if (list_of_min_rel_A.size() != list_of_target_A.size() && list_of_min_rel_A.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMinRel_A option with wrong number of arguments");

		if (list_of_max_abs_Beta.size() != list_of_target_Beta.size() && list_of_max_abs_Beta.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMaxAbs_Beta option with wrong number of arguments");
		if (list_of_min_abs_Beta.size() != list_of_target_Beta.size() && list_of_min_abs_Beta.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMinAbs_Beta option with wrong number of arguments");
		if (list_of_max_rel_Beta.size() != list_of_target_Beta.size() && list_of_max_rel_Beta.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMaxRel_Beta option with wrong number of arguments");
		if (list_of_min_rel_Beta.size() != list_of_target_Beta.size() && list_of_min_rel_Beta.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMinRel_Beta option with wrong number of arguments");

		if (list_of_max_abs_E_over_R.size() != list_of_target_E_over_R.size() && list_of_max_abs_E_over_R.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMaxAbs_E_over_R option with wrong number of arguments");
		if (list_of_min_abs_E_over_R.size() != list_of_target_E_over_R.size() && list_of_min_abs_E_over_R.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMinAbs_E_over_R option with wrong number of arguments");
		if (list_of_max_rel_E_over_R.size() != list_of_target_E_over_R.size() && list_of_max_rel_E_over_R.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMaxRel_E_over_R option with wrong number of arguments");
		if (list_of_min_rel_E_over_R.size() != list_of_target_E_over_R.size() && list_of_min_rel_E_over_R.size() != 0)
			OpenSMOKE::FatalErrorMessage("@ListOfMinRel_E_over_R option with wrong number of arguments");

		if (list_of_target_uncertainty_factors.size() != list_of_uncertainty_factors.size())
			OpenSMOKE::FatalErrorMessage("@ListOfTargetUncertaintyFactors and @ListOfUncertaintyFactors must have the same length");
	}

	// Initial values of kinetic constants
	violated_uncertainty = false;
	k0_500.resize(kineticsMapXML->NumberOfReactions());
	k0_1000.resize(kineticsMapXML->NumberOfReactions());
	k0_1500.resize(kineticsMapXML->NumberOfReactions());
	k0_2000.resize(kineticsMapXML->NumberOfReactions());
	for (unsigned int i = 0; i < kineticsMapXML->NumberOfReactions(); i++)
	{
		k0_500[i] = kineticsMapXML->A(i)*std::pow(500., kineticsMapXML->Beta(i))*std::exp(-kineticsMapXML->E_over_R(i) / 500.);
		k0_1000[i] = kineticsMapXML->A(i)*std::pow(1000., kineticsMapXML->Beta(i))*std::exp(-kineticsMapXML->E_over_R(i) / 1000.);
		k0_1500[i] = kineticsMapXML->A(i)*std::pow(1500., kineticsMapXML->Beta(i))*std::exp(-kineticsMapXML->E_over_R(i) / 1500.);
		k0_2000[i] = kineticsMapXML->A(i)*std::pow(2000., kineticsMapXML->Beta(i))*std::exp(-kineticsMapXML->E_over_R(i) / 2000.);
	}

	// Experiments
	{
		nExpBatch = list_of_batch_experiments.size();
		if (nExpBatch != 0)
		{
			batch_reactors = new OpenSMOKE::BatchReactorExperiment[nExpBatch];
			for (unsigned int k = 0; k < nExpBatch; k++)
			{
				batch_reactors[k].Setup(list_of_batch_experiments[k], thermodynamicsMapXML, kineticsMapXML);
				batch_reactors[k].Solve(true);
			}
		}

		nExpPlug = list_of_plugflow_experiments.size();
		if (nExpPlug != 0)
		{
			plugflow_reactors = new OpenSMOKE::PlugFlowReactorExperiment[nExpPlug];
			for (unsigned int k = 0; k < nExpPlug; k++)
			{
				plugflow_reactors[k].Setup(list_of_plugflow_experiments[k], thermodynamicsMapXML, kineticsMapXML);
				plugflow_reactors[k].Solve(true);
			}
		}

		nExpPremixed1D = list_of_premixed1d_experiments.size();
		if (nExpPremixed1D != 0)
		{
			premixed1d_flames = new OpenSMOKE::Premixed1DFlameExperiment[nExpPremixed1D];
			for (unsigned int k = 0; k < nExpPremixed1D; k++)
			{
				premixed1d_flames[k].Setup(list_of_premixed1d_experiments[k], thermodynamicsMapXML, kineticsMapXML, transportMapXML);
				premixed1d_flames[k].Solve(true);
			}
		}

		nExpCounterFlow1D = list_of_counterflow1d_experiments.size();
		if (nExpCounterFlow1D != 0)
		{
			counterflow1d_flames = new OpenSMOKE::CounterFlow1DFlameExperiment[nExpCounterFlow1D];
			for (unsigned int k = 0; k < nExpCounterFlow1D; k++)
			{
				counterflow1d_flames[k].Setup(list_of_counterflow1d_experiments[k], thermodynamicsMapXML, kineticsMapXML, transportMapXML);
				counterflow1d_flames[k].Solve(true);
			}
		}
	}

	// Minimization
	{
		number_parameters = list_of_target_A.size() + list_of_target_Beta.size() + list_of_target_E_over_R.size();

		// First guess values
		Eigen::VectorXd real_parameters_0(number_parameters);
		unsigned int count = 0;
		for (unsigned int k = 0; k < list_of_target_A.size(); k++)
		{
			const unsigned int index = list_of_target_A[k] - 1;
			real_parameters_0(count++) = kineticsMapXML->A(index);
		}
		for (unsigned int k = 0; k < list_of_target_Beta.size(); k++)
		{
			const unsigned int index = list_of_target_Beta[k] - 1;
			real_parameters_0(count++) = kineticsMapXML->Beta(index);
		}
		for (unsigned int k = 0; k < list_of_target_E_over_R.size(); k++)
		{
			const unsigned int index = list_of_target_E_over_R[k] - 1;
			real_parameters_0(count++) = kineticsMapXML->E_over_R(index);
		}

		// Get minimization parameters
		Eigen::VectorXd b0(number_parameters);
		Eigen::VectorXd bOpt(number_parameters);
		Eigen::VectorXd bMin(number_parameters);
		Eigen::VectorXd bMax(number_parameters);
		FromRealParametersToMinimizationParameters(real_parameters_0, b0, bMin, bMax);
		CorrectMinMaxValues(bMin, bMax);

		// Print data on screen
		std::cout << "Parameters: " << std::endl;
		for (unsigned int i = 0; i < number_parameters; i++)
			std::cout << i << " " << real_parameters_0(i) << " | " << b0(i) << ": " << bMin(i) << "-" << bMax(i) << std::endl;

		std::cout << std::endl;
		std::cout << "Performing additional checks on input data..." << std::endl;
	
		// Check
		for (unsigned int i = 0; i < number_parameters; i++)
		{
			if (bMin(i) >= bMax(i))
				OpenSMOKE::FatalErrorMessage("Error in min/max constraints: min > max");
			if (b0(i) <= bMin(i))
				OpenSMOKE::FatalErrorMessage("Error in min/max constraints: first guess value <= min");
			if (b0(i) <= bMin(i))
				OpenSMOKE::FatalErrorMessage("Error in min/max constraints: first guess value >= max");
		}

		// Initial objective function
		std::cout << "Calculating the initial value of objective function..." << std::endl;
		const double f0 = OptFunction(b0);

		// Start optimization procedure
		Eigen::VectorXd parameters_opt(number_parameters);
		fobj_rel_best = 1.e20;
		fobj_abs_best = 1.e20;

		// Counters
		numberOfFunctionEvaluations = 0;
		numberOfGradientEvaluations = 0;

		fMonitoring.open("log", std::ios::out);
		fMonitoring.setf(std::ios::scientific);
		fMonitoringBest.open("log.best", std::ios::out);
		fMonitoringBest.setf(std::ios::scientific);

		std::cout << "Starting calculation..." << std::endl;

		if (algorithm == "OpenSMOKEpp-Simplex")
		{
			// Initialize the optimization
			OpenSMOKE::MinimizationSimplex mr;
			mr(b0, OptFunction, bMin, bMax);

			// Solve the optimization
			mr();

			// Recover optimization parameters
			mr.GetSolution(bOpt);
		}

		#if NLOPT_LIB == 1
		else if (	algorithm == "DIRECT" || algorithm == "CRS" ||
					algorithm == "MLSL" || algorithm == "STOGO" ||
					algorithm == "ISRES" || algorithm == "ESCH" ||
					algorithm == "COBYLA" || algorithm == "BOBYQA" ||
					algorithm == "NEWUOA" || algorithm == "PRAXIS" ||
					algorithm == "NELDERMEAD" || algorithm == "SBPLX" ||
					algorithm == "SLSQP" || algorithm == "LBFGS" || 
					algorithm == "TNEWTON_PRECOND" || algorithm == "SLM_VAR" )
		{
			nlopt_opt opt;

			if (algorithm == "DIRECT")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GN_DIRECT, number_parameters);
				else if (variant == "L")
					opt = nlopt_create(NLOPT_GN_DIRECT_L, number_parameters);
				else if (variant == "L_RAND")
					opt = nlopt_create(NLOPT_GN_DIRECT_L_RAND, number_parameters);
				else if (variant == "NOSCAL")
					opt = nlopt_create(NLOPT_GN_DIRECT_NOSCAL, number_parameters);
				else if (variant == "L_NOSCAL")
					opt = nlopt_create(NLOPT_GN_DIRECT_L_NOSCAL, number_parameters);
				else if (variant == "L_RAND_NOSCAL")
					opt = nlopt_create(NLOPT_GN_DIRECT_L_RAND_NOSCAL, number_parameters);
				else if (variant == "ORIG")
					opt = nlopt_create(NLOPT_GN_ORIG_DIRECT, number_parameters);
				else if (variant == "L_ORIG")
					opt = nlopt_create(NLOPT_GN_ORIG_DIRECT_L, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for DIRECT algorithms: \
						NONE | L | L_RAND | NOSCAL | L_NOSCAL | L_RAND_NOSCAL | ORIG | L_ORIG");
			}
			else if (algorithm == "CRS")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GN_CRS2_LM, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for CRS algorithms: NONE ");
			}
			else if (algorithm == "MLSL")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_G_MLSL, number_parameters);
				else if (variant == "LDS")
					opt = nlopt_create(NLOPT_G_MLSL_LDS, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for MLSL algorithms: NONE | LDS");

				nlopt_opt local_opt;
				nlopt_create(NLOPT_LN_COBYLA, number_parameters);
				nlopt_set_min_objective(local_opt, NLOptFunction, NULL);
				nlopt_set_local_optimizer(opt, local_opt);
			}
			else if (algorithm == "STOGO")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GD_STOGO, number_parameters);
				else if (variant == "RAND")
					opt = nlopt_create(NLOPT_GD_STOGO_RAND, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for STOGO algorithms: NONE | RAND");
			}
			else if (algorithm == "ISRES")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GN_ISRES, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for ISRES algorithms: NONE");
			}
			else if (algorithm == "ESCH")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_GN_ESCH, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for ESCH algorithms: NONE");
			}
			else if (algorithm == "COBYLA")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_COBYLA, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for COBYLA algorithms: NONE");
			}
			else if (algorithm == "BOBYQA")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_BOBYQA, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for BOBYQA algorithms: NONE");
			}
			else if (algorithm == "NEWUOA")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_NEWUOA_BOUND, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for NEWUOA algorithms: NONE");
			}
			else if (algorithm == "PRAXIS")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_PRAXIS, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for PRAXIS algorithms: NONE");
			}
			else if (algorithm == "NELDERMEAD")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_NELDERMEAD, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for NELDERMEAD algorithms: NONE");
			}
			else if (algorithm == "SBPLX")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LN_SBPLX, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for SBPLX algorithms: NONE");
			}
			else if (algorithm == "SLSQP")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LD_SLSQP, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for SLSQP algorithms: NONE");
			}
			else if (algorithm == "LBFGS")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LD_LBFGS, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for LBFGS algorithms: NONE");
			}
			else if (algorithm == "TNEWTON_PRECOND")
			{
				if (variant == "NONE" || variant == "")
					opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND_RESTART, number_parameters);
				else if (variant == "NO_RESTART")
					opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND, number_parameters);
				else if (variant == "NO_PRECOND")
					opt = nlopt_create(NLOPT_LD_TNEWTON_RESTART, number_parameters);
				else if (variant == "NO_RESTART_NO_PRECOND")
					opt = nlopt_create(NLOPT_LD_TNEWTON, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for TNEWTON_PRECOND algorithms: NONE | NO_RESTART | NO_PRECOND | NO_RESTART_NO_PRECOND");
			}
			else if (algorithm == "SLM_VAR")
			{
				if (variant == "VAR2")
					opt = nlopt_create(NLOPT_LD_VAR2, number_parameters);
				else if (variant == "VAR1")
					opt = nlopt_create(NLOPT_LD_VAR1, number_parameters);
				else
					OpenSMOKE::FatalErrorMessage("Allowed variants for TNEWTON_PRECOND algorithms: VAR2 | VAR1");
			}


			double* lb = new double[number_parameters];// lower bounds
			double* ub = new double[number_parameters];// upper bounds
			double* x = new double[number_parameters];// first guess
			for (unsigned int i = 0; i < number_parameters; i++)
			{
				lb[i] = bMin(i);
				ub[i] = bMax(i);
				x[i] = b0(i);
			}

			nlopt_set_lower_bounds(opt, lb);
			nlopt_set_upper_bounds(opt, ub);
			nlopt_set_min_objective(opt, NLOptFunction, NULL);
			nlopt_set_maxeval(opt, max_eval);

			double fOpt;
			if (nlopt_optimize(opt, x, &fOpt) < 0)
			{
				std::cout << "NLopt failed!" << std::endl;
				getchar();
				exit(-1);
			}
			else
			{
				for (unsigned int i = 0; i < number_parameters; i++)
					bOpt(i) = x[i];
			}
		}
		#endif

		#if OPTIM_LIB == 1
		else if ( algorithm == "OPTIM-DE" || algorithm == "OPTIM-PSO" || algorithm == "OPTIM-NM")
		{
			arma::vec lb(number_parameters);
			arma::vec ub(number_parameters);
			arma::vec x(number_parameters);
			for (unsigned int i = 0; i < number_parameters; i++)
			{
				lb[i] = bMin(i);
				ub[i] = bMax(i);
				x[i] = b0(i);
			}

			optim::algo_settings_t settings;
			settings.vals_bound = true;
			settings.lower_bounds = lb;
			settings.upper_bounds = ub;

			std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

			bool success = false;
			if (algorithm == "OPTIM-DE")
				success = optim::de(x,OPTIMFunction,nullptr, settings);
			else if (algorithm == "OPTIM-PSO")
				success = optim::pso(x,OPTIMFunction,nullptr, settings);
			else if (algorithm == "OPTIM-NM")
				success = optim::nm(x,OPTIMFunction,nullptr, settings);			

			std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end-start;

			if (success == true) 
			{
				std::cout << "OPTIM: Optimization successfully completed" << std::endl;
				std::cout << "       Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

				for (unsigned int i = 0; i < number_parameters; i++)
					bOpt(i) = x[i];
			} 
			else 
			{
				std::cout << "OPTIM failed" << std::endl;
				getchar();
				exit(-1);
			}
		}
		#endif

		else
		{
			OpenSMOKE::FatalErrorMessage("Error @Algorithm option:	OpenSMOKEpp-Simplex | DIRECT | CRS | MLSL | STOGO | ISRES | ESCH | \
																	COBYLA | BOBYQA |NEWUOA | PRAXIS | NELDERMEAD | SBPLX | \
																	SLSQP | LBFGS | TNEWTON_PRECOND | SLM_VAR");
		}

		// Close files
		fMonitoring.close();
		fMonitoringBest.close();

		// Final objective function
		const double fOpt = OptFunction(bOpt);

		// Recover final solution
		FromMinimizationParametersToRealParameters(bOpt, parameters_opt);

		// Final message
		std::cout << "Objective function: " << f0 << " -> " << fOpt << std::endl;
		WriteFinalResult(real_parameters_0, parameters_opt);
		
		getchar();
	}

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_KineticOptimizer", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	return OPENSMOKE_SUCCESSFULL_EXIT;
}


double ReturnObjFunction(const Eigen::VectorXd parameters)
{
	double fobj_abs = 0.;
	double fobj_rel = 0.;
	{
		for (unsigned int k = 0; k < nExpBatch; k++)
		{
			batch_reactors[k].Solve();
			fobj_abs += batch_reactors[k].norm2_abs_error();
			fobj_rel += batch_reactors[k].norm2_rel_error();
		}

		for (unsigned int k = 0; k < nExpPlug; k++)
		{
			plugflow_reactors[k].Solve();
			fobj_abs += plugflow_reactors[k].norm2_abs_error();
			fobj_rel += plugflow_reactors[k].norm2_rel_error();
		}

		for (unsigned int k = 0; k < nExpPremixed1D; k++)
		{
			premixed1d_flames[k].Solve();
			fobj_abs += premixed1d_flames[k].norm2_abs_error();
			fobj_rel += premixed1d_flames[k].norm2_rel_error();
		}

		for (unsigned int k = 0; k < nExpCounterFlow1D; k++)
		{
			counterflow1d_flames[k].Solve();
			fobj_abs += counterflow1d_flames[k].norm2_abs_error();
			fobj_rel += counterflow1d_flames[k].norm2_rel_error();
		}
	}

	std::cout << "f(abs) = " << fobj_abs << "   -   f(rel) = " << fobj_rel << std::endl;

	// Objective function based on relative errors
	if (obj_function_relative_errors == true)
	{
		if (fobj_rel < fobj_rel_best)
		{
			fobj_rel_best = fobj_rel;

			fMonitoringBest << std::left << std::setw(16) << violated_uncertainty;
			fMonitoringBest << std::left << std::setw(16) << fobj_rel;
			for (unsigned int i = 0; i<number_parameters; i++)
				fMonitoringBest << std::left << std::setw(16) << parameters(i);
			fMonitoringBest << std::endl;
		}

		fMonitoring << std::left << std::setw(16) << violated_uncertainty;
		fMonitoring << std::left << std::setw(16) << fobj_rel;
		for (unsigned int i = 0; i<number_parameters; i++)
			fMonitoring << std::left << std::setw(16) << parameters(i);
		fMonitoring << std::endl;

		return fobj_rel;
	}
	// Objective function based on absolute errors
	else
	{
		if (fobj_abs < fobj_abs_best)
		{
			fobj_abs_best = fobj_abs;

			fMonitoringBest << std::left << std::setw(16) << violated_uncertainty;
			fMonitoringBest << std::left << std::setw(16) << fobj_abs;
			for (unsigned int i = 0; i<number_parameters; i++)
				fMonitoringBest << std::left << std::setw(16) << parameters(i);
			fMonitoringBest << std::endl;
		}

		fMonitoring << std::left << std::setw(16) << violated_uncertainty;
		fMonitoring << std::left << std::setw(16) << fobj_abs;
		for (unsigned int i = 0; i<number_parameters; i++)
			fMonitoring << std::left << std::setw(16) << parameters(i);
		fMonitoring << std::endl;

		return fobj_abs;
	}
}

bool CheckViolationUncertaintyFactors()
{
	for (unsigned int k = 0; k < list_of_target_uncertainty_factors.size(); k++)
	{
		const unsigned int index = list_of_target_uncertainty_factors[k]-1;
		const double gamma = std::pow(10., list_of_uncertainty_factors[k]);
		
		const double k500 = kineticsMapXML->A(index)*std::pow(500., kineticsMapXML->Beta(index))*std::exp(-kineticsMapXML->E_over_R(index) / 500.);
		if (k500 < k0_500[index] / gamma || k500 > k0_500[index] * gamma) 
		{
			std::cout << "Violation @500K, Reaction: " << index+1 << ", Current ratio: " << k500/k0_500[index] << ", Allowed ratio: " << 1./gamma << "-" << gamma << std::endl; 
			return true;
		}
		
		const double k1000 = kineticsMapXML->A(index)*std::pow(1000., kineticsMapXML->Beta(index))*std::exp(-kineticsMapXML->E_over_R(index) / 1000.);
		if (k1000 < k0_1000[index] / gamma || k1000 > k0_1000[index] * gamma)
		{
			std::cout << "Violation @1000K, Reaction: " << index + 1 << ", Current ratio: " << k1000 / k0_1000[index] << ", Allowed ratio: " << 1. / gamma << "-" << gamma << std::endl;
			return true;
		}
		
		const double k1500 = kineticsMapXML->A(index)*std::pow(1500., kineticsMapXML->Beta(index))*std::exp(-kineticsMapXML->E_over_R(index) / 1500.);
		if (k1500 < k0_1500[index] / gamma || k1500 > k0_1500[index] * gamma)
		{
			std::cout << "Violation @1500K, Reaction: " << index + 1 << ", Current ratio: " << k1500 / k0_1500[index] << ", Allowed ratio: " << 1. / gamma << "-" << gamma << std::endl;
			return true;
		}

		const double k2000 = kineticsMapXML->A(index)*std::pow(2000., kineticsMapXML->Beta(index))*std::exp(-kineticsMapXML->E_over_R(index) / 2000.);
		if (k2000 < k0_2000[index] / gamma || k2000 > k0_2000[index] * gamma)
		{
			std::cout << "Violation @2000K, Reaction: " << index + 1 << ", Current ratio: " << k2000 / k0_2000[index] << ", Allowed ratio: " << 1. / gamma << "-" << gamma << std::endl;
			return true;
		}
	}

	return false;
}

double OptFunction(const Eigen::VectorXd &b)
{
	Eigen::VectorXd real_parameters(b.size());
	FromMinimizationParametersToRealParameters(b, real_parameters);

	unsigned int count = 0;
	for (unsigned int k = 0; k < list_of_target_A.size(); k++)
	{
		const unsigned int index = list_of_target_A[k] - 1;
		kineticsMapXML->Set_A(index, real_parameters(count++));
	}
	for (unsigned int k = 0; k < list_of_target_Beta.size(); k++)
	{
		const unsigned int index = list_of_target_Beta[k] - 1;
		kineticsMapXML->Set_Beta(index, real_parameters(count++));
	}
	for (unsigned int k = 0; k < list_of_target_E_over_R.size(); k++)
	{
		const unsigned int index = list_of_target_E_over_R[k] - 1;
		kineticsMapXML->Set_E_over_R(index, real_parameters(count++));
	}

	if (list_of_target_uncertainty_factors.size() != 0)
	{
		violated_uncertainty = CheckViolationUncertaintyFactors();
		if (violated_uncertainty == true)
		{
			std::cout << "Violated constraints based on uncertainty factors..." << std::endl;
			OpenSMOKE::MinimizationSimplex::unfeasible = true;
			return UNFEASIBLE_BIG_NUMBER;
		}
	}

	return ReturnObjFunction(real_parameters);
}

#if NLOPT_LIB == 1
double NLOptFunction(unsigned n, const double *x, double *grad, void *my_func_data)
{
	Eigen::VectorXd b(number_parameters);
	for (unsigned int i = 0; i < number_parameters; i++)
		b(i) = x[i];

	const double f = OptFunction(b);

	if (grad)
	{
		std::cout << "Gradient evaluation..." << std::endl;

		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double ETA3 = std::pow(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE, 1./3.);
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_DOUBLE);
		
		// Dimensions of parameters vector
		Eigen::VectorXd b_dimensions(number_parameters);
		b_dimensions.setConstant(1.);

		// Choosing between forward and centrale difference
		double eta = ETA2;
		if (central_difference == true)
			eta = ETA3;
		
		// Estimate the increment for forward approximation
		Eigen::VectorXd deltab(number_parameters);
		for (unsigned int i = 0; i < number_parameters; i++)
		{
			if (b(i) < 0.)
				deltab(i) = -eta * std::max(std::fabs(b(i)), std::fabs(b_dimensions(i)));
			else
				deltab(i) =  eta * std::max(std::fabs(b(i)), std::fabs(b_dimensions(i)));
			
			if (deltab(i) == 0.)
				deltab(i) = ZERO_DER;
		}

		// Forward gradient
		if (central_difference == false)
		{
			Eigen::VectorXd b_plus = b;
			for (unsigned int j = 0; j < number_parameters; j++)
			{
				b_plus(j) = b(j) + deltab(j);
				const double f_plus = OptFunction(b_plus);
				
				grad[j] = (f_plus - f) / deltab(j);

				b_plus(j) = b(j);
			}
		}

		// Central gradient
		if (central_difference == true)
		{
			Eigen::VectorXd b_star = b;
			for (unsigned int j = 0; j < number_parameters; j++)
			{
				b_star(j) = b(j) + deltab(j);
				const double f_plus = OptFunction(b_star);
				b_star(j) = b(j) - deltab(j);
				const double f_minus = OptFunction(b_star);

				grad[j] = (f_plus - f_minus) / (2.*deltab(j));

				b_star(j) = b(j);
			}
		}

		numberOfGradientEvaluations++;
	}
	
	numberOfFunctionEvaluations++;

	return f;
}
#endif

#if OPTIM_LIB == 1
double OPTIMFunction(const arma::vec& x, arma::vec* grad_out, void* opt_data)
{
	Eigen::VectorXd b(number_parameters);
	for (unsigned int i = 0; i < number_parameters; i++)
		b(i) = x[i];

	const double f = OptFunction(b);

	return f;
}
#endif

void FromMinimizationParametersToRealParameters(const Eigen::VectorXd& b, Eigen::VectorXd& real_parameters)
{
	unsigned int count = 0;
	for (unsigned int k = 0; k < list_of_target_A.size(); k++)
	{
		real_parameters(count) = std::exp(b(count));
		count++;
	}
	for (unsigned int k = 0; k < list_of_target_Beta.size(); k++)
	{
		real_parameters(count) = b(count);
		count++;
	}
	for (unsigned int k = 0; k < list_of_target_E_over_R.size(); k++)
	{
		real_parameters(count) = b(count);
		count++;
	}
}

void FromRealParametersToMinimizationParameters(const Eigen::VectorXd& real_parameters, Eigen::VectorXd& b, Eigen::VectorXd& bMin, Eigen::VectorXd& bMax)
{
	const double psi_A = 2.;
	const double psi_Beta = 1.05;
	const double psi_E_over_R = 1.05;

	unsigned int count = 0;
	for (unsigned int k = 0; k < list_of_target_A.size(); k++)
	{
		b(count) = std::log(real_parameters(count));
		bMin(count) = std::log(real_parameters(count) / psi_A);
		bMax(count) = std::log(real_parameters(count) * psi_A);
		count++;
	}
	for (unsigned int k = 0; k < list_of_target_Beta.size(); k++)
	{
		b(count) = real_parameters(count);
		bMin(count) = real_parameters(count) / psi_Beta;
		bMax(count) = real_parameters(count) * psi_Beta;
		count++;
	}
	for (unsigned int k = 0; k < list_of_target_E_over_R.size(); k++)
	{
		b(count) = real_parameters(count);
		bMin(count) = real_parameters(count) / psi_E_over_R;
		bMax(count) = real_parameters(count) * psi_E_over_R;
		count++;
	}
}

void CorrectMinMaxValues(Eigen::VectorXd& bMin, Eigen::VectorXd& bMax)
{
	unsigned int count = 0;
	for (unsigned int k = 0; k < list_of_target_A.size(); k++)
	{
		if (list_of_min_abs_A.size() != 0)	bMin(count) = std::log(list_of_min_abs_A[k]);
		if (list_of_min_rel_A.size() != 0)	bMin(count) += std::log(list_of_min_rel_A[k]);

		if (list_of_max_abs_A.size() != 0)	bMax(count) = std::log(list_of_max_abs_A[k]);
		if (list_of_max_rel_A.size() != 0)	bMax(count) += std::log(list_of_max_rel_A[k]);

		count++;
	}
	for (unsigned int k = 0; k < list_of_target_Beta.size(); k++)
	{
		if (list_of_min_abs_Beta.size() != 0)	bMin(count) = list_of_min_abs_Beta[k];
		if (list_of_min_rel_Beta.size() != 0)	bMin(count) *= list_of_min_rel_Beta[k];

		if (list_of_max_abs_Beta.size() != 0)	bMax(count) = list_of_max_abs_Beta[k];
		if (list_of_max_rel_Beta.size() != 0)	bMax(count) *= list_of_max_rel_Beta[k];

		count++;
	}
	for (unsigned int k = 0; k < list_of_target_E_over_R.size(); k++)
	{
		if (list_of_min_abs_E_over_R.size() != 0)	bMin(count) = list_of_min_abs_E_over_R[k];
		if (list_of_min_rel_E_over_R.size() != 0)	bMin(count) *= list_of_min_rel_E_over_R[k];

		if (list_of_max_abs_E_over_R.size() != 0)	bMax(count)  = list_of_max_abs_E_over_R[k];
		if (list_of_max_rel_E_over_R.size() != 0)	bMax(count) *= list_of_max_rel_E_over_R[k];

		count++;
	}
}

void WriteFinalResult(const Eigen::VectorXd& real_parameters_0, const Eigen::VectorXd& real_parameters)
{
	std::ofstream fOut("Results.out", std::ios::out);
	fOut.setf(std::ios::scientific);

	unsigned int count = 0;
	for (unsigned int k = 0; k < list_of_target_A.size(); k++)
	{
		fOut << " * Reaction: " << list_of_target_A[k] << std::endl;
		fOut << "   A: " << real_parameters_0(count) << " -> " << real_parameters(count) << " (ratio: " << real_parameters(count) / real_parameters_0(count) << ")" << std::endl;
		fOut << std::endl;

		count++;
	}
	for (unsigned int k = 0; k < list_of_target_Beta.size(); k++)
	{
		fOut << " * Reaction: " << list_of_target_Beta[k] << std::endl;
		fOut << "   Beta: " << real_parameters_0(count) << " -> " << real_parameters(count) << " (ratio: " << real_parameters(count) / real_parameters_0(count) << ")" << std::endl;
		fOut << std::endl;
		count++;
	}
	for (unsigned int k = 0; k < list_of_target_E_over_R.size(); k++)
	{
		fOut << " * Reaction: " << list_of_target_E_over_R[k] << std::endl;
		fOut << "   E: " << real_parameters_0(count)*PhysicalConstants::R_cal_mol << " -> " << real_parameters(count)*PhysicalConstants::R_cal_mol << " (ratio: " << real_parameters(count) / real_parameters_0(count) << ")" << std::endl;
		fOut << std::endl;
		count++;
	}

	fOut.close();
}
