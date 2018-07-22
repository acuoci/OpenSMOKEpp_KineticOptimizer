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

#ifndef OpenSMOKE_OptimizationRules_PlugFlowReactorExperiment_H
#define OpenSMOKE_OptimizationRules_PlugFlowReactorExperiment_H

namespace OpenSMOKE
{
	class Grammar_OptimizationRules_PlugFlowReactorExperiment : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@IgnitionDelayTime",
				OpenSMOKE::SINGLE_MEASURE,
				"Ignition delay time",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Definition",
				OpenSMOKE::SINGLE_STRING,
				"Ignition delay time definition: slope | increase",
				true));
		}
	};

	class OptimizationRules_PlugFlowReactorExperiment
	{
	public:

		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, OpenSMOKE::OpenSMOKE_DictionaryManager& dictionaries);

		double tau() const { return tau_; }
		bool slope_definition() const { return slope_definition_; }

	private:

		double tau_;
		bool slope_definition_;
	};

	void OptimizationRules_PlugFlowReactorExperiment::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, OpenSMOKE::OpenSMOKE_DictionaryManager& dictionaries)
	{
		Grammar_OptimizationRules_PlugFlowReactorExperiment grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@IgnitionDelayTime") == true)
		{
			std::string units;
			double value;
			dictionary.ReadMeasure("@IgnitionDelayTime", value, units);
			if (units == "s")         tau_ = value;
			else if (units == "ms")   tau_ = value / 1000.;
			else if (units == "min")  tau_ = value * 60.;
			else if (units == "h")    tau_ = value * 3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown time units");
		}

		if (dictionary.CheckOption("@Definition") == true)
		{
			std::string value;
			dictionary.ReadString("@Definition", value);
			if (value == "slope")			slope_definition_ = true;
			else if (value == "increase")   slope_definition_ = false;
			else OpenSMOKE::FatalErrorMessage("Unknown ignition delay time definition. Available options: slope | increase");
		}

	}
}

#endif // OpenSMOKE_OptimizationRules_PlugFlowReactorExperiment_H