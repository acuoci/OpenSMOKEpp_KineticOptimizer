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

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_KineticOptimizer : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder",
				OpenSMOKE::SINGLE_PATH,
				"Name of the folder containing the kinetic scheme (XML Version)",
				true,
				"@KineticsPreProcessor",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary containing the list of kinetic files to be interpreted",
				true,
				"@KineticsFolder",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfBatchExperiments",
				OpenSMOKE::VECTOR_STRING,
				"List of input files for experiments (batch reactors)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfPlugFlowExperiments",
				OpenSMOKE::VECTOR_STRING,
				"List of input files for experiments (plug flow reactors)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfPremixed1DExperiments",
				OpenSMOKE::VECTOR_STRING,
				"List of input files for experiments (premixed 1d flames)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_A",
				OpenSMOKE::VECTOR_INT,
				"List of target reactions for frequency factors (indices starting from 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_Beta",
				OpenSMOKE::VECTOR_INT,
				"List of target reactions for temperature exponents (indices starting from 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_E_over_R",
				OpenSMOKE::VECTOR_INT,
				"List of target reactions for activation temperatures (indices starting from 1)",
				false));


			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_A",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of minimum absolute values for frequency factors (units: kmol, m, s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_A",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of maximum absolute values for frequency factors (units: kmol, m, s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_A",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of minimum relative values for frequency factors",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_A",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of maximum relative values for frequency factors",
				false));


			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_Beta",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of minimum absolute values for temperature exponent",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_Beta",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of maximum absolute values for temperature exponent",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_Beta",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of minimum relative values for temperature exponent",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_Beta",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of maximum relative values for temperature exponent",
				false));


			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_E_over_R",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of minimum absolute values for activation temperature",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_E_over_R",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of maximum absolute values for activation temperature",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_E_over_R",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of minimum relative values for activation temperature",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_E_over_R",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of maximum relative values for activation temperature",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTargetUncertaintyFactors",
				OpenSMOKE::VECTOR_INT,
				"List of reaction indices (starting from 1) for which uncertainty factors are defined",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors",
				OpenSMOKE::VECTOR_DOUBLE,
				"List of uncertainty factors",
				false));
			
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RelativeErrors",
				OpenSMOKE::SINGLE_BOOL,
				"Objective function of relative errors (default: true)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaxIterations",
				OpenSMOKE::SINGLE_INT,
				"Max number of evaluations (default: 100000)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Algorithm",
				OpenSMOKE::SINGLE_STRING,
				"Algorithm: OpenSMOKEpp-Simplex | DIRECT | CRS | MLSL | STOGO | ISRES | ESCH | \
							COBYLA | BOBYQA |NEWUOA | PRAXIS | NELDERMEAD | SBPLX | \
							SLSQP | LBFGS | TNEWTON_PRECOND | SLM_VAR (default: OpenSMOKEpp-Simplex)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Variant",
				OpenSMOKE::SINGLE_STRING,
				"Algorithm variant (default: none)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CentralGradient",
				OpenSMOKE::SINGLE_BOOL,
				"Central differences for calculating the gradient instead of forward differences (default: false)",
				false));
		}
	};
}

