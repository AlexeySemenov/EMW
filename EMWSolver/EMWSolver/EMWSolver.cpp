#include "Field3D.h"
#include "ParameterParser.h"
#include "ConservativeSolver3D.h"
#include "Log.h"
#include <iostream>

using namespace EMWSolver;

int main(int argc, char* argv[])
{
	Log::GetInstance().SetOutputStream(&std::cout);
	CParameterParser parser("default.xml");
	NumericalParameters dflt = parser.GetDefaultNumericalParameters();

	CField3D* field = new CField3D(dflt.sizeX, dflt.sizeY, dflt.sizeZ);

	CConservativeSolver3D* solver = new CConservativeSolver3D(field);

	solver->Create(parser.GetDefaultTaskParameters(), parser.GetDefaultNumericalParameters(), 
		parser.GetDefaultSourceParameters(), parser.GetDefaultBoundaryParameters());

	solver->Solve(200);
	

	//system("PAUSE");
	delete field;
	delete solver;
	system("PAUSE");

	Log::GetInstance().Close();
	return 0;
}