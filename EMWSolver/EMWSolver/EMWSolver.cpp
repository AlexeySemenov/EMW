#include "Field3D.h"
#include "ParameterParser.h"
#include "ConservativeSolver3D.h"
#include "Log.h"
#include <iostream>
#include <ctime>


using namespace EMWSolver;

int main(int argc, char* argv[])
{
	clock_t wrkTime;
	Log::GetInstance().SetOutputStream(&std::cout);
	CParameterParser parser("default.xml");
	NumericalParameters dflt = parser.GetDefaultNumericalParameters();

	for(int i = 4; i <= 7; i++)
	{
		CField3D* field = new CField3D(dflt.sizeX, dflt.sizeY, dflt.sizeZ);

		CConservativeSolver3D* solver = new CConservativeSolver3D(field);

		solver->Create(parser.GetDefaultTaskParameters(), parser.GetDefaultNumericalParameters(), 
			parser.GetDefaultSourceParameters(), parser.GetDefaultBoundaryParameters(), i);

		solver->Solve(3000);
	

		delete solver;
		delete field;
	}

	wrkTime = clock();
	std::cout << "Solve Time: "<<wrkTime/CLOCKS_PER_SEC << std::endl;
	system("PAUSE");

	Log::GetInstance().Close();
	return 0;
}