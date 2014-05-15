#ifndef CONSERVATIVESOLVER3D_H
#define CONSERVATIVESOLVER3D_H
#include "Solver.h"
#include "Field3D.h"

namespace EMWSolver
{
	class CConservativeSolver3D :
		public CSolver
	{
	public:
		CConservativeSolver3D(CField3D* _field);
		void Create(TaskParameters taskParam, NumericalParameters numParam, SourceParameters srcParam, 
			BoundaryParameters* bndParam);
		void Solve(int timeSteps);
		void SolveStep();
		~CConservativeSolver3D(void);

	private:
		double*** Eps;
		double*** Mju;

		double*** Sig;
		double*** SigS;

		double*** Ms;
		double*** Js;

		double dx, dy, dz;
		double dt;

		double S;
		int gridX, gridY, gridZ;
		int globalTimestep;

		void prepare();

		CField3D* field;
	};
}
#endif