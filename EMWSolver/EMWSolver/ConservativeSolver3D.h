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
			BoundaryParameters* bndParam, int eps);
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

		void setupWaveguideGeometry(int aSize, int bSize, int length, int bOff, double leps);
		void setupWaveguideMaterial(int length, int bOff, double leps);
		void storeRawEy(int i, int j);
		void calculateEyAmp(double** inputAmp, int startZ, int endZ);

		double** RawEyCenter;
		double** RawEyUp;
		double** RawEyDown;
		double** RawEyLeft;
		double** RawEyRight;

		double*** Ew;

		int startRawTime;
		int endRawTime;
		int rawIterator;

		double gamma0;
		double gamma1;

		double gamma0_num;
		double gamma1_num;

		double Ey0;
		int l0;

		int wgIndex;

		int _a;
		int _b;

		int _x1;
		int _x2;
		int _y1;
		int _y2;
		int _z1;
		int _z2;

		CField3D* field;
	};
}
#endif