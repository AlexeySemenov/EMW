#include "ConservativeSolver3D.h"
#include "PhysicConst.h"
#include "math.h"
#include "Helpers.h"
#include "Log.h"
#include <sstream>

namespace EMWSolver
{
	CConservativeSolver3D::CConservativeSolver3D(CField3D* _field)
	{
		field = _field;
	}

	void CConservativeSolver3D::Create(TaskParameters taskParam, NumericalParameters numParam, SourceParameters srcParam, 
			BoundaryParameters* bndParam)
	{
		std::stringstream logger;
		Log::GetInstance().WriteLine("Creating solver:");
		dx = LAMBDA / numParam.Nx;
		dy = LAMBDA / numParam.Ny;
		dz = LAMBDA / numParam.Nz;
		
		logger << "dx = " << dx << std::endl;
		logger << "dy = " << dy << std::endl;
		logger << "dz = " << dz << std::endl;
		

		dt = ((double)numParam.S) / (CC * sqrt((1.0 / dx) * (1.0 / dx) + (1.0 / dy) * (1.0 / dy) + (1.0 / dz) * (1.0 / dz)));

		logger << "dt = " << dt << std::endl;

		globalTimestep = 0;

		Log::GetInstance().Write(logger.str());
		logger.str("");

		Log::GetInstance().Write("Creating solver arrays...");

		Eps = CreateJaggedArray3D(numParam.sizeX, numParam.sizeY, numParam.sizeZ);
		Mju = CreateJaggedArray3D(numParam.sizeX, numParam.sizeY, numParam.sizeZ);
		Sig = CreateJaggedArray3D(numParam.sizeX, numParam.sizeY, numParam.sizeZ);
		SigS = CreateJaggedArray3D(numParam.sizeX, numParam.sizeY, numParam.sizeZ);
		Ms = CreateJaggedArray3D(numParam.sizeX, numParam.sizeY, numParam.sizeZ);
		Js = CreateJaggedArray3D(numParam.sizeX, numParam.sizeY, numParam.sizeZ);

		ValueJaggedArray3D(Eps, numParam.sizeX, numParam.sizeY, numParam.sizeZ, 1.0);
		ValueJaggedArray3D(Mju, numParam.sizeX, numParam.sizeY, numParam.sizeZ, 1.0);
		ZeroJaggedArray3D(Sig, numParam.sizeX, numParam.sizeY, numParam.sizeZ);
		ZeroJaggedArray3D(SigS, numParam.sizeX, numParam.sizeY, numParam.sizeZ);
		ZeroJaggedArray3D(Ms, numParam.sizeX, numParam.sizeY, numParam.sizeZ);
		ZeroJaggedArray3D(Js, numParam.sizeX, numParam.sizeY, numParam.sizeZ);

		

		Log::GetInstance().WriteLine("Finished.");

		gridX = numParam.sizeX;
		gridY = numParam.sizeY;
		gridZ = numParam.sizeZ;

		logger << "Task size: " << gridX << " " << gridY << " " << gridZ;
		Log::GetInstance().WriteLine(logger.str());

		

		prepare();

		//for(int i = 0; i < gridX; i++)
		//	for(int j = 0; j < gridY; j++)
		//		for(int k = gridZ-20; k < gridZ; k++)
		//		{
		//			Eps[i][j][k] = 4.0;
		//			Mju[i][j][k] = 4.0;
		//		}

	}

	void CConservativeSolver3D::SolveStep()
	{
		const double invdx = 1.0 / dx;
		const double invdy = 1.0 / dy;
		const double invdz = 1.0 / dz;

		Log::GetInstance().Write("Solving step: ");
		Log::GetInstance().Write(globalTimestep);
		Log::GetInstance().Write(" ...");
				//std::cout << "TimeStep:" << globalTimestep << std::endl;
		//Calculate E-field at boundaries
		//i = 0, gridX
		for(int j = 1; j < gridY - 1; j++)
		{
				#pragma omp parallel for
				for(int k = 1; k < gridZ - 1; k++)
				{
					double EpsX;
					double SigX;
					double dtEps, sigDt2Eps, plus1sigDt2Eps;
					int i = 0;

					EpsX = 0.25 * (Eps[i][j][k] + Eps[i][j][k - 1] + Eps[i][j - 1][k - 1] + Eps[i][j - 1][k]) * EPS_Z;
					SigX = 0.25 * (Sig[i][j][k] + Sig[i][j][k - 1] + Sig[i][j - 1][k - 1] + Sig[i][j - 1][k]);
					dtEps = dt / EpsX;
					sigDt2Eps = 0.5 * SigX * dtEps;
					plus1sigDt2Eps = 1.0 + sigDt2Eps;

					field->Ex[i][j][k] = ((1.0 - sigDt2Eps) / plus1sigDt2Eps) * field->Ex[i][j][k] + 
						(dtEps / plus1sigDt2Eps) * ((field->Hz[i][j][k] - field->Hz[i][j-1][k]) * invdy - 
						(field->Hy[i][j][k] - field->Hy[i][j][k-1]) * invdz - Js[i][j][k]);

					/*i = gridX - 1;

					EpsX = 0.25 * (Eps[i][j][k] + Eps[i][j][k - 1] + Eps[i][j - 1][k - 1] + Eps[i][j - 1][k]) * EPS_Z;
					SigX = 0.25 * (Sig[i][j][k] + Sig[i][j][k - 1] + Sig[i][j - 1][k - 1] + Sig[i][j - 1][k]);
					dtEps = dt / EpsX;
					sigDt2Eps = 0.5 * SigX * dtEps;
					plus1sigDt2Eps = 1.0 + sigDt2Eps;

					field->Ex[i][j][k] = ((1.0 - sigDt2Eps) / plus1sigDt2Eps) * field->Ex[i][j][k] + 
						(dtEps / plus1sigDt2Eps) * ((field->Hz[i][j][k] - field->Hz[i][j-1][k]) * invdy - 
						(field->Hy[i][j][k] - field->Hy[i][j][k-1]) * invdz - Js[i][j][k]);*/
				}
		}

		//j = 0, gridY -1
		for(int i = 1; i < gridX - 1; i++)
		{
				#pragma omp parallel for
				for(int k = 1; k < gridZ - 1; k++)
				{
					double EpsY;
					double SigY;
					double dtEps, sigDt2Eps, plus1sigDt2Eps;
					int j = 0;
					
					EpsY = 0.25 * (Eps[i][j][k] + Eps[i][j][k - 1] + Eps[i - 1][j][k - 1] + Eps[i - 1][j][k]) * EPS_Z;
					SigY = 0.25 * (Sig[i][j][k] + Sig[i][j][k - 1] + Sig[i - 1][j][k - 1] + Sig[i - 1][j][k]);
					dtEps = dt / EpsY;
					sigDt2Eps = 0.5 * SigY * dtEps;
					plus1sigDt2Eps = 1.0 + sigDt2Eps;

					field->Ey[i][j][k] = ((1.0 - sigDt2Eps) / plus1sigDt2Eps) * field->Ey[i][j][k] + 
						(dtEps / plus1sigDt2Eps) * ((field->Hx[i][j][k] - field->Hx[i][j][k - 1]) * invdz - 
						(field->Hz[i][j][k] - field->Hz[i - 1][j][k]) * invdx - Js[i][j][k]);

					/*j = gridY - 1;
					
					EpsY = 0.25 * (Eps[i][j][k] + Eps[i][j][k - 1] + Eps[i - 1][j][k - 1] + Eps[i - 1][j][k]) * EPS_Z;
					SigY = 0.25 * (Sig[i][j][k] + Sig[i][j][k - 1] + Sig[i - 1][j][k - 1] + Sig[i - 1][j][k]);
					dtEps = dt / EpsY;
					sigDt2Eps = 0.5 * SigY * dtEps;
					plus1sigDt2Eps = 1.0 + sigDt2Eps;

					field->Ey[i][j][k] = ((1.0 - sigDt2Eps) / plus1sigDt2Eps) * field->Ey[i][j][k] + 
						(dtEps / plus1sigDt2Eps) * ((field->Hx[i][j][k] - field->Hx[i][j][k - 1]) * invdz - 
						(field->Hz[i][j][k] - field->Hz[i - 1][j][k]) * invdx - Js[i][j][k]);*/

				}
		}
		
		//z = 0, gridZ -1
		for(int i = 1; i < gridX - 1; i++)
		{
				#pragma omp parallel for
				for(int j = 1; j < gridY - 1; j++)
				{
					double EpsZ;
					double SigZ;
					double dtEps, sigDt2Eps, plus1sigDt2Eps;
					int k = 0;
					
					EpsZ = 0.25 * (Eps[i][j][k] + Eps[i - 1][j][k] + Eps[i - 1][j - 1][k] + Eps[i][j - 1][k]) * EPS_Z;
					SigZ = 0.25 * (Sig[i][j][k] + Sig[i - 1][j][k] + Sig[i - 1][j - 1][k] + Sig[i][j - 1][k]);
					dtEps = dt / EpsZ;
					sigDt2Eps = 0.5 * SigZ * dtEps;
					plus1sigDt2Eps = 1.0 + sigDt2Eps;

					field->Ez[i][j][k] = ((1.0 - sigDt2Eps) / plus1sigDt2Eps) * field->Ez[i][j][k] + 
						(dtEps / plus1sigDt2Eps) * ((field->Hy[i][j][k] - field->Hy[i - 1][j][k]) * invdx - 
						(field->Hx[i][j][k] - field->Hx[i][j - 1][k]) * invdx - Js[i][j][k]);
					
					//k = gridZ - 1;
					//
					//EpsZ = 0.25 * (Eps[i][j][k] + Eps[i - 1][j][k] + Eps[i - 1][j - 1][k] + Eps[i][j - 1][k]) * EPS_Z;
					//SigZ = 0.25 * (Sig[i][j][k] + Sig[i - 1][j][k] + Sig[i - 1][j - 1][k] + Sig[i][j - 1][k]);
					//dtEps = dt / EpsZ;
					//sigDt2Eps = 0.5 * SigZ * dtEps;
					//plus1sigDt2Eps = 1.0 + sigDt2Eps;

					//field->Ez[i][j][k] = ((1.0 - sigDt2Eps) / plus1sigDt2Eps) * field->Ez[i][j][k] + 
					//	(dtEps / plus1sigDt2Eps) * ((field->Hy[i][j][k] - field->Hy[i - 1][j][k]) * invdx - 
					//	(field->Hx[i][j][k] - field->Hx[i][j - 1][k]) * invdx - Js[i][j][k]);
				}
		}
		//Calculate E - field at domain
		for(int i = 1; i < gridX - 1; i++)
			for(int j = 1; j < gridY - 1; j++)
			{
				#pragma omp parallel for
				for(int k = 1; k < gridZ - 1; k++)
				{
					double EpsX, EpsY, EpsZ;
					double SigX, SigY, SigZ;
					double dtEps, sigDt2Eps, plus1sigDt2Eps;

					EpsX = 0.25 * (Eps[i][j][k] + Eps[i][j][k - 1] + Eps[i][j - 1][k - 1] + Eps[i][j - 1][k]) * EPS_Z;
					SigX = 0.25 * (Sig[i][j][k] + Sig[i][j][k - 1] + Sig[i][j - 1][k - 1] + Sig[i][j - 1][k]);
					dtEps = dt / EpsX;
					sigDt2Eps = 0.5 * SigX * dtEps;
					plus1sigDt2Eps = 1.0 + sigDt2Eps;

					field->Ex[i][j][k] = ((1.0 - sigDt2Eps) / plus1sigDt2Eps) * field->Ex[i][j][k] + 
						(dtEps / plus1sigDt2Eps) * ((field->Hz[i][j][k] - field->Hz[i][j-1][k]) * invdy - 
						(field->Hy[i][j][k] - field->Hy[i][j][k-1]) * invdz - Js[i][j][k]);

					EpsY = 0.25 * (Eps[i][j][k] + Eps[i][j][k - 1] + Eps[i - 1][j][k - 1] + Eps[i - 1][j][k]) * EPS_Z;
					SigY = 0.25 * (Sig[i][j][k] + Sig[i][j][k - 1] + Sig[i - 1][j][k - 1] + Sig[i - 1][j][k]);
					dtEps = dt / EpsY;
					sigDt2Eps = 0.5 * SigY * dtEps;
					plus1sigDt2Eps = 1.0 + sigDt2Eps;

					field->Ey[i][j][k] = ((1.0 - sigDt2Eps) / plus1sigDt2Eps) * field->Ey[i][j][k] + 
						(dtEps / plus1sigDt2Eps) * ((field->Hx[i][j][k] - field->Hx[i][j][k - 1]) * invdz - 
						(field->Hz[i][j][k] - field->Hz[i - 1][j][k]) * invdx - Js[i][j][k]);

					EpsZ = 0.25 * (Eps[i][j][k] + Eps[i - 1][j][k] + Eps[i - 1][j - 1][k] + Eps[i][j - 1][k]) * EPS_Z;
					SigZ = 0.25 * (Sig[i][j][k] + Sig[i - 1][j][k] + Sig[i - 1][j - 1][k] + Sig[i][j - 1][k]);
					dtEps = dt / EpsZ;
					sigDt2Eps = 0.5 * SigZ * dtEps;
					plus1sigDt2Eps = 1.0 + sigDt2Eps;

					field->Ez[i][j][k] = ((1.0 - sigDt2Eps) / plus1sigDt2Eps) * field->Ez[i][j][k] + 
						(dtEps / plus1sigDt2Eps) * ((field->Hy[i][j][k] - field->Hy[i - 1][j][k]) * invdx - 
						(field->Hx[i][j][k] - field->Hx[i][j - 1][k]) * invdx - Js[i][j][k]);

				}
			}

			//Hard coded dipole source
				double rtau = 50.0e-12;
				double tau = rtau / dt;
				double ndelay = 3 * tau;
				double srcconst = -dt * 3.0e+11;
				for(int k = 2; k < gridZ - 1; k++)
				{

					field->Ez[gridX / 2][gridY / 2][k] = field->Ez[gridX / 2][gridY / 2][k] + 
						srcconst * (globalTimestep - ndelay) * exp(-((globalTimestep - ndelay) * (globalTimestep - ndelay)/(tau * tau)));
				}

			//Calculate H - field at boundaries
			//i = 0, gridX - 1
			#pragma omp parallel for
			for(int i = 1; i < gridX - 1; i++)
			{
				double MjuX;
				double SigSX;
				double dtMju, sigDt2Mju, plus1sigDt2Mju;
				int j = 0;
				int k = 0;

				MjuX = Mju[i][j][k] * MU_Z;
				SigSX = SigS[i][j][k];
				dtMju = dt / MjuX;
				sigDt2Mju = 0.5 * SigSX * dtMju;
				plus1sigDt2Mju = 1.0 + sigDt2Mju;

				field->Hx[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hx[i][j][k] + 
					(dtMju / plus1sigDt2Mju) * ((field->Ey[i][j][k+1] - field->Ey[i][j][k]) * invdz - 
					(field->Ez[i][j + 1][k] - field->Ez[i][j][k]) * invdy - Ms[i][j][k]);

				/*i = gridX - 1;

				MjuX = 0.5 * (Mju[i][j][k] + Mju[i - 1][j][k]) * MU_Z;
				SigSX = 0.5 * (SigS[i][j][k] + SigS[i - 1][j][k]);
				dtMju = dt / MjuX;
				sigDt2Mju = 0.5 * SigSX * dtMju;
				plus1sigDt2Mju = 1.0 + sigDt2Mju;

				field->Hx[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hx[i][j][k] + 
					(dtMju / plus1sigDt2Mju) * ((field->Ey[i][j][k+1] - field->Ey[i][j][k]) * invdz - 
					(field->Ez[i][j + 1][k] - field->Ez[i][j][k]) * invdy - Ms[i][j][k]);*/
					
			}
			//j = 0, gridY - 1
			#pragma omp parallel for
			for(int j = 1; j < gridY - 1; j++)
			{
				double MjuY;
				double SigSY;
				double dtMju, sigDt2Mju, plus1sigDt2Mju;
				int i = 0;
				int k = 0;

				MjuY = Mju[i][j][k] * MU_Z;
				SigSY = SigS[i][j][k];
				dtMju = dt / MjuY;
				sigDt2Mju = 0.5 * SigSY * dtMju;
				plus1sigDt2Mju = 1.0 + sigDt2Mju;

				field->Hy[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hy[i][j][k] + 
					(dtMju / plus1sigDt2Mju) * ((field->Ez[i + 1][j][k] - field->Ez[i][j][k]) * invdx - 
					(field->Ex[i][j][k + 1] - field->Ex[i][j][k]) * invdz - Ms[i][j][k]);
						
				/*j = gridY - 1;

				MjuY = 0.5 * (Mju[i][j][k] + Mju[i][j - 1][k]) * MU_Z;
				SigSY = 0.5 * (SigS[i][j][k] + SigS[i][j - 1][k]);
				dtMju = dt / MjuY;
				sigDt2Mju = 0.5 * SigSY * dtMju;
				plus1sigDt2Mju = 1.0 + sigDt2Mju;

				field->Hy[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hy[i][j][k] + 
					(dtMju / plus1sigDt2Mju) * ((field->Ez[i + 1][j][k] - field->Ez[i][j][k]) * invdx - 
					(field->Ex[i][j][k + 1] - field->Ex[i][j][k]) * invdz - Ms[i][j][k]);*/
					
			}
			//k = 0, gridZ - 1
			#pragma omp parallel for
			for(int k = 1; k < gridZ - 1; k++)
			{
				double MjuZ;
				double SigSZ;
				double dtMju, sigDt2Mju, plus1sigDt2Mju;
				int i = 0;
				int j = 0;

				MjuZ = Mju[i][j][k] * MU_Z;
				SigSZ = SigS[i][j][k];
				dtMju = dt / MjuZ;
				sigDt2Mju = 0.5 * SigSZ * dtMju;
				plus1sigDt2Mju = 1.0 + sigDt2Mju;

				field->Hz[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hz[i][j][k] + 
					(dtMju / plus1sigDt2Mju) * ((field->Ex[i][j + 1][k] - field->Ex[i][j][k]) * invdy - 
					(field->Ey[i + 1][j][k] - field->Ey[i][j][k]) * invdx - Ms[i][j][k]);

				/*k = gridZ - 1;

				MjuZ = 0.5 * (Mju[i][j][k] + Mju[i][j][k - 1]) * MU_Z;
				SigSZ = 0.5 * (SigS[i][j][k] + SigS[i][j][k - 1]);
				dtMju = dt / MjuZ;
				sigDt2Mju = 0.5 * SigSZ * dtMju;
				plus1sigDt2Mju = 1.0 + sigDt2Mju;

				field->Hz[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hz[i][j][k] + 
					(dtMju / plus1sigDt2Mju) * ((field->Ex[i][j + 1][k] - field->Ex[i][j][k]) * invdy - 
					(field->Ey[i + 1][j][k] - field->Ey[i][j][k]) * invdx - Ms[i][j][k]);*/
			
			}
			//Calculate H - field at domain
			for(int i = 1; i < gridX - 1; i++)
				for(int j = 1; j < gridY - 1; j++)
				{
					#pragma omp parallel for
					for(int k = 1; k < gridZ - 1; k++)
					{
						double MjuX, MjuY, MjuZ;
						double SigSX, SigSY, SigSZ;
						double dtMju, sigDt2Mju, plus1sigDt2Mju;

						MjuX = 0.5 * (Mju[i][j][k] + Mju[i - 1][j][k]) * MU_Z;
						SigSX = 0.5 * (SigS[i][j][k] + SigS[i - 1][j][k]);
						dtMju = dt / MjuX;
						sigDt2Mju = 0.5 * SigSX * dtMju;
						plus1sigDt2Mju = 1.0 + sigDt2Mju;

						field->Hx[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hx[i][j][k] + 
							(dtMju / plus1sigDt2Mju) * ((field->Ey[i][j][k+1] - field->Ey[i][j][k]) * invdz - 
							(field->Ez[i][j + 1][k] - field->Ez[i][j][k]) * invdy - Ms[i][j][k]);

						MjuY = 0.5 * (Mju[i][j][k] + Mju[i][j - 1][k]) * MU_Z;
						SigSY = 0.5 * (SigS[i][j][k] + SigS[i][j - 1][k]);
						dtMju = dt / MjuY;
						sigDt2Mju = 0.5 * SigSY * dtMju;
						plus1sigDt2Mju = 1.0 + sigDt2Mju;

						field->Hy[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hy[i][j][k] + 
							(dtMju / plus1sigDt2Mju) * ((field->Ez[i + 1][j][k] - field->Ez[i][j][k]) * invdx - 
							(field->Ex[i][j][k + 1] - field->Ex[i][j][k]) * invdz - Ms[i][j][k]);

						MjuZ = 0.5 * (Mju[i][j][k] + Mju[i][j][k - 1]) * MU_Z;
						SigSZ = 0.5 * (SigS[i][j][k] + SigS[i][j][k - 1]);
						dtMju = dt / MjuZ;
						sigDt2Mju = 0.5 * SigSZ * dtMju;
						plus1sigDt2Mju = 1.0 + sigDt2Mju;

						field->Hz[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hz[i][j][k] + 
							(dtMju / plus1sigDt2Mju) * ((field->Ex[i][j + 1][k] - field->Ex[i][j][k]) * invdy - 
							(field->Ey[i + 1][j][k] - field->Ey[i][j][k]) * invdx - Ms[i][j][k]);

							
					}
			}
		globalTimestep++;
		Log::GetInstance().WriteLine("done.");
	}

	void CConservativeSolver3D::Solve(int timeSteps)
	{
		EMCrop crop;
		crop.up = 0;
		crop.down = 0;
		crop.left = 0;
		crop.right = 0;
		crop.top = 0;
		crop.bottom = 0;
		for(int t = 0; t < timeSteps; t++)
		{
			SolveStep();
			
			field->WriteFieldToBinary(Ez, crop, t, "C:\\Ez");
		}
	}

	CConservativeSolver3D::~CConservativeSolver3D(void)
	{
		DeleteJaggedArray3D(Eps, gridX, gridY, gridZ);
		DeleteJaggedArray3D(Mju, gridX, gridY, gridZ);
		DeleteJaggedArray3D(Sig, gridX, gridY, gridZ);
		DeleteJaggedArray3D(SigS, gridX, gridY, gridZ);
		DeleteJaggedArray3D(Ms, gridX, gridY, gridZ);
		DeleteJaggedArray3D(Js, gridX, gridY, gridZ);
	}

	void CConservativeSolver3D::prepare()
	{

	}

}