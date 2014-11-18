#define _USE_MATH_DEFINES
#include "ConservativeSolver3D.h"
#include "PhysicConst.h"
#include "math.h"
#include "Helpers.h"
#include "Log.h"
#include <sstream>
#include <iostream>

namespace EMWSolver
{
	CConservativeSolver3D::CConservativeSolver3D(CField3D* _field)
	{
		field = _field;
	}

	void CConservativeSolver3D::Create(TaskParameters taskParam, NumericalParameters numParam, SourceParameters srcParam, 
			BoundaryParameters* bndParam, int eps)
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

		gridX = numParam.sizeX + 1;
		gridY = numParam.sizeY + 1;
		gridZ = numParam.sizeZ + 1;

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

		for(int i = 0; i < field->sizeX / 3; i++)
			for(int j = 0; j < field->sizeY; j++)
				for(int k = 0; k < field->sizeZ; k++)
				{
					Eps[i][j][k] = 4;
					Mju[i][j][k] = 4;
				}

		Log::GetInstance().WriteLine("Finished.");

		

		logger << "Task size: " << gridX << " " << gridY << " " << gridZ;
		Log::GetInstance().WriteLine(logger.str());

		

		//prepare();
		//int aWg = 50;//numParam.Nx * 0.75;
		//int bWg = aWg / 2;

		//wgIndex = eps;
		//setupWaveguideGeometry(aWg, bWg, aWg / 4, 0 , eps);


		//startRawTime = 2500;
		//endRawTime = 2700;
		//int rawSize = endRawTime - startRawTime;
		//rawIterator = 0;
		//RawEyCenter = CreateJaggedArray2D(rawSize, gridZ);
		//RawEyUp = CreateJaggedArray2D(rawSize, gridZ);
		//RawEyDown = CreateJaggedArray2D(rawSize, gridZ);
		//RawEyLeft = CreateJaggedArray2D(rawSize, gridZ);
		//RawEyRight = CreateJaggedArray2D(rawSize, gridZ);

		//Ew = CreateJaggedArray3D(numParam.sizeX, numParam.sizeY, numParam.sizeZ);
		//Ey0 = 1;//M_PI/(_a*dx);
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
			//double rtau = 1.0 / FREQ;
			//double tau = rtau / dt;
			//double ndelay = 3 * tau;
			//double taper = 1 - exp(-((globalTimestep) * (globalTimestep)/(tau * tau)));

			//for(int i = 1; i < gridX-1; i++)
			//	for(int j = 1; j < gridY-1; j++)
			//	{
			//		#pragma omp parallel for
			//		for(int k = 1; k < gridZ - 1; k++)
			//		{
			//			//Ey0 = 1.0 / sqrt(Eps[i][j][k]);
			//			double Eyinc1 = Ey0 * cos(OMEGA * (globalTimestep + 0.5) * dt - gamma0_num * (k - gridZ/2 - l0/2) * dz) * sin(M_PI * (i - _x1) / (_a ));
			//			double Eyinc0 = Ey0 * cos(OMEGA * (globalTimestep - 0.5) * dt - gamma0_num * (k- gridZ/2 - l0/2) * dz) * sin(M_PI * (i - _x1)  / (_a ));
			//			Ew[i][j][k] = Eyinc1;
			//			field->Ey[i][j][k] = field->Ey[i][j][k] - (1.0 - 1.0/(0.25 * (Eps[i][j][k] + Eps[i][j][k - 1] + Eps[i - 1][j][k - 1] + Eps[i - 1][j][k])))*
			//				(Eyinc1 - Eyinc0) - (dt * Sig[i][j][k] / (0.25 * (Eps[i][j][k] + Eps[i][j][k - 1] + Eps[i - 1][j][k - 1] + Eps[i - 1][j][k])*EPS_Z))*Eyinc1;
			//		}
			//	}

			//for(int i = _x1; i <_x2; i++)
			//	for(int j = _y1; j <_y2; j++)
			//	{
			//		int k = gridZ / 2 - 100;
			//		field->Ey[i][j][k] = taper * cos(OMEGA * (globalTimestep) * dt - gamma0 * k * dz) * sin(M_PI * (i - _x1) * dx / (_a * dx));
			//		//data->Ey[i][j][1] = data->Ey[i][j][0];
			//	}	
			
		for(int k = 0; k < field->gridZ; k++)
		{
			field->Ez[field->gridX/2][field->gridY/2][k] = cos(OMEGA* globalTimestep * dt);
		}
			//Calculate H - field at domain
			for(int i = 0; i < field->sizeX; i++)
				for(int j = 0; j < field->sizeY; j++)
				{
					#pragma omp parallel for
					for(int k = 0; k < field->sizeZ; k++)
					{
						double MjuX, MjuY, MjuZ;
						double SigSX, SigSY, SigSZ;
						double dtMju, sigDt2Mju, plus1sigDt2Mju;

						if(i != 0)
						{
							MjuX = 0.5 * (Mju[i][j][k] + Mju[i - 1][j][k]) * MU_Z;
							SigSX = 0.5 * (SigS[i][j][k] + SigS[i - 1][j][k]);
						}
						else
						{
							MjuX = Mju[i][j][k] * MU_Z;
							SigSX = SigS[i][j][k];
						}
						dtMju = dt / MjuX;
						sigDt2Mju = 0.5 * SigSX * dtMju;
						plus1sigDt2Mju = 1.0 + sigDt2Mju;

 						field->Hx[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hx[i][j][k] + 
							(dtMju / plus1sigDt2Mju) * ((field->Ey[i][j][k+1] - field->Ey[i][j][k]) * invdz - 
							(field->Ez[i][j + 1][k] - field->Ez[i][j][k]) * invdy - Ms[i][j][k]);


						if(j != 0)
						{
							MjuY = 0.5 * (Mju[i][j][k] + Mju[i][j - 1][k]) * MU_Z;
							SigSY = 0.5 * (SigS[i][j][k] + SigS[i][j - 1][k]);
						}
						else
						{
							MjuY = Mju[i][j][k] * MU_Z;
							SigSY = SigS[i][j][k];
						}
						dtMju = dt / MjuY;
						sigDt2Mju = 0.5 * SigSY * dtMju;
						plus1sigDt2Mju = 1.0 + sigDt2Mju;


						field->Hy[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hy[i][j][k] + 
							(dtMju / plus1sigDt2Mju) * ((field->Ez[i + 1][j][k] - field->Ez[i][j][k]) * invdx - 
							(field->Ex[i][j][k + 1] - field->Ex[i][j][k]) * invdz - Ms[i][j][k]);

						

						if(k != 0)
						{
							MjuZ = 0.5 * (Mju[i][j][k] + Mju[i][j][k - 1]) * MU_Z;
							SigSZ = 0.5 * (SigS[i][j][k] + SigS[i][j][k - 1]);
						}
						else
						{
							MjuZ = Mju[i][j][k - 1]* MU_Z;
							SigSZ = SigS[i][j][k];
						}
						dtMju = dt / MjuZ;
						sigDt2Mju = 0.5 * SigSZ * dtMju;
						plus1sigDt2Mju = 1.0 + sigDt2Mju;

						field->Hz[i][j][k] = ((1.0 - sigDt2Mju) / plus1sigDt2Mju) * field->Hz[i][j][k] + 
							(dtMju / plus1sigDt2Mju) * ((field->Ex[i][j + 1][k] - field->Ex[i][j][k]) * invdy - 
							(field->Ey[i + 1][j][k] - field->Ey[i][j][k]) * invdx - Ms[i][j][k]);

					}
			}

		//Calculate E - field at domain
		for(int i = 1; i < field->sizeX ; i++)
			for(int j = 1; j < field->sizeY; j++)
			{
				#pragma omp parallel for
				for(int k = 1; k < field->sizeZ; k++)
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

			//Boundary -X
			for(int j = 1; j < field->sizeY; j++)
			{
				#pragma omp parallel for
				for(int k = 1; k < field->sizeZ; k++)
				{
					double EpsX, EpsY, EpsZ;
					double SigX, SigY, SigZ;
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
				}
			}
			//Boundary -Y
			for(int i = 1; i < field->sizeX ; i++)
					#pragma omp parallel for
					for(int k = 1; k < field->sizeZ; k++)
					{
						double EpsX, EpsY, EpsZ;
						double SigX, SigY, SigZ;
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
				}

			//Boundary -Z
			for(int i = 1; i < field->sizeX ; i++)
			{
				#pragma omp parallel for
				for(int j = 1; j < field->sizeY; j++)
				{		
						int k = 0;
						double EpsX, EpsY, EpsZ;
						double SigX, SigY, SigZ;
						double dtEps, sigDt2Eps, plus1sigDt2Eps;
				
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
			//if(t >= 500)
				//if(t%10 == 0)
					//field->WriteFieldToBinary(Ey, crop, t, "D:\\Ez");

			//if( t>= startRawTime && t < endRawTime)
			//	storeRawEy(gridX / 2, gridY / 2);
		}
		//calculateEyAmp(RawEyCenter, 0, gridZ);
	}

	CConservativeSolver3D::~CConservativeSolver3D(void)
	{
		int sizeX = field->sizeX;
		int sizeY = field->sizeY;
		int sizeZ = field->sizeZ;
		field = NULL;
		DeleteJaggedArray3D(Eps, sizeX, sizeY, sizeZ);
		DeleteJaggedArray3D(Mju, sizeX, sizeY, sizeZ);
		DeleteJaggedArray3D(Sig, sizeX, sizeY, sizeZ);
		DeleteJaggedArray3D(SigS, sizeX, sizeY, sizeZ);
		DeleteJaggedArray3D(Ms, sizeX, sizeY, sizeZ);
		DeleteJaggedArray3D(Js, sizeX, sizeY, sizeZ);

		//int rawSize = endRawTime - startRawTime;
		//DeleteJaggedArray2D(RawEyCenter, rawSize, gridZ);
		//DeleteJaggedArray2D(RawEyUp, rawSize, gridZ);
		//DeleteJaggedArray2D(RawEyDown, rawSize, gridZ);
		//DeleteJaggedArray2D(RawEyLeft, rawSize, gridZ);
		//DeleteJaggedArray2D(RawEyRight, rawSize, gridZ);

		//DeleteJaggedArray3D(Ew, gridX, gridY, gridZ);
	}

	void CConservativeSolver3D::prepare()
	{

	}

	void CConservativeSolver3D::setupWaveguideGeometry(int aSize, int bSize, int length, int bOff, double leps)
	{
		int Nxh = gridX / 2;
		int Nyh = gridY / 2;
		int Nzh = gridZ / 2;
		_a = aSize;
		_b = bSize;
		double L0 = length;
		l0 = length;
		double blockSize = 12;
		double L1 = L0;
		double L2 = 2.5 * bSize;
		double _w = 2.5;
		double blockOffsetZ = bOff;
		_x1 = Nxh - _a/2;
		_x2 = Nxh + _a/2;
		_y1 = Nyh - _b/2;
		_y2 = Nyh + _b/2;

		_z1 = Nzh - length / 2;
		_z2 = Nzh + length / 2;

		int startZ = Nzh - L0/2;

		int aThick = _a/2 - blockSize/2;
		int bThick = _b/2 - blockSize/2;
		int cThick = 0;//5;
		int ofs = 0;
		double gamma0 = sqrt(OMEGA*OMEGA*EPS_Z*MU_Z - M_PI*M_PI/(_a*_a * dx * dx));
		double gamma1 = sqrt(OMEGA*OMEGA*EPS_Z*MU_Z*leps - M_PI*M_PI/(_a *_a * dx * dx));
		double k0 = OMEGA*OMEGA*EPS_Z*MU_Z;

		double Vf = CC / sqrt(1.0 - (M_PI * M_PI * CC * CC / (_a * _a * dx * dx * OMEGA * OMEGA)));
		double Krlambda = 2 * M_PI * Vf / OMEGA;
		double WGlambda = LAMBDA / sqrt(1 - LAMBDA*LAMBDA/(2 * _a * 2 * _a * dx * dx));
		double SPlambda = LAMBDA;
		int perL = LAMBDA / dx;
		int perW = WGlambda / dx;
		if ((WGlambda < 2*_a*dx) || (WGlambda > 4*_a*dx))
		{
			std::cout << "ERROR! Lambda = " << WGlambda << std::endl;
			std::cout << "2a = " << 2*_a*dx << std::endl;
			system("PAUSE");
		}
		//double Nlz = Nl;
		double Na = _a;
		double Sz = CC * dt / dz;
		double g0A = (1.0 / Sz) * sin(Sz * M_PI * k0 / (L0 * gamma0)) ;
		//g0A *= g0A;
		double g0B = (2.0 * Na / (gamma0 * Na * dx * L0)) * sin(M_PI / (2.0 * Na));
		g0B *= g0B;
		double g1A = sin(OMEGA*dt/2.0) / dt;
		double g1B = sin(M_PI * dx / (2.0 * _a * dx)) / dx;
		g1A *= g1A;
		g1B *= g1B;
		gamma0_num = gamma0 * (L0 / M_PI) * asin(sqrt(g0A - g0B)); 
		gamma0_num = (2.0/dz) * asin(dz * sqrt(EPS_Z*MU_Z * g1A - g1B));
		gamma1_num = (2.0/dz) * asin(dz * sqrt(leps*EPS_Z*MU_Z * g1A - g1B));

		std::cout << "g0 = " << gamma0 << std::endl;
		std::cout << "g1 = " << gamma1 << std::endl;
		std::cout << "g0num = " << gamma0_num << std::endl;
		std::cout << "g1num = " << gamma1_num << std::endl;

		double LambdaK = 2*_a*dz;
		double A1 = 1;
		//double gamma1 = sqrt(OMEGA*OMEGA*EPS_Z*MU_Z*leps - M_PI*M_PI/(_a*_a));
		double regz = cos(gamma1 * L0 * dz);
		double imgz = sin(gamma1 * L0 * dz * (1 + (gamma1/gamma0)*(gamma1/gamma0))/2.0);
		double modF = cos(gamma1*L0*dx)*cos(gamma1*L0*dx) + (gamma1*0.5/gamma0 + gamma0*0.5/gamma1)*(gamma1*0.5/gamma0 + gamma0*0.5/gamma1);
		double ReF = (cos(gamma0 * L0 * dz) * regz + sin(gamma0 * L0 * dz) * imgz)/(regz*regz + imgz*imgz);
		double ImF = (-cos(gamma0 * L0 * dz) * imgz + sin(gamma0 * L0 * dz) * regz)/(regz*regz + imgz*imgz);

		setupWaveguideMaterial(length, bOff, leps);
		//system("pause");
	}

	void CConservativeSolver3D::setupWaveguideMaterial(int length, int bOff, double leps)
	{
		int Nzh = gridZ / 2;
		int _z1 = Nzh - length / 2;
		int _z2 = Nzh + length / 2;

		double L0 = length;
		int startZ = Nzh - length/2;
		int endZ = Nzh + length/2;
		for(int i = _x1; i < _x2; i++)
			for(int j = _y1; j < _y2; j++)
				for(int k = startZ; k < endZ; k++)
				{
					Eps[i][j][k] = leps;
				}
	}

	void CConservativeSolver3D::storeRawEy(int i, int  j)
	{
		for(int k = 0; k < gridZ; k++)
		{
			RawEyCenter[rawIterator][k] = field->Ey[i][j][k];// +  Ew[i][j][k];
			RawEyUp[rawIterator][k] = field->Ey[i + 20][j][k];//  + Ew[i + 20][j][k];
			RawEyDown[rawIterator][k] = field->Ey[i - 20][j][k] ;// + Ew[i - 20][j][k];
			RawEyLeft[rawIterator][k] = field->Ey[i][j - 10][k] ;// + Ew[i][j - 10][k];
			RawEyRight[rawIterator][k] = field->Ey[i][j + 10][k] ;// + Ew[i][j + 10][k];
		}

		rawIterator++;
	}

	void CConservativeSolver3D::calculateEyAmp(double** inputAmp, int startZ, int endZ)
	{
	   double* amp;
	   double* ampIM;
	   double* ReF;
	   double* ImF;
	   double* modA;

	   int zdim = endZ - startZ;
	   int nmax = endRawTime - startRawTime;

	   int period = 1.0 / FREQ / dt;
	   amp = new double[zdim];
	   ampIM = new double[zdim];
	   ReF = new double[zdim];
	   ImF = new double[zdim];
	   modA = new double[zdim];

	   //fftw_plan p,s;
	   //cout << "Creating FFTW plan." << endl;
	   //outP = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nmax/2 + 1));
	   //outS = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nmax/2 + 1));
	   //s =  fftw_plan_dft_r2c_1d(nmax, inputAmpS, outS, FFTW_ESTIMATE);
	   //p =  fftw_plan_dft_r2c_1d(nmax, inputAmpP, outP, FFTW_ESTIMATE);
	   //cout << "FFTW Transform...";
	   //fftw_execute(s);
	   //fftw_execute(p); /* repeat as needed */
	   //cout << "done." << endl;
	   //fftw_complex* tem1;
	   //fftw_complex* tem2;
	   int kharm=1;
	   double beta = 1.535;
	   double dnmax = period;

	   for(int i = 0; i <zdim; i++)
	   {
		   amp[i] = 0;
		   ampIM[i] = 0;

		   if(period > nmax)
			   std::cout << "Not enough time data" << std::endl;
		   for(int rec = 1; rec <= period; rec++)
		   {
			   double drec = startRawTime + rec;
			   //inputAmp[rec-1][i]+=cos(OMEGA*dt*start -gamma0 * i * dx);
			   //inputAmp[rec][i]+=sin(OMEGA*dt*start -gamma0 * i * dx);
	//			amp[i] += inputAmp[rec][i]*(cos(2*M_PI*drec/dnmax));//+sin(2*M_PI*kharm*rec/nmax));
			   //Trapezoidal method requires previous time step for second order accurace fourier transform
				amp[i] += 0.5*(inputAmp[rec][i+startZ]*cos(2*M_PI*drec/dnmax)+inputAmp[rec-1][i+startZ]*cos(2*M_PI*(drec-1.0)/dnmax));//+sin(2*M_PI*kharm*rec/nmax));
				ampIM[i] += 0.5*(inputAmp[rec][i+startZ]*sin(2*M_PI*drec/dnmax)+inputAmp[rec-1][i+startZ]*sin(2*M_PI*(drec-1.0)/dnmax));//+sin(2*M_PI*kharm*rec/nmax));
		   }

		    amp[i]=(2.0 * (amp[i])/ (dnmax));+ cos( - gamma0_num * (i - gridZ/2 - l0/2) * dz);
			ampIM[i]=(2.0 * (ampIM[i])/ (dnmax)); + sin( - gamma0_num * (i - gridZ/2 - l0/2) * dz);;

			 ReF[i] = -(sin(- gamma0_num * (i - gridZ/2 - l0/2) * dz) * ampIM[i] - 
					cos(- gamma0_num * (i - gridZ/2 - l0/2) * dz) * amp[i])/cos(2*(- gamma0_num * (i - gridZ/2 - l0/2) * dz));
			   if(fabs(ReF[i]) > 20) ReF[i] = 0;
			   //else
				   //ReF[i] = 0;
			   //cout << ReF[i] <<endl;
			   //if(cos(2*gamma0*i*dx) <= 0.000001)
			  ImF[i] = (cos(- gamma0_num * (i - gridZ/2 - l0/2) * dz) * ampIM[i] - 
					sin(- gamma0_num * (i - gridZ/2 - l0/2) * dz) * amp[i])/cos(2*(- gamma0_num * (i - gridZ/2 - l0/2) * dz));

			modA[i] = sqrt(ReF[i]*ReF[i] + ImF[i]*ImF[i]);
		}
	   WriteJaggedArray1DToFile(modA, zdim, "D:\\Ez", "modF", wgIndex);
	   WriteJaggedArray1DToFile(ReF, zdim, "D:\\Ez", "ReF", wgIndex);
	   WriteJaggedArray1DToFile(ImF, zdim, "D:\\Ez", "ImF", wgIndex);
	}

}