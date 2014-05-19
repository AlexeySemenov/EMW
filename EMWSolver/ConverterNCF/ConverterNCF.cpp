#define _USE_MATH_DEFINES
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include "netcdf.h"
#include "fftw3.h"
#include "PhysicConst.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

const int gridX = 100;
const int gridY = 50;
const int gridZ = 1500;

double ez[gridX][gridY][gridZ];
double ex[gridX][gridY][gridZ];
double ey[gridX][gridY][gridZ];


double xDim[gridX];
double yDim[gridY];
double zDim[gridZ];
double* data;

FILE* dFile;

const double dx = LAMBDA/20.0;
const double dt = 0.99*dx/(sqrt(3.0)*CC);
using namespace std;
void readFile(int timestp, string path)
{   
 //   string fileName;
	//stringstream s;
	//s << "K:\\Ez\\Ez[" << timestp << "].dat";
	//fileName = s.str();
	// ifstream file(fileName.c_str(), ios::in | ios::binary);

	// string fileName2;
	//stringstream s2;
	// s2 << path<<"\\Ez[" << timestp << "].dat";
	//fileName2 = s2.str();
	// ifstream file2(fileName2.c_str(), ios::in | ios::binary);

	// string fileName3;
	//stringstream s3;
	// s3 << path<<"\\Ex[" << timestp << "].dat";
	//fileName3 = s3.str();
	// ifstream file3(fileName3.c_str(), ios::in | ios::binary);

	 string fileName4;
	stringstream s4;
	 s4 << path<<"\\Ey[" << timestp << "]" << ".dat";
	fileName4 = s4.str();
	 ifstream file4(fileName4.c_str(), ios::in | ios::binary);


	 double tmp1,tmp2, tmp3, tmp4;
	 //if(!file.is_open())
		//cout << "ERROR::Cannot open file for read! " << fileName << endl ;
	//if(!file2.is_open())
	//	cout << "ERROR::Cannot open file for read! " << fileName2 << endl ;

	//if(!file3.is_open())
	//	cout << "ERROR::Cannot open file for read! " << fileName3 << endl ;

	if(!file4.is_open())
		cout << "ERROR::Cannot open file for read! " << fileName4 << endl ;

	cout << "Reading file... " << timestp <<"...";

	for(int i = 0; i < gridX; i++)
		for(int j = 0; j < gridY; j++)
			for(int k = 0; k < gridZ; k++)
			{
				//file2.read((char *) &tmp2, sizeof(double));
				//file3.read((char *) &tmp3, sizeof(double));
				file4.read((char *) &tmp4, sizeof(double));
				//ez[i][j][k] = tmp2;
				//ex[i][j][k] = tmp3;
				ez[i][j][k] = tmp4;

			}
	cout << gridX*gridY*gridZ*sizeof(double) << " BYTE read succesfully" << endl;
	//file2.close();
	//file3.close();
	file4.close();

	     
};

int writeNCF(int from, int to)
{
	 /* IDs for the netCDF file, dimensions, and variables. */
   int ncid, z_dimid, y_dimid, x_dimid, t_dimid;
   int x_varid, y_varid, z_varid, ez_varid, ex_varid, ey_varid, t_varid;
   int dimids[4];

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
   size_t start[4], count[4];
   size_t* start2;
   size_t* count2;


   /* Loop indexes. */
   int lvl, lat, lon, rec, i = 0;
   
   /* Error handling. */
   int retval;

   string fileName;
     stringstream s;
	 s << "C:\\ncf\\testEMW.nc";
	 fileName = s.str();

   /* Create the file. */
	 if ((retval = nc_create(fileName.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid)))
      ERR(retval);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is

    * the time dimension.*/
	    if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &t_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "xdim", gridX, &x_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "ydim", gridY, &y_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "zdim", gridZ, &z_dimid)))
      ERR(retval);


   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&lat_dimid) and
      similarly for (&lon_dimid). */
   if ((retval = nc_def_var(ncid, "xdim", NC_DOUBLE, 1, &x_dimid, 
			    &x_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "ydim", NC_DOUBLE, 1, &y_dimid, 
			    &y_varid)))
      ERR(retval);

   if ((retval = nc_def_var(ncid, "zdim", NC_DOUBLE, 1, &z_dimid, 
			    &z_varid)))
      ERR(retval);




   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = t_dimid;
   dimids[1] = x_dimid;
   dimids[2] = y_dimid;
   dimids[3] = z_dimid;

   /* Define the netCDF variables for the pressure and temperature
    * data. */
   //if ((retval = nc_def_var(ncid, "Ez", NC_DOUBLE, 4, 
			//    dimids, &ez_varid)))
   //   ERR(retval);

   //if ((retval = nc_def_var(ncid, "Ex", NC_DOUBLE, 4, 
			//    dimids, &ex_varid)))
   //   ERR(retval);

   if ((retval = nc_def_var(ncid, "Ez", NC_DOUBLE, 4, 
			    dimids, &ey_varid)))
      ERR(retval);

   //if ((retval = nc_def_var(ncid, "Ez", NC_DOUBLE, 3, 
			//    dimids, &ez_varid)))
   //   ERR(retval);



   /* End define mode. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

   /* Write the coordinate variable data. This will put the latitudes
      and longitudes of our data grid into the netCDF file. */
   if ((retval = nc_put_var_double(ncid, x_varid, &xDim[0])))
      ERR(retval);
   if ((retval = nc_put_var_double(ncid, y_varid, &yDim[0])))
      ERR(retval);
   if ((retval = nc_put_var_double(ncid, z_varid, &zDim[0])))
      ERR(retval);

   /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
   count[0] = 1;
   count[1] = gridX;
   count[2] = gridY;
   count[3] = gridZ;
   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   start[3] = 0;

   start2 = new size_t;
   count2 = new size_t;
   *start2 = 1;
   *count2 = 1;
   double* re = new double;
   *re = 1;

   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays only hold one timestep worth
      of data. We will just rewrite the same data for each timestep. In
      a real application, the data would change between timesteps. */
   int nmax = to - from;
   for (rec = 0; rec < nmax; rec++)
   {
      start[0] = rec;
	  readFile(rec+from, "C:\\Ez");
  //    if ((retval = nc_put_vara_double(ncid, ez_varid, start, count, 
		//		      &ez[0][0][0])))
		//ERR(retval);
	 // if ((retval = nc_put_vara_double(ncid, ex_varid, start, count, 
		//		      &ex[0][0][0])))
		//ERR(retval);
	  if ((retval = nc_put_vara_double(ncid, ey_varid, start, count, 
				      &ez[0][0][0])))
		ERR(retval);
	  *start2 += 1;
	  *re += 1;
	  //if ((retval = nc_put_vara_double(ncid, t_varid, start2, count2, &(*re))))
   //   ERR(retval);

	 // if ((retval = nc_put_vara_double(ncid, ez_varid, start, count, 
		//		      &ez[0][0])))
	 //ERR(retval);

   }

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);
   
   printf("*** SUCCESS writing  file %s!\n", "Hz.nc");
   return 0;
}

int countAmp(int start, int end, int arstart, int arend, char* name, int id)
{
	int nmax = end - start;
	int zdim = (arend - arstart)+1;
	int adim = 80;
	double** inputAmp;
	inputAmp = new double*[nmax+1];
	for (int rec = 0; rec <= nmax; rec++)
	{
		inputAmp[rec] = new double[zdim];
	}
	int px,py,pz;
	int sx,sy,sz;

   double gamma0 = sqrt(OMEGA*OMEGA*EPS_Z*MU_Z - M_PI*M_PI/(adim*adim));
   double g1A = sin(OMEGA*dt/2.0) / dt;
	double g1B = sin(M_PI * dx / (2.0 * adim * dx)) / dx;
	g1A *= g1A;
	g1B *= g1B; 
	gamma0 = (2.0/dx) * asin(dx * sqrt(EPS_Z*MU_Z * g1A - g1B));

   cout <<" GAMMA0 = " << gamma0 << endl;
   //system("PAUSE");

   for (int rec = 0; rec <= nmax; rec++)
   {
	  readFile(rec+start, "G:\\Ez");
	  
	  for (int k = 0; k < zdim; k++)
		inputAmp[rec][k] =ey[gridX/2][gridY/2][arstart+k];//cos(sqrt(OMEGA*OMEGA*EPS_Z*MU_Z - M_PI*M_PI/(adim*adim))*k*dx - OMEGA*rec*dt);//ey[50][50][k];
   }

   //fftw_complex *outP, *outS;
   double* amp;
   double* ampIM;
   double* ReF;
   double* ImF;
   double* ExactCos;
   double* ExactSin;
   double* NumerCos;
   double* NumerSin;
   amp = new double[zdim];
   ampIM = new double[zdim];
   ReF = new double[zdim];
   ImF = new double[zdim];
   ExactCos = new double[zdim];
   ExactSin = new double[zdim];
   NumerCos = new double[zdim];
   NumerSin = new double[zdim];
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
   double dnmax = nmax;

   for(int i = 0; i <zdim; i++)
   {
	   amp[i] = 0;
	   ampIM[i] = 0;
	   ReF[i] = 0;
	   ImF[i] = 0;
	   ExactCos[i] = 0;
	   ExactSin[i] = 0;
	   NumerCos[i] = 0;
	   NumerSin[i] = 0;
	   for(int rec = 1; rec <= nmax; rec++)
	   {
		   double drec = rec;
		   //inputAmp[rec-1][i]+=cos(OMEGA*dt*start -gamma0 * i * dx);
		   //inputAmp[rec][i]+=sin(OMEGA*dt*start -gamma0 * i * dx);
//			amp[i] += inputAmp[rec][i]*(cos(2*M_PI*drec/dnmax));//+sin(2*M_PI*kharm*rec/nmax));
		   //Trapezoidal method requires previous time step for second order accurace fourier transform
		    amp[i] += 0.5*(inputAmp[rec][i]*cos(2*M_PI*drec/dnmax)+inputAmp[rec-1][i]*cos(2*M_PI*(drec-1.0)/dnmax));//+sin(2*M_PI*kharm*rec/nmax));
			ampIM[i] += 0.5*(inputAmp[rec][i]*sin(2*M_PI*drec/dnmax)+inputAmp[rec-1][i]*sin(2*M_PI*(drec-1.0)/dnmax));//+sin(2*M_PI*kharm*rec/nmax));
	   }
	   amp[i]=(2.0 * (amp[i])/ (dnmax));
	   //cout << amp[i] <<endl;
	   ampIM[i]=(2.0 * (ampIM[i])/ (dnmax));
	   //if(cos(2*gamma0*i*dx) <= 0.000001)

	  // amp[i] += cos(OMEGA*dt*start -gamma0 * i * dx);
	   //ampIM[i] += sin(OMEGA*dt*start -gamma0 * i * dx);
	   ReF[i] = -(sin(OMEGA*dt*start -gamma0 * (i) * dx) * ampIM[i] - 
					cos(OMEGA*dt*start -gamma0 * (i-4) * dx) * amp[i])/cos(2*(OMEGA*dt*start -gamma0*i*dx));
	   if(fabs(ReF[i]) > 20) ReF[i] = 0;
	   //else
		   //ReF[i] = 0;
	   //cout << ReF[i] <<endl;
	   //if(cos(2*gamma0*i*dx) <= 0.000001)
	  ImF[i] = (cos(OMEGA*dt*start -gamma0 * i * dx) * ampIM[i] - 
					sin(OMEGA*dt*start -gamma0 * i * dx) * amp[i])/cos(2*(OMEGA*dt*start -gamma0*i*dx));
	  if(fabs(ImF[i]) > 20) ImF[i] = 0;
	   //else
		   //ImF[i] = 0;

	   NumerCos[i] = -(sin(OMEGA*dt*start -gamma0 * (i+4+2*M_PI) * dx) * ampIM[i] - 
					cos(OMEGA*dt*start -gamma0 * (i+4) * dx) * amp[i]);//cos(gamma0 * (i-1) * dx) * amp[i];
	   NumerSin[i] = -(cos(OMEGA*dt*start -gamma0 * (i+4) * dx) * ampIM[i] - 
					sin(OMEGA*dt*start -gamma0 * (i+4+2*M_PI) * dx) * amp[i]);//sin(gamma0 * (i-1) * dx) * ampIM[i];
	   ExactCos[i] = -(sin(OMEGA*dt*start -gamma0 * (i+4+2*M_PI) * dx) * sin(OMEGA*dt*start - gamma0 * (i+4+2*M_PI) * dx) - 
					cos(OMEGA*dt*start -gamma0 * (i+4) * dx) * cos(OMEGA*dt*start -gamma0 * (i+4) * dx));///cos(2*gamma0*i*dx);//cos(-gamma0 * (i-1) * dx + OMEGA * start * dt) * cos(gamma0 * (i-1) * dx);
	   ExactSin[i] = -(cos(OMEGA*dt*start - gamma0 * (i+4) * dx) * sin(OMEGA*dt*start -gamma0 * (i+4+2*M_PI) * dx) - 
					sin(OMEGA*dt*start -gamma0 * (i+4+2*M_PI) * dx) * cos(OMEGA*dt*start -gamma0 * (i+4) * dx));//cos(-gamma0 * (i-1) * dx + OMEGA * start * dt) * sin(gamma0 * (i-1) * dx);
	   //cout << ampIM[i] <<endl;
	   //cout << ImF[i] <<endl;
	  //// cout << (double)((outP[i])[1]) <<endl;
   }
    double Ff1 = 0;
	double Ff2 = 0;
   for(int i = 11; i <= 42; i++)
	   Ff1 += NumerCos[i];

   for(int i = 11; i <= 42; i++)
	   Ff2+= NumerSin[i];

   Ff1 /= 32.0;
   Ff2 /= 32.0;
   cout << "FinalF1 = " << Ff1 << endl;
   cout << "FinalF2 = " << Ff2 << endl;
   double maxF = 0;
   int maxFind1[4] = {0, 0, 0, 0};
   //int maxFind1[2] = {0, 0};
   bool increase = amp[0] < amp[1];
   int findIterator = 0;
   for(int i = 0; i < zdim-1; i++)
   {

	   if(amp[i] < amp[i+1])
	   {
		   if(!increase) {maxFind1[findIterator] = i; findIterator++;}
			   increase = true;
	   }
	   else 
	   {
		   if(increase) {maxFind1[findIterator] = i; findIterator++;}
		   increase = false;
	   }

   }
	   //if(fabs(amp[i]) > maxF) { maxF = fabs(amp[i]); maxFind = i; }

   double mreF = 0.25*(fabs(amp[maxFind1[0]]) + fabs(amp[maxFind1[1]]) + fabs(amp[maxFind1[2]]) + fabs(amp[maxFind1[3]]))*2.274/M_PI;//(fabs(amp[212])+fabs(amp[276]))/2.0;//(fabs(amp[17])+fabs(amp[38])+fabs(amp[58])+fabs(amp[79]))/4.0;
   double mimF = 0.25*(ampIM[maxFind1[0]]+ampIM[maxFind1[1]]+ampIM[maxFind1[2]]+ampIM[maxFind1[3]])*2.274/M_PI;//(fabs(ampIM[180])+fabs(ampIM[246]))/2.0;//(fabs(ampIM[7])+fabs(ampIM[27])+fabs(ampIM[48])+fabs(ampIM[68]))/4.0;
   printf("%lf %lf\n", mreF, mimF);
   fprintf(dFile, "%lf %lf\n", mreF, mimF);
   int ncid, omega_dimid;
   int ey_varid, omega_varid, ampim_varid, ReF_varid, ImF_varid, 
	   ecos_varid, esin_varid, ncos_varid, nsin_varid;
   int dimids[1];

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
   size_t start[1], count[1];


   /* Loop indexes. */
   int lvl, lat, lon, rec, i = 0;
   
   /* Error handling. */
   int retval;

   string fileName;
     stringstream str;
	 str << "C:\\ncf\\WGAmpsBlockVarEps" << id << ".nc";
	 fileName = str.str();

   /* Create the file. */
	 if ((retval = nc_create(fileName.c_str(), NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/

   if ((retval = nc_def_dim(ncid, "omega",zdim, &omega_dimid)))
      ERR(retval);

   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = omega_dimid;

   /* Define the netCDF variables for the pressure and temperature
    * data. */
   if ((retval = nc_def_var(ncid, "Amp", NC_DOUBLE, 1, 
			    dimids, &ey_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "AmpIM", NC_DOUBLE, 1, 
			    dimids, &ampim_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "ReF", NC_DOUBLE, 1, 
			    dimids, &ReF_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, "ImF", NC_DOUBLE, 1, 
			    dimids, &ImF_varid)))
      ERR(retval);

   if ((retval = nc_def_var(ncid, "ExactCos", NC_DOUBLE, 1, 
			    dimids, &ecos_varid)))
      ERR(retval);

   if ((retval = nc_def_var(ncid, "ExactSin", NC_DOUBLE, 1, 
			    dimids, &esin_varid)))
      ERR(retval);

   if ((retval = nc_def_var(ncid, "NumerCos", NC_DOUBLE, 1, 
			    dimids, &ncos_varid)))
      ERR(retval);

   if ((retval = nc_def_var(ncid, "NumerSin", NC_DOUBLE, 1, 
			    dimids, &nsin_varid)))
      ERR(retval);


   if ((retval = nc_def_var(ncid, "Omega", NC_DOUBLE, 1, 
			    dimids, &omega_varid)))
      ERR(retval);



   /* End define mode. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);


      if ((retval = nc_put_var_double(ncid, ey_varid,  &amp[0])))
	 ERR(retval);

	  if ((retval = nc_put_var_double(ncid, ampim_varid,  &ampIM[0])))
	 ERR(retval);

	  if ((retval = nc_put_var_double(ncid, ReF_varid,  &ReF[0])))
	 ERR(retval);

	  if ((retval = nc_put_var_double(ncid, ImF_varid,  &ImF[0])))
	 ERR(retval);

	  if ((retval = nc_put_var_double(ncid, ecos_varid,  &ExactCos[0])))
	 ERR(retval);

	  if ((retval = nc_put_var_double(ncid, esin_varid,  &ExactSin[0])))
	 ERR(retval);

	  if ((retval = nc_put_var_double(ncid, ncos_varid,  &NumerCos[0])))
	 ERR(retval);

	  if ((retval = nc_put_var_double(ncid, nsin_varid,  &NumerSin[0])))
	 ERR(retval);

	  if ((retval = nc_close(ncid)))
      ERR(retval);
   
   printf("*** SUCCESS writing  file!");



   //cin >> px;
}

int main(int argc, char** argv)
{
	for(int i = 0; i < gridX; i++)
		xDim[i] = i/80.0;

	for(int i = 0; i < gridY; i++)
		yDim[i] = i/80.0;

	for(int i = 0; i < gridZ; i++)
		zDim[i] = i/80.0;

	//dFile = fopen ("ampfileSizeVarEps1602.txt","w");
	writeNCF(0, 98);
	//for(int i = 2; i <= 20; i++)
		//countAmp(760,970, 620, 840,"ey", i+1600);

	//fclose(dFile);
	//cout << 1.0/(FREQ*dt);
	//CountQ(2000);
	//double maxQ = 0;
	//double avrQ = 0;

	//for(int t = 960; t < 999; t++)
	//{
	//readFile(t);
	//for(int i = 0; i < gridX; i++)
	//	for(int j = 0; j < gridY; j++)
	//	{
	//		expr[i][j] = ex[i][j];
	//		eypr[i][j] = ey[i][j];
	//		hzpr[i][j] = hz[i][j];
	//	}

	//readFile(t+1);

	//double energy = 0;
	//double p =0;
	//double dx = LAMBDA/120;
	//double dt = dx / (2* CC);
	//for(int j = 15; j < 375; j++)
	//	for(int i = 0; i < gridX; i++)
	//	{
	//		double en = (ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j])/2.0 + hz[i][j]*hz[i][j]/2.0;
	//		double e1 = expr[i][j]*expr[i][j] + eypr[i][j]*eypr[i][j];
	//		double e2 = ex[i][j]*ex[i][j] + ey[i][j]*ey[i][j];
	//		//double en = hz[i][j]*hz[i][j]/2.0;
	//		//double e1 = hzpr[i][j];
	//		//double e2 = hz[i][j];
	//		energy+= en;
	//		p+= (e2-e1)/(5.0*dt);
	//	}
	////cout << "Energy: " << energy << endl;
	////cout << "P: " << p << endl;
	//	p = abs(p);
	//double Q = 2.0*PI*FREQ*energy/p;
	//cout << "Q-factor: " << Q << endl;
	//avrQ+=Q;
	//if (maxQ < Q) maxQ = Q;
	//}
	//avrQ = avrQ/1000.0;
	//cout << "MAX Q-factor: " << maxQ << endl;
	//cout << "Aver Q-factor: " << avrQ << endl;
	return 1;
}