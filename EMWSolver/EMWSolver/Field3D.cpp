#include "Field3D.h"
#include "Helpers.h"
#include "Log.h"
#include <sstream>

namespace EMWSolver
{
	CField3D::CField3D(int _sizeX, int _sizeY, int _sizeZ)
	{
		sizeX = _sizeX;
		sizeY = _sizeY;
		sizeZ = _sizeZ;

		gridX = sizeX + 1;
		gridY = sizeY + 1;
		gridZ = sizeZ + 1;

		Log::GetInstance().Write("Creating em-fields 3D arrays...");

		Ex = CreateJaggedArray3D(gridX, gridY, gridZ);
		Ey = CreateJaggedArray3D(gridX, gridY, gridZ);
		Ez = CreateJaggedArray3D(gridX, gridY, gridZ);

		//H field arrays size is equals to number of cells in each dimension
		Hx = CreateJaggedArray3D(sizeX, sizeY, sizeZ);
		Hy = CreateJaggedArray3D(sizeX, sizeY, sizeZ);
		Hz = CreateJaggedArray3D(sizeX, sizeY, sizeZ);

		Log::GetInstance().WriteLine("Finished.");

		Log::GetInstance().Write("Zeroing em-fields 3D arrays...");

		ZeroJaggedArray3D(Ex, gridX, gridY, gridZ);
		ZeroJaggedArray3D(Ey, gridX, gridY, gridZ);
		ZeroJaggedArray3D(Ez, gridX, gridY, gridZ);

		ZeroJaggedArray3D(Hx, sizeX, sizeY, sizeZ);
		ZeroJaggedArray3D(Hy, sizeX, sizeY, sizeZ);
		ZeroJaggedArray3D(Hz, sizeX, sizeY, sizeZ);

		Log::GetInstance().WriteLine("Finished.");
	}


	CField3D::~CField3D(void)
	{
		DeleteJaggedArray3D(this->Ex, gridX, gridY, gridZ);
		DeleteJaggedArray3D(this->Ey, gridX, gridY, gridZ);
		DeleteJaggedArray3D(this->Ez, gridX, gridY, gridZ);

		DeleteJaggedArray3D(this->Hx, sizeX, sizeY, sizeZ);
		DeleteJaggedArray3D(this->Hy, sizeX, sizeY, sizeZ);
		DeleteJaggedArray3D(this->Hz, sizeX, sizeY, sizeZ);
	}

	void CField3D::WriteFieldToBinary(EMField component, EMCrop crop, int timestep, const std::string& path)
	{
		double*** outputArray;
		std::string outputName;
		switch (component)
		{
			case EMField::Ex: outputArray = Ex; outputName = "Ex"; break;
			case EMField::Ey: outputArray = Ey; outputName = "Ey"; break;
			case EMField::Ez: outputArray = Ez; outputName = "Ez"; break;

			case EMField::Hx: outputArray = Hx; outputName = "Hx"; break;
			case EMField::Hy: outputArray = Hy; outputName = "Hy"; break;
			case EMField::Hz: outputArray = Hz; outputName = "Hz"; break;
		}

		std::string fileName;
		std::stringstream s;
		s << path << "\\" << outputName << "[" << timestep << "].dat";
		fileName = s.str();
		std::ofstream file(fileName.c_str(), std::ios::out | std::ios::binary);
		double tmp;

		if(!file.is_open())
			Log::GetInstance().WriteLine("ERROR::Cannot open file for write!");

		Log::GetInstance().Write("Writing to file... ");
		for(int i = crop.left; i < sizeX - crop.right; i++)
			for(int j = crop.down; j < sizeY - crop.up; j++)
				for(int k = crop.bottom; k < sizeZ - crop.top; k++)
			{
				tmp = outputArray[i][j][k];
				file.write((char *) &tmp, sizeof(double));
			}
		Log::GetInstance().Write(sizeX * sizeY * sizeZ * sizeof(double));
		Log::GetInstance().WriteLine(" BYTE written succesfully");
		file.close();
	}
}
