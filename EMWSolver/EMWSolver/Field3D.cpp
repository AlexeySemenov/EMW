#include "Field3D.h"
#include "Helpers.h"
#include "Log.h"
#include <sstream>

namespace EMWSolver
{
	CField3D::CField3D(int sizeX, int sizeY, int sizeZ)
	{
		gridX = sizeX;
		gridY = sizeY;
		gridZ = sizeZ;

		Log::GetInstance().Write("Creating em-fields 3D arrays...");

		Ex = CreateJaggedArray3D(gridX, gridY, gridZ);
		Ey = CreateJaggedArray3D(gridX, gridY, gridZ);
		Ez = CreateJaggedArray3D(gridX, gridY, gridZ);

		Hx = CreateJaggedArray3D(gridX, gridY, gridZ);
		Hy = CreateJaggedArray3D(gridX, gridY, gridZ);
		Hz = CreateJaggedArray3D(gridX, gridY, gridZ);

		Log::GetInstance().WriteLine("Finished.");

		Log::GetInstance().Write("Zeroing em-fields 3D arrays...");

		ZeroJaggedArray3D(Ex, gridX, gridY, gridZ);
		ZeroJaggedArray3D(Ey, gridX, gridY, gridZ);
		ZeroJaggedArray3D(Ez, gridX, gridY, gridZ);

		ZeroJaggedArray3D(Hx, gridX, gridY, gridZ);
		ZeroJaggedArray3D(Hy, gridX, gridY, gridZ);
		ZeroJaggedArray3D(Hz, gridX, gridY, gridZ);

		Log::GetInstance().WriteLine("Finished.");
	}


	CField3D::~CField3D(void)
	{
		DeleteJaggedArray3D(this->Ex, gridX, gridY, gridZ);
		DeleteJaggedArray3D(this->Ey, gridX, gridY, gridZ);
		DeleteJaggedArray3D(this->Ez, gridX, gridY, gridZ);

		DeleteJaggedArray3D(this->Hx, gridX, gridY, gridZ);
		DeleteJaggedArray3D(this->Hy, gridX, gridY, gridZ);
		DeleteJaggedArray3D(this->Hz, gridX, gridY, gridZ);
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
		for(int i = crop.left; i < gridX - crop.right; i++)
			for(int j = crop.down; j < gridY - crop.up; j++)
				for(int k = crop.bottom; k < gridZ - crop.top; k++)
			{
				tmp = outputArray[i][j][k];
				file.write((char *) &tmp, sizeof(double));
			}
		Log::GetInstance().Write(gridX * gridY * gridZ * sizeof(double));
		Log::GetInstance().WriteLine(" BYTE written succesfully");
		file.close();
	}
}
