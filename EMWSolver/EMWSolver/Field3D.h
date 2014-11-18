#ifndef FIELD3D_H
#define FIELD3D_H

#include <string>

namespace EMWSolver
{
	enum EMField
	{
		Ex = 0,
		Ey,
		Ez,

		Hx,
		Hy,
		Hz
	};

	struct EMCrop
	{
		int left;
		int right;
		int up;
		int down;
		int top;
		int bottom;
	};

	class CField3D
	{
	public:
		double*** Ex;
		double*** Ey;
		double*** Ez;

		double*** Hx;
		double*** Hy;
		double*** Hz;

		int gridX;
		int gridY;
		int gridZ;

		int sizeX;
		int sizeY;
		int sizeZ;

		//Number of cells per each dimension
		CField3D(int _sizeX, int _sizeY, int _sizeZ);
		void WriteFieldToBinary(EMField component, EMCrop crop, int timestep, const std::string& path);
		~CField3D(void);

	};
}

#endif