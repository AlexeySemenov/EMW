#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <string>
namespace EMWSolver
{
	enum BoundaryType
	{
		Dirichlet = 0,
		PEC,
		PMC,
		UPML,
		Custom
	};

	enum BoundaryPosition
	{
		Left = 0,
		Right,
		Up,
		Down,
		Top,
		Bottom,
		Inside
	};

	struct NumericalParameters
	{
		//Task size
		int sizeX;
		int sizeY;
		int sizeZ;

		//Courant stability factor: default 1
		double S;

		//Number of cells per wavelength
		int Nx;
		int Ny;
		int Nz;
	};

	struct SourceParameters
	{
		//Is total field - scattered field enabled
		bool isTFSF;

		//TF/SF region offset from boundaries
		int leftTFCellOffset;
		int rightTFCellOffset;
		int upTFCellOffset;
		int downTFCellOffset;
		int topTFCellOffset;
		int bottomCellOffset;

		//Arrays of point source X,Y,Z cell position
		int* pointSourceX;
		int* pointSourceY;
		int* pointSourceZ;

	};

	struct TaskParameters
	{
		//Task name
		std::string Name;
		
		//Number of OpenMP threads
		int ompThreads;

	};

	struct BoundaryParameters
	{
		//Boundary condition type
		BoundaryType boundaryType;

		//Boundary position
		BoundaryPosition boundaryPosition;
		
		//UPML cell offset
		int leftUPMLCells;
		int rightUPMLCells;
		int upUPMLCells;
		int downUPMLCells;
		int topUPMLCells;
		int bottomUPMLCells;

	};
}

#endif