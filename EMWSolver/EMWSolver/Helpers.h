#ifndef HELPERS_H
#define HELPERS_H
namespace EMWSolver
{
	double*		const	CreateJaggedArray1D(int size);
	double**	const	CreateJaggedArray2D(int sizeX, int sizeY);
	double***	const	CreateJaggedArray3D(int sizeX, int sizeY, int sizeZ);

	double*		const	CreateElongatedArray1D(int size);
	double*		const	CreateElongatedArray2D(int sizeX, int sizeY);
	double*		const	CreateElongatedArray3D(int sizeX, int sizeY, int sizeZ);

	void	DeleteJaggedArray1D(double*  const array1d, int size);
	void	DeleteJaggedArray2D(double** const array2d, int sizeX, int sizeY);
	void	DeleteJaggedArray3D(double***const array3d, int sizeX, int sizeY, int sizeZ);

	void	DeleteElongatedArray1D(double* const array1d, int size);
	void	DeleteElongatedArray2D(double* const array2d, int sizeX, int sizeY);
	void	DeleteElongatedArray3D(double* const array3d, int sizeX, int sizeY, int sizeZ);

	void	ZeroJaggedArray1D(double*  const array1d, int size);
	void	ZeroJaggedArray2D(double** const array2d, int sizeX, int sizeY);
	void	ZeroJaggedArray3D(double***const array3d, int sizeX, int sizeY, int sizeZ);

	void	ZeroElongatedArray1D(double* const array1d, int size);
	void	ZeroElongatedArray2D(double* const array2d, int sizeX, int sizeY);
	void	ZeroElongatedArray3D(double* const array3d, int sizeX, int sizeY, int sizeZ);

	void	ValueJaggedArray1D(double*  const array1d, int size, double value);
	void	ValueJaggedArray2D(double** const array2d, int sizeX, int sizeY, double value);
	void	ValueJaggedArray3D(double***const array3d, int sizeX, int sizeY, int sizeZ, double value);
}
#endif