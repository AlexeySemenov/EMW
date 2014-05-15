namespace EMWSolver
{
	double* const	CreateJaggedArray1D(int size)
	{
		return new double[size];
	}

	double** const	CreateJaggedArray2D(int sizeX, int sizeY)
	{
		double** const tempArray = new double*[sizeX];
		for(int i = 0; i < sizeX; i++)
		{
			tempArray[i] = new double[sizeY];
		}
		return tempArray;
	}

	double***const	CreateJaggedArray3D(int sizeX, int sizeY, int sizeZ)
	{
		double*** const tempArray = new double**[sizeX];
		for(int i = 0; i < sizeX; i++)
		{
			tempArray[i] = new double*[sizeY];
		}

		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++)
			{
				tempArray[i][j] = new double[sizeZ];
			}

		return tempArray;
	}

	double* const	CreateElongatedArray1D(int size)
	{
		return CreateJaggedArray1D(size);
	}

	double* const	CreateElongatedArray2D(int sizeX, int sizeY)
	{
		return CreateJaggedArray1D(sizeX * sizeY);
	}

	double* const	CreateElongatedArray3D(int sizeX, int sizeY, int sizeZ)
	{
		return CreateJaggedArray1D(sizeX * sizeY * sizeZ);
	}

	const  void	DeleteJaggedArray1D(double* const array1d, int size)
	{
		delete[] array1d;
	}

	const  void	 DeleteJaggedArray2D(double** const array2d, int sizeX, int sizeY)
	{	
		for(int i = 0; i < sizeX; i++)
		{
			delete[] array2d[i];
		}

		delete[] array2d;
	}

	const  void	DeleteJaggedArray3D(double*** const array3d, int sizeX, int sizeY, int sizeZ)
	{
		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++)
			{
				delete[] array3d[i][j];
			}
	
		for(int i = 0; i < sizeX; i++)
		{
			delete[] array3d[i];
		}

		delete[] array3d;;
	}

	const  void	DeleteElongatedArray1D(double* const array1d, int size)
	{
		DeleteJaggedArray1D(array1d, size);
	}

	const  void	DeleteElongatedArray2D(double* const array2d, int sizeX, int sizeY)
	{
		DeleteJaggedArray1D(array2d, sizeX * sizeY);
	}

	const  void	DeleteElongatedArray3D(double* const array3d, int sizeX, int sizeY, int sizeZ)
	{
		DeleteJaggedArray1D(array3d, sizeX * sizeY * sizeZ);
	}

	const  void	ZeroJaggedArray1D(double* const array1d, int size)
	{
		for(int i = 0; i < size; i++)
			array1d[i] = 0.0;
	}

	const  void	 ZeroJaggedArray2D(double** const array2d, int sizeX, int sizeY)
	{	
		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++)
				array2d[i][j] = 0.0;
	}

	const  void	ZeroJaggedArray3D(double*** const array3d, int sizeX, int sizeY, int sizeZ)
	{
		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++)
				for(int k = 0; k < sizeZ; k++)
					array3d[i][j][k] = 0.0;
	}

	const  void	ZeroElongatedArray1D(double* const array1d, int size)
	{
		ZeroJaggedArray1D(array1d, size);
	}

	const  void	ZeroElongatedArray2D(double* const array2d, int sizeX, int sizeY)
	{
		ZeroJaggedArray1D(array2d, sizeX * sizeY);
	}

	const  void	ZeroElongatedArray3D(double* const array3d, int sizeX, int sizeY, int sizeZ)
	{
		ZeroJaggedArray1D(array3d, sizeX * sizeY * sizeZ);
	}

	const  void	ValueJaggedArray1D(double* const array1d, int size, double value)
	{
		for(int i = 0; i < size; i++)
			array1d[i] = value;
	}

	const  void	 ValueJaggedArray2D(double** const array2d, int sizeX, int sizeY, double value)
	{	
		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++)
				array2d[i][j] = value;
	}

	const  void	ValueJaggedArray3D(double*** const array3d, int sizeX, int sizeY, int sizeZ, double value)
	{
		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++)
				for(int k = 0; k < sizeZ; k++)
					array3d[i][j][k] = value;
	}
}