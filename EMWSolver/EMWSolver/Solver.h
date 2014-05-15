#pragma once
#include "Parameters.h"

namespace EMWSolver
{
	class CSolver
	{
	public:
		virtual void Create(TaskParameters taskParam, NumericalParameters numParam, SourceParameters srcParam, 
				BoundaryParameters* bndParam) = 0;
		virtual void Solve(int timeStep) = 0;
	private:
		virtual void prepare() = 0;
	};
}
