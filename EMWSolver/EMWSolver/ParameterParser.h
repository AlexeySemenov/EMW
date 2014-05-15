#pragma once
#include <string>
#include <vector>
#include "Parameters.h"

namespace EMWSolver
{

	class CParameterParser
	{
	public:
		CParameterParser(std::string fileName);
		void Parse(std:: string fileName);

		TaskParameters GetTaskParameters();
		SourceParameters GetSourceParameters();
		BoundaryParameters* GetBoundaryParameters();
		NumericalParameters GetNumericalParameters();

		TaskParameters GetDefaultTaskParameters();
		SourceParameters GetDefaultSourceParameters();
		BoundaryParameters* GetDefaultBoundaryParameters();
		NumericalParameters GetDefaultNumericalParameters();



		int NextTask();

		~CParameterParser(void);
	private:
		int numberTasks;
		int currentTask;
		TaskParameters defaultTaskParameters;
		SourceParameters defaultSourceParameters;
		BoundaryParameters defaultBoundaryParameters[6];
		NumericalParameters defaultNumericalParameters;

		std::vector<TaskParameters> taskParameters;
		std::vector<SourceParameters> sourceParameters;
		std::vector<BoundaryParameters*> boundaryParameters;
		std::vector<NumericalParameters> numericalParameters;

		int inline getTaskNumber();
	};
}
