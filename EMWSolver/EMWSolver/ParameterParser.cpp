#include "ParameterParser.h"
#include <assert.h>

namespace EMWSolver
{
	CParameterParser::CParameterParser(std::string fileName)
	{
		numberTasks = 0;
		currentTask = 0;
		defaultTaskParameters.Name = "Untitled";
		defaultTaskParameters.ompThreads = 4;

		defaultSourceParameters.isTFSF = false;
		defaultSourceParameters.pointSourceX = NULL;
		defaultSourceParameters.pointSourceY = NULL;
		defaultSourceParameters.pointSourceZ = NULL;

		defaultNumericalParameters.S = 1;
		defaultNumericalParameters.sizeX = 100;
		defaultNumericalParameters.sizeY = 50;
		defaultNumericalParameters.sizeZ = 1500;

		defaultNumericalParameters.Nx = 80;
		defaultNumericalParameters.Ny = 80;
		defaultNumericalParameters.Nz = 40;

		defaultBoundaryParameters[0].boundaryPosition = BoundaryPosition::Left;
		defaultBoundaryParameters[0].boundaryType = BoundaryType::Dirichlet;

		defaultBoundaryParameters[1].boundaryPosition = BoundaryPosition::Right;
		defaultBoundaryParameters[1].boundaryType = BoundaryType::Dirichlet;

		defaultBoundaryParameters[2].boundaryPosition = BoundaryPosition::Up;
		defaultBoundaryParameters[2].boundaryType = BoundaryType::Dirichlet;

		defaultBoundaryParameters[3].boundaryPosition = BoundaryPosition::Down;
		defaultBoundaryParameters[3].boundaryType = BoundaryType::Dirichlet;

		defaultBoundaryParameters[4].boundaryPosition = BoundaryPosition::Top;
		defaultBoundaryParameters[4].boundaryType = BoundaryType::Dirichlet;

		defaultBoundaryParameters[5].boundaryPosition = BoundaryPosition::Bottom;
		defaultBoundaryParameters[5].boundaryType = BoundaryType::Dirichlet;
	}

	void CParameterParser::Parse(std::string fileName)
	{

	}

	//return 0 if no more tasks in the queue
	int CParameterParser::NextTask()
	{
		currentTask++;
		assert(currentTask < numberTasks);
		return 1;
	}

	TaskParameters CParameterParser::GetTaskParameters()
	{
		assert(currentTask <= numberTasks);

		return taskParameters[getTaskNumber()];
	}

	SourceParameters CParameterParser::GetSourceParameters()
	{
		assert(currentTask <= numberTasks);

		return sourceParameters[getTaskNumber()];
	}

	BoundaryParameters* CParameterParser::GetBoundaryParameters()
	{
		assert(currentTask <= numberTasks);

		return boundaryParameters[getTaskNumber()];
	}

	NumericalParameters CParameterParser::GetNumericalParameters()
	{
		assert(currentTask <= numberTasks);

		return numericalParameters[getTaskNumber()];
	}

	TaskParameters CParameterParser::GetDefaultTaskParameters()
	{
		return defaultTaskParameters;
	}
	SourceParameters CParameterParser::GetDefaultSourceParameters()
	{
		return defaultSourceParameters;
	}

	BoundaryParameters* CParameterParser::GetDefaultBoundaryParameters()
	{
		return defaultBoundaryParameters;
	}

	NumericalParameters CParameterParser::GetDefaultNumericalParameters()
	{
		return defaultNumericalParameters;
	}

	int inline CParameterParser::getTaskNumber()
	{
		return currentTask;
	}

	CParameterParser::~CParameterParser(void)
	{
	}
}