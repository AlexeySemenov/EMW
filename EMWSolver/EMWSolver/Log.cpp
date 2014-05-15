#include "Log.h"
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>

namespace EMWSolver
{
	Log::Log(void)
	{
		 std::time_t rawtime;
		 struct tm * timeinfo;
		 char buffer[80];

		 time (&rawtime);
		 timeinfo = localtime(&rawtime);

		 strftime(buffer,80,"%d-%m-%Y %I_%M_%S",timeinfo);
		 std::string str(buffer);
	 
		 std::stringstream fileNameStream;
		 fileNameStream << "emwsolver" << str << ".log";
		 std::string fileName =  fileNameStream.str();
		 logFile = new std::fstream();
		 logFile->open( fileName, std::ios::out | std::ios::trunc);
		 logStream = logFile;
		 if(!logFile->is_open())
		 {
			 logStream = &std::cout;
			 std::cerr << "ERROR::Can't create log file. Using: standard out instead." << std::endl;
		 }

	}

	void Log::SetOutputStream(std::ostream* stream)
	{
		if(stream != NULL)
		{
			if(stream != static_cast<std::ostream*>(logFile))
			{
				logFile->flush();
				logFile->close();
				logStream  = stream;
				delete logFile;
				logFile = 0;
			}
		}
	}

	void Log::WriteLine(const std::string& str)
	{
		if(logStream)
		{
			*logStream << str << std::endl;

		}
	}

	void Log::Write(const std::string& str)
	{
		if(logStream)
		{
			*logStream << str;

		}
	}

	void Log::Write(const int val)
	{
		if(logStream)
		{
			*logStream << val;

		}
	}

	void Log::Close()
	{
		if(logFile)
			if(logFile->is_open())
				logFile->close();
	}
}