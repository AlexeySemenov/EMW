#ifndef LOG_H
#define LOG_H
#include <ostream>
#include <fstream>
#include <string>

namespace EMWSolver
{
	class Log
	{
	public:
		static Log& GetInstance()
			{
				static Log    instance;                              
				return instance;
			}
		void WriteLine(const std::string& str);
		void Write(const std::string& str);
		void Write(const int val);
		void SetOutputStream(std::ostream* stream);
		void Close();
	private:
		Log(void);
		Log(Log const&);             
        void operator=(Log const&);
		std::ostream* logStream;
		std::fstream* logFile;
	};
}
#endif
