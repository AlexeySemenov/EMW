#include "Field3D.h"
#include <iostream>

int FieldTest();

int main(int argc, char* argv[])
{
	std::cerr << "Memory allocation test...";
	if(FieldTest()) std::cerr << "FAILED."; 
	else std::cerr << "PASSED.";
	std::cerr << std::endl; 

	return 0;
}

int FieldTest()
{
	try
	{
		EMWSolver::CField3D* field = new EMWSolver::CField3D(100, 100, 100);
		delete field;
	}
	catch(...)
	{
		return 1;
	}
	return 0;
}