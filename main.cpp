#include "Energy_Post_Processing.h"
#include <ctime>

#define Energy_Error std::cerr << std::endl << __PRETTY_FUNCTION__ << " - ERROR: \n"
#define Energy_Warning std::cout << std::endl <<  __PRETTY_FUNCTION__ << " - Warning: \n"
#define Energy_Out std::cout << std::endl << __PRETTY_FUNCTION__ << " : \n"

int main(int argc, char *argv[])
{
	std::clock_t start;
    double duration;

    start = std::clock();
 
	if (argc != 2) 
	{
        Energy_Error << "Usage: " << argv[0] << " ESSI OUTPUT FILE \n";
        return 1;
    }

	Energy_Post_Processing* test1 = new Energy_Post_Processing();
	//test1.HDF5filename = "1D_Damping_numerical_dynamic.h5.feioutput";
	test1->HDF5filename = argv[1];

	//cout << test1.is_initialized << endl;

	test1->initialize();
	//cout << test1.is_initialized << endl;
	test1->computeAndWriteEnergyDensity();

	delete test1;

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	std::cout << "\n=========================================\n";
	std::cout << "Program end normally! \n";
	std::cout << "Time used: " << duration << "\n";
	std::cout << "Time to go fxxk yourself :)";
	std::cout << "\n=========================================\n";



	return 0;
}