#ifndef TEST_PRIMORDIA
#define TEST_PRIMORDIA

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "../include/common.h"


/**
 * @class test_p
 * @author barden
 * @date 17/04/18
 * @file test_p.h
 * @brief Class to instantiate a object with a collection of memeber functions for PRIMoRDiA-libs unit testing.
 */
class test_p {
	public:	
	
	std::vector<std::string> mopac_aux;
	std::vector<std::string> mopac_out;
	std::vector<std::string> mopac_mgf;
	std::vector<std::string> orca_out;
	std::vector<std::string> gauss_fchk;
	std::vector<std::string> gauss_log;
	std::vector<std::string> gamess_log;
	std::vector<std::string> molden;
	std::vector<std::string> pdb_files;
	
	
	//-----------------------------------------------------------------------
	/**
	 * @brief Default constructor for the class.
	 * Initializes the data for the classes tests.
	 */	
	test_p();
	
	//-----------------------------------------------------------------------
	/**
	 * @brief Destructor for the class. 
	 */
	~test_p();
	//-----------------------------------------------------------------------
	void test_primordia_1();
	void test_primordia_2();
	void test_primordia_3();
};

#endif