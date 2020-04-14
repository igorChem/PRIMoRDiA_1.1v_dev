#ifndef INTERFACE 
#define INTERFACE


//--------------------------------------------
#include <iostream>
#include <memory>
#include <cstring>
#include <string>
#include <ctime>


#include "../include/log_class.h"
#include "../include/common.h"
#include "../include/Imolecule.h"
#include "../include/QMparser.h"
#include "../include/Icube.h"
#include "../include/Ibuffer.h"
#include "../include/Iline.h"
#include "../include/gridgen.h"
#include "../include/primordia.h"
#include "../include/autoprimordia.h"
#include "../include/test_p.h"

//------------------------------------------
using std::unique_ptr;
using std::cout;
using std::endl;
using std::stoi; 
using std::string;

//==================================================
/**
 * @class interface
 * @author igor
 * @date 14/03/20
 * @file interface.h
 * @brief 
 */
class interface{
	public:
	/**
	 *
	 */
		int m_argc;

		std::vector<std::string> m_argv;
	
		std::string runtyp;
	
		interface();
	
		interface(int argc, char* argv[] );
		
		interface(const interface& rhs) = delete;
		
		interface& operator=(const interface& rhs) = delete;
	
		~interface();
	
		void run();
	
		void write_help();
	
		void test_run();
		
		void write_input();
		
		void print_options();
	
};

#endif