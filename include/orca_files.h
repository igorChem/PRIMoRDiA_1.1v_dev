//header file for ORCA QM output file parser
// gamess_files.h

#ifndef  ORCAFILES
#define ORCAFILES

//Including c++ headers
#include <string>
#include <sstream>
#include <cstring>
#include <memory>

#include "../include/common.h"

class Ibuffer;
class Imolecule;

//====================================================
/**
 * @class orca_files
 * @author igor
 * @date 10/03/20
 * @file orca_files.h
 * @brief 
 */
class orca_files{
	public:
		//-------------------------------------
		/**
		 *
		 */
		const char* name_f;
		std::unique_ptr<Ibuffer> buffer;
		std::unique_ptr<Imolecule> molecule;
		orca_files();
		~orca_files();
		orca_files(const char* file_name);
		orca_files(const orca_files& rhs) = delete;
		orca_files& operator=(const orca_files& rhs) = delete;
		void parse_out();
		void get_overlap(int ov_in, int ov_fin);
		void parse_molden();
		
};
//===================================================

/**
 * @class basis_orca
 * @author igor
 * @date 10/03/20
 * @file orca_files.h
 * @brief 
 */
 class basis_orca{
	public:
		std::string element_type;
		std::vector<double> shell_size;
		std::vector<std::string> shell_sym;	
		std::vector<double> coefficients;
		std::vector<double> exp;
		basis_orca(){};
 };

#endif