//header file for GAUSSIAN QM output file parser
// gamess_files.h

#ifndef  GAUSSFILES
#define GAUSSFILES

//Including c++ headers
#include <string>
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
class gaussian_files{
	public:
		
		//--------------------------------
		/**
		 *
		 */
		const char* name_f;
		
		//--------------------------------
		/**
		 *
		 */
		std::unique_ptr<Ibuffer> buffer;
		
		//---------------------------------
		/**
		 *
		 */
		std::unique_ptr<Imolecule> molecule;
		
		//---------------------------------
		/**
		 *
		 */
		gaussian_files();
		
		//---------------------------------
		/**
		 *
		 */
		~gaussian_files();
		
		//---------------------------------
		/**
		 *
		 */
		gaussian_files(const char* file_name);
		
		//---------------------------------
		/**
		 *
		 */
		gaussian_files(const gaussian_files& rhs) = delete;
		
		//---------------------------------
		/**
		 *
		 */
		gaussian_files& operator=(const gaussian_files& rhs) = delete;
		
		//---------------------------------
		/**
		 *
		 */
		void parse_fchk();
		
		//---------------------------------
		/**
		 *
		 */
		void get_overlap_m();
};

#endif