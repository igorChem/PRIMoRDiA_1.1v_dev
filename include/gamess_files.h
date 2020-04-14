//header file for GAMESS QM output file parser
// gamess_files.h

#ifndef  GAMESSFILES
#define GAMESSFILES

//Including c++ headers
#include <string>
#include <sstream>
#include <cstring>
#include <memory>

#include "../include/common.h"

class Iline; //fowards declarations 
class Ibuffer;
class Imolecule;

/**
 * @class gamess_files
 * @author igor
 * @date 10/03/20
 * @file gamess_files.h
 * @brief 
 */
class gamess_files{
	public:
		//-------------------------------------------------------
		/**
		* 
		*/		
		const char* name_f;
		
		//-------------------------------------------------------
		/**
		* 
		*/
		std::unique_ptr<Ibuffer> buffer;
		
		//-------------------------------------------------------
		/**
		* 
		*/
		std::unique_ptr<Imolecule> molecule;
		
		//-------------------------------------------------------
		/**
		* 
		*/
		gamess_files();
		
		//-------------------------------------------------------
		/**
		* 
		*/
		gamess_files(const char* file_name);
		
		//-------------------------------------------------------
		/**
		* 
		*/
		gamess_files( const gamess_files& rhs) = delete;
		
		//-------------------------------------------------------
		/**
		* 
		*/
		gamess_files& operator=( const gamess_files& rhs) = delete;
		
		//-------------------------------------------------------
		/**
		* 
		*/
		~gamess_files();
		//-------------------------------------------------------
		/**
		* 
		*/
		void parse_log();
		
		
};

#endif