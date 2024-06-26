//header file for QM output file parser
// currently done for mopac and GAMESS
// mopac_files.h

#ifndef  MOPACFILES
#define MOPACFILES

//Including c++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <memory>

#include "../include/common.h"

#endif

class mopac_files{
	public:
		//-----------------------------------------------------------------------------
		/**
		* @brief Conditional if the file is from MOZYME calculations
		*/
		bool LMO;
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		bool RHF;
		
		double f_chg;
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		unsigned int band;
 		//-----------------------------------------------------------------------------
		/**
		*
		*/
		const char* name_f;
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		std::string type;
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		std::unique_ptr<Ibuffer> buffer;
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		std::unique_ptr<Imolecule> molecule;
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		mopac_files();
		
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		mopac_files(const char* file_name,unsigned int MOband);
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		mopac_files(const mopac_files& rhs_mop) = delete;
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		mopac_files& operator=(const mopac_files& rhs_mop) = delete;
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		~mopac_files();
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		void parse_aux();
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		void parse_out();
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		void parse_mgf();
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		void get_overlap_m();
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		void get_mo(bool beta);
		
		//-----------------------------------------------------------------------------
		/**
		*
		*/
		void get_mo_energies(bool beta);
		
		
};
