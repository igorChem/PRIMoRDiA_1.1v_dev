//Including c++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <algorithm>  

//Including PRIMoRDiA headers
//-------------------------------------------------------
#include "../include/log_class.h"
#include "../include/common.h"
#include "../include/Imolecule.h"
#include "../include/QMparser.h"
#include "../include/gamess_files.h"
#include "../include/orca_files.h"
#include "../include/gaussian_files.h"
#include "../include/mopac_files.h"
#include "../include/terachem_files.h"
//-------------------------------------------------------
// Aliases for standard c++ scope functions
using std::cout;
using std::endl;
using std::string;
/************************************************************************************/
QMparser::QMparser()	:
	program("none")		,
	name_f("none")			,
	parsed(false)				{
}
/************************************************************************************/
QMparser::QMparser(const char* file_name,
								string Program			):
	program(Program)									,
	name_f(file_name)									{
}
/************************************************************************************/
Imolecule QMparser::get_molecule(unsigned int mob){
	if (program == "mopac"){
		mopac_files file_obj(name_f,mob);
		if ( file_obj.type == "AUX" ){
			file_obj.parse_aux();
		}else if ( file_obj.type == "OUT"){
			file_obj.parse_out();
			file_obj.parse_mgf();
		}
		else if (file_obj.type == "MGF"){
			file_obj.parse_mgf();
		}
		return *file_obj.molecule;
	}else if (program == "gamess"){
		gamess_files file_obj(name_f);
		file_obj.parse_log();
		return *file_obj.molecule;
	}else if (program == "orca"){
		orca_files file_obj(name_f);
		file_obj.parse_out();
		return *file_obj.molecule;
	}else if (program == "gaussian"){
		gaussian_files file_obj(name_f);
		file_obj.parse_fchk();
		file_obj.get_overlap_m();
		return *file_obj.molecule;
	}else if (program == "terachem"){
		terachem_files file_obj(name_f);
		file_obj.parse_molden();
		return *file_obj.molecule;
	}else{
		cout << "Program keyword not recognized\n Exiting Program! " << endl;
		exit(-1);
	}
}
/**********************************************************************************/
QMparser::~QMparser(){}
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////END OF CLASS//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
