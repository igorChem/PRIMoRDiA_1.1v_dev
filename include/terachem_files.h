#ifndef TERACHEM_FILES
#define TERACHEM_FILES

#include <iostream>
#include <string>


#include "../include/common.h"


class Ibuffer;
class Imolecule;

class terachem_files{
	public:
		const char* name_f;
		std::unique_ptr<Ibuffer> buffer;
		std::unique_ptr<Imolecule> molecule;
		terachem_files();
		~terachem_files();
		terachem_files(const char* file_name);
		terachem_files(const terachem_files& rhs) = delete;
		terachem_files& operator=(const terachem_files& rhs) = delete;
		void parse_molden();
};

#endif