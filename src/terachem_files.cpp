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
#include "../include/Iaorbital.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/Iline.h"
#include "../include/Ibuffer.h"
#include "../include/terachem_files.h"
#include "../include/overlap_int.h"

using std::vector;
using std::string;
using std::stoi;
using std::stod;
using std::endl;
using std::cout;

//-------------------------------------------------------
/********************************************************/
terachem_files::terachem_files(){}
/********************************************************/
terachem_files::terachem_files(const char* file_name)	:
	molecule( new Imolecule() )							,
	buffer( new Ibuffer() )								{
	if ( !check_file_ext(".molden",file_name) ) {
		cout << "Warning! The file has wrong etension name!" << endl;
		m_log->input_message("Warning! The file has wrong etension name!");
		exit(-1);
	}
	name_f = file_name; 
}
/********************************************************/
void terachem_files::parse_molden(){
	unsigned int i,j,k;
	int at_in		= 0;
	int at_fin		= 0;
	int gto_in		= 0;
	int gto_fin		= 0;
	int mo_in		= 0;
	int n_gto		= 0;
	int n_at_gto	= 0;
	
	
	
	buffer.reset( new Ibuffer (name_f,true) );
	for (i=0;i<buffer->nLines;i++){
		if ( buffer->lines[i].IF_word("[Atoms]",0,7) ){
			at_in = i;
		}
		else if ( buffer->lines[i].IF_word("[GTO]",0,5) ){
			at_fin	= i;
			gto_in	= i;
		}
		else if ( buffer->lines[i].IF_word("[MO]",0,4) ){
			gto_fin	=i;
			mo_in	=i;
			break;
		}
	}
	for(i=0;i<buffer->nLines-1;i++){
		if ( i > at_in && i < at_fin ){
			double xx, yy, zz;
			xx = buffer->lines[i].get_double(3);
			yy = buffer->lines[i].get_double(4);
			zz = buffer->lines[i].get_double(5);
			string type_= buffer->lines[i].get_string(0);
			molecule->add_atom(xx,yy,zz,type_);
		}
		else if ( i > gto_in && i < gto_fin ){
			if ( buffer->lines[i].words[0] == "s" ){
				Iaorbital orb;
				orb.gto			= true;
				orb.symmetry	= "S";
				n_gto			= buffer->lines[i].get_int(1);
				if ( buffer->lines[i-1].words[0].size() == 1 ){
					n_at_gto	= buffer->lines[i-1].get_int(0);
				}
				for(j=0;j<n_gto;j++){
					orb.add_primitive(buffer->lines[i+j+1].get_double(0),buffer->lines[i+j+1].get_double(1));
				}
				molecule->atoms[n_at_gto-1].add_orbital(orb);
				
			}
			else if ( buffer->lines[i].words[0] == "p"){
				Iaorbital orbX,orbY,orbZ;
				orbX.gto			= true;
				orbX.symmetry	= "PX";
				n_gto				= buffer->lines[i].get_int(1);
				for(j=0;j<n_gto;j++){
					orbX.add_primitive(buffer->lines[j+i+1].get_double(0),buffer->lines[j+i+1].get_double(1));
				}
				orbY = orbZ = orbX;
				orbX.powx		= 1;
				orbY.powy		= 1;
				orbZ.powz		= 1;
				orbY.symmetry	= "PY";
				orbZ.symmetry	= "PZ";
				molecule->atoms[n_at_gto-1].add_orbital(orbX);
				molecule->atoms[n_at_gto-1].add_orbital(orbY);
				molecule->atoms[n_at_gto-1].add_orbital(orbZ);
			}
		}
		else if( i  > mo_in ){
			if ( buffer->lines[i].words[0] == "Ene=" ){
				molecule->orb_energies.push_back( buffer->lines[i].get_double(1) );
				molecule->MOnmb++;
			}
			else if (buffer->lines[i].words[0] == "Occup=" ) {
				molecule->occupied.push_back( buffer->lines[i].get_int(1) );
				molecule->num_of_electrons += buffer->lines[i].get_int(1);
			}
			else if ( buffer->lines[i].words[0] == "Spin=" ) {continue;}
			else	molecule->coeff_MO.push_back( buffer->lines[i].get_double(1) ); 
		}
		cout << n_at_gto << endl;
	}
	molecule->ang_to_bohr();
	cout << molecule->num_of_atoms << endl;
	overlap_int overlap_integrals(molecule->atoms);
	overlap_integrals.calculate_overlap(molecule->m_overlap);
}

/********************************************************/
terachem_files::~terachem_files(){}
//////////////////////////////////////////////////////////
