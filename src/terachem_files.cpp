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
	
	int nLines		= 2;
	
	
	molecule->name = remove_extension(name_f);
	m_log->input_message("Starting to parse molden file!");
	buffer.reset( new Ibuffer (name_f,"[Atoms]","[GTO]") );
	nLines++;
	for (i=1;i<buffer->nLines;i++){
		if ( buffer->lines[i].words.size() > 1 ){
			double xx, yy, zz;
			xx = buffer->lines[i].get_double(3);
			yy = buffer->lines[i].get_double(4);
			zz = buffer->lines[i].get_double(5);
			string type_= buffer->lines[i].get_string(0);
			molecule->add_atom(xx,yy,zz,type_);
		}
		nLines++;
	}
	
	buffer.reset( new Ibuffer (name_f,"[GTO","[MO]") );
	nLines++;
	for ( i=1;i<buffer->nLines;i++) {
		if (  buffer->lines[i].words.size() > 1){
			if ( buffer->lines[i].words[0] == "s" ){
				Iaorbital orb;
				orb.gto			= true;
				orb.symmetry	= "S";
				n_gto			= buffer->lines[i].get_int(1);
				if ( buffer->lines[i-1].words.size() == 2 && buffer->lines[i-1].words[0].size() < 8 ){
					n_at_gto	= buffer->lines[i-1].get_int(0);
				}
				for(j=0;j<n_gto;j++){
					orb.add_primitive(buffer->lines[i+j+1].get_double(0),buffer->lines[i+j+1].get_double(1));
				}
				molecule->atoms[n_at_gto-1].add_orbital(orb);
			}
			else if ( buffer->lines[i].words[0] == "p"){
				Iaorbital orbX,orbY,orbZ;
				orbX.gto		= true;
				orbX.symmetry	= "PX";
				n_gto			= buffer->lines[i].get_int(1);
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
			else if ( buffer->lines[i].words[0] == "d"){
				Iaorbital orbXX,orbYY,orbZZ,orbXY,orbXZ,orbYZ;
				orbXX.gto		= true;
				orbXX.symmetry	= "X2";
				n_gto			= buffer->lines[i].get_int(1);
				for(j=0;j<n_gto;j++){
					orbXX.add_primitive(buffer->lines[j+i+1].get_double(0),buffer->lines[j+i+1].get_double(1));
				}
				orbYY = orbZZ = orbXY = orbXZ = orbYZ = orbXX;
				orbXX.powx		= 2;
				orbYY.powy		= 2;
				orbZZ.powz		= 2;
				orbXY.powx		= 1;
				orbXY.powy		= 1;
				orbXZ.powx		= 1;
				orbXZ.powz		= 1;
				orbYZ.powy		= 1;
				orbYZ.powz		= 1;
				orbYY.symmetry	= "YY";
				orbZZ.symmetry	= "ZZ";
				orbXZ.symmetry	= "XZ";
				orbXY.symmetry	= "XY";
				orbYZ.symmetry	= "YZ";
				
				molecule->atoms[n_at_gto-1].add_orbital(orbXX);
				molecule->atoms[n_at_gto-1].add_orbital(orbYY);
				molecule->atoms[n_at_gto-1].add_orbital(orbZZ);
				molecule->atoms[n_at_gto-1].add_orbital(orbXZ);
				molecule->atoms[n_at_gto-1].add_orbital(orbXY);
				molecule->atoms[n_at_gto-1].add_orbital(orbYZ);
			}
			
		}
		nLines++;
	}
	
	mo_in = nLines;
	char tmp_line[50];
	buffer.reset(nullptr);
	int in_indx = 0;
	std::ifstream buf(name_f);
	while( !buf.eof() ){
		buf.getline(tmp_line,50);
		if ( in_indx >= mo_in  ){
			Iline Line(tmp_line);
			if ( Line.IF_line("Ene=",0,2) ) {
				molecule->orb_energies.push_back( Line.get_double(1) );
				molecule->MOnmb++;
			}else if ( Line.IF_line("Occup=",0,2) ) {
				molecule->occupied.push_back( Line.get_int(1) );
				molecule->num_of_electrons += Line.get_int(1);
			}
			else if ( Line.IF_line("Spin=",0,2 ) ) {continue;}
			else	molecule->coeff_MO.push_back( Line.get_double(1) ); 
		}
		in_indx++;
	}
	buf.close();

	//molecule->print_basis();
	m_log->input_message("# of atoms: ");
	m_log->input_message( int(molecule->num_of_atoms));
	m_log->input_message("# of atomic orbitals: ");
	m_log->input_message(molecule->get_ao_number());
	m_log->input_message("# of electrons: ");
	m_log->input_message( int(molecule->num_of_electrons));
	m_log->input_message("# of energy levels: ");
	m_log->input_message( int(molecule->MOnmb) );
	m_log->input_message("# of molecular orbital coefficients: ");
	molecule->ang_to_bohr();
	m_log->input_message("Starting do calculate the overlap 1e matrix: ");
	overlap_int overlap_integrals(molecule->atoms);
	overlap_integrals.calculate_overlap(molecule->m_overlap);
	molecule->bohr_to_ang();
	molecule->update();
	if ( !molecule->check() ){ 
		cout << "Problems in reading the .molden file: " << name_f << endl;
	}
	molecule->norm_orbs();
	m_log->input_message("Electronic energy: ");
	m_log->input_message( double_to_string(molecule->energy_tot) );
	m_log->input_message("HOMO energy: ");
	m_log->input_message( double_to_string(molecule->homo_energy) );
	m_log->input_message("HOMO number: ");
	m_log->input_message( double_to_string(molecule->homoN) );
	m_log->input_message("LUMO energy: ");
	m_log->input_message( double_to_string(molecule->lumo_energy) );
	m_log->input_message("LUMO number: ");
	m_log->input_message( double_to_string(molecule->lumoN) );
}

/********************************************************/
terachem_files::~terachem_files(){}
//////////////////////////////////////////////////////////
