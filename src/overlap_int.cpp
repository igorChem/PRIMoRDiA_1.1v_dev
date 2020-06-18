#include <iostream>
#include <string>
#include <vector>

#include "../include/log_class.h"
#include "../include/common.h"
#include "../include/Iaorbital.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/overlap_int.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
// Libint Gaussian integrals library
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#  include <libint2/statics_definition.h>
#endif
using libint2::Shell;
using std::vector;
using std::cout;
using std::endl;
/************************************************/
overlap_int::overlap_int(){}
/************************************************/
overlap_int::overlap_int(const vector<Iatom>& atoms):
	nAO(0)											{
	int spd = 0;
	for(int i=0;i<atoms.size();i++){
		for(int j=0;j<atoms[i].orbitals.size();j++){
			spd=atoms[i].orbitals[j].powx +
				atoms[i].orbitals[j].powy +
				atoms[i].orbitals[j].powz;
			//std::cout << i << " " << j << " " << atoms[i].orbitals[j].symmetry << " " << atoms[i].element << std::endl;
			if ( atoms[i].orbitals[j].powy == 0 && atoms[i].orbitals[j].powz == 0 ){
				if ( spd == 0 ) nAO += 1;
				if ( spd == 1 ) nAO += 3;
				if ( spd == 2 ) nAO += 6;
				if		( atoms[i].orbitals[j].gtos.size() == 1 ){
					shells.push_back( 
						{
							{atoms[i].orbitals[j].gtos[0].exponent},
							{
								{spd,false,{atoms[i].orbitals[j].gtos[0].c_coef}}
							},
							{{atoms[i].xcoord,atoms[i].ycoord,atoms[i].zcoord}}
						}
					);
				}
				else if	( atoms[i].orbitals[j].gtos.size() == 3 ){
					shells.push_back( 
						{
							{atoms[i].orbitals[j].gtos[0].exponent,atoms[i].orbitals[j].gtos[1].exponent,atoms[i].orbitals[j].gtos[2].exponent},
							{
								{spd,false,{atoms[i].orbitals[j].gtos[0].c_coef,atoms[i].orbitals[j].gtos[1].c_coef,atoms[i].orbitals[j].gtos[2].c_coef}}
							},
							{{atoms[i].xcoord,atoms[i].ycoord,atoms[i].zcoord}}
						}
					);
				}
				else if ( atoms[i].orbitals[j].gtos.size() == 6 ){
					shells.push_back( 
						{
							{atoms[i].orbitals[j].gtos[0].exponent,atoms[i].orbitals[j].gtos[1].exponent,atoms[i].orbitals[j].gtos[2].exponent,atoms[i].orbitals[j].gtos[3].exponent,atoms[i].orbitals[j].gtos[4].exponent,atoms[i].orbitals[j].gtos[5].exponent},
							{
								{spd,false,{atoms[i].orbitals[j].gtos[0].c_coef,atoms[i].orbitals[j].gtos[1].c_coef,atoms[i].orbitals[j].gtos[2].c_coef,atoms[i].orbitals[j].gtos[3].c_coef,atoms[i].orbitals[j].gtos[4].c_coef,atoms[i].orbitals[j].gtos[5].c_coef}}
							},
							{{atoms[i].xcoord,atoms[i].ycoord,atoms[i].zcoord}}
						}
					);
				}
			}
		}
	}
	
	for(unsigned int i=0;i<shells.size();i++ ){
		shell_size.push_back( shells[i].size() );
	}
	
	unsigned int n = 0;
	for (auto shell: shells) {
		shell_map.push_back(n);
		n += shell.size();
	}
	
	m_log->input_message("# of contracted shells: ");
	m_log->input_message( int( shells.size() ) );
	m_log->input_message("# of basis functions:	" );
	m_log->input_message( int(nAO) );
}
/************************************************/
void overlap_int::calculate_overlap(vector<double>& overlap_matrix){
	libint2::initialize();
	
	size_t	n = 0;
	int 	l = 0;
	for (unsigned int i=0;i<shells.size();i++){
		n = std::max(shells[i].nprim(), n);
		//std::cout << shells[i].nprim() << std::endl;
	}
	
	vector< vector <double> > matrix_holder(nAO);
	for(int i=0;i<nAO;i++){
		matrix_holder[i].resize(nAO);
	}
	
	for (auto shell: shells)
		for (auto c: shell.contr)
			l = std::max(c.l, l);
	
	libint2::Engine engine(libint2::Operator::overlap, n, l, 0);
	const auto& buf_vec = engine.results();
	
	auto ints_shellset = buf_vec[0];
	for( unsigned int i=0;i<shells.size();i++){
		for ( unsigned int j=0;j<shells.size();j++){
			engine.compute(shells[i],shells[j]);
			auto ints_shellset = buf_vec[0];
			auto bf1 = shell_map[i];
			auto n1  = shell_size[i];
			auto bf2 = shell_map[j];
			auto n2  = shell_size[j];
			
			for(int k=0;k<n1;k++){
				for(int m=0;m<n2;m++){
					//std::cout << "  " << bf1+k << " " << bf2+m << " " << ints_shellset[k*n2+m] << std::endl;
					matrix_holder[bf1+k][bf2+m] = ints_shellset[k*n2+m];
				}
			}
		}
	}
	
	for(int i=0;i<nAO;i++){
		for(int j=0;j<=i;j++){
			//std::cout << i <<  " " << j << " " << matrix_holder[i][j] << std::endl;
			overlap_matrix.push_back( matrix_holder[i][j]);
		}
	}
	
	libint2::finalize();
	
}
/************************************************/
overlap_int::~overlap_int(){}
