// source file for the local_rd class 
// local_rd.cpp

// include statements from c++ library
#include <iostream>
#include <string> 
#include <vector>
#include <cmath>
#include <vector>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <cstdlib>
#include <experimental/filesystem>
// include statements from PRIMORDiA-libs
#include "../include/common.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/QMdriver.h"
#include "../include/global_rd.h"
#include "../include/local_rd_cnd.h"
#include "../include/Iprotein.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::move;
using std::abs;
using std::unique_ptr;
namespace fs = std::experimental::filesystem;

const double precision = 1e-08;

/***********************************************************************************/
local_rd_cnd::local_rd_cnd()	:
	name("nonamed")				,
	finite_diff(true)					,
	charge(0)							,
	ed(false)							,
	mep(false)							{		
}
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(Imolecule&& mol) noexcept:
	name(mol.name)													,
	finite_diff(false)														,
	charge(1)																,
	ed(false)																,
	mep(false)																,
	molecule( move(mol) )											{
	
	EAS.resize(molecule.num_of_atoms);
	NAS.resize(molecule.num_of_atoms);
	RAS.resize(molecule.num_of_atoms);
	dual.resize(molecule.num_of_atoms);
	Softness_Dual.resize(molecule.num_of_atoms);
	Hyper_softness.resize(molecule.num_of_atoms);
	hardness_A.resize(molecule.num_of_atoms);
	hardness_B.resize(molecule.num_of_atoms);
	hardness_C.resize(molecule.num_of_atoms);
	hardness_D.resize(molecule.num_of_atoms);
	multifilic.resize(molecule.num_of_atoms);
	electrophilicity.resize(molecule.num_of_atoms);
	relative_EAS.resize(molecule.num_of_atoms);
	relative_NAS.resize(molecule.num_of_atoms);
	electron_density.resize(molecule.num_of_atoms);
	fukushima.resize(molecule.num_of_atoms);
}
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(const Imolecule& mol_neut	,
										const Imolecule& mol_cation	, 
										const Imolecule& mol_anion)	:
	name(mol_neut.name)												, 
	finite_diff(true)															,
	charge(1)																	,
	ed(false)																	,
	mep(false)																	,
	molecule( Imolecule() )												{
	
	molecule.atoms					= mol_neut.atoms;
	molecule.num_of_atoms	= mol_neut.num_of_atoms;
	molecule.num_of_electrons= mol_neut.num_of_electrons;
	//molecule.update();
		
	EAS.resize(molecule.num_of_atoms);
	NAS.resize(molecule.num_of_atoms);
	RAS.resize(molecule.num_of_atoms);
	dual.resize(molecule.num_of_atoms);
	Softness_Dual.resize(molecule.num_of_atoms);
	Hyper_softness.resize(molecule.num_of_atoms);
	hardness_A.resize(molecule.num_of_atoms);
	hardness_B.resize(molecule.num_of_atoms);
	hardness_C.resize(molecule.num_of_atoms);
	hardness_D.resize(molecule.num_of_atoms);
	multifilic.resize(molecule.num_of_atoms);
	electrophilicity.resize(molecule.num_of_atoms);
	relative_EAS.resize(molecule.num_of_atoms);
	relative_NAS.resize(molecule.num_of_atoms);
	electron_density.resize(molecule.num_of_atoms);
	fukushima.resize(molecule.num_of_atoms);	
	
	for (unsigned int i=0;i<molecule.num_of_atoms;i++){
		EAS[i]				= ( (-molecule.atoms[i].charge)		- (-mol_cation.atoms[i].charge) )/charge;
		NAS[i]				= ( (-mol_anion.atoms[i].charge)		- (-molecule.atoms[i].charge) ) /charge;
		RAS[i]				=  ( ( (-mol_anion.atoms[i].charge) 	- (-mol_cation.atoms[i].charge) )/2 ) /charge;
		dual[i]				= NAS[i] - EAS[i];
		relative_EAS[i]	= EAS[i]/(NAS[i]+precision);
		relative_NAS[i] 	= NAS[i]/(EAS[i]+precision);
	}
}         
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(const local_rd_cnd& lrd_rhs)	:
	name(lrd_rhs.name)													,
	ed(lrd_rhs.ed)															,
	mep(lrd_rhs.mep)														,
	finite_diff(lrd_rhs.finite_diff)										,
	charge(lrd_rhs.charge)												,
	molecule( lrd_rhs.molecule )										,
	EAS(lrd_rhs.EAS)														, 
	NAS(lrd_rhs.NAS)														,
	RAS(lrd_rhs.RAS)														,
	dual(lrd_rhs.dual)														,
	Softness_Dual(lrd_rhs.Softness_Dual)						,
	Hyper_softness(lrd_rhs.Hyper_softness)						,
	hardness_A(lrd_rhs.hardness_A)								,
	hardness_B(lrd_rhs.hardness_B)								,
	hardness_C(lrd_rhs.hardness_C)								,
	hardness_D(lrd_rhs.hardness_D)								,
	multifilic(lrd_rhs.multifilic)											,
	electrophilicity(lrd_rhs.electrophilicity)						,
	electron_density(lrd_rhs.electron_density)					,
	relative_EAS(lrd_rhs.relative_EAS)								,
	relative_NAS(lrd_rhs.relative_NAS)								,
	fukushima(lrd_rhs.fukushima)									{
}
/***********************************************************************************/
local_rd_cnd& local_rd_cnd::operator=(const local_rd_cnd& lrd_rhs){
	if( this!=&lrd_rhs){
		name				= lrd_rhs.name;
		ed 					= lrd_rhs.ed;
		mep					= lrd_rhs.mep;
		finite_diff			= lrd_rhs.finite_diff;                 
		charge				= lrd_rhs.charge;                          
		molecule			= lrd_rhs.molecule;     
		EAS					= lrd_rhs.EAS;              
		NAS					= lrd_rhs.NAS;             
		RAS					= lrd_rhs.RAS;       
		dual					= lrd_rhs.dual;            
		Softness_Dual	= lrd_rhs.Softness_Dual;
		Hyper_softness	= lrd_rhs.Hyper_softness;
		hardness_A		= lrd_rhs.hardness_A;   
		hardness_B		= lrd_rhs.hardness_B;   
		hardness_C		= lrd_rhs.hardness_C;   
		hardness_D		= lrd_rhs.hardness_D;   
		multifilic			= lrd_rhs.multifilic;         
		electrophilicity	= lrd_rhs.electrophilicity;         
		relative_EAS		= lrd_rhs.relative_EAS;      
		relative_NAS		= lrd_rhs.relative_NAS;       
		fukushima		= lrd_rhs.fukushima;
	}       	
	return *this;
}                
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(local_rd_cnd&& lrd_rhs) noexcept	:
	name( move(lrd_rhs.name) )												,
	ed( move(lrd_rhs.ed) )															,
	mep( move(lrd_rhs.mep) )													,
	finite_diff( move(lrd_rhs.finite_diff) )									,
	charge( move(lrd_rhs.charge) )											,
	molecule( move(lrd_rhs.molecule) )										,
	EAS( move(lrd_rhs.EAS) )													, 
	NAS( move(lrd_rhs.NAS) )													,
	RAS( move(lrd_rhs.RAS) )													,
	dual( move(lrd_rhs.dual) )													,
	Softness_Dual( move(lrd_rhs.Softness_Dual) )						,
	Hyper_softness( move(lrd_rhs.Hyper_softness) )					,
	hardness_A( move(lrd_rhs.hardness_A) )								,
	hardness_B( move(lrd_rhs.hardness_B) )								,
	hardness_C( move(lrd_rhs.hardness_C) )								,
	hardness_D( move(lrd_rhs.hardness_D) )							,
	multifilic( move(lrd_rhs.multifilic) )										,
	electrophilicity( move(lrd_rhs.electrophilicity) )						,
	relative_EAS( move(lrd_rhs.relative_EAS) )							,
	relative_NAS( move(lrd_rhs.relative_NAS) )							,
	fukushima( move(lrd_rhs.fukushima) )								{
}
/***********************************************************************************/
local_rd_cnd& local_rd_cnd::operator=(local_rd_cnd&& lrd_rhs) noexcept {
	if( this!=&lrd_rhs){
		name					= move(lrd_rhs.name);  
		ed 						= move(lrd_rhs.ed);
		mep 						= move(lrd_rhs.mep);
		finite_diff				= move(lrd_rhs.finite_diff);                 
		charge					= move(lrd_rhs.charge);                          
		molecule				= move(lrd_rhs.molecule);     
		EAS						= move(lrd_rhs.EAS);              
		NAS						= move(lrd_rhs.NAS);             
		RAS						= move(lrd_rhs.RAS);       
		dual						= move(lrd_rhs.dual);            
		Softness_Dual		= move(lrd_rhs.Softness_Dual);
		Hyper_softness		= move(lrd_rhs.Hyper_softness);
		electron_density 	= move(lrd_rhs.electron_density);
		hardness_A			= move(lrd_rhs.hardness_A);   
		hardness_B			= move(lrd_rhs.hardness_B);   
		hardness_C			= move(lrd_rhs.hardness_C);   
		hardness_D			= move(lrd_rhs.hardness_D);   
		multifilic				= move(lrd_rhs.multifilic);         
		electrophilicity		= move(lrd_rhs.electrophilicity);         
		relative_EAS			= move(lrd_rhs.relative_EAS);      
		relative_NAS			= move(lrd_rhs.relative_NAS);  
		fukushima			= move(lrd_rhs.fukushima);
	}       	
	return *this;
}
/***********************************************************************************/
local_rd_cnd operator-(const local_rd_cnd& lrd_lhs,const local_rd_cnd& lrd_rhs){
	local_rd_cnd Result = lrd_lhs;
	if ( lrd_lhs.molecule.num_of_atoms == lrd_rhs.molecule.num_of_atoms) {
			for (unsigned int i=0;i<lrd_lhs.molecule.num_of_atoms;i++){
				Result.EAS[i]						= lrd_lhs.EAS[i]						- lrd_rhs.EAS[i];
				Result.NAS[i]						= lrd_lhs.NAS[i]						- lrd_rhs.NAS[i];
				Result.RAS[i]						= lrd_lhs.RAS[i]						- lrd_rhs.RAS[i];
				Result.dual[i]						= lrd_lhs.dual[i]						- lrd_rhs.dual[i];
				Result.Softness_Dual[i]		= lrd_lhs.Softness_Dual[i] 		- lrd_rhs.Softness_Dual[i];
				Result.Hyper_softness[i]		= lrd_lhs.Hyper_softness[i] 	- lrd_rhs.Hyper_softness[i];
				Result.hardness_A[i]			= lrd_lhs.hardness_A[i]			- lrd_rhs.hardness_A[i];
				Result.hardness_B[i]			= lrd_lhs.hardness_B[i]			- lrd_rhs.hardness_B[i];
				Result.hardness_C[i]			= lrd_lhs.hardness_C[i]			- lrd_rhs.hardness_C[i];
				Result.hardness_D[i]			= lrd_lhs.hardness_D[i]			- lrd_rhs.hardness_D[i];
				Result.electron_density[i] 	= lrd_lhs.electron_density[i] 	- lrd_rhs.electron_density[i];
				Result.multifilic[i]				= lrd_lhs.multifilic[i]				- lrd_rhs.multifilic[i];				
				Result.electrophilicity[i]		= lrd_lhs.electrophilicity[i]		- lrd_rhs.electrophilicity[i];
				Result.relative_EAS[i]			= lrd_lhs.relative_EAS[i]			- lrd_rhs.relative_EAS[i];
				Result.relative_NAS[i]			= lrd_lhs.relative_NAS[i]			- lrd_rhs.relative_NAS[i];
				Result.fukushima[i]			= lrd_lhs.fukushima[i]			- lrd_rhs.fukushima[i];
			}
	}
	return Result;
}		
/***********************************************************************************/	
void local_rd_cnd::calculate_Fukui(){
	int i =0;
	if ( !finite_diff ) {
		unique_ptr<QMdriver> qm_calc ( new QMdriver( move(molecule) ) );
		#pragma omp parallel for default(shared) private(i) 
		for (i=0;i<molecule.num_of_atoms;i++){
			EAS[i]				= qm_calc->EAS_in_atoms(i,0);
			NAS[i]				= qm_calc->NAS_in_atoms(i,0);
			fukushima[i]		= qm_calc->fukushima(i,0);
			RAS[i]				= (EAS[i] + NAS[i])/2;
			dual[i]				= NAS[i] - EAS[i];
			relative_EAS[i]	= EAS[i]/(NAS[i]+precision);
			relative_NAS[i]	= NAS[i]/(EAS[i]+precision);
		}
		//molecule.reset( new Imolecule(*qm_calc->molecule) );
		//molecule.reset( new Imolecule() );
		molecule = move( qm_calc->molecule );
	}
}
/***********************************************************************************/
void local_rd_cnd::calculate_Fukui_band(int band){
	unique_ptr<QMdriver> qm_calc ( new QMdriver( move(molecule) ) );
	unsigned int nof,i;
	nof  = molecule.num_of_atoms;
	omp_set_num_threads(NP);
	#pragma omp parallel for default(shared) private(i)	
		for (i=0;i<nof;i++){
			EAS[i] 				= qm_calc->EAS_in_atoms(i,band);
			NAS[i] 				= qm_calc->NAS_in_atoms(i,band);
			fukushima[i]		= qm_calc->fukushima(i,band);
			RAS[i] 				= (EAS[i] + NAS[i])/2;
			dual[i] 				= NAS[i] - EAS[i];
		}
		//molecule.reset( new Imolecule(*qm_calc->molecule) );
		//molecule.reset( new Imolecule() );
		molecule = move( qm_calc->molecule );
}
/***********************************************************************************/
void local_rd_cnd::calculate_Fukui_EW(int band){
	unsigned int nof ,i;
	nof  = molecule.num_of_atoms;
	unique_ptr<QMdriver> qm_calc ( new QMdriver( move(molecule) ) );
	omp_set_num_threads(NP);
	#pragma omp parallel for default(shared) private(i) 
	for (i=0;i<nof;i++){
		EAS[i]				= qm_calc->EAS_in_atoms_EW(i);
		NAS[i] 				= qm_calc->NAS_in_atoms_EW(i);
		fukushima[i]		= qm_calc->fukushima(i,band);
		RAS[i] 				= (EAS[i] + NAS[i])/2;
		dual[i]				= NAS[i] - EAS[i];
	}
	molecule = move(qm_calc->molecule);
}
/***********************************************************************************/
void local_rd_cnd::calculate_RD(const global_rd& grd){
	for(unsigned int i=0;i<molecule.num_of_atoms;i++){
		Softness_Dual[i]	= dual[i]*grd.softness;
		Hyper_softness[i]	= RAS[i]*grd.softness;
		multifilic[i]				= dual[i]*grd.Electrophilicity;
		electrophilicity[i]	= NAS[i]*grd.Electrophilicity;
	}
}
/***********************************************************************************/
void local_rd_cnd::calculate_Hardness(const global_rd& grd){
	unsigned int nof = molecule.num_of_atoms;
	unsigned int i,j;
	
	unique_ptr<QMdriver> qm_calc ( new QMdriver( move(molecule) ) );
	
	if ( !finite_diff ){
	#pragma omp parallel for default(shared) private(i) 
		for(i=0;i<nof;i++){
			electron_density[i] = qm_calc->density_in_atoms(i);
		}
	}else{
		vector<Iatom> atoms = qm_calc->molecule.atoms;
		#pragma omp parallel for default(shared) private(i) 
		for(i=0;i<nof;i++){
			electron_density[i] = atoms[i].atomicN - atoms[i].charge;
		}
	}
	
	#pragma omp parallel for default(shared) private(i) 
	for(i=0;i<nof;i++){
		hardness_A[i] = (EAS[i] - NAS[i]/nof)
							*(grd.chemical_pot/nof)
							+ (electron_density[i]/nof)*grd.hardness*2;
	}
	molecule = move(qm_calc->molecule);
	double r  = 0;
	double xi, yi,zi;
	
	#pragma omp parallel for default(shared) private(i)
	for (i=0;i<nof;i++){
		for (j=0;j<nof;j++){
			if ( i != j){
				xi = molecule.atoms[i].xcoord - molecule.atoms[j].xcoord;
				xi *= xi; 
				yi = molecule.atoms[i].ycoord - molecule.atoms[j].ycoord;
				yi *= yi;
				zi = molecule.atoms[i].zcoord - molecule.atoms[j].zcoord;
				zi *= zi;
				r  = sqrt(xi+yi+zi);
				hardness_B[i] += electron_density[j]/r;
			}
		}
		hardness_B[i] /= qm_calc->molecule.num_of_electrons;
	}
	
	
	r  = 0;
	for (i=0;i<nof;i++){
		for (j=0;j<nof;j++){
			if ( i != j){
				xi = molecule.atoms[i].xcoord - molecule.atoms[j].xcoord;
				xi *= xi; 
				yi = molecule.atoms[i].ycoord - molecule.atoms[j].ycoord;
				yi *= yi;
				zi = molecule.atoms[i].zcoord - molecule.atoms[j].zcoord;
				zi *= zi;
				r  = sqrt(xi+yi+zi);
				hardness_C[i] += EAS[j]/r;
			}
		}
	}
	for(i=0;i<nof;i++){	hardness_D[i] = NAS[i]*grd.lumo_en - EAS[i]*grd.homo_en; }
}
/*************************************************************************************/
void local_rd_cnd::write_rd_protein(const Iprotein& prot){
	std::string temps,temps1,temps2;
	vector< vector< vector<double> > > res_rd;
	res_rd.resize( prot.residues.size() );
	unsigned int i,nof=0;
	nof = molecule.atoms.size();
	
	vector<double> lig_interaction(14);
	if ( prot.ligand ){
		int lig_n = 0;
		
		double xc, yc,zc,r,xi,yi,zi = 0.0;
		xc = molecule.atoms[prot.residues[prot.lig_num].atoms_index[0]].xcoord;
		yc = molecule.atoms[prot.residues[prot.lig_num].atoms_index[0]].ycoord;
		zc = molecule.atoms[prot.residues[prot.lig_num].atoms_index[0]].zcoord;
		for( int j =0;j<nof;j++){
			xi = molecule.atoms[j].xcoord-xc;
			yi = molecule.atoms[j].ycoord-yc;
			zi = molecule.atoms[j].zcoord-zc;
			r = sqrt(xi*xi + yi*yi + zi*zi);
			if ( r != 0.0 ) {
				lig_interaction[0]	+= EAS[j]/r;
				lig_interaction[1]	+= NAS[j]/r;
				lig_interaction[2]	+= RAS[j]/r;
				lig_interaction[3]	+= dual[j]/r;
				lig_interaction[4]	+= Hyper_softness[j]/r;
				lig_interaction[5]	+= hardness_A[j]/r;
				lig_interaction[6]	+= hardness_B[j]/r;
				lig_interaction[7]	+= hardness_C[j]/r;
				lig_interaction[8]	+= hardness_D[j]/r;
				lig_interaction[9]	+= multifilic[j]/r;
				lig_interaction[10]	+= electrophilicity[j]/r;
				lig_interaction[11]	+= fukushima[j]/r;
				lig_interaction[12] 	+= electron_density[j]/r;
				lig_interaction[13]	+= Softness_Dual[j]/r;
			}
		}
	}
	
	for(unsigned int i=0;i<prot.residues.size();i++) {
		res_rd[i].resize(14);
		for(int j=0;j<14;j++)
			res_rd[i][j].resize(3);
	}

	int cnt = 0;
	for(unsigned int i=0;i<prot.residues.size();i++){
		for(unsigned int j=0;j<prot.residues[i].atom_s;j++){
			res_rd[i][0][0] += EAS[cnt];
			res_rd[i][1][0] += NAS[cnt];
			res_rd[i][2][0] += RAS[cnt];
			res_rd[i][3][0] += dual[cnt];
			res_rd[i][4][0] += Hyper_softness[cnt];
			res_rd[i][5][0] += hardness_A[cnt];
			res_rd[i][6][0] += hardness_B[cnt];
			res_rd[i][7][0] += hardness_C[cnt];
			res_rd[i][8][0] += hardness_D[cnt];
			res_rd[i][9][0] += multifilic[cnt];
			res_rd[i][10][0] += electrophilicity[cnt];
			res_rd[i][11][0] += fukushima[cnt];
			res_rd[i][12][0] += electron_density[cnt];
			res_rd[i][13][0] += Softness_Dual[cnt];
			cnt++;
		}
	}
	for(unsigned int i=0;i<prot.residues.size();i++){		
		res_rd[i][0][1] = res_rd[i][0][0]/prot.residues[i].atom_s;
		res_rd[i][1][1] = res_rd[i][1][0]/prot.residues[i].atom_s;
		res_rd[i][2][1] = res_rd[i][2][0]/prot.residues[i].atom_s;
		res_rd[i][3][1] = res_rd[i][3][0]/prot.residues[i].atom_s;
		res_rd[i][4][1] = res_rd[i][4][0]/prot.residues[i].atom_s;
		res_rd[i][5][1] = res_rd[i][5][0]/prot.residues[i].atom_s;
		res_rd[i][6][1] = res_rd[i][6][0]/prot.residues[i].atom_s;
		res_rd[i][7][1] = res_rd[i][7][0]/prot.residues[i].atom_s;
		res_rd[i][8][1] = res_rd[i][8][0]/prot.residues[i].atom_s;
		res_rd[i][9][1] = res_rd[i][9][0]/prot.residues[i].atom_s;
		res_rd[i][10][1] = res_rd[i][10][0]/prot.residues[i].atom_s;
		res_rd[i][11][1] = res_rd[i][11][0]/prot.residues[i].atom_s;
		res_rd[i][12][1] = res_rd[i][12][0]/prot.residues[i].atom_s;
		res_rd[i][13][1] = res_rd[i][13][0]/prot.residues[i].atom_s;
	}
	cnt = 0;
	for(unsigned int i=0;i<prot.residues.size();i++){
		for(unsigned int j=0;j<prot.residues[i].atom_s;j++){
			if ( j == 0 ) res_rd[i][0][2] = EAS[cnt];
			else{
				if ( EAS[cnt] > EAS[cnt-1] ) 
					res_rd[i][0][2] = EAS[cnt];
			}
			if ( j == 0 ) res_rd[i][1][2] = NAS[cnt];
			else{
				if ( NAS[cnt] > NAS[cnt-1] ) 
					res_rd[i][1][2] = NAS[cnt];
			}
			if ( j == 0 ) res_rd[i][2][2] = RAS[cnt];
			else{
				if ( RAS[cnt] > RAS[cnt-1] ) 
					res_rd[i][2][2] = RAS[cnt];
			}
			if ( j == 0 ) res_rd[i][3][2] = dual[cnt];
			else{
				if ( abs(dual[cnt]) > abs(dual[cnt-1]) ) 
					res_rd[i][3][2] = dual[cnt];
			}
			if ( j == 0 ) res_rd[i][4][2] = Hyper_softness[cnt];
			else{
				if ( abs(Hyper_softness[cnt]) > abs(Hyper_softness[cnt-1]) ) 
					res_rd[i][4][2] = Hyper_softness[cnt];
			}
			if ( j == 0 ) res_rd[i][5][2] = hardness_A[cnt];
			else{
				if ( abs(hardness_A[cnt]) > abs(hardness_A[cnt-1]) ) 
					res_rd[i][5][2] = hardness_A[cnt];
			}
			if ( j == 0 ) res_rd[i][6][2] = hardness_B[cnt];
			else{
				if ( abs(hardness_B[cnt]) > abs(hardness_B[cnt-1]) ) 
					res_rd[i][6][2] = hardness_B[cnt];
			}
			if ( j == 0 ) res_rd[i][7][2] = hardness_C[cnt];
			else{
				if ( abs(hardness_C[cnt]) > abs(hardness_C[cnt-1]) ) 
					res_rd[i][7][2] = hardness_C[cnt];
			}
			if ( j == 0 ) res_rd[i][8][2] = hardness_D[cnt];
			else{
				if ( abs(hardness_D[cnt]) > abs(hardness_D[cnt-1]) ) 
					res_rd[i][8][2] = hardness_D[cnt];
			}
			if ( j == 0 ) res_rd[i][9][2] = multifilic[cnt];
			else{
				if ( abs(multifilic[cnt]) > abs(multifilic[cnt-1]) ) 
					res_rd[i][9][2] = multifilic[cnt];
			}
			if ( j == 0 ) res_rd[i][10][2] = electrophilicity[cnt];
			else{
				if ( abs(electrophilicity[cnt]) > abs(electrophilicity[cnt-1]) ) 
					res_rd[i][10][2] = electrophilicity[cnt];
			}
			if ( j == 0 ) res_rd[i][11][2] = fukushima[cnt];
			else{
				if ( fukushima[cnt] > fukushima[cnt-1] ) 
					res_rd[i][11][2] = fukushima[cnt];
			}
			if ( j == 0 ) res_rd[i][12][2] = electron_density[cnt];
			else{
				if ( fukushima[cnt] > electron_density[cnt-1] ) 
					res_rd[i][12][2] = electron_density[cnt];
			}
			if ( j == 0 ) res_rd[i][13][2] = Softness_Dual[cnt];
			else{
				if ( abs(Softness_Dual[cnt]) > abs(Softness_Dual[cnt-1]) ) 
					res_rd[i][13][2] = Softness_Dual[cnt];
			}
			cnt++;
		}
	}
	
	if ( finite_diff ) {
		temps  = name+"pro.lrd";
		temps1 = name+"pro_mean.lrd";
		temps2 = name+"pro_max.lrd";
	}else{
		temps  = name+"residues.lrd";
		temps1 = name+"residues_mean.lrd";
		temps2 = name+"residues_max.lrd";
	}
	const char* names  = temps.c_str();
	const char* names1 = temps1.c_str();
	const char* names2 = temps2.c_str();
	std::ofstream lrd_file(names);
	std::ofstream lrd_file1(names1);
	std::ofstream lrd_file2(names2);
	std::ofstream r_script_lrd;

	lrd_file.precision(8);
	lrd_file  << std::fixed;
	lrd_file1.precision(8);
	lrd_file1 << std::fixed;
	lrd_file2.precision(8);
	lrd_file2 << std::fixed;
	
	lrd_file << "res EAS NAS RAS Dual Softness Hardness_A Hardness_B Hardness_C Hardness_D  Multiphilic Electrophilic Fukushima Electron_Density Softness_dual\n";
	for(unsigned int i=0;i<prot.residues.size();i++){
		lrd_file << (i+1)                        << ""
				 << prot.residues[i].type        << " "
				 << res_rd[i][0][0]              << " " 
				 << res_rd[i][1][0]              << " " 
				 << res_rd[i][2][0]              << " " 
				 << res_rd[i][3][0]              << " " 
				 << res_rd[i][4][0]              << " " 
				 << res_rd[i][5][0]              << " " 
				 << res_rd[i][6][0]              << " " 
				 << res_rd[i][7][0]              << " " 
				 << res_rd[i][8][0]              << " "
				 << res_rd[i][9][0]              << " "
				 << res_rd[i][10][0]            << " "
				 << res_rd[i][11][0]            << " "
				 << res_rd[i][12][0]            << " "
				 << res_rd[i][13][0]            << " "
				 << "\n";
	}
	if ( prot.ligand ){
		lrd_file << (i+1)										<< ""
					 << "ligand_interaction_sum"		<< " "
					<<  lig_interaction[0]					<< " " 
					<<  lig_interaction[1]					<< " " 
					<<  lig_interaction[2]					<< " " 
					<<  lig_interaction[3]					<< " " 
					<<  lig_interaction[4]					<< " " 
					<<  lig_interaction[5]					<< " " 
					<<  lig_interaction[6]					<< " " 
					<<  lig_interaction[7]					<< " " 
					<<  lig_interaction[8]					<< " "
					<<  lig_interaction[9]					<< " "
					<<  lig_interaction[10]				<< " "
					<< "\n";
	}
	lrd_file.close();
	/*
	lrd_file1 << "res EAS NAS RAS Dual Softness Hardness Multiphilic Electrophilic Fukushima Electron_Density Softness_dual\n";
	for(unsigned int i=0;i<prot.residues.size();i++){
		lrd_file1<< (i+1)							<< ""
				 << prot.residues[i].type		<< " "
				 << res_rd[i][0][1]				<< " " 
				 << res_rd[i][1][1]				<< " " 
				 << res_rd[i][2][1]				<< " " 
				 << res_rd[i][3][1]				<< " " 
				 << res_rd[i][4][1] 				<< " " 
				 << res_rd[i][5][1]				<< " " 
				 << res_rd[i][6][1]				<< " " 
				 << res_rd[i][7][1]				<< " " 
				 << res_rd[i][8][1]				<< " "
				 << res_rd[i][9][1]				<< " "
				 << res_rd[i][10][1]				<< " "
				 << "\n";
	}
	lrd_file1.close();
	lrd_file2 << "res EAS NAS RAS Dual Softness Hardness Multiphilic Electrophilic Fukushima Electron_density Softness_dual\n";
	for(unsigned int i=0;i<prot.residues.size();i++){
		lrd_file2<< (i+1)                         << ""
				 << prot.residues[i].type  << " "
				 << res_rd[i][0][2]              << " " 
				 << res_rd[i][1][2]              << " " 
				 << res_rd[i][2][2]              << " " 
				 << res_rd[i][3][2]              << " " 
				 << res_rd[i][4][2]              << " " 
				 << res_rd[i][5][2]              << " " 
				 << res_rd[i][6][2]              << " " 
				 << res_rd[i][7][2]              << " " 
				 << res_rd[i][8][2]              << " " 
				 << res_rd[i][9][2]              << " " 
				 << res_rd[i][10][2]            << " " 
				 << "\n";
	}
	lrd_file2.close();
	*/
}
/*******************************************************************************************/
void local_rd_cnd::write_rd_protein_pdb(const Iprotein& protein){
	
	Iprotein prot = protein;
	pdb rd_results;
	rd_results.name = this->name;
	prot.load_b_column(EAS);
	prot.title 		= "EAS";
	prot.remark	= "Electrophilic Attack Suscptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(NAS);
	prot.title 		= "NAS";
	prot.remark	= "Nucleophilic Attack Susceptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(RAS);
	prot.title 		= "RAS";
	prot.remark	= "Radical Attack Susceptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(dual);
	prot.title 		= "Dual";
	prot.remark	= "Local Dual Descriptor Attack Susceptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(hardness_A);
	prot.title 		= "hardness_A";
	prot.remark	= "Local Hardness (local chemical potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(hardness_B);
	prot.title 		= "hardness_B";
	prot.remark	= "Local Hardness (Electron-electron potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(hardness_C);
	prot.title 		= "hardness_C";
	prot.remark	= "Local Hardness (Fukui potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(hardness_D);
	prot.title 		= "hardness_D";
	prot.remark	= "Local Hardness (Fukui distribution based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(Softness_Dual);
	prot.title 		= "softness";
	prot.remark	= "Local Softness Dual at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(Hyper_softness);
	prot.title 		= "hypersoftness";
	prot.remark	= "Local Hyper Softness at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(fukushima);
	prot.title 		= "fukushima";
	prot.remark	= "Localization of frontier band orbitals at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(multifilic);
	prot.title 		= "multiphilic";
	prot.remark	= "Local Multiphilic descriptor  at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(electrophilicity);
	prot.title 		= "electrophilicity";
	prot.remark	= "Local Electrophilicity descriptor  at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);	
	prot.load_b_column(electron_density);
	prot.title 		= "electron_density";
	prot.remark	= "Electronic density per atom at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	m_log->input_message("Finishing the writting of the condensed local reactivity descriptors in PDBs.");

	vector<double> pcharges(molecule.atoms.size());
	for(unsigned int i=0; i<molecule.atoms.size();i++) {
		pcharges[i] = molecule.atoms[i].charge;
	}
	prot.load_b_column(pcharges);
	prot.title 		= "mep";
	prot.remark	= "Partial charge per atom at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	rd_results.write_pdb(rd_results.name +"_RD.pdb");
	fs::create_directory(name+"_PDB_RD");
	rd_results.write_models(name+"_PDB_RD");
}
/***************************************************************************************/
void local_rd_cnd::write_rd_protein_reaction(const Iprotein& prot){
	if ( prot.p1.size() > 0 ){
		string pair_n = name;
		if ( finite_diff ){
			pair_n = name + "_FD.prd";
		}else {
			pair_n = name + "_FOA.prd";
		}
		std::ofstream pair( pair_n);
		vector<double> values;	
		double temp = 0;
		double x,y,z,r,ri =0;
		
		if ( prot.p2.size() > 0 ){
			x = molecule.atoms[ prot.p1[0] ] .xcoord - molecule.atoms[ prot.p1[1] ] .xcoord;
			y = molecule.atoms[ prot.p1[0] ] .ycoord - molecule.atoms[ prot.p1[1] ] .ycoord;
			z = molecule.atoms[ prot.p1[0] ] .zcoord - molecule.atoms[ prot.p1[1] ] .zcoord;
			r = sqrt(x*x + y*y + z*z );
			x = molecule.atoms[ prot.p2[0] ] .xcoord - molecule.atoms[ prot.p2[1] ] .xcoord;
			y = molecule.atoms[ prot.p2[0] ] .ycoord - molecule.atoms[ prot.p2[1] ] .ycoord;
			z = molecule.atoms[ prot.p2[0] ] .zcoord - molecule.atoms[ prot.p2[1] ] .zcoord;
			ri = sqrt(x*x + y*y + z*z );			
			if (  dual[ prot.p1[0] ] <=  dual[ prot.p1[1] ] ){
				temp = EAS[ prot.p1[0] ]*NAS[ prot.p1[0] ]/r;
				values.push_back(temp);
			}else if( dual[ prot.p1[0] ] >  dual[ prot.p1[1] ] ){
				temp =- NAS[ prot.p1[0] ]*EAS[ prot.p1[0] ]/r;
				values.push_back(temp);
			}
			if (  dual[ prot.p2[0] ] <= dual[ prot.p2[1] ] ){
				temp = EAS[ prot.p2[0] ]*NAS[ prot.p2[0] ]/ri;
				values.push_back(temp);
			}else if( dual[ prot.p2[0] ] >  dual[ prot.p2[1] ] ){
				temp = NAS[ prot.p2[0] ]*EAS[ prot.p2[0] ]/ri;
				values.push_back(temp);
			}
			temp = Hyper_softness[ prot.p1[0] ]*Hyper_softness[ prot.p1[1] ]/r;
			values.push_back(temp);
			temp = Hyper_softness[ prot.p2[0] ]*Hyper_softness[ prot.p2[1] ]/ri;
			values.push_back(temp);
			temp = hardness_A[ prot.p1[0] ]*hardness_A[ prot.p1[1] ]/r;
			values.push_back(temp);
			temp = hardness_A[ prot.p2[0] ]*hardness_A[ prot.p2[1] ]/ri;
			values.push_back(temp);
			
			pair << "RD atom1 atom2 atom3 atom4 pair12 pair34 \n"
					<< "dist " << r << " " << ri << "\n "
					<< "EAS " << EAS[ prot.p1[0] ] << "  " << EAS[ prot.p1[1]]  << " " << EAS[ prot.p2[0] ] << " " << EAS[ prot.p2[1] ] << "  - - \n"
					<< "NAS " << NAS[ prot.p1[0] ] << "  " << NAS[ prot.p1[1]]  << " " << NAS[ prot.p2[0] ] << " " << NAS[ prot.p2[1] ] << "  - - \n"
					<< "RAS " << RAS[ prot.p1[0] ] << "  " << RAS[ prot.p1[1]]  << " " << RAS[ prot.p2[0] ] << " " << RAS[ prot.p2[1] ] << "  - - \n"
					<< "dual " << dual[ prot.p1[0] ] << "  " << dual[ prot.p1[1]]  << " " << dual[ prot.p2[0] ] << " " << dual[ prot.p2[1] ] << "  - - \n"
					<< "Hardness " << hardness_A[ prot.p1[0] ] << "  " << hardness_A[ prot.p1[1]]  << " " << hardness_A[ prot.p2[0] ] << " " << hardness_A[ prot.p2[1] ] << "  - - \n"
					<< "Charge " << molecule.atoms[ prot.p1[0] ].charge << "  " << molecule.atoms[ prot.p1[1] ].charge  << " " << molecule.atoms[ prot.p2[0] ].charge << " " << molecule.atoms[ prot.p2[1] ].charge << "  - - \n"
					<< "Elecdens " << electron_density[ prot.p1[0] ] << "  " << electron_density[ prot.p1[1]]  << " " << electron_density[ prot.p2[0] ] << " " << electron_density[ prot.p2[1] ] << "  - - \n";
			pair << "CT  - - - - " << values[0]  << " " << values[1] << "\n"
					<< "SPI  - - - - " << values[2] << " " << values[3] << "\n"
					<< "HPI  - - - - " << values[4] << " " << values[5] << endl;
					
		}else{
			x = molecule.atoms[ prot.p1[0] ] .xcoord - molecule.atoms[ prot.p1[1] ] .xcoord;
			y = molecule.atoms[ prot.p1[0] ] .ycoord - molecule.atoms[ prot.p1[1] ] .ycoord;
			z = molecule.atoms[ prot.p1[0] ] .zcoord - molecule.atoms[ prot.p1[1] ] .zcoord;
			r = sqrt(x*x + y*y + z*z );
			if (  dual[ prot.p1[0] ] <=  dual[ prot.p1[1] ] ){
				temp = EAS[ prot.p1[0] ]*NAS[ prot.p1[0] ]/r;
				values.push_back(temp);
			}else if( dual[ prot.p1[0] ] >  dual[ prot.p1[1] ] ){
				temp = NAS[ prot.p1[0] ]*EAS[ prot.p1[0] ]/r;
				values.push_back(temp);
			}
			temp = Hyper_softness[ prot.p1[0] ]*Hyper_softness[ prot.p1[1] ]/r;
			values.push_back(temp);			
			temp = hardness_A[ prot.p1[0] ]*hardness_A[ prot.p1[1] ]/r;
			values.push_back(temp);
			pair << "RD atom1 atom2 pair12 \n"
					<< "dist " << r <<  " \n"
					<< "EAS " << EAS[ prot.p1[0] ] << "  " << EAS[ prot.p1[1]]  << "   - - \n"
					<< "NAS " << NAS[ prot.p1[0] ] << "  " << NAS[ prot.p1[1]]  << "  - - \n"
					<< "RAS " << RAS[ prot.p1[0] ] << "  " << RAS[ prot.p1[1]]  << "  - - \n"
					<< "dual " << dual[ prot.p1[0] ] << "  " << dual[ prot.p1[1]]  <<  "  - - \n"
					<< "Hardness " << hardness_A[ prot.p1[0] ] << "  " << hardness_A[ prot.p1[1]]  <<"  - - \n"
					<< "Charge " << molecule.atoms[ prot.p1[0] ].charge << "  " << molecule.atoms[ prot.p1[1]].charge  << "  - - \n";
			pair << "CT  - - - - " << values[0]  << "\n"
					<< "SPI  - - - - " << values[1] << "\n"
					<< "HPI  - - - - " << values[2] <<  endl;
		}
		pair.close();
	}
}
/*************************************************************************************/
void local_rd_cnd::write_LRD(){
	std::string temps;
	if ( finite_diff ) { 	temps = name+"FD.lrd";
	}else{	temps = name+"FOA.lrd";	}
	const char* names = temps.c_str();
	std::ofstream lrd_file;
	lrd_file.open(names);
	
	lrd_file.precision(8);
	lrd_file << std::fixed;
	lrd_file << name << " " << "\n" <<  std::left;
	lrd_file << "n atom charge EAS NAS RAS dual Softness Hardness_A Hardness_B Hardness_C Hardness_D Multiphilic Electrohilicity Fukushima Electron_density Softness_dual\n";
	
	for(int i=0;i<molecule.num_of_atoms;i++){
		lrd_file	<< (i+1) 
					<< " "
					<< molecule.atoms[i].element	<< " "
					<< molecule.atoms[i].charge	<< " "
					<< EAS[i]								<< " "
					<< NAS[i]								<< " "
					<< RAS[i]								<< " "
					<< dual[i]								<< " "
					<< Hyper_softness[i] 				<< " "
					<< hardness_A[i]					<< " "
					<< hardness_B[i]					<< " "
					<< hardness_C[i]					<< " "
					<< hardness_D[i]					<< " "
					<< multifilic[i]						<< " "
					<< electrophilicity[i]				<< " "
					<< fukushima[i] 					<< " "
					<< electron_density[i]			 	<< " "
					<< Softness_Dual[i]				<< "\n";
	}
	lrd_file.close();
	m_log->input_message("Finishing the writting of the condensed local reactivity descriptors.");
}
/*************************************************************************************/	

//==============================================================================