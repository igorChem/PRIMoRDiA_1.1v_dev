#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <omp.h>

#include "../include/log_class.h"
#include "../include/Itimer.h"
#include "../include/common.h"
#include "../include/Iline.h"
#include "../include/Ibuffer.h"
#include "../include/Iatom.h"
#include "../include/Iaorbital.h"
#include "../include/Imolecule.h"
#include "../include/Icube.h"
#include "../include/global_rd.h"
#include "../include/local_rd.h"
#include "../include/QMdriver.h"
#include "../include/gridgen.h"
#include "../include/QMparser.h"
#include "../include/primordia.h"
#include "../include/test_p.h"
#include "../include/compPrimordia.h"
#include "../include/autoprimordia.h"
#include "../include/Iprotein.h"

//#include <gtest/gtest.h>

using std::cout;
using std::endl;
using std::vector;
using std::unique_ptr;
using std::string;
using std::move;

/************************************************************************************/
test_p::test_p(){
	M_logfile = true;
	
	M_verbose = true;
	M_logfile 	=	true;
	string teste_folder1 = "/home/igorchem/Documents/primordia_data_test/mopac/";
	string teste_folder2 = "/home/igorchem/Documents/primordia_data_test/gamess/";
	string teste_folder3 = "/home/igorchem/Documents/primordia_data_test/gaussian/";
	string teste_folder4 = "/home/igorchem/Documents/primordia_data_test/orca/";
	string teste_folder5 = "/home/igorchem/Documents/primordia_data_test/others/";
	
	pdb_files.push_back(teste_folder1+"1l2y.pdb"); // # 0 
	pdb_files.push_back(teste_folder1+"1l2yPM7neutro.pdb"); // # 1 
	
	mopac_aux.push_back(teste_folder1+"1l2y_pm7.aux");		//#0
	mopac_aux.push_back(teste_folder1+"1l2y_pm7_an.aux"); 	//#1
	mopac_aux.push_back(teste_folder1+"1l2y_pm7_cat.aux"); 	//#2
	mopac_aux.push_back(teste_folder1+"1l2y_pm7_lmo.aux");	//#3
	mopac_aux.push_back(teste_folder1+"acrolein.aux");			//#4
	mopac_aux.push_back(teste_folder1+"acrolein_cat.aux");	//#5
	mopac_aux.push_back(teste_folder1+"acrolein_an.aux"); 	//#6
	mopac_aux.push_back(teste_folder1+"benzene.aux");			//#7
	mopac_aux.push_back(teste_folder1+"benzene_cat.aux");	//#8
	mopac_aux.push_back(teste_folder1+"benzene_an.aux");		//#9
	
	mopac_mgf.push_back(teste_folder1+"1l2yPM7neutro.mgf");//#0
	mopac_mgf.push_back(teste_folder1+"trp_lmo.mgf"); 			//#1
	mopac_mgf.push_back(teste_folder1+"benzene.mgf"); 			//#2
	mopac_mgf.push_back(teste_folder1+"benzene_cat.mgf"); 	//#3
	mopac_mgf.push_back(teste_folder1+"benzene_an.mgf"); 	//#4
	mopac_mgf.push_back(teste_folder1+"trp_pm6_neuto.mgf");//#5

	mopac_out.push_back(teste_folder1+"acrolein.out");		//#0
	mopac_out.push_back(teste_folder1+"acrolein_cat.out");	//#1
	mopac_out.push_back(teste_folder1+"acrolein_an.out"); 	//#2
	mopac_out.push_back(teste_folder1+"benzene.out");		//#3
	mopac_out.push_back(teste_folder1+"benzene_cat.out");	//#4
	mopac_out.push_back(teste_folder1+"benzene_an.out");	//#5

	orca_out.push_back(teste_folder4+"acrolein_orca.out");		//#0	
	orca_out.push_back(teste_folder4+"acrolein_orca_cat.out"); //#1
	orca_out.push_back(teste_folder4+"acrolein_orca_an.out");  //#2
	
	gauss_fchk.push_back(teste_folder3+"acrolein_gauss.fchk");//#0
	gauss_log.push_back(teste_folder3+"acrolein_gauss.log");     //#0
	
	gamess_log.push_back(teste_folder2+"trp_dft.log"); //#0
	gamess_log.push_back(teste_folder2+"trp_convHF.log"); //#1
	gamess_log.push_back(teste_folder2+"trp_cat_dft.log"); //#2
	gamess_log.push_back(teste_folder2+"trp_an_dft.log"); //#3
	gamess_log.push_back(teste_folder2+"acrolein_gam.log"); //#4
	gamess_log.push_back(teste_folder2+"acrolein_gam_cat.log"); //#5
	gamess_log.push_back(teste_folder2+"acrolein_gam_an.log"); //#6
	gamess_log.push_back(teste_folder2+"trp_an_dft.log"); //#7
	gamess_log.push_back(teste_folder2+"trp_an_dft.log"); //#8
}
/**********************************************************************************/
void test_p::test_primordia_1(){
	//testing option 1 of primordia class for mopac
	dos = true;
	primordia mop_a;
	mop_a.init_FOA(mopac_aux[4].c_str(),0,"mepEE",true,"mopac");
	//mop_a.init_FOA(orca_out[1].c_str(),0,"LCP",true,"orca");
	//mop_a.init_FOA(gamess_log[5].c_str(),40,"mepFukui",true,"gamess");
	//mop_a.init_FOA(gauss_fchk[0].c_str(),0,"fukui",true,"gaussian");
}
/**********************************************************************************/
void test_p::test_primordia_2(){
	//testing option 2 of primordia class for mopac
	dos = true;
	primordia mop_a;
	//mop_a.init_FD(mopac_out[3].c_str(),mopac_out[4].c_str(), mopac_out[5].c_str(),30,1,true,"mepEE","mopac");
	//mop_a.init_FD(gamess_log[4].c_str(),gamess_log[5].c_str(),gamess_log[6].c_str(),30,1,true,"mepFukui","gamess");
	mop_a.init_FD(orca_out[0].c_str(),orca_out[1].c_str(),orca_out[2].c_str(),30,1,true,"mepFukui","orca");
	
}
/**********************************************************************************/
void test_p::test_primordia_3(){
	//testing option 2 of primordia class for mopac
	dos = true;
	double ref_coord[3];
	ref_coord[0] = 0.779;
	ref_coord[1] = 0.524;
	ref_coord[2] = 2.812;
	primordia mop_a;
	//mop_a.init_protein_RD(mopac_aux[0].c_str(),"mepEE",40,10,ref_coord,0,pdb_files[0].c_str(),true,"EW","mopac"); // mopac rhf
	//mop_a.init_protein_RD(mopac_aux[3].c_str(),"mepEE",40,10,ref_coord,0,pdb_files[0].c_str(),true,"BD","mopac"); // mopac lmo
	//mop_a.init_protein_RD(mopac_aux[3].c_str(),"mepFukui",40,10,ref_coord,10,pdb_files[0].c_str(),true,"BD","mopac"); // mopac lmo
	//mop_a.init_protein_RD(mopac_aux[1].c_str(),"mepFukui",40,10,ref_coord,0,pdb_files[0].c_str(),true,"BD","mopac"); // mopac lmo
	//mop_a.init_protein_RD(mopac_aux[2].c_str(),"mepFukui",40,10,ref_coord,0,pdb_files[0].c_str(),true,"BD","mopac"); // mopac lmo
	//mop_a.init_protein_RD(gamess_log[0].c_str(),"mepFukui",40,10,ref_coord,0,pdb_files[0].c_str(),true,"BD","gamess"); // mopac lmo
	mop_a.init_protein_RD(gamess_log[2].c_str(),"mepFukui",40,10,ref_coord,0,pdb_files[0].c_str(),true,"BD","gamess"); // mopac lmo
	//mop_a.init_protein_RD() // mopac uhf
	//mop_a.init_protein_RD() // mopac gamess rhf
	//mop_a.init_protein_RD() // mopac gamess uhf
	//mop_a.init_protein_RD() // mopac gamess dft
	//mop_a.init_protein_RD() // mopac gamess dftb ( sem grid )1140
	
}

test_p::~test_p(){}
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////END OF CLASS//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////1105  