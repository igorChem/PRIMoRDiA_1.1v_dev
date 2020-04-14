#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

#include "../include/log_class.h"
#include "../include/Imolecule.h"
#include "../include/Icube.h"
#include "../include/global_rd.h"
#include "../include/local_rd.h"
#include "../include/local_rd_cnd.h"
#include "../include/primordia.h"
#include "../include/compPrimordia.h"

using std::abs;
/*************************************************************/
compPrimordia::compPrimordia():
	result_pr( new primordia() )  {
}
/*************************************************************/
compPrimordia::compPrimordia(const primordia& pr1, 
							 const primordia& pr2):
	result_pr( new primordia( pr1 - pr2 ) )       {
	/*	
	if ( result_pr->glob ) {
		result_pr->grd->name = pr1.name + "_diff_grd_" + pr2.name;
		result_pr->grd->write_rd();
	}
	if ( result_pr->local ){
		if ( result_pr->condensed ){
			result_pr->lrdCnd->name = pr1.name + "_diff_cnd_" + pr2.name;
			result_pr->lrdCnd->write_LRD();
			unsigned int NOF = result_pr->lrdCnd->molecule->num_of_atoms;
			
			double RMSE = 0.0;
			for(unsigned int i=0;i<NOF;i++)	{ RMSE += abs(result_pr->lrdCnd->fukui_elec_cnd[i]); }
			condensedRMSE.push_back(RMSE);
			
			RMSE = 0.0;
			for(unsigned int i=0;i<NOF;i++){ RMSE  += abs(result_pr->lrdCnd->fukui_nuc_cnd[i]); }
			condensedRMSE.push_back(RMSE);
			
			RMSE = 0.0;			
			for(unsigned int i=0;i<NOF;i++){ RMSE  += abs(result_pr->lrdCnd->fukui_radical_cnd[i]); }
			condensedRMSE.push_back(RMSE);
			
			RMSE = 0.0;
			for(unsigned int i=0;i<NOF;i++){ RMSE  += abs(result_pr->lrdCnd->deltFukui_cnd[i]); }
			condensedRMSE.push_back(RMSE);
			
			RMSE = 0.0;
			for(unsigned int i=0;i<NOF;i++)	{ RMSE += abs(result_pr->lrdCnd->loc_softness_dual_cnd[i]); }
			condensedRMSE.push_back(RMSE);
			
			RMSE = 0.0;
			for(unsigned int i=0;i<NOF;i++)	{ RMSE += abs(result_pr->lrdCnd->loc_hardness_cnd[i]); }
			condensedRMSE.push_back(RMSE);
			
			RMSE = 0.0;
			for(unsigned int i=0;i<NOF;i++)	{ RMSE += abs(result_pr->lrdCnd->multifilic_cnd[i]); }
			condensedRMSE.push_back(RMSE);
			
			RMSE = 0.0;
			for(unsigned int i=0;i<NOF;i++)	{ RMSE += abs(result_pr->lrdCnd->electrophilicity_cnd[i]); }
			condensedRMSE.push_back(RMSE);
			
			RMSE = 0.0;
			for(unsigned int i=0;i<NOF;i++)	{ RMSE += abs(result_pr->lrdCnd->rel_suc_elec_cnd[i]); }
			condensedRMSE.push_back(RMSE);
			
			RMSE = 0.0;
			for(unsigned int i=0;i<NOF;i++)	{ RMSE += abs(result_pr->lrdCnd->rel_suc_nuc_cnd[i]); } 
			condensedRMSE.push_back(RMSE);
		}		
		else if ( result_pr->volumetric ){ 
			double integralErr = 0.0;
			result_pr->lrdVol->name = pr1.lrdVol->name + "_diff_Vol_" + pr2.lrdVol->name;
			result_pr->lrdVol->write_LRD();
			integralErr = pr1.lrdVol->fukui_elec.similarity_index(pr2.lrdVol->fukui_elec,"default");
			integralError.push_back(integralErr);
			integralErr = pr1.lrdVol->fukui_nuc.similarity_index(pr2.lrdVol->fukui_nuc,"default");
			integralError.push_back(integralErr);
			integralErr = pr1.lrdVol->fukui_radical.similarity_index(pr2.lrdVol->fukui_radical,"default");
			integralError.push_back(integralErr);
			integralErr = pr1.lrdVol->deltFukui.similarity_index(pr2.lrdVol->deltFukui,"default");
			integralError.push_back(integralErr);
			integralErr = pr1.lrdVol->loc_hardness.similarity_index(pr2.lrdVol->loc_hardness,"default");
			integralError.push_back(integralErr);
			integralErr = pr1.lrdVol->loc_softness_dual.similarity_index(pr2.lrdVol->loc_softness_dual,"default");
			integralError.push_back(integralErr);
			integralErr = pr1.lrdVol->multifilic.similarity_index(pr2.lrdVol->multifilic,"default");
			integralError.push_back(integralErr);
		}
	}
	 */
}
/*************************************************************/
compPrimordia::compPrimordia(const compPrimordia& rhs_compPri):
	result_pr( new primordia(*rhs_compPri.result_pr) )        {
}
/*************************************************************/
compPrimordia& compPrimordia::operator=(const compPrimordia& rhs_compPri){	
	if ( this != &rhs_compPri ) {
		*result_pr = *rhs_compPri.result_pr;
	}
	return *this;
}
/*************************************************************/
compPrimordia::compPrimordia(compPrimordia&& rhs_compPri) noexcept:
	result_pr( std::move(rhs_compPri.result_pr) )        {
}
/*************************************************************/
compPrimordia& compPrimordia::operator=(compPrimordia&& rhs_compPri) noexcept{	
	if ( this != &rhs_compPri ) {
		result_pr = std::move(rhs_compPri.result_pr);
	}
	return *this;	
}
/*************************************************************/
void compPrimordia::write_report(){
	/*
	std::ofstream diff_repfile;
	if ( result_pr->condensed ) diff_repfile.open( result_pr->lrdCnd->name.c_str() );
	else if ( result_pr->volumetric ) diff_repfile.open( result_pr->lrdVol->name.c_str() );
	diff_repfile << " Report with the values of the calculated indices for differences in local reactivity descriptors \n";
	diff_repfile.precision(6);
	diff_repfile << std::scientific;
	
	if ( result_pr->condensed ){
		diff_repfile << "The local reactivity descriptors are condensed to atoms \n" 
					 << "The following values are the Sum of Root Mean Square difference betewen the atoms \n"
					 << "Fukui Electrophilic Suc: "     << condensedRMSE[0] << " \n"
					 << "Fukui Nucleophilic Suc: "      << condensedRMSE[1] << " \n"
					 << "Fukui Radical Suc: "           << condensedRMSE[2] << " \n"
					 << "Dual Fukui Suc: "              << condensedRMSE[3] << " \n"
					 << "Local Softness dual: "         << condensedRMSE[4] << " \n"
					 << "Multiphilic dual descriptor: " << condensedRMSE[5] << " \n"
					 << "Electrophilicity descriptor: " << condensedRMSE[6] << " \n"
					 << "Relatie Fukui Elec Suc: "      << condensedRMSE[7] << " \n"
					 << "Relatie Fukui Nuc Suc: "       << condensedRMSE[8] << " \n";
	}
	else if ( result_pr->volumetric ){
		diff_repfile << "The local reactivity descriptors 3D grid representation \n"
					 << "The following values are the difference of integral of the grids in all calculated space \n"
					 << "Fukui Electrophilic Suc: "     << integralError[0] << " \n"
					 << "Fukui Nucleophilic Suc: "      << integralError[1] << " \n"
					 << "Fukui Radical Suc: "           << integralError[2] << " \n"
					 << "Dual Fukui Suc: "              << integralError[3] << " \n"
					 << "Local Softness dual: "         << integralError[4] << " \n"
					 << "Multiphilic dual descriptor: " << integralError[5] << " \n"
					 << "Relatie Fukui Elec Suc: "      << integralError[6] << " \n"
					 << "Relatie Fukui Nuc Suc: "       << integralError[7] << " \n"
					 << "Local Hardness: "              << integralError[8] << " \n";
	}
	diff_repfile.close();
	 */
}
/*************************************************************/
compPrimordia::~compPrimordia(){}
