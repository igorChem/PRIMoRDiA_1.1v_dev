// header file for local_rd_cnd class
// local_rd_cnd.h
//-------------------------------------------------------------------------------
#ifndef LOCAL_RD_CND
#define LOCAL_RD_CND
//----------------------------------------------------------------------------------
#include <iostream>
#include <string> 
#include <vector>
#include <cmath>
//include statements from PRIMORDiA-libs
#include "../include/common.h"
//--------------------------------------------------------------------------------------
class Iatom;
class Imolecule;
class QMdriver;
class global_rd;
class Iprotein;

//===============================================================
class local_rd_cnd{
	public:
		/**
		 * @brief String with the system name.
		 */
		std::string name;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Bool to indicate if the condensed electron density must be calculated; 
		 */
		bool ed;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Bool to indicate if the condensed mep will be made for pdb vizualization
		 */
		bool mep;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Bool to indicate if the calculations will be done with finite differences method.
		 */		
		bool finite_diff;
			
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Integer value with the charge difference between the anion and cations state.
		 */
		int charge;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Imolecule smart pointer for molecular information.
		 */	
		Imolecule molecule;
		
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Icube pointer to store the calculated scalar values of the electrophilicic succetibility
		 * attack local descriptor.
		 */			
		std::vector<double> EAS;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the  nucleophilic succetibility
		 * attack condensed to atoms local descriptor.
		 */			
		std::vector<double> NAS;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the radical succetibility
		 * attack condensed to atoms local descriptor.
		 */			
		std::vector<double> RAS;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the radical succetibility
		 * condensed to atoms local descriptor.
		 */			
		std::vector<double> dual;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the local dual softness 
		 * condensed to atoms local descriptor.
		 */					 
		std::vector<double> Softness_Dual;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the local  hyper dual softness 
		 * condensed to atoms local descriptor.
		 */
		std::vector<double> Hyper_softness;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the local dual softness 
		 * condensed to atoms local descriptor.
		 */
		std::vector<double> hardness_A;
		std::vector<double> hardness_B;
		std::vector<double> hardness_C;
		std::vector<double> hardness_D;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the multifilic condensed 
		 * to atoms local descriptor.
		 */	
		std::vector<double> multifilic;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the electrophilicity condensed 
		 * to atoms local descriptor.
		 */	
		std::vector<double> electrophilicity;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the relative electrophilic 
		 * succetibility condensed to atoms local descriptor.
		 */
		std::vector<double> relative_EAS;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles to store the calculated scalar values of the relative nucleophilic 
		 * succetibility condensed to atoms local descriptor.
		 */
		std::vector<double> relative_NAS;
		
		//---------------------------------------------------------------------------------------
		/**
		 *
		 */
		std::vector<double> fukushima;
		
		//---------------------------------------------------------------------------------------
		/**
		 *
		 */
		std::vector<double> electron_density;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Default constructor.
		 */	
		local_rd_cnd();
				
		//---------------------------------------------------------------------------------------
		/**
		 * @brief One-arg constructor.
		 * @param Imolecule constant reference 
		 */	
		local_rd_cnd(Imolecule&& mol) noexcept;
		
		//----------------------------------------------------------------------------------------
		/**
		 * Instantiates the object with the necessary info to calculate condensed to atoms local 
		 * reactivity descriptors using the finite differences methods from atom charges from three different charged
		 * states loaded in Imolecules objects passed in the argument list. 
		 * @brief Instantiates local_rd object to calculated the RD with finite differences method from Imolecules
		 * objects.
		 * @param Imolecule constant reference to object with the info of neutro charge state of the system.
		 * @param Imolecule constant reference to object with the info of cation charge state of the system.
		 * @param Imolecule constant reference to object with the info of anion charge state of the system.  
		 */
		local_rd_cnd(const Imolecule& mol_neut,const Imolecule& mol_cation,const Imolecule& mol_anion);
	
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Copy constructor.
		 * @param local_rd_cnd object to be copied 
		 */	
		local_rd_cnd(const local_rd_cnd& lrd_rhs);
		
		//----------------------------------------------------------------------------------------
		/**
		 * Overwrites this memeber data with the member data of the object at right hand side of 
		 * assigment operator. 
		 * @brief Assigment operator overloading for the class.
		 * @param local_rd constant reference to a constant object.
		 * @return local_rd reference.
		 */
		local_rd_cnd& operator=(const local_rd_cnd& lrd_rhs);
		
		//----------------------------------------------------------------------------------------
		/**
		 * Instantiates a object with data stolen ownership from another local_rd instance
		 * to moved in the current object.
		 * @brief Move Constructor for the local_rd class.
		 * @param local_rd reference object to be moved. 
		 */
		local_rd_cnd(local_rd_cnd&& lrd_rhs) noexcept;
		
		//----------------------------------------------------------------------------------------
		/**
		 * Overwrites this memeber data with the member data of the object at right hand side of 
		 * assigment operator. 
		 * @brief Move assigment operator overloading for the class.
		 * @param local_rd constant reference to object to be moved from.
		 * @return local_rd reference.
		 */
		local_rd_cnd& operator=(local_rd_cnd&& lrd_rhs) noexcept;
		
		//----------------------------------------------------------------------------------------
		/**
		 * Friend function of overloaded subtraction operator to return a object with the numeric 
		 * member data as a difference to the two others objects. This method is meant to be a way 
		 * to estimate the difference between two quantum chemical methods to get the local reactivity
		 * descriptors.
		 * @brief Arithmetic subtraction operator overloading for the class. 
		 * @param local_rd constant reference at left hand side of minus sign operator.
		 * @param local_rd constant reference at right hand side of minus sign operator.
		 * @return local_rd reference object with the subtracted local RD
		 */
		friend local_rd_cnd operator-(const local_rd_cnd& lrd_lhs,const local_rd_cnd& lrd_rhs);
		
		//----------------------------------------------------------------------------------------
		/**
		 * When the volumetric representation is used ( condensed = false ): perform the Icube 
		 * objects arithmetics to calculate the Icube objects of the reactivity descriptors 
		 * stored as meber in this object. 
		 * @brief Calculates the basic local descriptors that are the the fukui function 
		 * and stores in the meber varibales of the object.
		 * @return None.
		 */	
		void calculate_Fukui();
		
		/**
		 * @brief 
		 * @param band
		 */
		void calculate_Fukui_band(int band);
		
		/**
		 * @brief 
		 * @param band
		 */
		void calculate_Fukui_EW(int band);
				
		
		//----------------------------------------------------------------------------------------
		/**
		 * Member function to calculate more specific local descriptor.
		 * @brief Calculates local reactivity descriptor that are a combination of fukui functions 
		 * with global reactivity descriptor
		 * @param global_rd reference to object with global reactivity descriptors to combine with 
		 * Fukui functions.
		 * @return none.
		 */
		void calculate_RD(const global_rd& grd);
		
		//----------------------------------------------------------------------------------------
		/**
		 * Member function to calculate local Hardness.
		 * @brief Calculates the approximation to local hardness. 
		 * @param global_rd reference to object with global reactivity descriptors to combine with 
		 * Fukui functions.
		 * @return none.
		 */
		void calculate_Hardness(const global_rd& grd);
		//----------------------------------------------------------------------------------------
		/**
		 * @brief 
		 * @param prot
		 */
		//----------------------------------------------------------------------------------------
		void write_rd_protein(const Iprotein& prot);
		/**
		 * @brief 
		 * @param prot
		 */
		 //----------------------------------------------------------------------------------------
		void write_rd_protein_pdb(const Iprotein& protein);
		/**
		 * @brief 
		 * @param prot
		 */
		//----------------------------------------------------------------------------------------
		/**
		 * @brief 
		 * @param prot
		 */
		void write_rd_protein_reaction(const Iprotein& prot);

		//----------------------------------------------------------------------------------------
		/**
		 * Member function to write the report local reactivtiy descriptors.
		 * @brief Creates and writes the files to cube files. 
		 * @return none.
		 */
		void write_LRD();
		
		//------------------------------------------------------------------------------------------
		/**
		 * 
		 */
		~local_rd_cnd(){};
	
};



#endif 