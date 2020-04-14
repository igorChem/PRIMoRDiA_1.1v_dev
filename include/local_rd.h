// header file for local_rd class
// local_rd.h

#ifndef LOCAL_RD
#define LOCAL_RD

//include statements from c++ library
#include <iostream>
#include <string> 
#include <vector>
//include statements from PRIMORDiA-libs
#include "../include/common.h"

//foward declarations 
class Iatom;
class Imolecule;
class Icube;
class global_rd;
//=================================================================================================
/**
 * This class is meant to represent a collection of local reactivity descriptors for a given molecular system,
 * in volumetric representation hold by Icube objects or in condensed form hold by STL vectors of doubles. 
 * The constructors of the current class instantiates objects passing necessary information to make the calculation
 * ans reporting of the local RD using the Koopman approximation or finite differences method.
 * @class local_rd
 * @author barden
 * @date 20/03/18
 * @file local_rd.h
 * @brief Class to calculate local reactivity descriptors from scalar fields and atomic charges from 
 * quantum chemistry info provided in the programs outputs. 
 */
class local_rd {
	public:
		/**
		 * @brief String with the system name.
		 */
		std::string name;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Bool to indicate if the calculations will be done with finite differences method.
		 */
		bool finite_diff;
		
		/**
		 * @brief Bool to indicate if local hardness was calculated.
		 */
		bool locHardness;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Integer value with the charge difference between the anion and cations state.
		 */	
		int charge;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Icube pointer to store scalar values of HOMO density of the system.
		 */			
		Icube homo;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Icube pointer to store scalar values of LUMO density of the system.
		 */			
		Icube lumo; 
		
		//----------------------------------------------------------------------------------------
		/**
		 *	@brief Icube pointer to store scalar values of total density of the system.
		 */			
		Icube elec_dens;
		
		//----------------------------------------------------------------------------------------
		/**
		 *	@brief Icube pointer to store scalar values of cation charge state density of the system.
		 */		
		Icube cation;
		
		//----------------------------------------------------------------------------------------
		/**
		 *	@brief Icube pointer to store scalar values of anion charge state density of the system.
		 */			
		Icube anion;
		
		//----------------------------------------------------------------------------------------
		/**
		 *	@brief Icube pointer to store the calculated scalar values of the electrophilicic succetibility 
		 * 	attack local descriptor. 
		 */			
		Icube EAS;
		
		//----------------------------------------------------------------------------------------
		/**
		 * 	@brief Icube pointer to store the calculated scalar values of the nucleophilic succetibility 
		 * 	attack local descriptor.
		 */			
		Icube NAS;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Icube pointer to store the calculated scalar values of the radical succetibility 
		 * attack local descriptor.
		 */			
		Icube RAS;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Icube pointer to store the calculated scalar values of the dual Fukui local descriptor.
		 */	 
		Icube Dual; 
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Icube pointer to store the calculated scalar values of the local chemical hardness 
		 * local descriptor.
		 */			
		Icube Hardness;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Icube pointer to store the calculated scalar values of the local chemical hardness 
		 * local descriptor.
		 */					
		Icube Softness_Dual;
		
		Icube Hyper_Softness;
		
		//----------------------------------------------------------------------------------------
		/**
		 * @brief Icube pointer to store the calculated scalar values of the local softness dual 
		 * local descriptor.
		 */	
		Icube multifilic;
		
		//---------------------------------------------------------------------------------------
		/**
		 * Instantiates a empty object to be filled later and initializes the member data to default values.
		 * @brief Deafult constructor for local_rd class.
		 */
		local_rd();
		
		
		//----------------------------------------------------------------------------------------
		/**
		 * Instantiates local_rd object with the necessary info to calculate volumetric representation for local 
		 * reactivity descriptors using the koopman approximation method from frontier molecular orbital densities
		 * loaded in Icube objects. 
		 * @brief Instantiates local_rd object with the necessary info to calculate local reactivity descriptors
		 * in its volumetric representation version using the Koopman approximation.
		 * @param Icube constant reference to object with volumetric HOMO scalar values.
		 * @param Icube constant reference to object with volumetric LUMO scalar values.	
		 */
		local_rd(const Icube& HOmo, const Icube& LUmo);
		
		/**
		*
		*/
		local_rd(const Icube& elec_dens, const Icube& HOmo, const Icube& LUmo);
		
		//----------------------------------------------------------------------------------------
		/**
		 * Instantiates the a object of the class passing the Icube objects with the scalar field values for the
		 * electron density of three different charged states to make the calculations in finite differences aproximation
		 * and return it in its volumetric form writting cube files with the reactivity maps. 
		 * @brief Instantiates the object with the necessary info to calculate local reactivity descriptors
		 * in its volumetric representation version using the finite difference approximation.
		 * @param Icube constant reference to object with electronic density of the neutro charge state for the system.
		 * @param Icube constant reference to object with electronic density of the cation charge state for the system.
		 * @param Icube constant reference to object with electronic density of the anion charge state for the system.
		 * @param Integer to value with the charge difference between the anion and cation charge states.
		 */
		local_rd(const Icube& elecDens, const Icube& cationDens, const Icube& anionDens,int chg);
		
		//----------------------------------------------------------------------------------------
		/**
		 * Instantiates a object as a copy of another local_rd instance.
		 * @brief Copy Constructor for the local_rd class.
		 * @param local_rd reference constant object. 
		 */
		local_rd(const local_rd& lrd_rhs);
		
		//----------------------------------------------------------------------------------------
		/**
		 * Overwrites this memeber data with the member data of the object at right hand side of 
		 * assigment operator. 
		 * @brief Assigment operator overloading for the class.
		 * @param local_rd constant reference to a constant object.
		 * @return local_rd reference.
		 */
		local_rd& operator=(const local_rd& lrd_rhs);
		
		//----------------------------------------------------------------------------------------
		/**
		 * Instantiates a object with data stolen ownership from another local_rd instance
		 * to moved in the current object.
		 * @brief Move Constructor for the local_rd class.
		 * @param local_rd reference object to be moved. 
		 */
		local_rd(local_rd&& lrd_rhs) noexcept;		
		
		//----------------------------------------------------------------------------------------
		/**
		 * Overwrites this memeber data with the member data of the object at right hand side of 
		 * assigment operator. 
		 * @brief Move assigment operator overloading for the class.
		 * @param local_rd constant reference to object to be moved from.
		 * @return local_rd reference.
		 */
		local_rd& operator=(local_rd&& lrd_rhs) noexcept;
		
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
		friend local_rd operator-(const local_rd& lrd_lhs,const local_rd& lrd_rhs);
		
		//----------------------------------------------------------------------------------------
		/**
		 * When the volumetric representation is used ( condensed = false ): perform the Icube 
		 * objects arithmetics to calculate the Icube objects of the reactivity descriptors 
		 * stored as meber in this object. 
		 * @brief Calculates the basic local descriptors that are the the fukui function 
		 * and stores in the meber varibales of the object.
		 * @return None.
		 */	
		void calculate_fukui();
		
		/**
		 * @brief 
		 * @param mol
		 */
		void calculate_fukui_EW(const Imolecule& mol);
		
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
		
		/**
		 * @brief 
		 * @param grd
		 * @param method
		 */
		void calculate_hardness(const global_rd& grd, std::string method);
		
		//----------------------------------------------------------------------------------------
		/**
		 * Member function to write the report local reactivtiy descriptors.
		 * @brief Creates and writes the files to cube files. 
		 * @return none.
		 */
		void write_LRD();
		
		//----------------------------------------------------------------------------------------
		/**
		 * Deallocate the memomory dinamically allocated in the constructors.
		 * @brief Destructor of the class.
		 */
		~local_rd();
};

#endif
