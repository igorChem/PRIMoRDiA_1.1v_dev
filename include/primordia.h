#ifndef PRIMORDIA 
#define PRIMORDIA

//C++ headers files
#include <iostream>
#include <string>
#include <fstream>
#include <vector>


//foward declarations
class Imolecule;
class Icube;
class global_rd;
class local_rd;
class local_rd_cnd;
//============================================================================================================
/**
 * These class wraps all the functionalities of this libraries and make available in different ways of initializations
 * in respect of the type of the data passed, the type of approxiamtion calculus and representation of the reactivity
 * descriptors.
 * @class primordia
 * @author Igor Barden Grillo
 * @date 11/04/18
 * @file primordia.h
 * @brief Interface class for calculate the reactivity descriptors based on conceptual DFT.
 */
class primordia {
	public:
		/**
		 * String indicating the name of object.
		 */
		std::string name;
		
		//------------------------------------------------------------------
		/**
		 * Unsigned integer with the band size 
		 */
		unsigned int band; 
		
		//-------------------------------------------------------------------
		/**
		 *
		 */
		unsigned int gridn;
		
		//------------------------------------------------------------------------------
		/**
		 * std::string with name of the program
		 */	
		std::string program;
		
		//-------------------------------------------------------------------------------
		/**
		 * global_rd object smart pointer to store the global reactivity descriptors calculated. 
		 */
		std::unique_ptr<global_rd> grd;
		
		//--------------------------------------------------------------------------------
		/**
		 * local_rd object smart pointer to store the global reactivity descriptors calculated
		 * with volumetric representation.
		 */
		std::unique_ptr<local_rd> lrdVol;
		
		//--------------------------------------------------------------------------------
		/**
		 * local_rd object smart pointer to store the global reactivity descriptors calculated
		 * with condensed to atoms values representation.
		 */
		std::unique_ptr<local_rd_cnd> lrdCnd;
		
		//--------------------------------------------------------------------------------
		/**
		 * Instantiates an emptry object and initializes the memeber data to default values.
		 * @brief Default constructor of the class.
		 */
		primordia();
		
		//---------------------------------------------------------------------
		/**
		 * @brief Copy constructor for the class.
		 * @param primordia constant reference to the object to be copied.
		 */
		primordia(const primordia& pr_rhs);
		
		//---------------------------------------------------------------------
		/**
		 * @brief Assigment operator overloading for the class. 
		 * @param primordia constant reference to the object to be copied.
		 * @return primordia reference object.
		 */
		primordia& operator=(const primordia& pr_rhs);
		
		//---------------------------------------------------------------------
		/**
		 * @brief Move constructor for the class.
		 * @param primordia constant reference to the object to be moved.
		 */		
		primordia(primordia&& pr_rhs) noexcept;
		
		//---------------------------------------------------------------------
		/**
		 * @brief Move assigment operator overloading for the class. 
		 * @param primordia constant reference to the object to be moved.
		 * @return primordia reference object. 
		 */
		primordia& operator=(primordia&& pr_rhs) noexcept;
		
		//---------------------------------------------------------------------
		/**
		 * @brief Subtraction operator overloading for comparsion between methods for
		 * RD calculations.
		 * @param primordia constant reference to the object to the current object to be subtracted.
		 */
		friend primordia operator-(const primordia& pr_lhs,const primordia& pr_rhs);
		
		//-----------------------------------------------------------------------------------
		/**
		 * @brief Member function for Koopman approximation to get local RD from output of QM program.
		 * @param Char pointer to the name of the file for neutro charge state of QM computation out.
		 * @param Integer value with the grid resolution value.
		 * @param String indicating the QM program to output reading.
		 */
		void init_FOA(const char* file_neutro,int gridN,std::string loc_hard,bool mep, std::string Program);
		
		//------------------------------------------------------------------------------------
		/**
		 * @brief Member function for finite difference approximation to get local RD from output of QM program.
		 * @param Char pointer to the name of the file for neutro charge state of QM computation out.
		 * @param Char pointer to the name of the file for cation charge state of QM computation out.
		 * @param Char pointer to the name of the file for anion charge state of QM computation out.
		 * @param Integer value with the grid resolution value.
		 * @param String indicating the QM program to output reading.
		 */
		void init_FD(const char* file_neutro,const char* file_cation,const char* file_anion, int grdN, int charge,bool mep,std::string loc_hard, std::string Program);

		//-------------------------------------------------------------------------------------
		/**
		 * @brief Member function to get local RD for proteins/biomolecules
		 * @param Bool to indicate if the local hardness to be indicated 
		 * @param Integer with the grid number to get volumetric local RD
		 * @param Integer with the number of frontier molecular orbitals to be considered in the local RD calcualtions
		 * @param String indicating the QM program to output reading.
		 */				
		void init_protein_RD(const char* file_neutro,std::string locHardness,int grdN,int bandgap,double* ref_atom,int size,const char* pdb, bool mep , std::string bt, std::string Program);
		
		//---------------------------------------------------------------------
		/**
		 * Member function made to be used internally in the devolpment of another programs
		 * that will use this library.
		 * @brief Instatiates object from quantum chemical information loaded from Imolecule to
		 * calculate the reactivity descriptors using Koopamn approximation.
		 * @param Imolecule reference to object with molecular information.
		 * @param Integer with grid finess to calculate the reactivity maps.
		 */
		void init_QS_KA(Imolecule& mol, int gridN);
		
		//---------------------------------------------------------------------
		/**
		 * Member function made to be used internally in the devolpment of another programs
		 * that will use this library. Calculates the reactivity descriptors using finite difference
		 * approximations. Output global and local reactivity descriptors with both volumetric
		 * and condensed representations.
		 * @brief Instatiates object from quantum chemical information loaded from Imolecules to 
		 * calculate the reactivity descriptors using finite differences approximation.
		 * @param Imolecule reference object with neutro charged state molecular information.
		 * @param Imolecule reference object with cation charged state molecular information.
		 * @param Imolecule reference object with anion charged state molecular information.
		 * @param Integer value indicating the charge difference between the charged states.
		 * @param Integer with grid finess to calculate the reactivity maps.
		 */
		void init_QS_FD( Imolecule& mol1, Imolecule& mol2, Imolecule& mol3, int charge ,int gridN);
		
		//-------------------------------------------------------------------------------------
		/**
		 * @brief Function to write pymol scripts to automatize the descriptors vizualizations
		 */
		void write_pymol(const std::string& pdb, bool fixed);
		
		
		//----------------------------------------------------------------------
		/**
		 * @brief Destructor for the class.
		 */
		~primordia();
};

#endif