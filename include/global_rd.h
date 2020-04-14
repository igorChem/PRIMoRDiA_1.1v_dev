// header file for global_rd.cpp source file 
// class for global reactivity descriptors 

#ifndef GLOBAL_RD
#define GLOBAL_RD
//include C++ headers
#include <iostream>
#include <string> 
#include <fstream>
//include our library header 
#include "../include/common.h"

//===================================================================================
/**
 * Class to represent objects that hold necessary information to and calculate the global reactivity
 * descriptors following the Conceptual Density Functional Theory.  
 * Contatc: barden.igor@gmail.com
 * @class global_rd
 * @author barden
 * @date 20/03/18
 * @file global_rd.h
 * @brief Instantiate objects with quantum chemistry info extracted from other programs calculation to get the
 * global reactivity descriptors for the system. * 
 */
class global_rd{
	public:
		/**
		* @brief String with the system's name.
		*/
		std::string name;
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with energy of Highest energy Occupied Molecular Orbital extracted from a 
		 * quantum chemistry calculation output.
		 */		
		double homo_en; 
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with energy of Lowest energy Unoccupied Molecular Orbital extracted from a 
		 * quantum chemistry calculation output.
		 */	
		double lumo_en;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the total energy extracted from a 
		 * quantum chemistry calculation output.
		 */	
		double energ_tot;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the total energy from the cation charge state extracted from a 
		 * quantum chemistry calculation output
		 */	
		double energ_cat;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the total energy from the anion charge state extracted from a 
		 * quantum chemistry calculation output.
		 */	
		double energ_an;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the chemical hardness calculated by this class methods.
		 */	
		double hardness;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the electronic chemical potential calculated by this class methods.
		 */	
		double chemical_pot;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the chemical softness calculated by this class methods.
		 */	
		double softness;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the value of Ionization Potential.
		 */	
		double IP;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the value of Electron Affinity.
		 */	
		double EA;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the value of Electrophilicity calculated by the mehotds of the class.
		 */	
		double Electrophilicity;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the value of the difference between LUMO and HOMO energies.
		 */	
		double gap;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the tranferred electrons to one system to another 
		 */	
		double deltaN;
		
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Double value with the maximum number of electron that can be transferred to the system.
		 */	
		double nMax;
		
		double HOF; 
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Bool indicating if the calculation is done with the koopman approximation.
		 */
		bool KA;
		 
		//---------------------------------------------------------------------------------------
		/**
		 * @brief Bool indicating if the calculation is done with the finite differences approximation.
		 */
		bool DF;
		
		
		//---------------------------------------------------------------------------------------
		/**
		 * Instantiate a empty object to be filled later and initializes the member data to default values.
		 * @brief Default constructor. 
		 */
		global_rd();
		
		//---------------------------------------------------------------------------------------
		/**
		 * Instantiate the object with necessary info from a Imolecule object with quamtum chemical information
		 * extracted from program output to calculate the global reactivity descriptors using the 
		 * Koopman approximation.
		 * @brief Instatiates global_rd object from Imolecule object to use the Koopman approximaion. 
		 * @param Imolecule const reference with quantum molecular information loaded.
		 */
		global_rd(const Imolecule& mol);
		
		//---------------------------------------------------------------------------------------
		/**
		 * Instantiates global_rd from three Imolecule objects for the three different charged
		 * states with quantum chemical information extracted from program outputs to calculate 
		 * the RD with finite difference approximation.
		 * @brief Constructor to calculate global RD with finite differences approx or koopman.
		 * @param Imolecule const reference object to the neutro state molecule.
		 * @param Imolecule const reference object to the cation state molecule.
		 * @param Imolecule const reference object to the anion state molecule.
		 */
		global_rd(const Imolecule& mol_neutro,const Imolecule& mol_cation,const Imolecule& mol_anion);
		
		//---------------------------------------------------------------------------------------
		/**
		 * Instantiates a global_rd object with member varibales copied from another object.
		 * @brief Copy constructor of the global_rd class.
		 * @param global_rd object constant reference to be copied.
		 */
		global_rd(const global_rd& rd_rhs);
		
		//---------------------------------------------------------------------------------------
		/**
		 * Return a class object from a copy to another object provided in the right hand side of
		 * the assigment opérator.
		 * @brief Assigment operator overloading.
		 * @param globa_rd object constant reference.
		 * @return global_rd reference to this object. 
		 */
		global_rd& operator=(const global_rd& rd_rhs);
		
		//---------------------------------------------------------------------------------------
		/**
		 * Instantiates a global_rd object with member data stolen ownership from another object.
		 * @brief Move constructor of the global_rd class.
		 * @param global_rd object constant reference to be copied.
		 */
		global_rd(global_rd&& rd_rhs) noexcept;
		
		//---------------------------------------------------------------------------------------
		/**
		 * Return a class object from a memeber data stolen ownership to another object provided 
		 * in the right hand side of the assigment opérator.
		 * @brief Assigment operator overloading.
		 * @param globa_rd object constant reference.
		 * @return global_rd reference to this object. 
		 */
		global_rd& operator=(global_rd&& rd_rhs) noexcept;
		
		//---------------------------------------------------------------------------------------
		/** 
		 * Friend function for subtraction operator overloading for global reactivity descriptors
		 * objects member data. The goal with this method is to return a object with the difference of global
		 * descriptors between to methods used to calculated for a given molecular system.
		 * @brief Arithmetic subtraction operator overloadingg for the class member data.
		 * @param global_rd constant reference where the RD values will be subtracted 
		 */
		friend global_rd operator-(const global_rd& lhs_grd,const global_rd& rhs_grd);
		
		//---------------------------------------------------------------------------------------
		/**
		 * Member function to calculate the global reactivitiy descriptors using a given type of
		 * approximation that were triggered by the construction used to instatiate this object.
		 * The present method modify the other member data that were declare to hold the values of 
		 * global reactivity descriptors.
		 * @brief Used the member data to calculate the global reactivity descriptors.
		 * @return None.
		 */
		void calculate_rd();
		
		//---------------------------------------------------------------------------------------
		/** 
		 * Prints in the console the values calulated by calculate_rd that are stored in the
		 * memeber varibales.
		 * @brief Member function to print in the screen the values of the descriptors.
		 * @return None.
		 */
		void print_rd();
		
		//---------------------------------------------------------------------------------------
		/**
		 * Creates and writes in a file the global reactivity descriptors calulated by the class. 
		 * @brief Member function to output the global reactivity descriptors to a file.
		 * @return none.
		 */
		void write_rd();
		
		//---------------------------------------------------------------------------------------
		/** 
		 * Destructor of the class.
		 */
		~global_rd();
};


#endif 