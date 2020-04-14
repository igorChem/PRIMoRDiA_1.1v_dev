#ifndef COMP_PRIMORDIA
#define COMP_PRIMORDIA

// C++ include files 
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <memory>
//===================================================================================
class primordia; // foward declaration
//===================================================================================
/**
 * Class to calculate the reactivity descriptors difference between two quantum chemical different methods.
 * Holds a primordia object with the difference between two primorida objects. 
 * Calculates the differences depending the type of local representations and write the report with the
 * indices summary and/or the differences maps and/the local condensed descriptors differences for each atom.
 * @class compPrimordia
 * @author Igor Barden Grillo
 * @date 28/04/18
 * @file compPrimordia.h
 * @brief Class to calculate the reactivity descriptors difference between two quantum chemical different methods
 */
class compPrimordia {
	public:
		/**
		 * @brief primordia smart pointer pointing the object with the reactivity descriptors information from the
		 * differences between other two.
		 */
		std::unique_ptr<primordia> result_pr;
		
		//---------------------------------------------------------------------
		/**
		 * @brief STL vector of dobles with the integrated value of the differences between scalar fields of the 
		 * reactivity descriptors.
		 */
		std::vector<double> integralError;
		
		//---------------------------------------------------------------------
		/**
		 * @brief STL vector of doubles with the Root mean square error of condensed to atom local reactivity descriptors 
		 */
		std::vector<double> condensedRMSE;
			
		//---------------------------------------------------------------------
		/**
		 * Instatiates an empty object and intilalizes the memeber data do default values.
		 * @brief Default constructor for the class. 
		 */
		compPrimordia();
		
		//---------------------------------------------------------------------
		/**
		 * @brief Instantiates an object and calculated the differences between two collections of reactivity
		 * descriptors.
		 * @param primordia constant reference to the first object in comparsion.
		 * @param primordia constant reference to the second object in comparsion.
		 */
		compPrimordia(const primordia& pr1,const primordia& pr2);
		
		//---------------------------------------------------------------------
		/**
		 * Instantiate an object of the class type initializing its memeber data with the info stored in another 
		 * object.
		 * @brief Copy constructor of the class.
		 * @param compPrimordia constante reference to the object to be copied.
		 */
		compPrimordia(const compPrimordia& rhs_compPri);
		
		//---------------------------------------------------------------------
		/**
		 * @brief Assigment operator overloading for the class. 
		 * @param compPrimordia object to be assigned.
		 * @return compPrimorida reference to this object.
		 */
		compPrimordia& operator=(const compPrimordia& rhs_compPri);
		
		//---------------------------------------------------------------------
		/**
		 * Instantiate an object of the class type initializing its memeber data with the info moved 
		 * from  another  object.
		 * @brief move constructor of the class.
		 * @param compPrimordia reference to the object to be moved.
		 */	
		compPrimordia( compPrimordia&& rhs_compPri) noexcept;
		
		//---------------------------------------------------------------------
		/**
		 * @brief Assigment operator overloading for the class. 
		 * @param rhs_compPri
		 * @return compPrimorida reference to this object.
		 */
		compPrimordia& operator=( compPrimordia&& rhs_compPri) noexcept;
		
		//---------------------------------------------------------------------
		/**
		 * @brief Output files with the member data with the summary of differences betrween the objects.
		 */
		void write_report();
		
		//---------------------------------------------------------------------
		/**
		 * @brief Destructor of the class.
		 */
		~compPrimordia();
	
};

#endif