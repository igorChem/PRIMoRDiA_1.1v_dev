// Icube.h 
// header file for the declarations ans documentation of the Icube class 
#ifndef CUBE
#define CUBE

#include <iostream>
#include <string>
#include <vector> 

#include "../include/common.h"

class Imolecule; // foward declaration
//===================================================================================
/**
 * This class is meant to represent an cube file object in the gaussian format that are accepted by 
 * the most of quantum chemical GUI programs. The Icube object holds a Imolecule object, thus molecular information,
 * with name and cooridnates of the atoms, molecular orbital number the the Icube object may represent, and an array 
 * with the scalar field data.
 * @class Icube
 * @author barden
 * @date 20/03/18
 * @file Icube.h
 * @brief A class to represent volumetric chemical data that are placed in grid files.
 * The present class methods works with these type of files, parsing and rewriting them.
 * @see Imolecule
 */ 
class Icube {
	public:
	    // member variables
		//--------------------------------------------------------------------
		/**
		 *  @brief String with cube's name.
		 */
		std::string name;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Bool to indicate if the scalar is of electron density of the a molecular system.
		 */
		bool elec_dens;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Integer with the number of molecular orbital represented in the scalar data.
		 * Is relevant when elecDens is setted false.
		 */
		unsigned int MOn;
		
		
		//--------------------------------------------------------------------
		/**
		 * @brief Number of voxels in the Icube object.
		 */
		unsigned int voxelN;
		
		//--------------------------------------------------------------------
		/**
		 * @brief String with text information that is placed in the beggining of 
		 * cube files.
		 */
		std::string header;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Array of double with the origin vertice (inferios coordinates) of the cube.
		 */
		double origin[3];
		
		//--------------------------------------------------------------------
		/**
		 * @brief Array od double with the lenght of the voxels sides.
		 */
		double gridsides[3];
		
		//--------------------------------------------------------------------
		/**
		 * @brief Array of integer with the number of voxels in each dimension.
		 */
		unsigned int grid[3];
		
		//--------------------------------------------------------------------
		/**
		 * @brief Imolecule object pointer with all important information about
		 * the molecular system. 
		 */
		Imolecule molecule;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief STL vector of doubles hold the values of the voxels. 
		 */
		std::vector<double> scalar;
		
		//-------------------------------------------------------------------
		/**
		 * Instantiates an empty object with default values initialized for the memebers.
		 * @brief Default constructor for the Icube class.
		 */
		Icube();
	
		//--------------------------------------------------------------------
		/**
		 * Instantiates an Icube object reading info from .cube file.
		 * Contructor with one-arg from .cube file.
		 * @brief Constructor that initializes the Icube object from a gaussian type cube file.
		 * @param name of the file to be readed.
		 */
		Icube(const char* file_nam);
		
		//--------------------------------------------------------------------
		/**
		 * Instatiate a Icube object from a copy of another object
		 * @brief Copy constructor for Icube class.
		 * @param Constant reference to a Icube object to be copied.  
		 */
		Icube(const Icube& rhs_cube); 
		
		//--------------------------------------------------------------------
		/**
		 * Creates and return a Icube object from a copy of another
		 * through the assgiment operator. 
		 * @brief Assigment operator overloading for Icube class.
		 * @param Icube const reference to Icube object to be copied.
		 * @return Icube reference to this.
		 */
		Icube& operator=(const Icube& rhs_cube);  
		
		//-------------------------------------------------------------------
		/**
		 * @brief Move constructor for the class.
		 * @param Icube rvalue reference to object to be moved. 
		 */
		Icube(Icube&& rhs_cube) noexcept;
		
		//-------------------------------------------------------------------
		/**
		 * @brief Move assigment operator overloading.
		 * @param Icube rvalue reference to object to be moved.
		 * @return Icube referemce to this.
		 */
		Icube& operator=(Icube&& rhs_cube) noexcept; 
		
		//--------------------------------------------------------------------
		/**
		 * Freind function to compare the dimensions and voxel number of two Icube objects to test if 
		 * mathematical operations can be performed between then.
		 * @brief Logical equal operator overloading. 
		 * @param Icube const reference object to compare at LHS side of the operator.
		 * @param Icube const reference object to compare at RHS side of the operator.
		 * @return Bool indicate if the dimensionas are equal. 
		 */
		friend bool operator==(const Icube& lhs_cube,const Icube& rhs_cube);
		
		//--------------------------------------------------------------------
		/**
		 * Friend subtraction operator overloading function to subtract the scalar values of an
		 * Icube object at left hand side of the minus operator from the 
		 * scalar values of an Icube to the right hand side of the operator.
		 * @brief Arithmetic subtraction operator overloading for the Icube class objects
		 * @param Icube const reference object to be subtracted at LHS side of the minus sign operator.
		 * @param Icube const reference object to subtract at RHS side of the minus operator.
		 * @return Icube object with the resultant scalar values.
		 */		 
		friend Icube operator-(const Icube& lhs_cube,const Icube& rhs_cube);
		
		//--------------------------------------------------------------------
		/**
		 * Friend addition operator overloading function to sum the scalar values of an
		 * Icube object at left hand side of the plus sign operator from the 
		 * scalar values of an Icube to the right hand side of the operator.
		 * @brief Arithmetic addition operator overloading for the class objects.
		 * @param Icube const reference object to be summed at LHS side of the plus sign operator.
		 * @param Icube const reference object to be summed at RHS side of the plus sign operator.
		 * @return Icube object with the resultant scalar values. 
		 */
		friend Icube operator+(const Icube& lhs_cube,const Icube& rhs_cube);
		
		//--------------------------------------------------------------------
		/**
		 * Friend product operator overloading function to multiply the scalar values of an
		 * Icube object at left hand side of the product operator from the 
		 * scalar values of an Icube to the right hand side of the operator.
		 * @brief Arithmetic product operator overloading for the class objects.
		 * @param Icube const reference object to be multiplied at LHS side of the product sign operator.
		 * @param Icube const reference object to be multiplied at RHS side of the product sign operator.
		 * @return Icube object with the resultant scalar values.
		 */
		friend Icube operator*(const Icube& lhs_cube, const Icube& rhs_cube);
		
		//--------------------------------------------------------------------
		/**
		 * Friend division operator overloading function to divide the scalar values of an
		 * Icube object at left hand side of the division operator from the 
		 * scalar values of an Icube to the right hand side of the operator.
		 * @brief Arithmetic division operator overloading for the class objects.
		 * @param Icube const reference object to be multiplied at LHS side of the division sign operator.
		 * @param Icube const reference object to be multiplied at RHS side of the division sign operator
		 * @return Icube object with the resultant scalar values.
		 */
		friend Icube operator/(const Icube& lhs_cube,const Icube& rhs_cube);
		
		//--------------------------------------------------------------------
		/**
		 * Friend subtraction operator overloading function to subtract the scalar values of an
		 * Icube object at left hand side of the minus sign operator from the double value passed to 
		 * function in the right hand side of the operator.  
		 * @brief Arithmetic subtraction operator overloading for the class objects with doubles values.
		 * @param Icube const reference object to be multiplied at LHS side of the minus sign operator.
		 * @param Double value to subtract the scalar values of the Icube object.
		 * @return Icube object with the resultant scalar values.
		 */
		friend Icube operator-(const Icube& lhs_cube,double value);
		
		//--------------------------------------------------------------------
		/**
		 * Friend addition operator overloading function to sum the scalar values of an
		 * Icube object at left hand side of the plus sign operator with the double value passed to 
		 * function in the right hand side of the operator. 
		 * @brief Arithmetic addition operator overloading for the class objects with doubles values. 
		 * @param Icube const reference object to be summed at LHS side of the plus sign operator.
		 * @param Double value to sum with the scalar values of the Icube object.
		 * @return Icube object with the resultant scalar values.
		 */
		friend Icube operator+(const Icube& lhs_cube, double value);
		
		//--------------------------------------------------------------------
		/**
		 * Friend product operator overloading function to multiply the scalar values of an
		 * Icube object at left hand side of the product sign operator with the double value passed to 
		 * function in the right hand side of the operator. 
		 * @brief Arithmetic product operator overloading for the class objects with doubles values.
		 * @param Icube const reference object to be multiplied at LHS side of the product sign operator.
		 * @param Double value to multiply with the scalar values of the Icube object.
		 * @return Icube object with the resultant scalar values.
		 */
		friend Icube operator*(const Icube& lhs_cube,double value);
		
		//--------------------------------------------------------------------
		/**
		 * Friend division operator overloading function to divide the scalar values of an
		 * Icube object at left hand side of the division sign operator with the double value passed to 
		 * function in the right hand side of the operator. 
		 * @brief Arithmetic division operator overloading for the class objects with doubles values.
		 * @param Icube const reference object to be divided at LHS side of the division sign operator.
		 * @param Double value to divide with the scalar values of the Icube object.
		 * @return Icube object with the resultant scalar values.
		 */		 
		friend Icube operator/(const Icube& lhs_cube, double value);
		
		//--------------------------------------------------------------------
		/**
		 * Exponentiate each voxel values from scalar member with a double value
		 * provided. 
		 * @brief Method to scale the scalar values of the cube ogbject. 
		 * @param Double to a value to scale the voxels.
		 * @return Icube object with the calculated voxels.
		 */
		Icube scale_cube(double val);	

		//----------------------------------------------------------------------
		/**
		 * @brief 
		 * @return 
		 */
		Icube log_cube();
		
		//--------------------------------------------------------------------
		/**
		 * @brief Class method to calculate the square of the voxels values.
		 * @return Icube object with the calculated voxels.
		 */
		Icube SQ();
		
		//---------------------------------------------------------------------
		/**
		 * Calculates the numerical integral in all scalar field data spaces.
		 * @brief Class method to calculate the sum of the voxel values. 
		 * @return Double with the total value of the sum.
		 */
		double calc_cube_integral();
		
		//---------------------------------------------------------------------
		/**
		 * @brief Normalizes the cube integral value to one.
		 * @return None
		 */
		Icube normalize();
		
		//--------------------------------------------------------------------
		/**
		 * @brief calculates the integral of the difference between two cubes
		 * @return double with the value of the difference
		 */
		double diff_integral(const Icube& cube);
		
		//--------------------------------------------------------------------
		/**
		 * @brief Calculates the similarity index of the current object from another Icube.
		 * @param Icube constant object to  make the comparison.
		 * @param Sting with the method type.
		 * @return Double with the calculated similarity value.
		 */
		double similarity_index(Icube& cube, std::string type);
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculates the electron density complement for the scalar values.
		 * @param Bool indicating if the complement calculation will be done using log method.
		 * @return Icube to the complement of the density.
		 */
		Icube calculate_complement(bool clog);
		
		//---------------------------------------------------------------------
		/**
		 * @brief Class method to get the value of the voxels with the three-dimensional indices.
		 * @param Integer value of inddex for first grid  dimension.
		 * @param Integer value of inddex for second grid  dimension.
		 * @param Integer value of inddex for third grid  dimension.
		 * @return Double value of the pointed voxel.
		 */
		double get_scalar(int x, int y, int z);
		
		//---------------------------------------------------------------------
		/**
		 * @brief  Class method to get the value of the indicated voxel in the vicinities of given coordinates.
		 * @param  Double with x coordinate.
		 * @param  Double with y coordinate.
		 * @param  Double with z coordinate.
		 * @return Double with the value of the voxel in these cooridnates.
		 */
		double get_scalar(double x, double y, double z);
		 
		//----------------------------------------------------------------------
		/** 
		 *  Made to fill earlier instantiate objects with basic data needed to work
		 *  with Icube objects. 
		 *  @brief Fill this object with fundamental data from another class from these libs.
		 *  @param Pointer to double with array of scalar data.
		 *  @param Pointer to double with array origin vertice.
		 *  @param Pointer to double with array with size.
		 *  @param Reference to a Imolecule object with molecular info.
		 *  @return None.
		 */		 
		void add_data(std::vector < std::vector < std::vector<double> > >& data);
		
		//--------------------------------------------------------------------
		/**
		 * Write a cube file in gaussian formatted program style to be readed by vizualization programs,
		 * or simply to store in a file the object data.* 
		 * @brief Method to write a cube file with cube information stored in the object.
		 * @param String with the name of the file to be written. 
		 * @return None.
		 */
		void write_cube(std::string cubeName);
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculates the main distribution stat of the scalar field values.
		 * @parameter Double reference to receive the average mean value.
		 * @parameter Double reference to receive the sum value.
		 * @parameter Double reference to receive the minimum  value.
		 * @parameter Double reference to receive the maximum  value.
		 * @return None.
		 */
		void get_cube_stats(double& mean, double& sum, double& min, double& max);
		
		//---------------------------------------------------------------------
		/**
		 * @brief Class method to print to the console basic informations about the Icube object.
		 * @return None.
		 */
		void print();
		
		//--------------------------------------------------------------------
		/**
		 * Deallocate memory for molecule and scalar member data pointers.
		 * @brief Destructor for the class. 
		 */
		~Icube();
};
//==============================================================================
//===============================END OF CLASS===================================
//==============================================================================
//----------------------------------------------------------------------------------------------------------------------------------------------
/**
 * @class cube_diffs
 * @author igorchem
 * @date 22/05/19
 * @file Icube.h
 * @brief Class to calculate the similarity index between cube files given by a inpu files with the names of cubes.
 */
class cube_diffs {
	public:
		//--------------------------------------------------------------------------------------------
		/**
		* @brief STL vector to hold the calculated values
		*/
		std::vector<double> diffs;
		
		//--------------------------------------------------------------------------------------------
		/**
		* @brief STL vector to hold the names of the cube files compared.
		*/
		std::vector<std::string> labels;
		
		//--------------------------------------------------------------------------------------------
		/**
		* Reads the input file with the cube file names and calculates the similarity index
		* @brief Default constructor. 
		*/
		cube_diffs(const char* name);
		
		//--------------------------------------------------------------------------------------------
		/**
		* @brief STL vector to hold the calculated values
		*/
		void write();
};

#endif 
