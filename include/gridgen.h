//claas header for grigen.cpp to generating electron density cubes from MO 
// grigen.h

#ifndef GRIDGEN_H
#define GRIDGEN_H
// C++ header files
#include <iostream>
#include <string>
//Our library includes
#include "../include/common.h"

class Iaorbital;
class Imolecule;
class Icube; //foward declaration
//======================================================================================
/**
 * Class to automate the generation of cube file and Icube objects from the generation of
 * three dimensional scalar fields. 
 * @class gridgen
 * @author barden
 * @date 20/03/18
 * @file gridgen.h
 * @brief Class to generate grid scalar values from quantum computational calculations.
 */
class gridgen {
	public:
		//-----------------------------------------------------------------------
		/**
		 * @brief String with the system name.
		 */
		std::string name;
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Bool to indicate if the grid to be generated is about a molecular orbitals.
		 */
		bool orbital;
			
		//-----------------------------------------------------------------------
		/**
		 * @brief Integer with the number of molecular orbital that is to be generated.
		 */
		int Norb;
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Integer unsigned value with the number of voxels in density Icube
		 */
		unsigned int points;
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Imolecule smart pointer with molecular information.
		 */
		Imolecule molecule;

		//-----------------------------------------------------------------------
		/**
		 * @brief Double array to the origin values of the grid. 
		 */
		double origin[3];
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Double array with the coordinates of the top vertice in the grid.
		 */
		double top_corner[3];
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Double array of grid lenghts. 
		 */
		double grid_sides[3];
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Integer value with number of voxels in one of the three dimensions.
		 */
		unsigned int grid_len[3];
		
		//----------------------------------------------------------------------
		/**
		 * @brief STL double vector with the x dimensions coordinates for the atomic orbitals. 
		 */
		std::vector<double> AOxcoords;
		
		//----------------------------------------------------------------------
		/**
		 * @brief STL double vector  with the y dimensions coordinates for the atomic orbitals. 
		 */		
		std::vector<double> AOycoords;
		
		//----------------------------------------------------------------------
		/**
		 * @brief STL double vector with the z dimensions coordinates for the atomic orbitals. 
		 */	
		std::vector<double> AOzcoords; 
		
		//----------------------------------------------------------------------
		/**
		 * @brief STL Iaorbital vector  of Iaorbital objects in the molecule for the grid calculation.
		 */
		std::vector<Iaorbital> orbs;
		
		//-----------------------------------------------------------------------
		/**
		 * @brief STL vector of vectors of double vectors to hold the voxel values of the three 
		 * dimensional grid;
		 */
		std::vector < std::vector < std::vector<double> > > psi;
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Icube obect with grid hold the grid calculated information
		 */
		Icube density;
		
		//----------------------------------------------------------------------
		/**
		 * Instatiates empty object with default initialization values for its member data.
		 * @brief Default constructor of the class.
		 */
		gridgen();
		
		//-----------------------------------------------------------------------
		/**
		 * Instantiates the object with necessary information to calculate 
		 * the grid with quantum chemical info scalar values to write a cube file.
		 * @brief Instatiates the object from molecular information and grid finess required.
		 * @param Integer reference to the grid lenght of one of its dimensions.
		 * @param Imolecule reference to object with molecular geometrical and qm 
		 * info. 
		 */
		gridgen(int grd, Imolecule&& mol) noexcept;
		
		//-----------------------------------------------------------------------
		/**
		 * Instatiates a class object from the copied member data from another.
		 * @brief Copy constructor.
		 * @param gridgen const reference to object to be copied.
		 */
		gridgen(const gridgen& rhs_grd) = delete;
		 
		//-----------------------------------------------------------------------
		/**
		 * Overwrites the memeber data of this object with a copy from another class 
		 * object in the right hand side of the assigment operator. 
		 * @brief Assigment operator overloading.
		 * @param gridgen const reference to object to be overwritted. 
		 */
		gridgen& operator=(const gridgen& rhs_grd) = delete;
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Calculate the value of atomic orbital for the index of grid passed in the function
		 * using the Slater atomic orbital in the cartesiona scheme.
		 * @param Integer value to the orbital index  stored in the current object.
		 * @param Integer value with the x dimension index of the grid to be calculated by the function.
		 * @param Integer value with the y dimension index of the grid to be calculated by the function.
		 * @param Integer value with the z dimension index of the grid to be calculated by the function.
		 * @return Double value of the calculated orbital.
		 */
		double calc_slater_orb(int i, int x, int y, int z);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Calculate the value of atomic orbital for the index of grid passed in the function
		 * using the gaussian atomic orbital in the cartesiona scheme. 
		 * @param Integer value to the orbital index  stored in the current object.
		 * @param Integer value with the x dimension index of the grid to be calculated by the function.
		 * @param Integer value with the y dimension index of the grid to be calculated by the function.
		 * @param Integer value with the z dimension index of the grid to be calculated by the function.
		 * @return Double value of the calculated orbital.
		 */
		double calc_gauss_orb(int i, int x, int y, int z);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Calculate the value of atomic orbital for the index of grid passed in the function
		 * using the gaussian atomic orbital in the spherical scheme for orca data. 
		 * @param Integer value to the orbital index  stored in the current object.
		 * @param Integer value with the x dimension index of the grid to be calculated by the function.
		 * @param Integer value with the y dimension index of the grid to be calculated by the function.
		 * @param Integer value with the z dimension index of the grid to be calculated by the function.
		 * @return Double value of the calculated orbital.
		 */
		double calc_orca_sphe(int i, int x, int y, int z);	
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Calculates atomic orbital for a given point in space.
		 * @param Integer value to the orbital index  stored in the current object.
		 * @param Integer value with the x dimension index of the grid to be calculated by the function.
		 * @param Integer value with the y dimension index of the grid to be calculated by the function.
		 * @param Integer value with the z dimension index of the grid to be calculated by the function.
		 * @return Double value of the calculated orbital.
		 */ 
		double calc_aorb(int i, int x, int y, int z);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Memeber function made to calculate the value of a given molecular orbital to
		 * a specific voxel of te grid.
		 * @param Integer with molecular orbital number to be calculated.
		 * @param Integer with x index of the desired voxel.
		 * @param Integer with y index of the desired voxel.
		 * @param Integer with z index of the desired voxel.
		 * @return Double with orbital value in the voxel.
		 */
		double calc_orb_voxel(int nm,int x,int y,int z, bool beta);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Memeber function made to calculate the value of a given molecular orbital to
		 * a specific voxel of te grid.
		 * @param Integer with molecular orbital number to be calculated.
		 * @param Integer with x index of the desired voxel.
		 * @param Integer with y index of the desired voxel.
		 * @param Integer with z index of the desired voxel.
		 * @return Double with orbital value in the voxel.
		 */
		double calc_orb_voxel_orca(int nm,int x,int y,int z,bool beta);
		
		//-----------------------------------------------------------------------
		/**
		 * Member function to calculate the grid values to a specific orbital.
		 * @param Integer reference with orbital indice to be calculated.
		 * @param Bool indicanting if the slater type orbital must be calculated in cartesian form.
		 * @return None.
		 */
		void calculate_orb(int Nmo,bool beta);
		
		//-----------------------------------------------------------------------
		/**
		 * Member function to calculate the grid values to a specific orbital for orca program.
		 * @param Integer reference with orbital indice to be calculated.
		 * @param Bool indicanting if the slater type orbital must be calculated in cartesian form.
		 * @return None.
		 */
		void calculate_orb_orca(int Nmo,bool beta);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Class function to calculate electron density from molecular orbitals to specific 
		 * voxel using.
		 * @param Integer with grid x index. 
		 * @param Integer with grid y index. 
		 * @param Integer with grid z index.
		 * @return Double with electron density from the voxel value for the indexes passed.
		 */
		double electron_density_mo(int x, int y, int z);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Class function to calculate electron density from molecular orbitals to specific 
		 * voxel using for orca.
		 * @param Integer with grid x index. 
		 * @param Integer with grid y index. 
		 * @param Integer with grid z index.
		 * @return Double with electron density from the voxel value for the indexes passed.
		 */
		double electron_density_mo_orca(int x, int y, int z);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Class function to calculate electron density from density matrix to specific 
		 * voxel using density matrix in Imolecule object member data.
		 * @param Integer with grid x index. 
		 * @param Integer with grid y index. 
		 * @param Integer with grid z index.
		 * @return Double with electron density from the voxel value for the indexes passed.
		 */
		double electron_density(int x, int y, int z);
	
		//-----------------------------------------------------------------------
		/**
		 * Member function to calculate the grid values of the electron density.
		 * @brief Takes the information loaded in the contructor and calculate the 
		 * scalar values of the gri. Creates a Icube object and writes a cube file.
		 * @param Bool indicanting if the slater type orbital must be calculated in cartesian form.
		 * @return None.
		 */
		 void calculate_density();
		 
		 //-----------------------------------------------------------------------
		/**
		 * Member function to calculate the grid values of the electron density for orca program information.
		 * @brief Takes the information loaded in the contructor and calculate the 
		 * scalar values of the gri. Creates a Icube object and writes a cube file.
		 * @param Bool indicanting if the slater type orbital must be calculated in cartesian form.
		 * @return None.
		 */
		 void calculate_density_orca();
		 
		//-----------------------------------------------------------------------
		/**
		 * Member function to calculate the scalar fields of molecular orbital overlap.
		 * @brief Calculate the "overlap density". 
		 * @return None.
		 */
		void orbs_overlap();
		
		//----------------------------------------------------------------------
		/**
		 * @brief Get value of electrostatic potential map from partial charges.
		 * @return None.
		 */
		void calculate_mep_from_charges();
		
		//-----------------------------------------------------------------------
		/**
		 * @brief member function ot return a Icube object with the generated grid information.
		 * @return Icube reference object. 
		 */
		Icube& get_cube();
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculate the grid for the HOMO.
		 * @return Icube reference to the calculated grid.
		 */
		Icube& calc_HOMO();
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculate the grid for the LUMO.
		 * @return Icube reference to the calculated grid
		 */
		Icube& calc_LUMO();
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculate the density from the values of Highest energy occuied molecular orbital 
		 * from Imolecule object member data.
		 * @return Icube reference object of grid from the calculated scalar values. 
		 */
		Icube& calc_HOMO_density();
		 
		 //---------------------------------------------------------------------
		 /**
		  * @brief Calculate the density from the values of Lowest energy unoccuied molecular orbital 
		 * from Imolecule object member data.
		  * @return Icube reference object of grid from the calculated scalar value. 
		  */
		Icube& calc_LUMO_density();
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculated the grid for HOMO band orbitals.
		 * @param Integer with the size of the band gap.
		 * @return Icube reference object of grid from the calculated scalar values.
		 */
		Icube calc_HOMO_band(int bandn);
		
		//----------------------------------------------------------------------
		/**
		 * @brief Integer with the size of the band gap.
		 * @param Integer with the size of the band gap.
		 * @return Icube reference object of grid from the calculated scalar values.
		 */
		Icube calc_LUMO_band(int bandn);
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculates the electrophilic attack susceptibility from a molecular orbital
		 * band near the HOMO. 
		 * @param Integer with the size of the band gap.
		 * @return Icube reference object of grid from the calculated scalar values.
		 */
		Icube calc_band_EAS(int bandn);
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculates the electrophilic attack susceptibility from a molecular orbital
		 * band near the LUMO.
		 * @param Integer with the size of the band gap.
		 * @return Icube reference object of grid from the calculated scalar values.
		 */
		Icube calc_band_NAS(int bandn);
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculate the electrophilic attack susceptibility from a linear
		 * combination pondered by energy using a band of molecular orbitals near 
		 * the HOMO.
		 * @param Integer with the size of the band gap.
		 * @return Icube reference object of grid from the calculated scalar values.
		 */
		Icube calc_EBLC_EAS();
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculate the electrophilic attack susceptibility from a linear
		 * combination pondered by energy using a band of molecular orbitals near 
		 * the HOMO.
		 * @param Integer with the size of the band gap.
		 * @return Icube reference object of grid from the calculated scalar values.
		 */
		Icube calc_EBLC_NAS();
		
		Icube grid_from_atoms(std::vector<double> values);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Class method  to wrtie the generated grid in a cube file
		 * @return None.
		 */
		void write_grid(); 
		
		//-------------------------------------------------------------------------
		/**
		 * @brief Set new values for limits for the scalar video.
		 * @param Double array with minimum values for the grid vertices.
		 * @param Double array with values of the grid lenghts
		 * @param Integer array with the number of points for each dimension of the scalar 
		 * @return None.
		 */
		void set_lim(double* Min, double* gridSides, int* gridSize);
		
		//-------------------------------------------------------------------------
		void redefine_lim(int atom,double size);
		void redefine_lim(double xc,double yc, double zc,int size);
		//-----------------------------------------------------------------------
		/**
		 * Destructor of the class. 
		 */
		~gridgen();
};

#endif