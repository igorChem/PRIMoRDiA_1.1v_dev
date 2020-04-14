//QMdriver.h

#ifndef QMDRIVER
#define QMDRIVER
//C++ headers
#include <iostream>
#include <memory>
//PRIMoRDiA headers
#include "../include/common.h"
//-------------------------
//foward declarations 
class Iaorbital;
class Imolecule;
//-------------------------
//======================================================================================
/**
 * @class QMdriver
 * @author Igor Barden Grillo
 * @date 21/09/18
 * @file QMdriver.h
 * @brief This class is meant to manipulate quantum chemical coefficients from output programs 
 * to calculate properties for each atoms
 */
class QMdriver{
	public:		
		//-------------------------------------------------------------
		/**
		 * @brief Imolecule smart pointer to hold molecular quantum chemical 
		 * information.
		 */
		Imolecule molecule;

		//-------------------------------------------------------------
		/**
		 * @brief Default constructor.
		 */
		QMdriver();
		
		//-------------------------------------------------------------
		/**
		 * @brief Constructor with one-arg, loading the object from molecular information.
		 * @param Imolecule rvalue reference to the object of molecular information.
		 */
		QMdriver(Imolecule&& mol) noexcept;
		
		//--------------------------------------------------------------
		/**
		 * @brief Copy constructor.
		 * @param QMdriver constant reference of the object to be copied.
		 */
		QMdriver(const QMdriver& qmd_rhs) = delete;
		
		//---------------------------------------------------------------
		/**
		 * @brief Assigment operator overloading.
		 * @param QMdriver constant reference of the object to be copied.
		 * @return QMdriver reference to the current object.
		 */
		QMdriver& operator=(const QMdriver& qmd_rhs) = delete ;
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Calculata the molecular orbital contribution for a given atoms.
		 * @param Integer with the index of the atom.
		 * @param Integer with the index of the molecular orbital.
		 * @return Double to value of the atom for the given molecular orbital.
		 */
		double MO_in_atoms(int atom, int MO,bool beta);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Calculate electrophilic attack susceptibility to the atom indicated in the molecule
		 * @param Integer with the index of the atom.
		 * @return Double with the reactivity descriptor assigned to the atom.
		 */
		double EAS_in_atoms(int atom,int band);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Calculate electrophilic attack susceptibility to the atom indicated in the molecule
		 * @param Integer with the index of the atom.
		 * @return Double with the reactivity descriptor assigned to the atom.
		 */
		double EAS_in_atoms_EW(int atom);
		
		//----------------------------------------------------------------------
		/**
		 * @brief Calculate nucleophilic attack susceptibility
		 * @param Integer with the index of the atom.
		 * @return Double with the reactivity descriptor assigned to the atom.
		 */
		double NAS_in_atoms(int atom,int band);
		
		//----------------------------------------------------------------------
		/**
		 * @brief
		 * @param Integer with the index of the atom.
		 * @return Double with the reactivity descriptor assigned to the atom.
		 */
		double NAS_in_atoms_EW(int atom);	
		
		//----------------------------------------------------------------------
		/**
		 * @brief
		 * @param Integer with the index of the atom.
		 * @return Double with the reactivity descriptor assigned to the atom.
		 */
		double fukushima(int atom,int band);
		//-----------------------------------------------------------------------
		/**
		 * @brief Calculate electron density to the atom indicated in the molecule.
		 * @param Integer with the index of the atom.
		 * @return Double with electron density calculated for the atom indicated.
		 */
		double density_in_atoms(int atom);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Destructor.
		 */
		~QMdriver();
};

#endif 