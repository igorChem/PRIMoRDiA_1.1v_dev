//header file for QM output file parser
// currently done for mopac and GAMESS
// QMparser.h

#ifndef  QMPARSER
#define QMPARSER
//---------------------------------------------
//Including c++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <memory>
//---------------------------------------------
#include "../include/common.h"
//---------------------------------------------
class Imolecule;
/**
 * @class QMparser
 * @author igor Barden Grillo
 * @date 24/03/18
 * @file QMparser.h
 * @brief Class to read and stores files from output of quantum chemical calculations programs
 * to generate objects with info to calculate reactivity descriptors.
 * 
 * This class instatiates a Ibuffer object to load all the lines of a given file from a quantum
 * chemistry output calculation and parse. For now the member functions to parse file of gamess and mopac
 * are implemented, but the class make available the needed tools to parse formated file and store the chemical
 * relevant information. 
 */
class QMparser{
	public:
	//---------------------------------------------------------------------------------------------------
	/**
	 * @brief std::string specifying the program type.
	 */
	std::string program;
	 
	//---------------------------------------------------------------------------------------------------
	/**
	 * @brief Indicating if the file was parsed. 
	 */
	bool parsed;
	
	//---------------------------------------------------------------------------------------------------
	const char* name_f;
	//-------------------------------------------------------------------------------------------------
	/**
	 * @brief Default constructor.
	 */
	QMparser();
	
	//---------------------------------------------------------------------------------------------------
	/**
	 *
	 */
	QMparser(const char* file_name, std::string Program);

	//--------------------------------------------------------------------------------------------------
	/**
	 * @brief Copy constructor. 
	 * @param QMparser reference object to be copied to the current object. 
	 */
	QMparser(const QMparser& rhs_QMp) = delete;
	 
	//-------------------------------------------------------------------------------------------------
	/**
	 * @brief Assigment operator overloading.
	 * @param QMparser object to be copied.
	 * @return QMparser reference to this object.
	 */
	QMparser& operator=(const QMparser& rhs_QMp) = delete;
	
	//------------------------------------------------------------------------------------------------
	Imolecule get_molecule(unsigned int mob);
	
	//------------------------------------------------------------------------------------------------
	/**
	 * @brief Destructor of the class.
	 */
	~QMparser();
	
};
 //====================================================================
 //=================END OF QMPARSER CLASS==============================
 //===================================================================
 
#endif
