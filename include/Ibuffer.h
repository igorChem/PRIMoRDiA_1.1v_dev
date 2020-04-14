#ifndef IBUFFER
#define IBUFFER
//------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
//------------------------------------------
#include "../include/common.h"
//------------------------------------------
class Iline; //foward declaration
//-------------------------------------------------------------------------------------------
/**
 * Class to hold and manipulate text files. Creates Ilines objects and sotre in a STL vector 
 * to hold each line of the file to be easily acessed and parsed.
 * @class Ibuffer
 * @author Igor Barden Grillo
 * @date 07/04/18
 * @file Ibuffer.h
 * @brief Class to open files and intatiate and hold Iline objects for each line in the file.
 */
class Ibuffer {
	public: 
		/**
		 * @brief Char pointer constant with file name to be open.
		 */
		const char* name;
		
		//------------------------------------------------------------
		/**
		 * @brief Integer value with the number of lines in the file.
		 */
		unsigned int nLines;
		
		//------------------------------------------------------------
		/**
		 * @brief bool indicating if the lines of the file were parsed and stored.
		 */
		bool parsed;
		 
		//------------------------------------------------------------
		/**
		 * @brief STL vector with the Iline objects to hold each file line.
		 */   
		std::vector<Iline> lines;
		  
		//------------------------------------------------------------
		/**
		 * Default costructor to instantiate a empty object to be fille later.
		 * initializes the member data to default values.
		 * @brief Default constructo of the class.
		 */ 
		Ibuffer();
		
		//------------------------------------------------------------
		/**
		 * Instantites the object from the file name to be loaded.
		 * @brief Constructor with file name information.
		 * @param Char pointer with file name. 
		 * @param Bool indicating if the lines will be parsed and stored.
		 */ 
		Ibuffer(const char* file_name,bool parse);
		
		//------------------------------------------------------------
		/**
		 * Instantites the object from the file name to be loaded, 
		 * storing only the lines between the indices given.
		 * @brief Constructor with file name information.
		 * @param Char pointer with file name. 
		 */
		Ibuffer(const char* file_name,int in, int fin);
		
		//-------------------------------------------------------------
		/**
		 * @brief Constructor to parser the lines of the file initiatin where it hits the first 
		 * string given until the second string in the params.
		 * @param Char pointer with file name.
		 * @param String object with starting line to begin the parsing;  
		 * @param String object with the line to stop parsing.
		 */
		Ibuffer(const char* file_name,std::string wrdin,std::string wrdfin);
		
		//------------------------------------------------------------
		/**
		 * Instantites the object from the file name to be loaded.
		 * @brief Constructor with file name information.
		 * @param Char pointer with file name. 
		 * @param STL vector of String object with starting possible lines to begin the parsing;  
		 * @param STL vector of String object with the possible lines to stop parsing.
		 */ 		
		Ibuffer(const char* file_name, std::vector<std::string>& wrds_in,std::vector<std::string>& wrds_fin);
		
		//------------------------------------------------------------
		/**
		 * Instatiates the object from using member data of another class type object 
		 * @brief Copy constructor of the class.
		 * @param Ibuffer constant reference to the object to be copied.
		 */
		Ibuffer(const Ibuffer& rhs_ibuf) = delete;
		
		//------------------------------------------------------------
		/**
		 * Return this object with member data overwritted from member data of another class type
		 * object passed in the right hand side of assigment operator.
		 * @brief Assigment operator overloading for the class.
		 * @param Ibuffer const reference to this object.
		 */
		Ibuffer& operator=(const Ibuffer& rhs_ibuf) = delete;

		//------------------------------------------------------------
		/**
		 * @brief Get a Ibuffer object with only the lines betwen the indices given.
		 * @param int with initial index.
		 * @param int with final index.
		 * @return Reference to the self object.
		 */
		Ibuffer& get_block(int in, int fin);
		
		//------------------------------------------------------------
		/**
		 * @brief class method to print Ibuffer info to the console.
		 */
		void print();
		
		//------------------------------------------------------------
		/**
		 * Destructor for the class.
		 */
		~Ibuffer();	
	
};


#endif