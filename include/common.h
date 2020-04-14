//commom.h
//Declarations for functions and source file with constant values. 

#ifndef COMMON
#define COMMON
//--------------------------------------------------------------------------------
//Including header from the c++
#include <iostream>
#include <string>
#include <vector> 
#include <memory>
//--------------------------------------------------------------------------------
//Including header from PRIMoRDiA library. 
#include "log_class.h"
#include "Itimer.h"

//===========================================================
// GLOBAL VARIABLES for internal usage: DECLARION
//===========================================================
//---------------------------------------------------------
/**
 * @brief Global variable to indicate if gnuplot scripts for radial distribution 
 * of scalar field data for gnuplot.
 */
extern bool dos;

//---------------------------------------------------------
/**
 * @brief bool global variable signilizing the outputting the extra reactivity descriptors.
 */ 
extern bool extra_RD;

//---------------------------------------------------------
/**
 * @brief bool global variable signilizing the outputting of pymol graphical software scripts
 */
extern bool pymol_script;

//---------------------------------------------------------
/**
 * @brief bool global variable signilizing the outputting of VMD graphical software scripts
 */
extern bool vmd_script;

//---------------------------------------------------------
/**
 * @brief bool global variable to indicate R scripts for descriptor heatmaps by residues
 */
extern bool M_R;

//---------------------------------------------------------
/**
 * @brief  integer global variable holding the maximum number of threads to execute the program
 */
extern unsigned int NP; 

//---------------------------------------------------------
/**
 * @brief double global variable indicating the energy criteria to account molecular orbitals in FOA calculations
 */
extern double energy_crit; 

//---------------------------------------------------------
/**
 * @brief Itimer global object that returns total wall time of the program execution.  
 */
extern Itimer chronometer;

//---------------------------------------------------------
/**
 * @brief Ilog object that writes to a log file and/or outputs messages to the console.
 */
extern std::unique_ptr<Ilog> m_log;

//===========================================================
// GLOBAL VARIABLES for main function usage DECLARION
//===========================================================
//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern bool M_verbose;

//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern bool M_logfile;

//------------------------------------------------------------------------------------------
/**
 * @fn get_atom_mass
 * @brief Function to get the atomic mass from array initialized in source file.
 * @param String indicating the element symbol.
 * @return Float with atomic mass value.
 */
float get_atom_mass(std::string sym);

//-----------------------------------------------------------------------------------------
/**
 * @fn get_atomic_number
 * @brief Function to get the atomic number from array initialized in source file.
 * @param String indicating the element symbol.
 * @return Integer with the atomic number.
 */
int get_atomic_number(std::string sym);

//-----------------------------------------------------------------------------------------
/**
 * @fn get_atomic_symbol
 * @brief Function to get the atomic symbol from array initialized in source file.
 * @param Integer indicating the atomic nuber.
 * @return String with the atomic number.
 */
std::string get_atomic_symbol(int i);

//---------------------------------------------------------------------------------------
/**
 * @fn factorial
 * @brief Calculate the factorial of an integer value.
 * @param Integer with the value to get the factorial from.
 * @return Integer with of the factorial of the integer passed. 
 */
int factorial(int n);

//---------------------------------------------------------------------------------------
/**
 * @fn double factorial
 * @brief Calculate the semi-factorial/double factorial of an integer value.
 * @param Integer with the value to get the factorial from.
 * @return Integer with of the factorial of the integer passed. 
 */
int doublefactorial(int n);

//---------------------------------------------------------------------------------------
/**
 * @fn int_to_string
 * @brief Convert integer value to an std string.
 * @param Integer to be converted to string.
 * @return String from integer passed.
 */
std::string int_to_string(int val);

//---------------------------------------------------------------------------------------
/**
 * @fn int_to_string
 * @brief Convert integer value to an std string.
 * @param Integer to be converted to string.
 * @return String from integer passed.
 */
std::string double_to_string(double val);

//---------------------------------------------------------------------------------------
/**
 * @fn IF_file
 * @brief Test if file can be open.
 * @param Const char with the name of the file.
 * @return Bool if the file can be open.
 */
bool IF_file(const char* name); 

//---------------------------------------------------------------------------------------
/**
 * @fn check_file_ext
 * @brief Test if file the file extension is the one that is required.
 * @param Const char with the name of the file.
 * @return Bool if the file can be open.
 */
bool check_file_ext(std::string extension,const char* file_name);

//---------------------------------------------------------------------------------------
/**
 * @fn remove_extension
 * @brief return a string with the filename without the extension part 
 * @param constant char* containing the file name.
 * @return string with the filename without the extension.
 */
std::string remove_extension(const char* file_name);

//---------------------------------------------------------------------------------------
/**
 * @fn change_extension
 * @brief Return a string with a file name with a new extension
 * @param char pointer constant with the file name
 * @param string with the new extension
 * @return string of the file name without the extension
 */
std::string change_extension(const char* file_name,std::string new_ext);

//---------------------------------------------------------------------------------------
/**
 * @brief Convert double in string format from file with "D" to "E"
 * @param string with the double to be converted
 * @return double with the converted values
 */
double D_E_conv(std::string sc_not);

//---------------------------------------------------------------------------------------
/**
 * @brief  Get the average mean from a vector.
 * @param STL vector of doubles
 * @return double with the average mean
 */
double mean_dvec(std::vector<double>& vec);
//---------------------------------------------------------------------------------------
/**
 * @brief  Get the maximum value from a vector.
 * @param STL vector of doubles
 * @return double with the maximum value
 */
double max_dvec(std::vector<double>& vec);
//---------------------------------------------------------------------------------------
/**
 * @brief  Get the minimum from a vector.
 * @param STL vector of doubles
 * @return double with the minimum value
 */
double min_dvec(std::vector<double>& vec);
//---------------------------------------------------------------------------------------
/**
 * @brief  Get the sum from a vector.
 * @param STL vector of doubles
 * @return double with the sum value
 */
double sum_dvec(std::vector<double>& vec);

//---------------------------------------------------------------------------------------
/**
 * @brief 
 * @param line
 * @param in
 * @param fin
 * @return 
 */
std::string str_array(std::string& line, int in, int fin);


#endif