#ifndef AUTO_PRIMORDIA 
#define AUTO_PRIMORDIA

#include<vector>

class primordia;
class compPrimordia;

//==========================================
/**
 * Class to automate the use of primordia class to calculate the RD from a collection of files in a folder.
 * The Files to be used in the calculations is indicated in a input-like file with structure 
 * @class AutoPrimordia
 * @author Igor Barden Grillo
 * @date 05/05/18
 * @file primordia.h
 * @brief Automate the reactivity descriptor calculations for a given list of converged QM calculations.
 */
class AutoPrimordia{
	public:
		//--------------------------------------------------------------------
		/**
		*
		*/
		const char* m_file_list;
		
		//--------------------------------------------------------------------
		/**
		*
		*/
		std::vector< std::vector< std::vector<double> > > traj_atoms_rd; 
		
		//--------------------------------------------------------------------
		/**
		*
		*/
		std::vector< int > atoms;
		
		//--------------------------------------------------------------------
		/**
		 * @brief STL vector of primordia objects to store reactivity descriptors calculated by the class.
		 */
		std::vector<primordia> RDs;
		
		//---------------------------------------------------------------------
		/**
		 * Default constructor to initialize empty class object.
		 * @brief Default constructor.
		 */
		AutoPrimordia();
		
		//---------------------------------------------------------------------
		/**
		 * Reads the file list and executes the primordia class workflow to calculate the reactivity 
		 * descriptors for a bunch of quantum chemical calculated systems, Sets and calculation the method
		 * comparsions if required in the input list.
		 * @brief Instantiates the object and calculate the reactivity descriptors from a give list.
		 * @param Char pointer constant with the name of the file with the list of RD to be calculated.
		 */
		AutoPrimordia(const char* file_list);
		
		//---------------------------------------------------------------------
		/**
		 * @brief Copy constructor of the class.
		 */
		AutoPrimordia(const AutoPrimordia& Apr_rhs) = delete;
		
		//---------------------------------------------------------------------
		/**
		 * @brief Assigment operator overloading for class type objects.
		 */
		AutoPrimordia& operator=(const AutoPrimordia& Apr_rhs) = delete;	
		
		//---------------------------------------------------------------------
		/**
		 * @brief Destrcutor of the class.
		 */
		~AutoPrimordia();
		
		//----------------------------------------------------------------------
		/**
		 * @brief 
		 */
		void write_global();
		
		//----------------------------------------------------------------------
		void traj_atoms();
	
};

#endif