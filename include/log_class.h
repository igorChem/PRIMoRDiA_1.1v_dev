//log.h
//------------------------------
#ifndef LOGCLASS
#define LOGCLASS
//------------------------------
#include <string>
#include <fstream>
//=================================================
/**
 * @class Ilog
 * @author igor Barden Grillo
 * @date 10/09/18
 * @file log_class.h
 * @brief Class with tools to output log of the program.
 * 
 * Class to store the messages for each library class to output 
 * in the scrren or/and in a file.
 */
class Ilog{
	public:
		/**
		* @brief Name used to write the log file.
		*/
		std::string name;
	
		//-----------------------------------------------------------------------
		/**
		* @brief Output stream to be written in log file or outputed in the screen.
		*/
		std::ofstream log_file;
		
		//-----------------------------------------------------------------------
		/**
		* @brief Bool variable to indicate if the messages will be outputed on the screen.
		*/
		bool screen_output;
		
		//-----------------------------------------------------------------------
		/**
		* @brief Bool variable to indicate if the messages will be written on a file.
		*/
		bool file_output;
		
		//-----------------------------------------------------------------------
		/**
		* @brief Default class constructor.
		*/
		Ilog();
		
		//-----------------------------------------------------------------------
		/**
		 *
		 */
		Ilog(const Ilog& rhs_log) = delete;
		
		//-----------------------------------------------------------------------
		/**
		 *
		 */
		Ilog& operator=(const Ilog& rhs_log) = delete;
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Initializes the object information.
		 * @param String with the name of the file to being writting.
		 * @param Bool indicating if the messages will be printed on the screen.
		 * @param Bool indicating if the messages will be written in the log file.
		 * @return None.
		 */
		void initialize(std::string naming,bool sout, bool fout);
		
		//-----------------------------------------------------------------------
		/**
		* @brief Function to receive the string with the message.
		* @param String with the message.
		* @return None.
		*/
		void input_message(std::string message);
		
		//-----------------------------------------------------------------------
		/**
		* @brief Function to receive the string with the message.
		* @param String with the message.
		* @return None.
		*/
		void input_message(int message );
		
		//-----------------------------------------------------------------------
		/**
		* @brief Function to receive the string with the message.
		* @param String with the message.
		* @return None.
		*/
		void input_message(double message);
		
		//------------------------------------------------------------------------
		/**
		 * @brief Get the time count and returns to the console or/and to the log file.
		 * @return None.
		 */
		void timer();
		
		//------------------------------------------------------------------------
		/**
		 * @brief 
		 * @param message
		 */
		void abort(std::string message);
		
		//-----------------------------------------------------------------------
		/**
		 * @brief Default destructor. Here the text is written in the output stream
		 * when the object goes out of scope.
		 */
		~Ilog();
		
};

#endif