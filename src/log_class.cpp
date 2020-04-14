//log.cpp

#include <iostream>
#include <string>
#include <fstream>
//----------------------------------------------------
#include "../include/common.h"
#include "../include/log_class.h"
//-----------------------------------------------------
using std::cout;
using std::endl;

/************************************************************/
Ilog::Ilog()						:
	name("noname")		,
	screen_output(false)	,
	file_output(false)		{
}
/************************************************************/
void Ilog::initialize(std::string naming	,
							bool sout				, 
							bool fout				){
								
	name				= naming;
	screen_output	= sout;
	file_output		= fout;
	if ( sout || fout ) log_file.open(name.c_str());
}
/************************************************************/
void Ilog::input_message(std::string message){
	if ( screen_output )	cout		<< message << endl;
	if ( file_output )			log_file	<< message << endl;
}
/************************************************************/
void Ilog::input_message(int message){
	if ( screen_output ) cout		<< message << endl; 
	if ( file_output )		log_file	<< message << endl;
}
/************************************************************/
void Ilog::input_message(double message){
	if ( screen_output ) cout		<< message << endl; 
	if ( file_output )		log_file	<< message << endl;
}
/************************************************************/
void Ilog::timer(){
	if ( screen_output )	cout		<< chronometer.return_wall_time() << " seconds" << endl;
	if ( file_output )			log_file	<< chronometer.return_wall_time() << " seconds" << endl;
}
/************************************************************/
void Ilog::abort(std::string message){
	cout << message << endl;
	exit(-1);
}
/************************************************************/
Ilog::~Ilog(){
	if ( file_output ) log_file << "Exiting PRIMoRDiA after " << chronometer.return_wall_time() << " seconds" << endl;
	if ( screen_output || file_output ) log_file.close();
}
/************************************************************/