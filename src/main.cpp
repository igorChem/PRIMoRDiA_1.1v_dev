#include <iostream>

#include "../include/interface.h"

/*********************************************************/
int main(int argc, char **argv){
	
	if ( argc <= 1 ){
		cout << "There are no provided arguments!" << endl;
		exit(-1);
	}
	interface PRIMoRDiA(argc,argv);
	PRIMoRDiA.run();
	return 0;
}
/*********************************************************/