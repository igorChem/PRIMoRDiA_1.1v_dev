#ifndef OVERLAP
#define OVERLAP

#include <iostream>
#include <string>
#include <vector> 
// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
// Libint Gaussian integrals library
#include <libint2.hpp>


#include "../include/common.h"

/****************************************/
class Iatom;
//--------------------------------------
class overlap_int{
	public:
		std::vector<libint2::Shell> shells;
		std::vector<int> shell_map;
		std::vector<int> shell_size;
		overlap_int();
		~overlap_int();
		overlap_int(const std::vector<Iatom>& atoms);
		overlap_int(const overlap_int& rhs) = delete;
		overlap_int& operator=(const overlap_int& rhs) = delete;
		void calculate_overlap(std::vector<double>& overlap_matrix);
};


#endif