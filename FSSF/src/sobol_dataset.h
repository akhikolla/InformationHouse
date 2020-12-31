#ifndef sobol_dataset_hpp
#define sobol_dataset_hpp

# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>
#include <vector>


using namespace std;

typedef std::vector<std::vector<double> > my_vector_of_vectors_t;

//int main ( int argc, char *argv[] );
int i8_bit_hi1 ( long long int n );
int i8_bit_lo0 ( long long int n );
void i8_sobol ( int dim_num, long long int *seed, double quasi[] );
void i8_sobol_generate ( int m, int n, int skip, my_vector_of_vectors_t r );
void r8mat_write ( string output_filename, int m, int n, my_vector_of_vectors_t &table );
void timestamp ( );

#endif /* Sobol_dataset_hpp */
