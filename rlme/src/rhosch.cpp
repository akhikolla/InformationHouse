#include "Rcpp.h"
#include "pairup.h"
using namespace Rcpp;

NumericMatrix expand_grid(NumericVector a, NumericVector b);
double cor(NumericVector x, NumericVector y);

// [[Rcpp::export]]
List rhoschC(NumericVector ahat, NumericVector section_counts, NumericMatrix student_counts) {
  // Count up the number of schools
  int num_schools = section_counts.size();
  
  NumericMatrix rho(num_schools, 4);
  
  // The design of the data can be unbalanced. This means
  // that schools may have differing number of sections, and
  // similarly sections may have different numbers of students.
  // this makes things more complicated.
  //
  // To account for this as we go through the data, we keep three index
  // counters. The first one, index_i, is the index to the first ahat value
  // for the current school.
  // The second, index_j, is the index to the first ahat value in the current
  // section. (section j)
  // index_k is the index to the first ahat value in section_k.
  int index_i = 0;
  int index_j = 0;
  int index_k = 0;
  
  // Count up the number of students in each school
  NumericVector school_student_counts(num_schools);
  
  for(int i = 0; i < num_schools; i++) {
    school_student_counts[i] = 0;
    
    for(int j = 0; j < section_counts[i]; j++) {
      school_student_counts[i] += student_counts(i, j);
    }
  }
  
  // Count up the total pair length for each school
  NumericVector total_pair_length(num_schools);
  
  // Loop through each school
  for(int i = 0; i < num_schools; i++) {
    int num_sections = section_counts[i];
    
    // Extract the ahats for this school
    NumericVector a_i(school_student_counts[i]);
    for(int a = 0; a < a_i.size(); a++) {
      a_i[a] = ahat[index_i + a];
    }
    
    total_pair_length[i] = 0;
    
    // Precalculate the mean of the ahats for this school
    double a_i_mean = mean(a_i);
    
    rho(i, 0) = 0;
    rho(i, 1) = 0;
    rho(i, 2) = 0;
    rho(i, 3) = 0;
    
    // Loop through the sections within the school
    for(int j = 0; j < num_sections - 1; j++) {
      
      NumericVector a_ij(student_counts(i, j));
      for(int a = 0; a < student_counts(i, j); a++) {
        a_ij[a] = ahat[index_j + a];
      }
      
      // Precalculate a_ij mean
      double a_ij_mean = mean(a_ij);
        
      index_k = index_j + student_counts(i, j);
      
      for(int k = j + 1; k < num_sections; k++) {
        // Extract ahat values for a_ik
        
        NumericVector a_ik(student_counts(i, k));
        for(int a = 0; a < student_counts(i, k); a++) {
          a_ik[a] = ahat[index_k + a];
        }
        
        NumericMatrix pairs = expand_grid(a_ij, a_ik);
        
        total_pair_length[i] += pairs.nrow();
        
        //std::cout << "Comparing school " << i << " section " << j << " to school " << i << " section " << k << std::endl;
        
        // Calculate pairs[,1] * pairs[,2]
        double pair_mult = sum(pairs.column(0) * pairs.column(1));
        
        // Precalculate a_ik mean
        double a_ik_mean = mean(a_ik);
        
        // Rho 1
        rho(i, 0) += cor(pairs.column(0), pairs.column(1));
        
        // Rho 2
        rho(i, 1) += pair_mult;
        
        // Rho 3
        rho(i, 2) += (pair_mult - pairs.nrow() * a_ij_mean * a_ik_mean) / (sd(pairs.column(0) - a_ij_mean) * sd(pairs.column(1) - a_ik_mean));
        
        // Rho 4
        rho(i, 3) += (pair_mult - pairs.nrow() * a_i_mean * pow(a_i_mean, 2)) / (sd(pairs.column(0) - a_i_mean) * sd(pairs.column(1) - a_i_mean));
        
        index_k += student_counts(i, k);
      }
      
      index_j += student_counts(i, j);
    }
    
    index_j += student_counts(i, num_sections - 1);
    
    index_i += school_student_counts[i];
  }
  
  std::map<std::string, NumericVector> result;
  
  result["rho1"] = rho.column(0);
  result["rho2"] = rho.column(1);
  result["rho3"] = rho.column(2);
  result["rho4"] = rho.column(3);
  
  result["nvec"] = school_student_counts;
  result["npair"] = total_pair_length;
  
  return wrap(result);
}

// [[Rcpp::export]]
List rhosectC(NumericVector ahat, NumericVector section_counts, NumericMatrix student_counts) {
  
  int num_schools = section_counts.size();
  
  NumericVector school_student_counts(num_schools);
  
  for(int i = 0; i < num_schools; i++) {
    school_student_counts[i] = 0;
    
    for(int j = 0; j < section_counts[i]; j++) {
      school_student_counts[i] += student_counts(i, j);
    }
  }
  
  // 3d matrix to hold 4 rho values for each section (column) nested within schools (rows)
  NumericMatrix rho1(num_schools, max(section_counts));
  NumericMatrix rho2(num_schools, max(section_counts));
  NumericMatrix rho3(num_schools, max(section_counts));
  NumericMatrix rho4(num_schools, max(section_counts));
  
  for(int i = 0; i < rho1.nrow(); i++) {
    for(int j = 0; j < rho1.ncol(); j++) {
      rho1(i, j) = 0;
      rho2(i, j) = 0;
      rho3(i, j) = 0;
      rho4(i, j) = 0;
    }
  }
  
  
  int index = 0;
  
  for(int i = 0; i < num_schools; i++) {
    int num_sections = section_counts[i];
    
    // Extract the ahats for this school
    
    NumericVector a_i(school_student_counts[i]);
    for(int a = 0; a < a_i.size(); a++) {
      a_i[a] = ahat[index + a];
    }
    
    double a_i_mean = mean(a_i);
    
    int index_j = index;
    
    // Loop through the sections within the school
    for(int j = 0; j < num_sections; j++) {
      NumericVector a_ij(student_counts(i, j));
        
      for(int a = 0; a < student_counts(i, j); a++) {
        a_ij[a] = ahat[index_j + a];
      }
      
      NumericMatrix pairs = pairup(a_ij);
        
      //std::cout << "Comparing school " << i << " section " << j << " to school " << i << " section " << k << std::endl;
        
      // Calculate pairs[,1] * pairs[,2]
      double pair_mult = sum(pairs.column(0) * pairs.column(1));
        
      // Precalculate mean
      double a_ij_mean = mean(a_ij);
        
      
      // Rho 1
      rho1(i, j) = cor(pairs.column(0), pairs.column(1));
      
      // Rho 2
      rho2(i, j) = pair_mult;
      
      
      // Rho 3
      rho3(i, j) = (1.0 / (double)pairs.nrow())
                   * ( sum(pairs.column(0) * pairs.column(1)) - pairs.nrow() * pow(a_ij_mean, 2) ) 
                   / ( sd(pairs.column(0)  - a_ij_mean) * sd(pairs.column(1) - a_ij_mean) );
        
      // Rho 4
      rho4(i, j) = (1.0 / (double)pairs.nrow())
                   * ( pair_mult - pairs.nrow() * pow(a_i_mean, 2) ) 
                   / ( sd(pairs.column(0) - a_i_mean) * sd(pairs.column(1) - a_i_mean) );
      
      index_j += student_counts(i, j);
    }
    
    index += school_student_counts[i];
  }
  
  
  // Set up a named list to return as the result
  // we do this by setting up a "map" datatype, which is like an
  // array but uses strings as keys
  std::map<std::string, NumericMatrix> result;
  
  result["rho1"] = rho1;
  result["rho2"] = rho2;
  result["rho3"] = rho3;
  result["rho4"] = rho4;
  
  return wrap(result);
}

double cor(NumericVector x, NumericVector y) {
  return 1.0 / ((double)x.size() - 1.0) * sum((x - mean(x)) / sd(x) * (y - mean(y)) / sd(y));
}

NumericMatrix expand_grid(NumericVector a, NumericVector b) {
  int as = a.size();
  int bs = b.size();
  
  NumericMatrix m(as * bs, 2);
    
  int index = 0;
    
  for(int i = 0; i < as; i++) {
    for(int j = 0; j < bs; j++) {
      m(index, 0) = a[i];
      m(index, 1) = b[j];
          
      index++;
    }
  }
    
  return(m);
}
