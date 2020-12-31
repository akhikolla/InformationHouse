/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Tree_properties.h
 * Author: valentinianmihailungu
 *
 * Created on 24 August 2020, 15:07
 */

#ifndef TREE_PROPERTIES_H
#define TREE_PROPERTIES_H
#include <vector>
#include <string>

// this class is used to store informarion needed for the bct/kbct function
class Tree_properties{
    public:
        double prior;  
        double log_prior;
        double posterior;
        double log_posterior;
        double odd_posterior;
        double bic;
        double aic;
        double mle;
        int n_leaves;
        int max_depth;
        
        vector <string> context;
};

#endif /* TREE_PROPERTIES_H */

