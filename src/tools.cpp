#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    Eigen::MatrixXd rmse(4);
    
    rmse << 0, 0, 0, 0;
    
  
    if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
        cout << "Error and Problems. James wrote this" << endl;
        return rmse;
    }
    
    for (unsigned int i = 0; i < estimations.size(); ++i) {
        Eigen::VectorXd residual = estimations[i] - ground_truth[i];
        
        residual = residual.array() * residual.array();
        
        rmse += residual;
    }
    // Calculating Array
    rmse = rmse / estimations.size();
    
    //Calculating Root Mean Squared
    rmse = rmse.array().sqrt();
    
    return rmse;
    
    
}
