#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 10;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 10;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Parameters above this line are scaffolding, do not modify
  
  // State Dimension
  n_x_ = 5;
    
  
    
    //Augmented state
    
  n_aug_ = 7;
    
  // Lambda
    lambda_ = 3 - n_aug_;
    
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    
    if (!is_initialized_) {
        cout << "UKF:" << endl;
        x_ << 1, 1, 1, 1, 0.1;
        
        // init covariance matrix
        P_ << 0.15, 0, 0, 0, 0,
        0, 0.15, 0, 0, 0,
        0,    0, 1, 0, 0,
        0,    0, 0, 1, 0,
        0,    0, 0, 0, 1;
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
            float rho = meas_package.raw_measurements_(0);
            float phi = meas_package.raw_measurements_(1);
            float rho_dot = meas_package.raw_measurements_(2);
            x_(0) = rho*cos(phi);
            x_(1) = rho*sin(phi);
            x_(2) = rho_dot*cos(phi);
            x_(3) = rho_dot*sin(phi);
        } else if (meas_package.sensor_type_ == MeasurementPackage:: LASER) {
           
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
        }
        
        is_initialized_ = true;
        
        return;
    }
    
    float dt = (meas_package.timestamp_  - time_us_) / 100000.0;
    time_us_ = meas_package.timestamp_;
    
    Prediction(dt);
    
    if(meas_package.MeasurementPackage:: RADAR) {
        UpdateRadar(meas_package);
    }
    else if (meas_package.MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    }
    
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    //Sigma point matrix
    Eigen::MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
    //Calculate Square Root of P
    Eigen::MatrixXd A = P_.llt().matrixL();
    
    
    

    
    
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
