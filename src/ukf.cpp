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
    
    is_initialized_ = false;
    
    //Augmented state
    
  n_aug_ = 7;
    
  // Lambda
    lambda_ = 3 - n_aug_;
    
    NIS_lidar = 0.0;
    NIS_radar = 0.0;
    
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
   */// skip processing if the both sensors are ignored
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
        (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_))
    {
   
        if (!is_initialized_) {
            cout << "UKF:" << endl;
            x_ << 1,1,1,1,1;
            
            P_ << 1,0,0,0,0,
                0,1,0,0,0,
                0,0,1,0,0,
                0,0,0,1,0,
                0,0,0,0,1;
            
            time_us_ = meas_package.timestamp_;

            if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
                cout << "UKF 2" << endl;

                float rho = meas_package.raw_measurements_(0);
                float phi = meas_package.raw_measurements_(1);
                x_(0) = rho*cos(phi);
                x_(1) = rho*sin(phi);
            
            } else if (meas_package.sensor_type_ == MeasurementPackage:: LASER) {
               cout << "UKF 3" << endl;
                x_(0) = meas_package.raw_measurements_(0);
                x_(1) = meas_package.raw_measurements_(1);

            }

            is_initialized_ = true;
            
            return;
        }
        
        float dt = (meas_package.timestamp_  - time_us_) / 100000.0;
        time_us_ = meas_package.timestamp_;
        
        Prediction(dt);
        
        if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            cout << "UKF 4" << endl;

            UpdateRadar(meas_package);

        }
        else if (meas_package.sensor_type_ == MeasurementPackage:: LASER) {
            cout << "UKF 5" << endl;

            UpdateLidar(meas_package);

        }
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
    Eigen::MatrixXd Xsig = MatrixXd(n_aug_, 2 * n_x_ + 1);
    //Calculate Square Root of P

    Eigen::VectorXd x_aug = Eigen::VectorXd(n_aug_);
    Eigen::MatrixXd P_aug = Eigen::MatrixXd(n_aug_, n_aug_);
    
    // Matrix with predicted sigma points
    Eigen::MatrixXd Xsig_pred = Eigen::MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5, 5) = P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;
    
    // Create square root matrix
    Eigen::MatrixXd L = P_aug.llt().matrixL();
    
    //Create sigma points
    Xsig.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++) {
        Xsig.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig.col(i + n_aug_ + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        
    }
    
    //predict sigma points
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        double px = Xsig(0, i);
        double py = Xsig(1, i);
        double v = Xsig(2, i);
        double yaw = Xsig(3, i);
        double yawd = Xsig(4, i);
        double nu_a = Xsig(5, i);
        double nu_yawdd = Xsig(6, i);
        
        double px_p, py_p;
        
        // avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = px + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = py + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = px + v*delta_t*cos(yaw);
            py_p = py + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;
        
        // Add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
        
    }
    
    // Predict Mean and Covariance
    
    //         Sigma point weights
    Eigen::VectorXd weights = Eigen::VectorXd(2 * n_aug_ + 1);
    
    
    weights(0) = lambda_ / (lambda_ + n_aug_);
    
    for (int i = 1; i < 2 * n_aug_ + 1; i++) {
        weights(i) =  0.5 * (lambda_ + n_aug_);
    }
    
    //         Predict State Mean
    x_.fill(0.0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        x_ = x_ + weights(i) * Xsig_pred.col(i);
    }
    
    P_.fill(0.0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        Eigen::VectorXd x_diff = Xsig_pred.col(i) - x_;
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
    }
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
    
    Eigen::VectorXd z = meas_package.raw_measurements_;
    
    //Measurement dimenstion
    int n_z = 2;
    
    Eigen::MatrixXd Zsig = Xsig_pred_.block(0,0, n_z, 2 * n_aug_ + 1);

    
    //predicted measurement mean
    Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
    
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    //Measurement covariance S
    
    Eigen::MatrixXd S = Eigen::MatrixXd(n_z, n_z);
    
    S.fill(0.0);
    for(int i = 0; i < 2 * n_aug_ + 1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    // Add noise
    
    Eigen::MatrixXd R = Eigen::MatrixXd(n_z, n_z);
    R << std_laspx_*std_laspx_, 0, 0, std_laspy_*std_laspy_;
    
    S = S + R;
    
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    Eigen::MatrixXd K = Tc * S.inverse();
    
    float p_x = meas_package.raw_measurements_(0);
    float p_y = meas_package.raw_measurements_(1);
    
    z << p_x, p_y;
    VectorXd z_diff = z - z_pred;
    
    //NIS
    NIS_lidar = z_diff.transpose() * S.inverse() * z_diff;
    
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
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
    
    //Measurement dimension
    int n_z = 3;
    
    Eigen::VectorXd z = meas_package.raw_measurements_;
    
    
    Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z, 2 * n_aug_ + 1);
    
    Eigen::MatrixXd Tc = Eigen::MatrixXd(n_z, n_z);
    
    Eigen::MatrixXd xsig = Eigen::MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    VectorXd z_pred = VectorXd(n_z);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        // extract values for better readibility
        double px = Xsig_pred_(0, i);
        double py = Xsig_pred_(1, i);
        double v  = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        Zsig(0, i) = sqrt(px*px + py*py);
        Zsig(1, i) = atan2(py, px);
        Zsig(2, i) = (px*v1 + py*v2) / sqrt(px*px + py*py);
    }
    
    //Mean predicted measurement
    Eigen::MatrixXd R = Eigen::MatrixXd(n_z, n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    //  Covariance matrix S
    Eigen::MatrixXd S = Eigen::MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
        Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;
        
        while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    Eigen::MatrixXd K = Tc * S.inverse();
    
    float rho = meas_package.raw_measurements_(0);
    float theta = meas_package.raw_measurements_(1);
    float rho_dot = meas_package.raw_measurements_(2);
    
 
    z << rho, theta, rho_dot;
    
    VectorXd z_diff = z - z_pred;

    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    
    //calculate NIS
    NIS_radar = z_diff.transpose() * S.inverse() * z_diff;
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    
    
    
}
