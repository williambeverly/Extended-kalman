#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

	//Firstly, initialise the private variables
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
							0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
							0, 0.0009, 0,
							0, 0, 0.09;

	H_laser_ << 1, 0, 0, 0,
							0, 1, 0, 0;

	//initialise the Jacobian matrix to zero
	Hj_ <<	0, 0, 0, 0,
					0, 0, 0, 0,
					0, 0, 0, 0;

	//secondly, initialise the inherented variables
	//create a 4D state vector, we don't know yet the values of the x state
	//ekf_.x_ = VectorXd(4);

	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ <<	1, 0, 0, 0,
							0, 1, 0, 0,
							0, 0, 1000, 0,
							0, 0, 0, 1000;

	//the initial transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ <<	1, 0, 1, 0,
							0, 1, 0, 1,
							0, 0, 1, 0,
							0, 0, 0, 1;

	noise_ax = 9;
	noise_ay = 9;

	zero_replacement = 0.001;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
  /**
    * Initialize the state ekf_.x_ with the first measurement.
    * Create the covariance matrix.
    * Remember: you'll need to convert radar from polar to cartesian coordinates.
  */
  // first measurement
  cout << "EKF: " << endl;
  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 1, 1, 1, 1;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /*
    Convert radar from polar to cartesian coordinates and initialize state.
		This is the x and y pos, as well as the x and y velocity, as radar provides these.
    */
		double rho = measurement_pack.raw_measurements_[0];
		double phi = measurement_pack.raw_measurements_[1];
		double dphi = measurement_pack.raw_measurements_[2];
		
		//now express the measurements in cart as follows:
		// px and py are the cosine and sine respectively
		// and dpx and dpy are the components of dphi

		//deal with the case that can make px and py both equal to zero
		if (rho == 0) {
			rho = zero_replacement;
		}
		ekf_.x_ << rho*cos(phi), rho*sin(phi), dphi*cos(phi), dphi*sin(phi);
		
  }

  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    /*
    Initialize state.
		This is simply the x and y pos, as lidar does not provide velocity.
		Be aware of the potential to initialise with px & py both equal to zero
    */
		double px = measurement_pack.raw_measurements_[0];
		double py = measurement_pack.raw_measurements_[1];
		if (px == 0 && py == 0)
		{
			px = zero_replacement;
			py = zero_replacement;
		}
		ekf_.x_ << px, py, 0, 0;
  }

	previous_timestamp_ = measurement_pack.timestamp_;
  // done initializing, no need to predict or update
  is_initialized_ = true;
  return;
  }

	/*****************************************************************************
		*  Prediction
		****************************************************************************/

	/**
		* Update the state transition matrix F according to the new elapsed time.
		- Time is measured in seconds.
		* Update the process noise covariance matrix.
		* Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
		*/

	//compute the time elapsed between the current and previous measurements
	double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

	previous_timestamp_ = measurement_pack.timestamp_;

	double dt_2 = dt * dt;
	double dt_3 = dt_2 * dt;
	double dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ <<	dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
							0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
							dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
							0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
		ekf_.H_ = MatrixXd(3, 4); //initialise the holder - - QUERY: Is this dynamic in allocation?
		Hj_ = tools.CalculateJacobian(ekf_.x_); //calculate the Jacob
		ekf_.H_ = Hj_; //set the ekf to the appropriate value
		ekf_.R_ = MatrixXd(3, 3); //initialise the holder - - QUERY: Is this dynamic in allocation?
		ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
	else {
    // Laser updates
		ekf_.H_ = MatrixXd(2, 4); //initialise the holder - QUERY: Is this dynamic in allocation?
		ekf_.H_ = H_laser_;
		ekf_.R_ = MatrixXd(2, 2); //initialise the holder - QUERY: Is this dynamic in allocation?
		ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
