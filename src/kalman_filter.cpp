#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */

	//Set z_predicted to polar coordinates using h(x')
	VectorXd z_pred = Cart2Pol(x_);
	//Ensure that the angle phi is between -pi and pi
	VectorXd y = NormaliseError(z, z_pred);

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

VectorXd KalmanFilter::Cart2Pol(const VectorXd &_x)
{
	/**
	*	This function transforms the known state variables from cartesian
		to polar coordinate system
	*/
	float px = _x[0];
	float py = _x[1];
	float vx = _x[2];
	float vy = _x[3];

	float c1 = sqrt(px * px + py * py);
	float phi = std::atan2(py,px);
	float drho = (px * vx + py * vy) / c1;

	VectorXd PolarState = VectorXd(3);
	PolarState << c1, phi, drho;
	return PolarState;

}

VectorXd KalmanFilter::NormaliseError(const VectorXd &z, const VectorXd &z_pred)
{
	/**
	*	This function ensures that the phi angle is inbetween -pi and pi
	*/
	VectorXd NormError = VectorXd(3);
	NormError = z - z_pred;
	while (NormError[1] < -M_PI) {
		NormError[1] += M_PI;
	}

	while (NormError[1] > M_PI) {
		NormError[1] -= M_PI;
	}
	return NormError;
}
