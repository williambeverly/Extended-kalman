#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */

	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	if (estimations.size() == 0 || (estimations.size() != ground_truth.size())) {
		std::cout << "Error in estimations or ground_truth vector sizes";
		return rmse;
	}


	for (unsigned int i = 0; i < estimations.size(); i++)
	{
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse / estimations.size();

	//calculate the sqrt
	rmse = rmse.array().sqrt();

	return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3, 4);
	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	double c1 = px*px + py*py;
	double c2 = sqrt(c1);
	double c3 = c1*c2;

	//check division by zero
	if (fabs(c1) < 0.001)
	{
		Hj << 0, 0, 0, 0,
					0, 0, 0, 0,
					0, 0, 0, 0;

		std::cout << "CalculateJacobian() - Error division by zero\n";
	}
	else //compute the Jacobian matrix
	{
		Hj << px / c2, py / c2, 0, 0,
					-py / c1, px / c1, 0, 0,
					py*(vx*py - vy*px) / c3, px*(vy*px - vx*py) / c3, px / c2, py / c2;
		/*
		Writing out in full to check
		Hj(0, 0) = px/c2;
		Hj(0, 1) = py/c2;
		Hj(0, 2) = 0;
		Hj(0, 3) = 0;
		Hj(1, 0) = -py/c1;
		Hj(1, 1) = px/c1;
		Hj(1, 2) = 0;
		Hj(1, 3) = 0;
		Hj(2, 0) = py*(vx*py-vy*px)/c3;
		Hj(2, 1) = px*(vy*px-vx*py)/c3;
		Hj(2, 2) = Hj(0, 0);
		Hj(2, 3) = Hj(0, 1);
		*/
	}

	return Hj;
}
