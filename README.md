# Extended-kalman
Udacity - Term 2 - Project 1

This project involved the construction of a Kalman and Extended Kalman filter, to read in RADAR and LIDAR data, to determine the relative position and speed of a tracked object. At each time step, the estimated state is compared to the ground-truth value, in order to produce a Root Mean Squared Error (RMSE) of the positive (in x and y) and the velocity (in x and y).

The project was built in VS 2015, however, the VS project files have been removed, and the files have been tested in Windows 10 using cmake, minGW and make versions as follows:
* cmake v 3.7.2
* make v 3.81
* minGW c++ (GCC) v5.3.0

The user may enter the following commands in cmd to build, in Windows:
* `mkdir build && cd build`
* `cmake .. -G "Unix Makefiles" && make`
* `ExtendedKF.exe ../data/sample-laser-radar-measurement-data-1.txt output.txt`

The files were also tested on Ubuntu 16.04. The user may enter the following similar commands:
* `mkdir build && cd build`
* `cmake .. -G "Unix Makefiles" && make`
* `./ExtendedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`


