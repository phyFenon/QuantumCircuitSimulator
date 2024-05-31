#ifndef SINGLEGATE_MAT_H
#define SINGLEGATE_MAT_H
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

class SingleGate {
 private:
     std::string gateType;
     std::vector<double> angles;
 public:
     SingleGate(std::string gateType);
     SingleGate(std::string gateType, double angle);
     SingleGate(std::string gateType, std::vector<double> angles);
     Eigen::Matrix2cd hGate();
     Eigen::Matrix2cd xGate();
     Eigen::Matrix2cd yGate();
     Eigen::Matrix2cd zGate();
     Eigen::Matrix2cd rxGate();
     Eigen::Matrix2cd ryGate();
     Eigen::Matrix2cd rzGate();
     Eigen::Matrix2cd u3Gate();
     Eigen::Matrix2cd gateMat();
};
#endif // SINGLEGATE_MAT_H