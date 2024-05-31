#ifndef QUANTUM_SIMULATOR_H
#define QUANTUM_SIMULATOR_H
#include <string>
#include <complex>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <random>
#include <bitset>
#include <omp.h>



class QuantumSimulator
{
private:
    std::vector<std::string> singleGateList={"h", "x", "y", "z", "rx", "ry", "rz", "u3"};
    std::vector<std::string> doubleGateList={"cx", "swp", "iswp", "cz", "cp", "rxx", "ryy", "rzz"};
    std::vector<std::string> tripleGateList={"ccx"};
public:
    int num_qubits;
    std::vector<std::complex<double> > psi;
    QuantumSimulator(int num_qubits);
    void SingleGate_ActOn_State(Eigen::Matrix2cd gateMat, int gatePos);
    void DoubleGate_ActOn_State(std::string gateType, std::vector<int> gatePos, double angle=0.0);
    void TripleGate_ActOn_State(std::string gateType, std::vector<int> gatePos, double angles=0.0);
    void add_gate(std::string gateType, std::vector<int> gatePos, std::vector<double> angles={});
    std::unordered_map<std::string, int> measure(std::vector<int> measure_bits_list = {}, int shots =1024);
    };
#endif // QUANTUM_SIMULATOR_H