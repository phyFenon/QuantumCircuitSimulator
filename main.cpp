#include <iostream>
#include <chrono>
#include "quantum_circuit/quantum_simulator.h"
#include <fstream>
#include <sstream>
#include <omp.h>
//#include <yaml-cpp/yaml.h> 目前mac电脑安装yaml-cpp编译链接出现问题。

void printWavefunc(std::vector<std::complex<double>> vector_print)
{
    std::cout << "波函数:" << std::endl;
    std::cout << "[ ";

    for (const auto& wavefunc : vector_print)
    {
        std::cout << wavefunc << " ";
    }
    std::cout << "]";
    std::cout << std::endl;
}

void printFinalResult(std::unordered_map<std::string, int> res)
{
    std::cout << "测量结果:" << std::endl;
    for (const auto& [state, count] : res)
    {
        std::cout << " (" << state << ": " << count << ")"
                  << " ";
    }
    std::cout << std::endl;
}
 
void input_txtfile(std::string filename, int num_qubits, int n_shots) {
    QuantumSimulator simulator(num_qubits);
    std::unordered_map<std::string, int> res;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string gate_type;
        std::vector<int> gate_pos;
        std::vector<double> gate_angles;

        iss >> gate_type;
        if (gate_type == "m") {
            std::string pos_str;
            while (iss >> pos_str) {
                int pos = std::stoi(pos_str.substr(1)); // 提取出 'q' 后面的数字
                gate_pos.push_back(pos);
            }

            if (gate_pos.empty()) {
                res = simulator.measure({}, n_shots);
            } else {
                res = simulator.measure(gate_pos, n_shots);
            }
        } else {
            std::string pos_str;
            while (iss >> pos_str && pos_str[0] == 'q') {
                int pos = std::stoi(pos_str.substr(1)); // 提取出 'q' 后面的数字
                gate_pos.push_back(pos);
            }

            // 如果读取的不是 'q' 开头的字符串，那么它应该是角度参数
            if (!pos_str.empty() && pos_str[0] != 'q') {
                double angle = std::stod(pos_str);
                gate_angles.push_back(angle);
            }

            // 继续读取剩余的角度参数
            double angle;
            while (iss >> angle) {
                gate_angles.push_back(angle);
            }

            simulator.add_gate(gate_type, gate_pos, gate_angles);
        }
    }
    //printWavefunc(simulator.psi);  // Print the final wavefunction
    printFinalResult(res);          // Print the measurement result
      
}


 /*void runQuantumExperiment(const std::string& configFile, int num_qubits, int n_shots)
    {
        QuantumSimulator simulator(num_qubits); // Create a quantum simulator with the specified number of qubits

        YAML::Node config = YAML::LoadFile(configFile);
        const auto& steps = config["steps"];
        for (const auto& step : steps)
        {
            const auto& gates = step["gates"];
            for (const auto& gate : gates)
            {
                std::string gate_type = gate["name"].as<std::string>();
                std::vector<int> gate_pos = gate["targets"].as<std::vector<int>>();
                std::vector<double> gate_angles;
                if (gate["theta"])
                {
                    gate_angles.push_back(gate["theta"].as<double>());
                }

                simulator.add_gate(gate_type, gate_pos, gate_angles);
            }
        }

        std::vector<int> measure_positions = config["measure-position"].as<std::vector<int>>();
        std::unordered_map<std::string, int> res;
        if (measure_positions.empty())
        {
            res = simulator.measure({}, n_shots);
        }
        else
        {
            res = simulator.measure(measure_positions, n_shots);
        }

        printFinalResult(res);          // Print the measurement result
        printWavefunc(simulator.psi);   // Print the final wavefunction
    }*/

    int main()
    {
        auto start = std::chrono::high_resolution_clock::now();
        int num_qubits = 25;
        int n_shots = 10000;
        //std::string configFile = "h_test.txt";
        //input_txtfile(configFile, num_qubits, n_shots);//传入一个txt文件，文件中包含量子门信息

        QuantumSimulator simulator(num_qubits);//construct the quantum circuit.


        for (int i=0; i<num_qubits;i++)
        {
            simulator.add_gate("h",{i});
           // std::cout<<i<<std::endl;
       }

        /*for (int i=0; i<num_qubits-1;i++)
        {
            simulator.add_gate("cx",{i,i+1});
            simulator.add_gate("rxx",{i,i+1},{0.4});
            simulator.add_gate("cz",{i,i+1});
        }
        for (int i=0; i<num_qubits-2;i++)
        {
            simulator.add_gate("ccx",{i,i+1,i+2});
        }
        */
        std::unordered_map<std::string, int> res = simulator.measure({0,1,2,3},n_shots);
        printFinalResult(res);
        auto end_m = std::chrono::high_resolution_clock::now();
        auto duration_m = std::chrono::duration_cast<std::chrono::seconds>(end_m- start);
        std::cout << "用时:" << duration_m.count() << " seconds" << std::endl;
        return 0;
    }