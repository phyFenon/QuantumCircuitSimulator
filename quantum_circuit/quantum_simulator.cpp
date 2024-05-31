#include "singlegate_mat.h"
#include "quantum_simulator.h"

QuantumSimulator::QuantumSimulator(int num_qubits)
    : num_qubits(num_qubits)
{
    int size = 1 << num_qubits;
    psi.resize(size, {0.0, 0.0});

    // 使用OpenMP并行化循环
    #pragma omp parallel for
    for (int i = 1; i < size; i++) { // 从1开始，因为psi[0]已经被单独设置
        psi[i] = {0.0, 0.0};
    }

    psi[0] = {1.0, 0.0}; // 单独设置psi[0]
}
 void QuantumSimulator::SingleGate_ActOn_State(Eigen::Matrix2cd gateMat, int gatePos)
{
    int i = 1 << gatePos;
    int nstates = (1 << num_qubits) - (1 << (gatePos + 1));

    #pragma omp parallel for collapse(2)
    for (int k = 0; k <= nstates; k += i * 2)
    {
        for (int l = 0; l < i; l++)
        {
            int i0 = l | k; // 直接计算i0
            int i1 = i0 | i; // 直接计算i1
            // 直接在psi上操作，避免使用temp
            std::complex<double> psi_i0 = psi[i0];
            std::complex<double> psi_i1 = psi[i1];
            psi[i0] = gateMat(0, 0) * psi_i0 + gateMat(0, 1) * psi_i1;
            psi[i1] = gateMat(1, 0) * psi_i0 + gateMat(1, 1) * psi_i1;
        }
    }
}

void QuantumSimulator::DoubleGate_ActOn_State(std::string gateType, std::vector<int> gatePos, double angle)
{
    int j1 = std::max(gatePos[0], gatePos[1]);
    int j0 = std::min(gatePos[0], gatePos[1]);

    int i1_plus = 1 << j0;
    int i2_plus = 1 << j1;
    int i3_plus = i1_plus + i2_plus;
    int delta2 = 1 << (j1 + 1);
    int delta1 = 1 << (j0 + 1);
    int max2 = (1 << num_qubits) - delta2;
    int max1 = (1 << j1) - delta1;
    int max0 = (1 << j0) - 1;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k <= max2; k += delta2)
    {
        for (int l = 0; l <= max1; l += delta1)
        {
            for (int m = 0; m <= max0; m++)
            {
                int i0 = m | l | k;
                int i1 = i0 | i1_plus;
                int i2 = i0 | i2_plus;
                int i3 = i0 | i3_plus;

                std::complex<double> psi_i0 = psi[i0];
                std::complex<double> psi_i1 = psi[i1];
                std::complex<double> psi_i2 = psi[i2];
                std::complex<double> psi_i3 = psi[i3];

                if (gateType == "cx")
                {
                    if (gatePos[0] > gatePos[1])
                    {
                        std::swap(psi_i2, psi_i3);
                    }
                    else
                    {
                        std::swap(psi_i1, psi_i3);
                    }
                }
                else if (gateType == "swp")
                {
                    std::swap(psi_i1, psi_i2);
                }
                else if (gateType == "iswp")
                {
                    psi_i1 *= std::complex<double>(0, 1);
                    psi_i2 *= std::complex<double>(0, 1);
                    std::swap(psi_i1, psi_i2);
                }
                else if (gateType == "cz")
                {
                    psi_i3 *= -1.0;
                }
                else if (gateType == "cp")
                {
                    psi_i3 *= std::exp(std::complex<double>(0, angle));
                }
                else if (gateType == "rxx") 
                {
                    std::complex<double> temp_i0 = std::cos(angle/2.0)*psi_i0 - std::complex<double>(0, 1)*std::sin(angle/2.0)*psi_i3;
                    std::complex<double> temp_i3 = std::cos(angle/2.0)*psi_i3 - std::complex<double>(0, 1)*std::sin(angle/2.0)*psi_i0;
                    psi_i0 = temp_i0;
                    psi_i3 = temp_i3;

                    std::complex<double> temp_i1 = std::cos(angle/2.0)*psi_i1 - std::complex<double>(0, 1)*std::sin(angle/2.0)*psi_i2;
                    std::complex<double> temp_i2 = std::cos(angle/2.0)*psi_i2 - std::complex<double>(0, 1)*std::sin(angle/2.0)*psi_i1;
                    psi_i1 = temp_i1;
                    psi_i2 = temp_i2;
                }
                else if (gateType == "ryy") 
                {
                    std::complex<double> temp_i0 = std::cos(angle/2.0)*psi_i0 + std::complex<double>(0, 1)*std::sin(angle/2.0)*psi_i3;
                    std::complex<double> temp_i3 = std::cos(angle/2.0)*psi_i3 + std::complex<double>(0, 1)*std::sin(angle/2.0)*psi_i0;
                    psi_i0 = temp_i0;
                    psi_i3 = temp_i3;

                    std::complex<double> temp_i1 = std::cos(angle/2.0)*psi_i1 - std::complex<double>(0, 1)*std::sin(angle/2.0)*psi_i2;
                    std::complex<double> temp_i2 = std::cos(angle/2.0)*psi_i2 - std::complex<double>(0, 1)*std::sin(angle/2.0)*psi_i1;
                    psi_i1 = temp_i1;
                    psi_i2 = temp_i2;
                }
                else if (gateType == "rzz")
                {
                    psi_i0 *= std::exp(std::complex<double>(0, -0.5 * angle));
                    psi_i1 *= std::exp(std::complex<double>(0, 0.5 * angle));
                    psi_i2 *= std::exp(std::complex<double>(0, 0.5 * angle));
                    psi_i3 *= std::exp(std::complex<double>(0, -0.5 * angle));
                }

                psi[i0] = psi_i0;
                psi[i1] = psi_i1;
                psi[i2] = psi_i2;
                psi[i3] = psi_i3;
            }
        }
    }
}


void QuantumSimulator::TripleGate_ActOn_State(std::string gateType, std::vector<int> gatePos, double angles)
{
    std::vector<int> sortedGatePos = gatePos;
    std::sort(sortedGatePos.begin(), sortedGatePos.end());

    int j2 = sortedGatePos[2];
    int j1 = sortedGatePos[1];
    int j0 = sortedGatePos[0];

    int i1_plus = 1 << j0;
    int i2_plus = 1 << j1;
    int i3_plus = (1 << j0) + (1 << j1);
    int i4_plus = 1 << j2;
    int i5_plus = (1 << j2) + (1 << j0);
    int i6_plus = (1 << j2) + (1 << j1);
    int i7_plus = (1 << j2) + (1 << j1) + (1 << j0);
    int delta3 = 1 << (j2 + 1);
    int delta2 = 1 << (j1 + 1);
    int delta1 = 1 << (j0 + 1);
    int max3 = (1 << num_qubits) - (1 << (j2 + 1));
    int max2 = (1 << j2) - (1 << (j1 + 1));
    int max1 = (1 << j1) - (1 << (j0 + 1));
    int max0 = (1 << j0) - 1;

    std::complex<double> temp_psi[(1 << num_qubits)]; // Temporary array to store updated psi values

    #pragma omp parallel for collapse(4) schedule(static)
    for (int k = 0; k <= max3; k += delta3)
    {
        for (int l = 0; l <= max2; l += delta2)
        {
            for (int m = 0; m <= max1; m += delta1)
            {
                for (int n = 0; n <= max0; n++)
                {
                    int i0 = m | l | k | n;
                    int i1 = i0 | i1_plus;
                    int i2 = i0 | i2_plus;
                    int i3 = i0 | i3_plus;
                    int i4 = i0 | i4_plus;
                    int i5 = i0 | i5_plus;
                    int i6 = i0 | i6_plus;
                    int i7 = i0 | i7_plus;

                    temp_psi[i0] = psi[i0];
                    temp_psi[i1] = psi[i1];
                    temp_psi[i2] = psi[i2];
                    temp_psi[i3] = psi[i3];
                    temp_psi[i4] = psi[i4];
                    temp_psi[i5] = psi[i5];
                    temp_psi[i6] = psi[i6];
                    temp_psi[i7] = psi[i7];
                }
            }
        }
    }

    #pragma omp parallel for collapse(4) schedule(static)
    for (int k = 0; k <= max3; k += delta3)
    {
        for (int l = 0; l <= max2; l += delta2)
        {
            for (int m = 0; m <= max1; m += delta1)
            {
                for (int n = 0; n <= max0; n++)
                {
                    int i0 = m | l | k | n;
                    int i1 = i0 | i1_plus;
                    int i2 = i0 | i2_plus;
                    int i3 = i0 | i3_plus;
                    int i4 = i0 | i4_plus;
                    int i5 = i0 | i5_plus;
                    int i6 = i0 | i6_plus;
                    int i7 = i0 | i7_plus;

                    std::complex<double> psi_i0 = temp_psi[i0];
                    std::complex<double> psi_i1 = temp_psi[i1];
                    std::complex<double> psi_i2 = temp_psi[i2];
                    std::complex<double> psi_i3 = temp_psi[i3];
                    std::complex<double> psi_i4 = temp_psi[i4];
                    std::complex<double> psi_i5 = temp_psi[i5];
                    std::complex<double> psi_i6 = temp_psi[i6];
                    std::complex<double> psi_i7 = temp_psi[i7];

                    if (gatePos[0] < gatePos[2] && gatePos[1] < gatePos[2])
                    {
                        std::swap(psi_i3, psi_i7);
                    }
                    else if (gatePos[2] < gatePos[0] && gatePos[2] < gatePos[1])
                    {
                        std::swap(psi_i6, psi_i7);
                    }
                    else
                    {
                        std::swap(psi_i5, psi_i7);
                    }

                    psi[i0] = psi_i0;
                    psi[i1] = psi_i1;
                    psi[i2] = psi_i2;
                    psi[i3] = psi_i3;
                    psi[i4] = psi_i4;
                    psi[i5] = psi_i5;
                    psi[i6] = psi_i6;
                    psi[i7] = psi_i7;
                }
            }
        }
    }
}

    
    void QuantumSimulator::add_gate(std::string gateType, std::vector<int> gatePos, std::vector<double> angles)
    {    
        //将输入的字符串转换为小写。
        std::transform(gateType.begin(), gateType.end(), gateType.begin(),
                   [](unsigned char c) { return std::tolower(c); });

        // 单量子门类
        if (std::find(singleGateList.begin(), singleGateList.end(), gateType) != singleGateList.end())
        {
            if (gatePos.size() == 1)
            {
                SingleGate *singlegate = nullptr; // 在所有块之外声明 singlegate 指针
                if (gateType == "rx" || gateType == "ry" || gateType == "rz")
                {
                    if (angles.size() != 1)
                    {
                        throw std::invalid_argument("单比特旋转门需要一个角度参数");
                    }
                    else
                    {
                        singlegate = new SingleGate(gateType, angles[0]); // 初始化 singlegate
                    }
                }
                else if (gateType == "u3")
                {
                    if (angles.size() != 3)
                    {
                        throw std::invalid_argument("u3单比特门需要三个角度参数");
                    }
                    else
                    {
                        singlegate = new SingleGate(gateType, angles); // 初始化 singlegate
                    }
                }
                else
                {
                    angles.clear();                        // 对于 'x', 'y', 'z', 'h'，不需要 angles 参数
                    singlegate = new SingleGate(gateType); // 初始化 singlegate
                }
                Eigen::Matrix2cd gateMat = singlegate->gateMat();
                SingleGate_ActOn_State(gateMat, gatePos[0]);
                delete singlegate; // 释放 singlegate 所占用的内存
            }
            else
            {
                throw std::invalid_argument("单比特门作用位数为1");
            }
        }
        // 双比特量子门
        else if (std::find(doubleGateList.begin(), doubleGateList.end(), gateType) != doubleGateList.end())
        {
            if (gatePos.size() == 2 && gatePos[0] != gatePos[1])
            {
                //"cp", "rxx", "ryy", "rzz"
                if (gateType == "rxx" || gateType == "ryy" || gateType == "rzz" || gateType == "cp")
                {
                    if (angles.size() != 1)
                    {
                        throw std::invalid_argument("两比特门需要一个角度参数");
                    }
                    else
                    {
                        DoubleGate_ActOn_State(gateType, gatePos, angles[0]);
                    }  
                }  
                else
                {
                // 对于 'cx', 'cz'
                DoubleGate_ActOn_State(gateType, gatePos);
                } 
            }
            else
            {
            throw std::invalid_argument("两比特门作用位数为2且位数不相同。");
            }
        }
        // 三比特量子门
        else if (std::find(tripleGateList.begin(), tripleGateList.end(), gateType) != tripleGateList.end())
        {
            if (gatePos.size() == 3 && gatePos[0] != gatePos[1] && gatePos[1] != gatePos[2])
            {
                if (gateType == "ccx")
                {
                    if (angles.size() != 0)
                        throw std::invalid_argument("ccx不需要添加角度参数！");
                }
                else
                {
                    TripleGate_ActOn_State(gateType, gatePos); // 初始化 singlegate
                }
            }
            else
            {
                throw std::invalid_argument("三比特门作用位数为3且位数不相同。");
            }
        }
        else
        {
            // 其他的输出无效量子门
            std::cout << "Invalid gate:" << gateType << std::endl;
        }
    }
//测量函数
/*
std::unordered_map<std::string, int> QuantumSimulator::measure(std::vector<int> measure_bits_list, int shots) {
    std::cout << "Sampling begins" << std::endl;
    std::unordered_map<std::string, double> dictionary; // 用于存储所有字符串的测量结果
    std::unordered_map<std::string, int> result;

    // 并行计算概率分布
    #pragma omp parallel for reduction(+:dictionary[:])
    for (size_t i = 0; i < psi.size(); ++i) {
        double probability = std::norm(psi[i]);
        if (probability > 0.0) {
            std::string binaryString = std::bitset<32>(i).to_string().substr(32 - num_qubits);
            #pragma omp critical
            dictionary[binaryString] += probability;
        }
    }

    // 处理测量位
    if (!measure_bits_list.empty()) {
        std::unordered_map<std::string, double> new_dictionary;
        for (const auto& pair : dictionary) {
            std::string key;
            for (int pos : measure_bits_list) {
                key += pair.first[pos];
            }
            new_dictionary[key] += pair.second;
        }
        dictionary.swap(new_dictionary);
    }

    // 创建累积概率分布
    std::vector<std::pair<std::string, double>> cumulative;
    double cumulate = 0.0;
    for (const auto& pair : dictionary) {
        cumulate += pair.second;
        cumulative.emplace_back(pair.first, cumulate);
    }

    // 并行生成随机数并采样
    #pragma omp parallel
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        #pragma omp for reduction(+:result[:])
        for (int i = 0; i < shots; ++i) {
            double sample = dis(gen);
            for (const auto& pair : cumulative) {
                if (sample <= pair.second) {
                    #pragma omp critical
                    result[pair.first]++;
                    break;
                }
            }
        }
    }
    
    return result;
}
*/
std::unordered_map<std::string, int> QuantumSimulator::measure(std::vector<int> measure_bits_list, int shots) {
    std::cout << "Sampling begins" << std::endl;
    std::vector<std::unordered_map<std::string, double>> partialDictionaries(omp_get_max_threads());

    // Parallel computation of probability distribution
    #pragma omp parallel for
    for (size_t i = 0; i < psi.size(); ++i) {
        double probability = std::norm(psi[i]);
        if (probability > 0.0) {
            std::string binaryString = std::bitset<32>(i).to_string().substr(32 - num_qubits);
            int tid = omp_get_thread_num();
            #pragma omp critical
            partialDictionaries[tid][binaryString] += probability;
        }
    }

    // Merge partial dictionaries
    std::unordered_map<std::string, double> dictionary;
    for (const auto& pd : partialDictionaries) {
        for (const auto& pair : pd) {
            dictionary[pair.first] += pair.second;
        }
    }

    // Adjust dictionary according to measure_bits_list
    // This part remains unchanged
    // 处理测量位
    if (!measure_bits_list.empty()) {
        std::unordered_map<std::string, double> new_dictionary;
        for (const auto& pair : dictionary) {
            std::string key;
            for (int pos : measure_bits_list) {
                key += pair.first[pos];
            }
            new_dictionary[key] += pair.second;
        }
        dictionary.swap(new_dictionary);
    }
    // Prepare for sampling
    std::vector<std::pair<std::string, double>> cumulative;
    double cumulate = 0.0;
    for (const auto& pair : dictionary) {
        cumulate += pair.second;
        cumulative.emplace_back(pair.first, cumulate);
    }

    std::vector<std::unordered_map<std::string, int>> partialResults(omp_get_max_threads());

    // Parallel sampling
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::mt19937 gen(omp_get_thread_num());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        #pragma omp for nowait
        for (int i = 0; i < shots; ++i) {
            double sample = dis(gen);
            for (const auto& pair : cumulative) {
                if (sample <= pair.second) {
                    partialResults[tid][pair.first]++;
                    break;
                }
            }
        }
    }

    // Merge partial results
    std::unordered_map<std::string, int> result;
    for (const auto& pr : partialResults) {
        for (const auto& pair : pr) {
            result[pair.first] += pair.second;
        }
    }

    return result;
}