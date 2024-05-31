
#include "singlegate_mat.h"

SingleGate::SingleGate(std::string gateType)
{
    this->gateType = gateType;
}

SingleGate::SingleGate(std::string gateType, double angle)
{
    this->gateType = gateType;
    this->angles.push_back(angle);
}

SingleGate::SingleGate(std::string gateType, std::vector<double> angles)
{
    this->gateType = gateType;
    this->angles = angles;
}
Eigen::Matrix2cd SingleGate::hGate()
{
    Eigen::Matrix2cd hMat;
    hMat << 1, 1,
        1, -1;
    return hMat / std::sqrt(2.0);
}

Eigen::Matrix2cd SingleGate::xGate()
{
    Eigen::Matrix2cd xMat;
    xMat << 0, 1,
        1, 0;
    return xMat;
}

Eigen::Matrix2cd SingleGate::yGate()
{
    Eigen::Matrix2cd yMat;
    yMat << 0, -std::complex<double>(0, 1),
        std::complex<double>(0, 1), 0;
    return yMat;
}

Eigen::Matrix2cd SingleGate::zGate()
{
    Eigen::Matrix2cd zMat;
    zMat << 1, 0,
        0, -1;
    return zMat;
}

Eigen::Matrix2cd SingleGate::u3Gate()
{
    Eigen::Matrix2cd u3Mat;
    u3Mat << std::cos(angles[0] / 2.0), -std::exp(std::complex<double>(0, 1) * angles[2]) * std::sin(angles[0] / 2.0),
        std::exp(std::complex<double>(0, 1) * angles[1]) * std::sin(angles[0] / 2.0), std::exp(std::complex<double>(0, 1) * (angles[1] + angles[2])) * std::cos(angles[0] / 2.0);
    return u3Mat;
}

Eigen::Matrix2cd SingleGate::rxGate()
{
    Eigen::Matrix2cd rxMat;
    rxMat << std::cos(angles[0] / 2.0), -std::complex<double>(0, 1) * std::sin(angles[0] / 2.0),
        -std::complex<double>(0, 1) * std::sin(angles[0] / 2.0), std::cos(angles[0] / 2.0);
    return rxMat;
}

Eigen::Matrix2cd SingleGate::ryGate()
{
    Eigen::Matrix2cd ryMat;
    ryMat << std::cos(angles[0] / 2.0), -std::sin(angles[0] / 2.0),
        std::sin(angles[0] / 2.0), std::cos(angles[0] / 2.0);
    return ryMat;
}

Eigen::Matrix2cd SingleGate::rzGate()
{
    Eigen::Matrix2cd rzMat;
    rzMat << std::exp(-std::complex<double>(0, 1) * angles[0] / 2.0), 0,
        0, std::exp(std::complex<double>(0, 1) * angles[0] / 2.0);
    return rzMat;
}

Eigen::Matrix2cd SingleGate::gateMat()
{

    if (gateType == "h")
    {
        return hGate();
    }
    else if (gateType == "x")
    {
        return xGate();
    }
    else if (gateType == "y")
    {
        return yGate();
    }
    else if (gateType == "z")
    {
        return zGate();
    }
    else if (gateType == "rx")
    {
        return rxGate();
    }
    else if (gateType == "ry")
    {
        return ryGate();
    }
    else if (gateType == "rz")
    {
        return rzGate();
    }
    else if (gateType == "u3")
    {
        return u3Gate();
    }
    else
    {
        throw std::invalid_argument("Invalid gate type");
    }
}

/*  
{
    std::string str[] = {"x", "y", "z", "h"}; // Change char to std::string

    // 创建一个SingleGate对象
    SingleGate singleGate(str[2]); // Add a semicolon at the end

    // 调用gateMat方法并获取结果
    Eigen::Matrix2cd gateMatrix = singleGate.gateMat();

    // 打印结果
    std::cout << "Gate matrix:\n"
              << gateMatrix << std::endl;

    return 0;
}
*/