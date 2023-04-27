#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>


struct Data {
    double R, Z, b_phi, b_r, b_z, T_e, n_e, v_e, t_i, n_i, v_i;
};

std::vector<Data> readDataFromFile(const std::string& filename) {
    std::vector<Data> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        Data d;
        std::istringstream iss(line);
        iss >> d.R >> d.Z >> d.b_phi >> d.b_r >> d.b_z >> d.T_e >> d.n_e >> d.v_e >> d.t_i >> d.n_i >> d.v_i;

        data.push_back(d);
    }

    file.close();
    return data;
}

double lerp(double a, double b, double t) {
    return a + (b - a) * t;
}

double bilinearInterpolation(double x, double y,
                             double f00, double f10,
                             double f01, double f11) {
    double t1 = lerp(f00, f10, x);
    double t2 = lerp(f01, f11, x);
    return lerp(t1, t2, y);
}

Data bilinearInterpolation(const std::vector<Data>& data, double r, double z) {
    Data d00, d01, d10, d11;
    double minDist = std::numeric_limits<double>::max();

    for (const auto& d : data) {
        double dist = std::sqrt(std::pow(r - d.R, 2) + std::pow(z - d.Z, 2));
        if (dist < minDist) {
            d11 = d10;
            d10 = d01;
            d01 = d00;
            d00 = d;
            minDist = dist;
        }
    }

    double x = (r - d00.R) / (d10.R - d00.R);
    double y = (z - d00.Z) / (d01.Z - d00.Z);

    Data result;
    result.R = r;
    result.Z = z;
    result.b_phi = bilinearInterpolation(x, y, d00.b_phi, d10.b_phi, d01.b_phi, d11.b_phi);
    result.b_r = bilinearInterpolation(x, y, d00.b_r, d10.b_r, d01.b_r, d11.b_r);
    result.b_z = bilinearInterpolation(x, y, d00.b_z, d10.b_z, d01.b_z, d11.b_z);
    result.T_e = bilinearInterpolation(x, y, d00.T_e, d10.T_e, d01.T_e, d11.T_e);
    result.n_e = bilinearInterpolation(x, y, d00.n_e, d10.n_e, d01.n_e, d11.n_e);
    result.v_e = bilinearInterpolation(x, y, d00.v_e, d10.v_e, d01.v_e, d11.v_e);
    result.t_i = bilinearInterpolation(x, y, d00.t_i, d10.t_i, d01.t_i, d11.t_i);
    result.n_i = bilinearInterpolation(x, y, d00.n_i, d10.n_i, d01.n_i, d11.n_i);
    result.v_i = bilinearInterpolation(x, y, d00.v_i, d10.v_i, d01.v_i, d11.v_i);

    return result;
}

int main() {
    const std::string filename = "data.txt";
    double r = 2.8;
    double z = -0.5;

    std::vector<Data> data = readDataFromFile(filename);
    Data interpolatedData = bilinearInterpolation(data, r, z);

    std::cout << "Interpolated values at (R, Z) = (" << r << ", " << z << "):" << std::endl;
    std::cout << "b_phi: " << interpolatedData.b_phi << std::endl;
    std::cout << "b_r: " << interpolatedData.b_r << std::endl;
    std::cout << "b_z: " << interpolatedData.b_z << std::endl;
    std::cout << "T_e: " << interpolatedData.T_e << std::endl;
    std::cout << "n_e: " << interpolatedData.n_e << std::endl;
    std::cout << "v_e: " << interpolatedData.v_e << std::endl;
    std::cout << "t_i: " << interpolatedData.t_i << std::endl;
    std::cout << "n_i: " << interpolatedData.n_i << std::endl;
    std::cout << "v_i: " << interpolatedData.v_i << std::endl;

    return 0;
}

