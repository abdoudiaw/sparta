#ifndef ADASRATESINTERPOLATOR_H
#define ADASRATESINTERPOLATOR_H

#include <iostream>
#include <H5Cpp.h>
#include <vector>
#include <map>
#include <tuple>

struct CacheKey {
    int charge;
    double te;
    double ne;
    std::string material;

    bool operator==(const CacheKey& other) const {
        return charge == other.charge && te == other.te && ne == other.ne && material == other.material;
    }
};

namespace std {
    template <>
    struct hash<CacheKey> {
        size_t operator()(const CacheKey& key) const {
            return std::hash<int>()(key.charge) ^
                   std::hash<double>()(key.te) ^
                   std::hash<double>()(key.ne) ^
                     std::hash<std::string>()(key.material);
        }
    };
}

struct RateResults {
    double ionization;
    double recombination;
};

struct RateData {
    std::vector<float> Atomic_Number;
    std::vector<std::vector<std::vector<float>>> IonizationRateCoeff, RecombinationRateCoeff;
    std::vector<std::vector<float>> gridChargeState_Ionization, gridChargeState_Recombination;
    std::vector<float> gridDensity_Ionization, gridDensity_Recombination, gridTemperature_Ionization, gridTemperature_Recombination;

};

std::vector<float> read1DDataSet(const H5::H5File& file, const std::string& datasetPath);

RateData readRateData(const std::string& filePath);
RateData readRateDataOnce(const std::string& filename);

double computeRate(const std::vector<float>& rateGrid_Temp,
                       const std::vector<float>& rateGrid_Dens,
                       const std::vector<std::vector<float>>& Rates,
                       double te, double ne, std::string material);

RateResults interpolateRateData(int charge, double te, double ne, const RateData rateData, std::string material);

#endif // ADASRATESINTERPOLATOR_H
