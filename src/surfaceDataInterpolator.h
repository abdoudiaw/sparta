#ifndef SURFACEDATAINTERPOLATOR_H
#define SURFACEDATAINTERPOLATOR_H

#include <H5Cpp.h>
#include <vector>
#include <map>
#include <string>

// Surface Data
struct SurfaceData {
    std::vector<float> E;
    std::vector<float> A;
    std::vector<std::vector<float>> rpyld;
    std::vector<std::vector<float>> spyld;
};

// Surface Data Parameters
struct SurfaceDataParams {
    float E, A, reflectionCoeff, sputteringCoeff;
};

// Function declarations
SurfaceData readSurfaceData(const std::string& filePath);
SurfaceDataParams interpolateSurfaceData(const SurfaceData& data, float energy_val, float angle_val);

#endif // SURFACEDATAINTERPOLATOR_H
