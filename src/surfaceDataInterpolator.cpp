#include "surfaceDataInterpolator.h"


std::map<std::pair<float, float>, SurfaceDataParams> surfacedata_cache;
std::map<std::string, SurfaceData> surfaceDataCache;  // Define the variable here


SurfaceData readSurfaceData(const std::string& filePath) {
    // Check if data is already cached
    auto it = surfaceDataCache.find(filePath);
    if (it != surfaceDataCache.end()) {
        return it->second;
    }

    H5::H5File file(filePath, H5F_ACC_RDONLY);

    auto read1DDataSet = [&file](const std::string& datasetPath) {
        H5::DataSet ds = file.openDataSet(datasetPath);
        H5::DataSpace space = ds.getSpace();
        std::vector<hsize_t> dims(1);
        space.getSimpleExtentDims(dims.data(), NULL);
        
        std::vector<float> data(dims[0]);
        ds.read(data.data(), H5::PredType::NATIVE_FLOAT);
        return data;
    };

    auto read2DDataSet = [&file](const std::string& datasetPath) {
        H5::DataSet ds = file.openDataSet(datasetPath);
        H5::DataSpace space = ds.getSpace();
        std::vector<hsize_t> dims(2);
        space.getSimpleExtentDims(dims.data(), NULL);

        std::vector<std::vector<float>> data(dims[0], std::vector<float>(dims[1]));
        std::vector<float> rawData(dims[0] * dims[1]);
        ds.read(rawData.data(), H5::PredType::NATIVE_FLOAT);
        
        for (hsize_t i = 0; i < dims[0]; ++i)
            for (hsize_t j = 0; j < dims[1]; ++j)
                data[i][j] = rawData[i * dims[1] + j];

        return data;
    };

    SurfaceData data;
    data.E = read1DDataSet("E");
    data.A = read1DDataSet("A");
    data.rpyld = read2DDataSet("rfyld");
    data.spyld = read2DDataSet("spyld");

    file.close();

    // Cache and return the read data
    surfaceDataCache[filePath] = data;
    return data;
}

SurfaceDataParams interpolateSurfaceData(const SurfaceData& data, float energy_val, float angle_val) {

    std::pair<float, float> key = {energy_val, angle_val};
    if (surfacedata_cache.find(key) != surfacedata_cache.end()) {
        return surfacedata_cache[key];
    }

    int i_energy = -1, j_energy = -1, i_angle = -1, j_angle = -1;

    // Locate the i_energy, j_energy for energy_val
    for (int i = 0; i < data.E.size() - 1; i++) {
        if (data.E[i] <= energy_val && energy_val < data.E[i+1]) {
            i_energy = i;
            j_energy = i + 1;  // Next index
            break;
        }
    }

    // Locate the i_angle, j_angle for angle_val
    for (int i = 0; i < data.A.size() - 1; i++) {
        if (data.A[i] <= angle_val && angle_val < data.A[i+1]) {
            i_angle = i;
            j_angle = i + 1;  // Next index
            break;
        }
    }

    // Check if either coordinate was not found
    if (i_energy == -1 || i_angle == -1) {
        SurfaceDataParams params = {0, 0};  // Assuming a simple structure, adapt as needed
        return params;
    }

    // Calculate Interpolation Weights
    float alpha = (energy_val - data.E[i_energy]) / (data.E[j_energy] - data.E[i_energy]);
    float beta = (angle_val - data.A[i_angle]) / (data.A[j_angle] - data.A[i_angle]);

    // Interpolate Each Parameter (simplified for clarity)
    SurfaceDataParams params;
    params.reflectionCoeff = (1 - alpha) * (1 - beta) * data.rpyld[i_energy][i_angle] +
                              alpha * (1 - beta) * data.rpyld[j_energy][i_angle] +
                              (1 - alpha) * beta * data.rpyld[i_energy][j_angle] +
                              alpha * beta * data.rpyld[j_energy][j_angle];

    params.sputteringCoeff = (1 - alpha) * (1 - beta) * data.spyld[i_energy][i_angle] +
                              alpha * (1 - beta) * data.spyld[j_energy][i_angle] +
                              (1 - alpha) * beta * data.spyld[i_energy][j_angle] +
                              alpha * beta * data.spyld[j_energy][j_angle];

    surfacedata_cache[key] = params;
    return params;
}
