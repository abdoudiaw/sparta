
#include "adasRatesInterpolator.h"

static std::unordered_map<CacheKey, RateResults> rateCache;

std::vector<float> read1DDataSet(const H5::H5File& file, const std::string& datasetPath) {
    H5::DataSet ds = file.openDataSet(datasetPath);
    H5::DataSpace space = ds.getSpace();
    hsize_t dims[1];
    int ndims = space.getSimpleExtentDims(dims, NULL);
    if (ndims != 1) {
        throw std::runtime_error("Dataset " + datasetPath + " is not 1D.");
    }
    std::vector<float> data(dims[0]);
    ds.read(data.data(), H5::PredType::NATIVE_FLOAT);
    return data;
}

RateData readRateData(const std::string& filePath) {
    H5::H5File file(filePath, H5F_ACC_RDONLY);

    auto rawDataTo2DVector = [](const std::vector<float>& rawData, const std::vector<hsize_t>& dims) {
        std::vector<std::vector<float>> data(dims[0], std::vector<float>(dims[1]));
        for (hsize_t i = 0; i < dims[0]; ++i)
            for (hsize_t j = 0; j < dims[1]; ++j)
                data[i][j] = rawData[i * dims[1] + j];
        return data;
    };

    auto read2DDataSet = [&file, &rawDataTo2DVector](const std::string& datasetPath) {
        H5::DataSet ds = file.openDataSet(datasetPath);
        H5::DataSpace space = ds.getSpace();
        std::vector<hsize_t> dims(2);
        space.getSimpleExtentDims(dims.data(), NULL);

        int N = dims[0] * dims[1];
        std::vector<float> buffer(N);
        ds.read(buffer.data(), H5::PredType::NATIVE_FLOAT);

        return rawDataTo2DVector(buffer, dims); // Use buffer here instead of rawData

    };

    auto rawDataTo3DVector = [](const std::vector<float>& rawData, const std::vector<hsize_t>& dims) {
        std::vector<std::vector<std::vector<float>>> data(dims[0], std::vector<std::vector<float>>(dims[1], std::vector<float>(dims[2])));
        for (hsize_t i = 0; i < dims[0]; ++i)
            for (hsize_t j = 0; j < dims[1]; ++j)
                for (hsize_t k = 0; k < dims[2]; ++k)
                    data[i][j][k] = rawData[i * dims[1] * dims[2] + j * dims[2] + k];
        return data;
    };

    auto read3DDataSet = [&file, &rawDataTo3DVector](const std::string& datasetPath) {
        H5::DataSet ds = file.openDataSet(datasetPath);
        H5::DataSpace space = ds.getSpace();
        std::vector<hsize_t> dims(3);
        space.getSimpleExtentDims(dims.data(), NULL);

        std::vector<float> rawData(dims[0] * dims[1] * dims[2]);
        ds.read(rawData.data(), H5::PredType::NATIVE_FLOAT);

        return rawDataTo3DVector(rawData, dims);
    };

    RateData data;
    data.IonizationRateCoeff = read3DDataSet("IonizationRateCoeff");
    data.RecombinationRateCoeff = read3DDataSet("RecombinationRateCoeff");
    data.gridChargeState_Ionization = read2DDataSet("gridChargeState_Ionization");
    data.gridChargeState_Recombination = read2DDataSet("gridChargeState_Recombination");

    data.Atomic_Number = read1DDataSet(file, "Atomic_Number");
    data.gridDensity_Ionization = read1DDataSet(file, "gridDensity_Ionization");
    data.gridDensity_Recombination = read1DDataSet(file, "gridDensity_Recombination");
    data.gridTemperature_Ionization = read1DDataSet(file, "gridTemperature_Ionization");
    data.gridTemperature_Recombination = read1DDataSet(file, "gridTemperature_Recombination");

    return data;
}


double computeRate(const std::vector<float>& rateGrid_Temp,
                       const std::vector<float>& rateGrid_Dens,
                       const std::vector<std::vector<float>>& Rates,
                       double te, double ne, std::string material) {

    double logT = std::log10(te);
    double logn = std::log10(ne);

    double d_T = rateGrid_Temp[1] - rateGrid_Temp[0];
    double d_n = rateGrid_Dens[1] - rateGrid_Dens[0];

    int nT = rateGrid_Temp.size();
    int nD = rateGrid_Dens.size();

    int indT_temp = static_cast<int>(std::floor((logT - rateGrid_Temp[0]) / d_T));
    int indT = (indT_temp < 0) ? 0 : ((indT_temp > nT - 2) ? nT - 2 : indT_temp);

    int indN_temp = static_cast<int>(std::floor((logn - rateGrid_Dens[0]) / d_n));
    int indN = (indN_temp < 0) ? 0 : ((indN_temp > nD - 2) ? nD - 2 : indN_temp);

    double tempIndT1 = std::pow(10.0, rateGrid_Temp[indT + 1]);
    double tempIndT = std::pow(10.0, rateGrid_Temp[indT]);
    double densIndN1 = std::pow(10.0, rateGrid_Dens[indN + 1]);
    double densIndN = std::pow(10.0, rateGrid_Dens[indN]);

    double aT = tempIndT1 - te;
    double bT = te - tempIndT;
    double aN = densIndN1 - ne;
    double bN = ne - densIndN;

    double abT = aT + bT;
    double abN = aN + bN;

    double fx_z1 = (aN * std::pow(10.0, Rates[indT][indN]) +
                    bN * std::pow(10.0, Rates[indT][indN + 1])) / abN;
    double fx_z2 = (aN * std::pow(10.0, Rates[indT + 1][indN]) +
                    bN * std::pow(10.0, Rates[indT + 1][indN + 1])) / abN;

    return (aT * fx_z1 + bT * fx_z2) / abT;
}

RateResults interpolateRateData(int charge, double te, double ne, const RateData rateData, std::string material) {

    static std::unordered_map<CacheKey, RateResults> rateCache;

    CacheKey key{charge, te, ne, material};
    // printf("key: %d, %f, %f %s\n", key.charge, key.te, key.ne, key.material.c_str());

    // Check cache first
    auto it = rateCache.find(key);
    if (it != rateCache.end()) {
        return it->second;
    }

    // charge = (charge > 73) ? 0 : charge;

    double ionization_rates = computeRate(rateData.gridTemperature_Ionization,
                                          rateData.gridDensity_Ionization,
                                          rateData.IonizationRateCoeff[charge],
                                          te, ne, material);

    double recombination_rates = computeRate(rateData.gridTemperature_Recombination,
                                             rateData.gridDensity_Recombination,
                                             rateData.RecombinationRateCoeff[charge],
                                             te, ne, material);

    return RateResults{ionization_rates, recombination_rates};
}

static std::unordered_map<std::string, RateData> rateDataCache;

RateData readRateDataOnce(const std::string& filename) {
    auto it = rateDataCache.find(filename);
    if(it == rateDataCache.end()) {
        rateDataCache[filename] = readRateData(filename);
    }
    return rateDataCache[filename];
}