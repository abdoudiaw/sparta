
#include "adasRatesInterpolator.h"
#include <vector>
#include <cmath>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>



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


RateResults interpolateRateDataIonization(double charge, double te, double ne, const RateData& rateData, std::string material) {
    size_t charge_idx = static_cast<size_t>(charge);


    std::vector<double> tempDouble(rateData.gridTemperature_Ionization.begin(), rateData.gridTemperature_Ionization.end());
    std::vector<double> densityDouble(rateData.gridDensity_Ionization.begin(), rateData.gridDensity_Ionization.end());

    std::vector<double> interpRatesAtTe(rateData.gridDensity_Ionization.size());

 // Check if 'te' and 'ne' are within the range of data. If not, return zero.
    if (te < tempDouble.front() || te > tempDouble.back() ||
        ne < densityDouble.front() || ne > densityDouble.back()) {
        return RateResults{0.0, 0.0};
    }

    gsl_interp_accel *acc_te = gsl_interp_accel_alloc();
    gsl_spline *spline_te = gsl_spline_alloc(gsl_interp_cspline, tempDouble.size());

    for(size_t i = 0; i < rateData.gridDensity_Ionization.size(); ++i) {
        std::vector<double> ratesForDensity(rateData.IonizationRateCoeff[charge_idx][i].begin(), rateData.IonizationRateCoeff[charge_idx][i].end());
        
        gsl_spline_init(spline_te, tempDouble.data(), ratesForDensity.data(), tempDouble.size());
        interpRatesAtTe[i] = gsl_spline_eval(spline_te, te, acc_te);
    }

    gsl_spline_free(spline_te);
    gsl_interp_accel_free(acc_te);

    // Interpolate over densities
    gsl_interp_accel *acc_ne = gsl_interp_accel_alloc();
    gsl_spline *spline_ne = gsl_spline_alloc(gsl_interp_cspline, densityDouble.size());
    gsl_spline_init(spline_ne, densityDouble.data(), interpRatesAtTe.data(), densityDouble.size());

    double rate_final = gsl_spline_eval(spline_ne, ne, acc_ne);

    gsl_spline_free(spline_ne);
    gsl_interp_accel_free(acc_ne);

    return RateResults{rate_final, 0.0};  // Assuming recombination isn't being interpolated here
}

RateResults interpolateRateDataRecombination(double charge, double te, double ne, const RateData& rateData, std::string material) {
    size_t charge_idx = static_cast<size_t>(charge);

    // Convert temperature and density grids once for recombination
    std::vector<double> tempDouble(rateData.gridTemperature_Recombination.begin(), rateData.gridTemperature_Recombination.end());
    std::vector<double> densityDouble(rateData.gridDensity_Recombination.begin(), rateData.gridDensity_Recombination.end());

    std::vector<double> interpRatesAtTe(rateData.gridDensity_Recombination.size());

 // Check if 'te' and 'ne' are within the range of data. If not, return zero.
    if (te < tempDouble.front() || te > tempDouble.back() ||
        ne < densityDouble.front() || ne > densityDouble.back()) {
        return RateResults{0.0, 0.0};
    }

    gsl_interp_accel *acc_te = gsl_interp_accel_alloc();
    gsl_spline *spline_te = gsl_spline_alloc(gsl_interp_cspline, tempDouble.size());

    for(size_t i = 0; i < rateData.gridDensity_Recombination.size(); ++i) {
        std::vector<double> ratesForDensity(rateData.RecombinationRateCoeff[charge_idx][i].begin(), rateData.RecombinationRateCoeff[charge_idx][i].end());
        
        gsl_spline_init(spline_te, tempDouble.data(), ratesForDensity.data(), tempDouble.size());
        interpRatesAtTe[i] = gsl_spline_eval(spline_te, te, acc_te);
    }

    gsl_spline_free(spline_te);
    gsl_interp_accel_free(acc_te);

    // Interpolate over densities for recombination
    gsl_interp_accel *acc_ne = gsl_interp_accel_alloc();
    gsl_spline *spline_ne = gsl_spline_alloc(gsl_interp_cspline, densityDouble.size());
    gsl_spline_init(spline_ne, densityDouble.data(), interpRatesAtTe.data(), densityDouble.size());

    double rate_final = gsl_spline_eval(spline_ne, ne, acc_ne);

    gsl_spline_free(spline_ne);
    gsl_interp_accel_free(acc_ne);

    return RateResults{0.0, rate_final};  // Now the recombination rate is being interpolated and returned
}

RateResults interpolateRates(double charge, double te, double ne, const RateData& rateData, std::string material) {

    double ionization_rate = 0.0;
    double recombination_rate = 0.0;
    double nuclearCharge = rateData.Atomic_Number[0];

    // Do not ionize fully ionized and do not recombine neutral
    if (charge < nuclearCharge) {
        ionization_rate = interpolateRateDataIonization(charge, te, ne, rateData, material).ionization;
    }

    if (1 < charge && charge <= nuclearCharge) {
        double tableIndex = charge - 1; // Note: This assumes that charge is 1-based, which seems to be the case from your description.
        recombination_rate = interpolateRateDataRecombination(tableIndex, te, ne, rateData, material).recombination;
    }

    return RateResults{ionization_rate, recombination_rate};
}


static std::unordered_map<std::string, RateData> rateDataCache;

RateData readRateDataOnce(const std::string& filename) {
    auto it = rateDataCache.find(filename);
    if(it == rateDataCache.end()) {
        rateDataCache[filename] = readRateData(filename);
    }
    return rateDataCache[filename];
}