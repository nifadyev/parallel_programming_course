// Copyright Kiselev Denis 2019

#include <omp.h>

#include <cfloat>
#include <functional>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#define DEFAULT_NPOINTS 100000

// Ellipsoid options
#define ELLPS_A 1.0
#define ELLPS_B 2.0
#define ELLPS_C 3.0

// Limits of integration
#define X1 -1.0
#define X2 1.0
#define Y1 -2.0
#define Y2 2.0
#define Z1 -3.0
#define Z2 3.0

double ellipsoid(const std::vector<double>& args) {
    return -1.0 + pow(args[0] / ELLPS_A, 2)
                + pow(args[1] / ELLPS_B, 2)
                + pow(args[2] / ELLPS_C, 2);
}

double integrateByMonteCarlo(
        std::function<double(const std::vector<double>&)> func,
        std::vector<std::pair<double, double> > limits, int nPoints) {
    int dimension = limits.size();
    std::vector<std::uniform_real_distribution<> > distrs;
    distrs.reserve(dimension);

    double measure = 1.0;
    for (int i = 0; i < dimension; i++) {
        measure *= limits[i].second - limits[i].first;
        distrs.emplace_back(limits[i].first, limits[i].second);
    }

    int nPointsInEllipsoid = 0;
    std::vector<double> args(dimension);
    std::random_device r;
    std::default_random_engine generator(r());
    for (int i = 0; i < nPoints; i++) {
        for (int j = 0; j < dimension; j++)
            args[j] = distrs[j](generator);
        double value = func(args);
        if (value <= 0) nPointsInEllipsoid++;
    }

    double hitProbability = static_cast<double>(nPointsInEllipsoid) / nPoints;
    return measure * hitProbability;
}

double integrateByMonteCarloParallel(
        std::function<double(std::vector<double>&)> func,
        std::vector<std::pair<double, double> > limits, int nPoints) {
    int dimension = limits.size();
    std::vector<std::uniform_real_distribution<> > distrs;
    distrs.reserve(dimension);

    double measure = 1.0;
    for (int i = 0; i < dimension; i++) {
        measure *= limits[i].second - limits[i].first;
        distrs.emplace_back(limits[i].first, limits[i].second);
    }

    int nPointsInEllipsoid = 0;
    #pragma omp parallel reduction(+: nPointsInEllipsoid)
    {
        std::vector<double> args(dimension);
        std::default_random_engine generator(omp_get_thread_num());
        #pragma omp for schedule(static)
        for (int i = 0; i < nPoints; i++) {
            for (int j = 0; j < dimension; j++)
                args[j] = distrs[j](generator);
            double value = func(args);
            if (value <= 0) nPointsInEllipsoid++;
        }
    }

    double hitProbability = static_cast<double>(nPointsInEllipsoid) / nPoints;
    return measure * hitProbability;
}

int main(int argc, char *argv[]) {
    int nPoints = (argc > 1) ? atoi(argv[1]) : DEFAULT_NPOINTS;

    // Sequential
    double t1 = omp_get_wtime();
    double seqResult = integrateByMonteCarlo(
        ellipsoid, { { X1, X2 }, { Y1, Y2 }, { Z1, Z2 } }, nPoints);
    double seqTime = omp_get_wtime() - t1;

    // Parallel
    t1 = omp_get_wtime();
    double parResult = integrateByMonteCarloParallel(
        ellipsoid, { { X1, X2 }, { Y1, Y2 }, { Z1, Z2 } }, nPoints);
    double parTime = omp_get_wtime() - t1;

    double realRes = 4.0 / 3.0 * std::acos(-1) * ELLPS_A * ELLPS_B * ELLPS_C;
    double speedUp = seqTime / parTime;
    std::cout << "Sequential alg:\n"
                 "\tResult: " << seqResult << "\n"
                 "\tTime: " << seqTime << " sec\n"
                 "Parallel alg:\n"
                 "\tResult: " << parResult << "\n"
                 "\tTime: " << parTime << " sec\n"
                 "Real result: " << realRes << "\n"
                 "Speed up: " << speedUp << "\n";
    return 0;
}
