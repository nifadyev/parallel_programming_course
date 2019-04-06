// Copyright 2019 Pavel Vorobev

#include <omp.h>
#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <algorithm>
#include <ctime>

struct point {
    int x;
    int y;
};

bool moreOrEqualThan(point a, point c, point b) {
    return a.x * b.y - a.y * b.x
        + b.x * c.y - b.y * c.x
        + c.x * a.y - c.y * a.x <= 0;
}

bool notEquals(point a, point b) {
    if (a.x != b.x)
        return true;
    else if (a.y != b.y)
        return true;
    else
        return false;
}

int determinant(point a, point b, point c) {
    return (b.x - a.x) * (c.y - a.y)
        - (c.x - a.x) * (b.y - a.y);
}

int partition(std::vector<point> *a, int low, int high) {
    point pivot = a->at(high);
    int i = (low - 1);
    for (int j = low; j < high; j++) {
        if (moreOrEqualThan(pivot, a->at(j), a->at(0))) {
            i++;
            point tmp = a->at(i);
            a->at(i) = a->at(j);
            a->at(j) = tmp;
        }
    }
    point tmp = a->at(i + 1);
    a->at(i + 1) = a->at(high);
    a->at(high) = tmp;

    return i + 1;
}

void quickSort(std::vector<point> *a, int low, int high) {
    if (low < high) {
        int partIndex = partition(a, low, high);
        quickSort(a, low, partIndex - 1);
        quickSort(a, partIndex + 1, high);
    }
}

int grahamScan(std::vector<point> *a) {
    point c = a->at(0);
    int m = 0;
    for (size_t i = 1; i < a->size(); i++) {
        if (a->at(i).x < c.x) {
            c = a->at(i);
            m = i;
        } else if (a->at(i).x == c.x) {
            if (a->at(i).y < c.y) {
                c = a->at(i);
                m = i;
            }
        }
    }

    point tmp = a->at(0);
    a->at(0) = a->at(m);
    a->at(m) = tmp;
    m = 1;
    quickSort(a, 0, a->size() - 1);

    for (size_t i = 1; i < a->size(); i++) {
        if (notEquals(a->at(i), a->at(m))) {
             if (m >= 1) {
                 while (m >= 1 && determinant(a->at(m - 1), a->at(m),
                    a->at(i)) >= 0) {
                    m = m - 1;
               }
            }
            m = m + 1;
            point tmp = a->at(m);
            a->at(m) = a->at(i);
            a->at(i) = tmp;
        }
    }
    return m + 1;
}

char* getCmdOption(char **begin, char **end, const std::string& option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
        return *itr;
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

int main(int argc, char *argv[]) {
    std::srand((unsigned)time(NULL));
    double numOfPoints = 10000;
    double min = 1;
    double max = 900;
    int nthreads;
    double t1, t2, dt;
    std::vector<point>* data = new std::vector<point>;
    std::vector<point>* result = new std::vector<point>;
    std::vector<point>* singleData = new std::vector<point>;
    if (cmdOptionExists(argv, argv + argc, "-num")) {
        char *wcount = getCmdOption(argv, argv + argc, "-num");
        numOfPoints = atof(wcount);
    }

    if (cmdOptionExists(argv, argv + argc, "-min")) {
        char *wcount = getCmdOption(argv, argv + argc, "-min");
        min = atof(wcount);
    }

    if (cmdOptionExists(argv, argv + argc, "-max")) {
        char *wcount = getCmdOption(argv, argv + argc, "-max");
        max = atof(wcount);
    }

    for (size_t i = 0; i < numOfPoints; i++) {
        point p;
        p.x = static_cast<int> (min) + std::rand() %
               (static_cast<int> (max) - static_cast<int> (min) + 1);
        p.y = static_cast<int> (min) + std::rand() %
               (static_cast<int> (max) - static_cast<int> (min) + 1);
        data->push_back(p);
        singleData->push_back(p);
    }

    t1 = omp_get_wtime();
    int s = grahamScan(singleData);
    std::cout << "Hull!" << std::endl;

    for (int i = 0; i < s; i++) {
         std::cout << "(" << singleData->at(i).x << " ,"
             << singleData->at(i).y << " ) ";
    }

    t2 = omp_get_wtime();
    dt = t2 - t1;
    std::cout << std::endl << "time: " << dt;
    t1 = omp_get_wtime();

    #pragma omp parallel shared(data, result)
    {
        std::vector<point>* localData = new std::vector<point>;
        int tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
        int currIndex = tid * (data->size() / nthreads);
        int nextIndex = (tid + 1) * (data->size() / nthreads);
        if (tid == nthreads - 1) {
            for (size_t i = currIndex; i < data->size(); i++) {
                localData->push_back(data->at(i));
            }
        } else {
            for (int i = currIndex; i < nextIndex; i++) {
                localData->push_back(data->at(i));
            }
        }

        int hullSize = grahamScan(localData);
        for (int i = 0; i < hullSize; i++) {
            result->push_back(localData->at(i));
        }
        #pragma omp barrier
        #pragma omp single
        {
            std::cout << std::endl;
            int size = grahamScan(result);
            for (int i = 0; i < size; i++) {
                std::cout << "(" << result->at(i).x << " ,"
                     << result->at(i).y << " ) ";
            }
        }
        localData->clear();
    }
    t2 = omp_get_wtime();
    dt = t2 - t1;
    std::cout << std::endl << "parallel time: " << dt << std::endl;
    data->clear();
}

