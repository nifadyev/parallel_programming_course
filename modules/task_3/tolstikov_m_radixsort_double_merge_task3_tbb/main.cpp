// Copyright 2019 Tolstikov Maksim
// RadixSort with simple merge (double)
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <tbb/tbb.h>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <cstring>
#include <climits>

void CountingSort(double *inp, double *out, int byteNum, int size) {
    unsigned char *mas = (unsigned char *)inp;
    int counter[256];
    int tem;
    memset(counter, 0, sizeof(int) * 256);
    for (int i = 0; i < size; i++)
        counter[mas[8 * i + byteNum]]++;
    int j = 0;
    for (; j < 256; j++) {
        if (counter[j] != 0)
            break;
    }
    tem = counter[j];
    counter[j] = 0;
    j++;
    for (; j < 256; j++) {
        int b = counter[j];
        counter[j] = tem;
        tem += b;
    }
    for (int i = 0; i < size; i++) {
        out[counter[mas[8 * i + byteNum]]] = inp[i];
        counter[mas[8 * i + byteNum]]++;
    }
}

void LastCountingSort(double *inp, double *out, int size) {
    unsigned char *mas = (unsigned char *)inp;
    int counter[256];
    int tem;
    memset(counter, 0, sizeof(int) * 256);
    for (int i = 0; i < size; i++)
        counter[mas[8 * i + 7] + 128]++;
    int j = 0;
    for (; j < 256; j++) {
        if (counter[j] != 0)
            break;
    }
    tem = counter[j];
    counter[j] = 0;
    j++;
    for (; j < 256; j++) {
        int b = counter[j];
        counter[j] = tem;
        tem += b;
    }
    for (int i = 0; i < size; i++) {
        out[counter[mas[8 * i + 7] + 128]] = inp[i];
        counter[mas[8 * i + 7] + 128]++;
    }
}

void merge(double* mas, int sizel, int sizer) {
    int size = sizel + sizer;
    double* tempMas = new double[size];
    int i = 0, j = sizel, k = 0;

    while (i != sizel && j != size) {
        if (mas[i] <= mas[j]) {
            tempMas[k] = mas[i];
            ++i;
            ++k;
        } else {
            tempMas[k] = mas[j];
            ++j;
            ++k;
        }
    }

    if (i < sizel) {
        for (; i < sizel; ++i) {
            tempMas[k] = mas[i];
            ++k;
        }
    }

    if (j < size) {
        for (; j < size; ++j) {
            tempMas[k] = mas[j];
            ++k;
        }
    }

    for (i = 0; i < size; ++i) {
        mas[i] = tempMas[i];
    }

    delete[] tempMas;
}

void PrintArray(double *array, int size) {
    if (size < 15) {
        for (int i = 0; i < size; i++) {
            std::cout << array[i] << " ";
        }
        std::cout << std::endl;
    }
}

void LSDSortDouble(double *inp, int size) {
    double *out = new double[size];
    CountingSort(inp, out, 0, size);
    CountingSort(out, inp, 1, size);
    CountingSort(inp, out, 2, size);
    CountingSort(out, inp, 3, size);
    CountingSort(inp, out, 4, size);
    CountingSort(out, inp, 5, size);
    CountingSort(inp, out, 6, size);
    LastCountingSort(out, inp, size);
    delete[] out;
}

void CopyArray(double *mas, double* tmp, int size) {
    for (int i = 0; i < size; i++)
        tmp[i] = mas[i];
}

void CheckingSort(double *mas, double* tmp, int size) {
    for (int i = 0; i < size; i++) {
        if (mas[i] != tmp[i]) {
            std::cout << "Sort is incorrectly" << std::endl;
            break;
        }
    }
    std::cout << "Sort is correctly" << std::endl;
}

void GenerateArray(double *mas, int size) {
    int b = 100;
    int a = 0;
    for (int i = 0; i < size; i++) {
        mas[i] = static_cast<double>(std::rand())*(b - a + 1) / RAND_MAX + a;
    }
}

int main(int argc, char *argv[]) {
    double time_lsd = 0;
    double ptime_lsd = 0;
    double merge_time = 0;
    double time_plsd = 0;
    int size = 10;
    int n = 1;
    std::srand((unsigned)time(NULL));
    double *mas, *tmp, *lmas;
    if (argc == 4) {
        n = atoi(argv[1]);
        if (strcmp(argv[2], "-size") == 0)
            size = atoi(argv[3]);
    }
    mas = new double[size];
    tmp = new double[size];
    lmas = new double[size];
    int tail = size % n;
    if (size < 15)
        std::cout << "Array: ";
    GenerateArray(mas, size);
    CopyArray(mas, tmp, size);
    std::sort(tmp, tmp + size);
    CopyArray(mas, lmas, size);
    if (mas == NULL) {
        std::cout << "Error! Incorrect input data for array";
        return -1;
    }
    PrintArray(mas, size);
    tbb::task_scheduler_init init(n);
    int grain_size = 0;
    grain_size = size / n;
    tbb::tick_count t1 = tbb::tick_count::now();
    tbb::tick_count t3 = tbb::tick_count::now();
    tbb::parallel_for(tbb::blocked_range<int>(0, n), [=, &mas](const tbb::blocked_range<int> &r) {
        for (int f = r.begin(); f != r.end(); f++) {
            if (f == 0) {
                LSDSortDouble(mas, grain_size + tail);
            } else {
                LSDSortDouble(mas + tail + f * grain_size, grain_size);
            }
        }
    });
    tbb::tick_count t4 = tbb::tick_count::now();
    init.terminate();
    int j = 1;
    int k = n / 2;
    tbb::tick_count t5 = tbb::tick_count::now();
    while (k != 0) {
        init.initialize(k);
        tbb::parallel_for(tbb::blocked_range<int>(0, k), [=, &mas](const tbb::blocked_range<int> &r) {
            for (int f = r.begin(); f != r.end(); f++) {
                if (f == 0) {
                    merge(mas, grain_size * j + tail, grain_size * j);
                } else {
                    merge(mas + 2 * f * grain_size * j + tail, grain_size * j, grain_size * j);
                }
            }
        });
        k /= 2;
        j *= 2;
        init.terminate();
    }
    tbb::tick_count t6 = tbb::tick_count::now();
    tbb::tick_count t2 = tbb::tick_count::now();
    ptime_lsd = (t2 -t1).seconds();
    tbb::tick_count t11 = tbb::tick_count::now();
    LSDSortDouble(lmas, size);
    tbb::tick_count t22 = tbb::tick_count::now();
    time_lsd = (t22 - t11).seconds();
    merge_time = (t6 - t5).seconds();
    time_plsd = (t4 - t3).seconds();
    if (size < 15)
        std::cout << "Array after LSD sort: ";
    PrintArray(lmas, size);
    std::cout << "LSD ";
    CheckingSort(lmas, tmp, size);
    if (size < 15)
        std::cout << "Array after LSD sort with simle merge: ";
    PrintArray(mas, size);
    std::cout << "LSD with simple merge ";
    CheckingSort(mas, tmp, size);
    std::cout << "Execution time LSD sort: " << time_lsd << std::endl;
    std::cout << "Execution time parallel LSD sort: " << time_plsd << std::endl;
    std::cout << "Execution merge time: " << merge_time << std::endl;
    std::cout << "Execution time LSD sort with simple merge: " << ptime_lsd << std::endl;
    std::cout << "Boost is : " << time_lsd / ptime_lsd << std::endl;
    delete[] mas;
    delete[] tmp;
    delete[] lmas;
    return 0;
}
