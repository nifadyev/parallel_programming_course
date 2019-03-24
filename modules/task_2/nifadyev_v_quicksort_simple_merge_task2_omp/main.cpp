// Copyright 2019 Nifadyev Vadim
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>

void parallelQuickSort(int* array, const int SIZE, const int threads);
void quickSort(int* array, int low, int high);
void merge(int* p_first_array_start, int* p_first_array_end,
           int* p_second_array_start, int* p_second_array_end,
           int* p_destination_array);
bool isCorrectlySorted(int* customSortedArray, int* stdSortedArray,
                       const int size);
void printArray(int* array, const int size);

int main() {
    const int SIZE = 1000001;
    const int THREADS = omp_get_max_threads();
    // const int THREADS = 2;
    int* array = new int[SIZE];
    int* linear_sorted_array = new int[SIZE];
    double start_time = 0.0, end_time = 0.0;
    double linear_time = 0.0, parallel_time = 0.0;

    // srand((unsigned int)time(NULL));
#pragma omp master
    {
        srand(static_cast<unsigned int>(time(NULL)));
        std::generate(array, array + SIZE,
                      []() { return std::rand() % 100 + 1; });
        if (SIZE <= 50) {
            std::cout << "Initial array:  ";
            printArray(array, SIZE);
        }
        std::copy(array, array + SIZE, linear_sorted_array);
    }

    start_time = omp_get_wtime();
    parallelQuickSort(array, SIZE, THREADS);
    parallel_time = omp_get_wtime() - start_time;

#pragma omp master
    {
        start_time = omp_get_wtime();
        quickSort(linear_sorted_array, 0, SIZE - 1);
        linear_time = omp_get_wtime() - start_time;

        if (isCorrectlySorted(array, linear_sorted_array, SIZE)) {
            if (SIZE <= 50) {
                std::cout << "Sorted array:  ";
                printArray(array, SIZE);
            }
        } else {
            std::cout << "Quick sort failed. Array is not sorted" << std::endl;
        }
        std::cout << "Linear time: " << linear_time << std::endl;
        std::cout << "Parallel time: " << parallel_time << std::endl;
    }

    delete[] array;
    delete[] linear_sorted_array;

    return 0;
}

/*  Quick sort implementation.
    This function takes middle element of array as pivot, places
    the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
    to left of pivot and all greater elements to right
    of pivot.

    array --> Array to be sorted,
    low  --> Starting index,
    high  --> Ending index. */
void quickSort(int* array, int low, int high) {
    int i = low, j = high;
    int pivot = array[(low + high) / 2];

    do {
        while (array[i] < pivot) {
            i++;
        }
        while (array[j] > pivot) {
            j--;
        }

        if (i <= j) {
            if (i < j) {
                std::swap(array[i], array[j]);
            }
            i++;
            j--;
        }
    } while (i <= j);

    if (i < high) {
        quickSort(array, i, high);
    }
    if (low < j) {
        quickSort(array, low, j);
    }
}

// TODO: add docstring
void merge(int* p_first_array_start, int* p_first_array_end,
           int* p_second_array_start, int* p_second_array_end,
           int* p_destination_array) {
    while ((p_first_array_start <= p_first_array_end) &&
           (p_second_array_start <= p_second_array_end)) {
        auto elem1 = *p_first_array_start;
        auto elem2 = *p_second_array_start;

        if (elem1 < elem2) {
            *p_destination_array = elem1;
            p_first_array_start++;
        } else {
            *p_destination_array = elem2;
            p_second_array_start++;
        }
        p_destination_array++;
    }

    if (p_first_array_start > p_first_array_end) {
        std::copy(p_second_array_start, p_second_array_end + 1,
                  p_destination_array);
    } else {
        std::copy(p_first_array_start, p_first_array_end + 1,
                  p_destination_array);
    }
}

// void parallelQuickSort(int* a, int* sorted_array, int n) {
void parallelQuickSort(int* array, const int SIZE, const int threads) {
    int* sorted_array = new int[SIZE];
    omp_set_num_threads(threads);
    // int chunk_size = SIZE / threads ? SIZE / threads : 1;
    // By default SIZE cannot be less than number of threads
    int chunk_size = SIZE / threads;

#pragma omp parallel for schedule(static) shared(array)
    for (int i = 0; i < threads; i++) {
        auto low = i * chunk_size;
        auto high = low + chunk_size - 1;

        quickSort(array, low, high);
    }

    auto ret = SIZE % threads;
    if (ret) {
        quickSort(array, threads * chunk_size, SIZE - 1);
    }

    // reduction
    int array_size = chunk_size;  // ? current_chunk_size
    // int count = array_size << 1;  // ? step
    int j = 0;
    int merges = 0;
    // while (array_size < SIZE) {
    // #pragma omp parallel for
    // TODO(nifadyev): write shift and sizess of subarrays into 2 arrays
    // outside cycle to decrease complexity
    for (int count = array_size << 1; array_size < SIZE;
         count = array_size << 1) {
        // #pragma omp parallel for shared(array, sorted_array)
        for (j = 0; j <= SIZE - count; j += count) {
            int* first_array_start = array + j;
            int* first_array_end = first_array_start + array_size - 1;

            int* second_array_start = array + j + array_size;
            int* second_array_end = second_array_start + array_size - 1;

            merge(first_array_start, first_array_end, second_array_start,
                  second_array_end, sorted_array + j);
            merges++;
        }

        std::copy(sorted_array, sorted_array + j, array);

        array_size = (array_size << 1) > SIZE ? SIZE : array_size << 1;
        // count = array_size << 1;
    }

    if (ret) {
        int array_size = SIZE - ret;
        int* first_array_start = array;
        int* first_array_end = first_array_start + array_size - 1;

        int* second_array_start = array + array_size;
        int* second_array_end = second_array_start + ret - 1;

        merge(first_array_start, first_array_end, second_array_start,
              second_array_end, sorted_array);
        merges++;
    }
    std::cout << "Merge counter: " << merges << std::endl;

    std::copy(sorted_array, sorted_array + SIZE, array);
    delete[] sorted_array;
}

/* Compare array sorted by parallelQuickSort() to
array sorted by linear quickSort().

    customSortedArray[] --> Array sorted by QuickSort,
    stdSortedArray[]  --> Array sorted by qsort,
    size  --> Size of both arrays. */
bool isCorrectlySorted(int* customSortedArray, int* stdSortedArray,
                       const int size) {
    for (int i = 0; i < size; i++) {
        if (customSortedArray[i] != stdSortedArray[i]) {
            return false;
        }
    }

    return true;
}

void printArray(int* array, const int size) {
    for (int i = 0; i < size; i++) {
        std::cout << array[i] << " ";
    }
    std::cout << std::endl;
}
