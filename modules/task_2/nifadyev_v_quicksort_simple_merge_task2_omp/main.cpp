#define _SCL_SECURE_NO_WARNINGS
// Copyright 2019 Nifadyev Vadim
#include <omp.h>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <utility>

void quickSort(int* array, int low, int high);
void printArray(int* array, const int size);
void merge(int* array, const int first_subarray_size,
           const int second_subarray_size, const int first_subarray_start_index,
           const int second_subarray_start_index);
void parallelQuickSort(int* array, const int size, const int threads);
bool isCorrectlySorted(int* linear_sorted_array, int* parallel_sorted_array,
                       const int size);

int main(int argc, char** argv) {
    int SIZE = 10000;
    int THREADS = omp_get_max_threads();

    // Read `SIZE` and `THREADS` from console if they are presented
    if (argc == 2) {
        SIZE = atoi(argv[1]);
    } else if (argc == 3) {
        SIZE = atoi(argv[1]);
        THREADS = atoi(argv[2]);
    }

    int* array = new int[SIZE];
    int* linear_sorted_array = new int[SIZE];
    int* stl_sorted_array = new int[SIZE];

    srand(static_cast<unsigned int>(time(NULL)));
    std::generate(array, array + SIZE, []() { return std::rand() % 100 + 1; });
    std::copy(array, array + SIZE, linear_sorted_array);
    std::copy(array, array + SIZE, stl_sorted_array);

    double start_time = omp_get_wtime();
    parallelQuickSort(array, SIZE, THREADS);
    double finish_time = omp_get_wtime();
    double parallel_time = finish_time - start_time;

    start_time = omp_get_wtime();

    quickSort(linear_sorted_array, 0, SIZE - 1);
    finish_time = omp_get_wtime();
    double linear_time = finish_time - start_time;

    // Sort array using stl qsort to check custom quickSort algorithms
    qsort(stl_sorted_array, SIZE, sizeof(int),
          [](const void* a, const void* b) {
              return (*reinterpret_cast<const int*>(a) -
                      *reinterpret_cast<const int*>(b));
          });

    if (isCorrectlySorted(stl_sorted_array, array, SIZE) &&
        isCorrectlySorted(stl_sorted_array, linear_sorted_array, SIZE)) {
        if (SIZE < 50) {
            std::cout << "Initial array: ";
            printArray(array, SIZE);
            std::cout << "Sorted array:  ";
            printArray(array, SIZE);
        }
    } else {
        std::cout << "Quick sort failed. Array is not sorted" << std::endl;
    }
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Parallel time: " << parallel_time << std::endl;
    std::cout << "Linear time:   " << linear_time << std::endl;
    std::cout << "Boost: " << linear_time / parallel_time << std::endl;

    delete[] array;
    delete[] linear_sorted_array;
    delete[] stl_sorted_array;

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
            std::swap(array[i], array[j]);
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

void printArray(int* array, const int size) {
    for (int i = 0; i < size; i++) {
        std::cout << array[i] << " ";
    }
    std::cout << std::endl;
}

/*  Merge two subarrays into one. Compare elements of arrays
    and write element with lower value to merged subarray.
    Then copy merged subarrays into initial array.

    array --> Array which contains 2 subarrays to merge,
    first_subarray_size --> Length of first subarray,
    second_subarray_size --> Length of second subarray,
    first_subarray_start_index --> Starting index of first array,
    second_subarray_start_index --> Starting index of second array. */
void merge(int* array, const int first_subarray_size,
           const int second_subarray_size, const int first_subarray_start_index,
           const int second_subarray_start_index) {
    int i = 0, j = 0, k = 0;
    int* merged_subarray = new int[first_subarray_size + second_subarray_size];

    while (i < first_subarray_size && j < second_subarray_size) {
        if (array[first_subarray_start_index + i] <=
            array[second_subarray_start_index + j]) {
            merged_subarray[k] = array[first_subarray_start_index + i];
            i++;
        } else {
            merged_subarray[k] = array[second_subarray_start_index + j];
            j++;
        }
        k++;
    }

    // Remaining elements of second subarray is bigger than each element
    // of first subarray. They are written to the end of merged_subarray
    if (i < first_subarray_size) {
        for (int p = i; p < first_subarray_size; p++) {
            merged_subarray[k] = array[first_subarray_start_index + p];
            k++;
        }
    } else {
        // On the contrary
        for (int p = j; p < second_subarray_size; p++) {
            merged_subarray[k] = array[second_subarray_start_index + p];
            k++;
        }
    }

    // Copy merged array to initial array
    for (i = first_subarray_start_index;
         i < second_subarray_start_index + second_subarray_size; i++) {
        array[i] = merged_subarray[i - first_subarray_start_index];
    }

    delete[] merged_subarray;
}

/*  Sort array using parallel quick sort algorithm and openMP.
    Separately sort parts of initial array (subarrays) on threads, then
    merge sorted subarrays into initial array.

    array --> Array to be sorted,
    size --> Array size,
    threads --> Total number of threads involved in parallel algorithm. */
void parallelQuickSort(int* array, const int size, const int threads) {
    omp_set_num_threads(threads);
    int subarray_size = size / threads;

#pragma omp parallel for schedule(static) shared(array)
    for (int i = 0; i < threads; i++) {
        int low = i * subarray_size;
        int high = 0;
        // Array may contain remainder. Handle it on last thread
        if (i == threads - 1 && size % threads) {
            high = size - 1;
        } else {
            high = low + subarray_size - 1;
        }

        quickSort(array, low, high);
    }

    int step = 1;  // Distance between merging threads
    for (int i = threads / 2; i > 0; i /= 2) {
        int first_subarray_size = subarray_size * step;
        int second_subarray_size = first_subarray_size;

        // Handle remainded array elements during last iteration
        if (i / 2 <= 0) {
            second_subarray_size += size % threads;
        }

#pragma omp parallel for schedule(static) \
    shared(array, first_subarray_size, second_subarray_size, step)
        for (int j = 0; j < i; j++) {
            int thread_id = omp_get_thread_num();
            int first_subarray_start_index =
                thread_id * subarray_size * step * 2;
            int second_subarray_start_index =
                first_subarray_start_index + subarray_size * step;
            merge(array, first_subarray_size, second_subarray_size,
                  first_subarray_start_index, second_subarray_start_index);
        }
        step *= 2;
    }
}

/* Compare array sorted by parallelQuickSort() to
array sorted by linear quickSort().

    customSortedArray[] --> Array sorted by QuickSort,
    stdSortedArray[]  --> Array sorted by qsort,
    size  --> Size of both arrays. */
bool isCorrectlySorted(int* linear_sorted_array, int* parallel_sorted_array,
                       const int size) {
    for (int i = 0; i < size; i++) {
        if (linear_sorted_array[i] != parallel_sorted_array[i]) {
            return false;
        }
    }

    return true;
}
