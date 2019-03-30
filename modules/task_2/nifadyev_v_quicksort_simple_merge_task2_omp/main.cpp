#define _SCL_SECURE_NO_WARNINGS
// Copyright 2019 Nifadyev Vadim
#include <omp.h>
#include <algorithm>
#include <ctime>
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

int main() {
    double start_time = 0.0, finish_time = 0.0;
    const int SIZE = 1234567;
    const int THREADS = omp_get_max_threads();
    int* array = new int[SIZE];
    int* linear_sorted_array = new int[SIZE];

    srand(static_cast<unsigned int>(time(NULL)));
    std::generate(array, array + SIZE, []() { return std::rand() % 100 + 1; });
    std::copy(array, array + SIZE, linear_sorted_array);

    if (SIZE < 50) {
        std::cout << "Initial array: ";
        printArray(array, SIZE);
    }


    start_time = omp_get_wtime();
    parallelQuickSort(array, SIZE, THREADS);
    finish_time = omp_get_wtime();
    double parallel_time = finish_time - start_time;

    start_time = omp_get_wtime();
    quickSort(linear_sorted_array, 0, SIZE - 1);
    finish_time = omp_get_wtime();
    double linear_time = finish_time - start_time;

    if (isCorrectlySorted(linear_sorted_array, array, SIZE)) {
        if (SIZE < 50) {
            std::cout << "Sorted array:  ";
            printArray(array, SIZE);
        }
    } else {
        std::cout << "Quick sort failed. Array is not sorted" << std::endl;
    }

    std::cout << "Parallel time: " << parallel_time << std::endl;
    std::cout << "Linear time:   " << linear_time << std::endl;
    std::cout << "Boost: " << linear_time / parallel_time << std::endl;

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

    // Each element of second subarray is smaller than any element from first
    // subarray
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

    // Copy to initial array merged array
    for (i = first_subarray_start_index;
         i < second_subarray_start_index + second_subarray_size; i++) {
        array[i] = merged_subarray[i - first_subarray_start_index];
    }
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
    for (int i = 0; i < threads - 1; i++) {
        int low = i * subarray_size;
        int high = low + subarray_size - 1;

        quickSort(array, low, high);
    }

    // Last part of array may contain remainder
    // Thats why it should be handled outside of the cycle
    quickSort(array, (threads - 1) * subarray_size, size - 1);

    int step = 1;  // Distance between merging threads
    for (int i = threads / 2; i > 0; i /= 2) {
        int first_subarray_size = subarray_size * static_cast<int>(step);
        int second_subarray_size = first_subarray_size;

        // Handle remainded array elements during last iteration
        if (i / 2 <= 0) {
            second_subarray_size += size % threads;
        }
        omp_set_num_threads(i);

#pragma omp parallel for schedule(static) \
    shared(array, first_subarray_size, second_subarray_size, step)
        for (int f = 0; f < i; f++) {
            int thread_id = omp_get_thread_num();
            int first_subarray_start_index =
                thread_id * subarray_size * static_cast<int>(step * 2);
            int second_subarray_start_index =
                first_subarray_start_index +
                subarray_size * static_cast<int>(step);
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
