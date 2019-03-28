// Copyright 2019 Nifadyev Vadim
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <utility>

void parallelQuickSort(int* array, const int SIZE, const int threads);
void quickSort(int* array, int low, int high);
void merge(int* p_first_array_start, int* p_first_array_end,
           int* p_second_array_start, int* p_second_array_end,
           int* p_destination_array);
bool isCorrectlySorted(int* customSortedArray, int* stdSortedArray,
                       const int size);
void printArray(int* array, const int size);

int main() {
    // TODO(nifadyev): sort fails for little SIZE, ex 12
    const int SIZE = 12001;
    const int THREADS = omp_get_max_threads();
    int* array = new int[SIZE];
    int* linear_sorted_array = new int[SIZE];
    double start_time = 0.0, end_time = 0.0;
    double linear_time = 0.0, parallel_time = 0.0;

// Store initial array on root thread
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

// Perform linear quick sort and data output on root thread
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
            // if (i < j) {
            std::swap(array[i], array[j]);

            // }
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

/*  Merge two arrays into one using pointers to starting and ending elements of
   arrays. Compare elements of arrays and write element with lower value to
   destination array.

    p_first_array_start --> pointer to beginning of first array,
    p_first_array_end --> pointer to ending of first array,
    p_second_array_start --> pointer to beginning of second array,
    p_second_array_end --> pointer to beginning of second array,
    p_destination_array --> pointer to beginning of destination array. */
void merge(int* p_first_array_start, int* p_first_array_end,
           int* p_second_array_start, int* p_second_array_end,
           int* p_destination_array) {
    // ! WTF is going on in this function
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

    std::copy(p_first_array_start, p_first_array_end + 1, p_destination_array);
}

/*  Sort array using parallel quick sort algorithm and openMP.
    Separetly sort parts of initial array (subarrays) on threads, then
    merge sorted subarrays into initial array.

    array --> array to be sorted,
    size --> array size,
    threads --> total number of threads involved in parallel algorithm. */
void parallelQuickSort(int* array, const int SIZE, const int threads) {
    int* sorted_array = new int[SIZE];
    omp_set_num_threads(threads);
    // By default SIZE cannot be less than number of threads
    int chunk_size = SIZE / threads;

#pragma omp parallel for schedule(static) shared(array)
    for (int i = 0; i < threads - 1; i++) {
        auto low = i * chunk_size;
        auto high = low + chunk_size - 1;

        quickSort(array, low, high);
    }

    // Last part of array may contain remainder
    // Thats why it should be handled outside of the cycle
    quickSort(array, (threads - 1) * chunk_size, SIZE - 1);

    auto remainder_size = SIZE % threads;

    // reduction
    int current_chunk_size = chunk_size;  // ? current_chunk_size
    // TODO(nifadyev): write shift and sizes of subarrays into 2 arrays
    // outside cycle to decrease complexity
    for (int count = current_chunk_size << 1; current_chunk_size < SIZE;
         count = current_chunk_size << 1) {
        int j = 0;
        for (j = 0; j <= SIZE - count; j += count) {
            int* first_array_start = array + j;
            int* first_array_end = first_array_start + current_chunk_size - 1;

            int* second_array_start = array + j + current_chunk_size;
            int* second_array_end = second_array_start + current_chunk_size - 1;

            merge(first_array_start, first_array_end, second_array_start,
                  second_array_end, sorted_array + j);
        }

        std::copy(sorted_array, sorted_array + j, array);

        current_chunk_size =
            (current_chunk_size << 1) > SIZE ? SIZE : current_chunk_size << 1;
    }

    // FIXME: extra merge execution, should remove it
    if (remainder_size) {
        int current_chunk_size = SIZE - remainder_size;
        int* first_array_start = array;
        int* first_array_end = first_array_start + current_chunk_size - 1;

        int* second_array_start = array + current_chunk_size;
        int* second_array_end = second_array_start + remainder_size - 1;

        merge(first_array_start, first_array_end, second_array_start,
              second_array_end, sorted_array);
    }

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
