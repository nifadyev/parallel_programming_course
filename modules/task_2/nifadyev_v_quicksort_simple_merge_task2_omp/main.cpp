// Copyright 2019 Nifadyev Vadim
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>

// void parallelQuickSort(int* a, int* sorted_a, int n);
void parallelQuickSort(int* array, const int SIZE, const int threads);
void quickSort(int* array, int first, int last);
void merge(int* arr1_begin, int* arr1_end, int* arr2_begin, int* arr2_end,
           int* dest);
bool isCorrectlySorted(int* customSortedArray, int* stdSortedArray,
                       const int size);
void printArray(int* array, const int size);

int main() {
    const int SIZE = 12;
    const int THREADS = omp_get_max_threads();
    int array[SIZE] = {0};
    int buf[SIZE] = {0};
    int linear_sorted_array[SIZE] = {0};
    double tot_time = 0;

    // srand((unsigned int)time(NULL));
    srand(static_cast<unsigned int>(time(NULL)));
    std::generate(array, array + SIZE, []() { return std::rand() % 100 + 1; });
    std::cout << "Initial array: ";
    printArray(array, SIZE);
    std::copy(array, array + SIZE, linear_sorted_array);
    tot_time = omp_get_wtime();

    parallelQuickSort(array, SIZE, THREADS);

    tot_time = omp_get_wtime() - tot_time;
    quickSort(linear_sorted_array, 0, SIZE - 1);

    // std::cout << "linear sorted array: ";
    // printArray(linear_sorted_array, SIZE);

    if (isCorrectlySorted(array, linear_sorted_array, SIZE)) {
        std::cout << "Sorted array:  ";
        printArray(array, SIZE);
    } else {
        std::cout << "Quick sort failed. Array is not sorted" << std::endl;
    }
    std::cout << tot_time << std::endl;

    return 0;
}

void quickSort(int* array, int first, int last) {
    int i = first, j = last;
    int temp, mid = array[(first + last) / 2];

    do {
        while (array[i] < mid) i++;
        while (array[j] > mid) j--;

        if (i <= j) {
            if (i < j) {
                temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
            i++;
            j--;
        }
    } while (i <= j);

    if (i < last) quickSort(array, i, last);
    if (first < j) quickSort(array, first, j);
}

void merge(int* arr1_begin, int* arr1_end, int* arr2_begin, int* arr2_end,
           int* dest) {
    while ((arr1_begin <= arr1_end) && (arr2_begin <= arr2_end)) {
        auto elem1 = *arr1_begin;
        auto elem2 = *arr2_begin;

        if (elem1 < elem2) {
            *dest = elem1;
            ++arr1_begin;
        } else {
            *dest = elem2;
            ++arr2_begin;
        }
        ++dest;
    }
    if (arr1_begin > arr1_end) {
        // copyArray(arr2_begin, arr2_end + 1, dest);
        std::copy(arr2_begin, arr2_end + 1, dest);
    // /*arr2_begin > arr2_end*/
    } else {
        // copyArray(arr1_begin, arr1_end + 1, dest);
        std::copy(arr1_begin, arr1_end + 1, dest);
    }
}

// void parallelQuickSort(int* a, int* sorted_a, int n) {
void parallelQuickSort(int* array, const int SIZE, const int threads) {
    // int threads_count = omp_get_max_threads();
    int threads_count = threads;
    int sorted_a[SIZE] = {0};
    // int threads_count = 2;
    omp_set_num_threads(threads_count);
    // TODO(nifadyev): not working for 1 thread
    // std::cout << threads_count;
    int chunk_size = SIZE / threads_count ? SIZE / threads_count : 1;

#pragma omp parallel for schedule(static) shared(array)
    for (int i = 0; i < threads_count; ++i) {
        auto first = i * chunk_size;
        auto last = first + chunk_size - 1;
        // std::cout << last << std::endl;
        // std::cout << omp_get_thread_num() << std::endl;
        quickSort(array, first, last);
    }

    auto ret = SIZE % threads_count;
    if (ret) {
        quickSort(array, threads_count * chunk_size, SIZE - 1);
    }

    // reduction
    int array_size = chunk_size;
    int count = array_size << 1;
    int j = 0;
    while (array_size < SIZE) {
        // #pragma omp parallel for shared(array, sorted_a)
        for (j = 0; j <= SIZE - count; j += count) {
            int* first_array_start = array + j;
            int* first_array_end = first_array_start + array_size - 1;

            int* second_array_start = array + j + array_size;
            int* secon_array_end = second_array_start + array_size - 1;

            merge(first_array_start, first_array_end, second_array_start,
                  secon_array_end, sorted_a + j);
        }

        // copyArray(sorted_a, sorted_a + j, array);
        std::copy(sorted_a, sorted_a + j, array);

        array_size = (array_size << 1) > SIZE ? SIZE : array_size << 1;
        count = array_size << 1;
    }
    // std::cout << "j " << j << std::endl;

    if (ret) {
        int array_size = SIZE - ret;
        int* first_array_start = array;
        int* first_array_end = first_array_start + array_size - 1;

        int* second_array_start = array + array_size;
        int* secon_array_end = second_array_start + ret - 1;

        merge(first_array_start, first_array_end, second_array_start,
              secon_array_end, sorted_a);
    }

    // copyArray(sorted_a, sorted_a + SIZE, array);
    std::copy(sorted_a, sorted_a + SIZE, array);
}

/* Compare array sorted by QuickSort function to
array sorted by qsort() from standard library.

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
