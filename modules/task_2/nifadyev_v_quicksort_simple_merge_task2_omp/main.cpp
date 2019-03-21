// Copyright Nifadyev Vadim 2019
#include <omp.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <utility>

// QuickSort - OpenMP algorithm

int partition(int* array, int low, int high);
void quick_sort(int* array, int low, int high);
void print_array(int* array, const int size);
bool is_correctly_sorted(int* customSortedArray, int* stdSortedArray,
                         const int size);
int* merge(int* firstArray, const int firstArrayLength, int* secondArray,
           const int secondArrayLength);

int main() {
    const int SIZE = 50;
    int array[SIZE] = {0};
    int copied_array[SIZE] = {0};

    std::generate(array, array + SIZE, []() { return std::rand() % 100; });
    std::copy(array, array + SIZE, copied_array);

    std::cout << "Initial array: ";
    print_array(array, SIZE);

    quick_sort(array, 0, SIZE - 1);
    // Sort copy of initial array by default qsort to compare sorted arrays
    qsort(copied_array, SIZE, sizeof(int), [](const void* a, const void* b) {
        return (*reinterpret_cast<const int*>(a) -
                *reinterpret_cast<const int*>(b));
    });

    if (is_correctly_sorted(array, copied_array, SIZE)) {
        std::cout << "Sorted array:  ";
        print_array(array, SIZE);
    } else {
        std::cout << "Quick sort failed. Array is not sorted" << std::endl;
    }

    return 0;
}

/* Takes array[high] as pivot, places
   the pivot element at its correct position in sorted
   array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot.

   array[] --> Array to be sorted,
   low  --> Starting index,
   high  --> Ending index. */
int partition(int* array, int low, int high) {
    int pivot = array[high];
    int i = (low - 1);  // Index of smaller element

    for (int j = low; j <= high - 1; j++) {
        if (array[j] <= pivot) {
            i++;
            std::swap(array[i], array[j]);
        }
    }

    std::swap(array[i + 1], array[high]);

    return (i + 1);
}

/* quick_sort implementation.

   array[] --> Array to be sorted,
   low  --> Starting index,
   high  --> Ending index. */
void quick_sort(int* array, int low, int high) {
    if (low < high) {
        // array[partitioningIndex] is now at right place
        int partitioningIndex = partition(array, low, high);

        /* Separately sort elements before
        partition and after partition */
        quick_sort(array, low, partitioningIndex - 1);
        quick_sort(array, partitioningIndex + 1, high);
    }
}

void print_array(int* array, const int size) {
    for (int i = 0; i < size; i++) {
        std::cout << array[i] << " ";
    }
    std::cout << std::endl;
}

/* Compare array sorted by QuickSort function to
array sorted by qsort() from standard library.

   customSortedArray[] --> Array sorted by ,
   stdSortedArray[]  --> Array sorted by qsort,
   size  --> Size of both arrays. */
bool is_correctly_sorted(int* customSortedArray, int* stdSortedArray,
                         const int size) {
    for (int i = 0; i < size; i++) {
        if (customSortedArray[i] != stdSortedArray[i]) {
            return false;
        }
    }

    return true;
}

/* Merge two arrays and return merged array.*/
int* merge(int * firstArray, const int firstArrayLength, int* secondArray,
           const int secondArrayLength) {
    int i = 0, k = 0, j = 0;
    int* resultingArray = new int(firstArrayLength + secondArrayLength);

    for (i; i < firstArrayLength + secondArrayLength; i++) {
        if (j > firstArrayLength - 1) {
            resultingArray[i] = secondArray[k];
            k++;
        } else if (k > secondArrayLength - 1) {
            resultingArray[i] = firstArray[j];
            j++;
        } else if (firstArray[j] <= secondArray[k]) {
            resultingArray[i] = firstArray[j];
            j++;
        } else {
            resultingArray[i] = secondArray[k];
            k++;
        }
    }

    return resultingArray;
}
