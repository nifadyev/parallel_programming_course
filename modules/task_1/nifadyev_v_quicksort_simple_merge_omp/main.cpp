// Copyright Nifadyev Vadim 2019
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

// QuickSort - linear algorithm

void Swap(int *a, int *b) {
    int temp = *a;

    *a = *b;
    *b = temp;
}

/* This function takes array[high] as pivot, places
   the pivot element at its correct position in sorted
   array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot.

   array[] --> Array to be sorted,
   low  --> Starting index,
   high  --> Ending index. */
int Partition(int* array, int low, int high) {
    int pivot = array[high];
    int i = (low - 1);  // Index of smaller element

    for (int j = low; j <= high - 1; j++) {
        if (array[j] <= pivot) {
            i++;
            Swap(&array[i], &array[j]);
        }
    }

    Swap(&array[i + 1], &array[high]);

    return (i + 1);
}

/* QuickSort implementation.

   array[] --> Array to be sorted,
   low  --> Starting index,
   high  --> Ending index. */
void QuickSort(int* array, int low, int high) {
    if (low < high) {
        // array[partitioningIndex] is now at right place
        int partitioningIndex = Partition(array, low, high);

        // Separately sort elements before
        // Partition and after Partition
        QuickSort(array, low, partitioningIndex - 1);
        QuickSort(array, partitioningIndex + 1, high);
    }
}

void PrintArray(int* array, const int size) {
    int i;


    for (i = 0; i < size; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");
}

void GenerateArray(int *array, const int size) {
    int i;
    unsigned int seed = time(NULL);

    for (i = 0; i < size; i++) {
        array[i] = rand_r(&seed) % 100;
    }
}

int main() {
    const int size = 50;
    int array[size] = {0};
    int n = sizeof(array) / sizeof(array[0]);

    GenerateArray(array, size);
    printf("Unsorted array: ");
    PrintArray(array, size);

    QuickSort(array, 0, n - 1);

    printf("Sorted array: ");
    PrintArray(array, size);

    return 0;
}