//  Copyright 2019 Chvanov Leonid

#include <ctime>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <functional>

void shellSort(int* arr, int n) {
    int step;
    for (step = n / 2; step > 0; step /= 2)
        for (int i = step; i < n; ++i)
            for (int j = i - step; (j >= 0) &&
                (arr[j] > arr[j + step]); j -= step)
                std::swap(arr[j], arr[j + step]);
}

void printArray(int* arr, int size) {
    for (int i = 0; i < size; i++)
        std::cout << arr[i] << " ";
    std::cout << std::endl;
}

bool compareArrays(int* arr0, int* arr1, int size) {
    for (int i = 0; i < size; i++)
        if (arr0[i] != arr1[i])
            return false;
    return true;
}

void copyArray(int* src, int* dest, int size) {
    for (int i = 0; i < size; i++)
        dest[i] = src[i];
}

int* generateArray(int size) {
    int* arr = new int[size];
    for (int i = 0; i < size; i++)
        arr[i] = std::rand() % 100;
    return arr;
}

int main(int argc, char* argv[]) {
    std::srand((unsigned int)std::time(nullptr));
    int* arr1 = nullptr;
    int* arr2 = nullptr;
    int size = 20;

    if (argc > 1)   size = atoi(argv[1]);
    arr1 = generateArray(size);
    arr2 = new int[size];
    copyArray(arr1, arr2, size);

    std::cout << "Original array: ";
    printArray(arr1, size);

    shellSort(arr1, size);
    std::sort(arr2, arr2 + size, std::less<int>());
    std::cout << "Sorted array:  ";
    printArray(arr1, size);
    printArray(arr2, size);

    if (compareArrays(arr1, arr2, size))
        std::cout << "Shell sort is working\n";
    else
        std::cout << "Shell sort isn't working\n";

    return 0;
}
