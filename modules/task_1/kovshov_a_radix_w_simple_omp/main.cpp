// Copyright 2019 Kovshov Andrey
#include <iostream>
#include <ctime>
#include <algorithm>

int max(int arr[], int n) {
    int m = arr[0];
    for (int i = 1; i < n; i++) {
        if (arr[i] > m) {
            m = arr[i];
        }
    }
    return m;
}

void radixSort(int arr[], int n) {
    int exp = 1;
    int m = max(arr, n);

    while (m / exp > 0) {
        int *temp = new int[n];
        int count[10] = { 0 };

        for (int i = 0; i < n; i++) {
            count[(arr[i] / exp) % 10]++;
        }

        for (int i = 1; i < 10; i++) {
            count[i] += count[i - 1];
        }

        for (int i = n - 1; i >= 0; i--) {
            temp[--count[(arr[i] / exp) % 10]] = arr[i];
        }

        for (int i = 0; i < n; i++) {
            arr[i] = temp[i];
        }

        exp *= 10;
    }
}

int main(int argc, char* argv[]) {
    const int arrLength = 20;
    int *arr = new int[arrLength];

    std::cout << "Source array: ";
    srand((unsigned int)time(NULL));
    for (int i = 0; i < arrLength; i++) {
        arr[i] = std::rand() % 100;
        std::cout << arr[i] << " ";
    }

    radixSort(arr, arrLength);

    std::cout << "\n\nSorted array: ";
    for (int i = 0; i < arrLength; i++) {
        std::cout << arr[i] << " ";
    }
    return 0;
}
