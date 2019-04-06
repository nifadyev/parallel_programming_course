// Copyright 2019 Makarin Kirill
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

void print_arr(int* arr, int size) {
    for (int i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

void shellSort(int* A, int N) {
    for (int k = N / 2; k > 0; k /= 2) {
        for (int i = k; i < N; i++) {
            int j;
            int t = A[i];
            for (j = i; j >= k; j -= k) {
                if (t < A[j - k]) {
                    A[j] = A[j - k];
                } else {
                    break;
                }
            }
            A[j] = t;
        }
    }
}

int main() {
    const int SIZE = 50;
    int* a = new int[SIZE];

    std::generate(a, a + SIZE, []() { return std::rand() % 100 + 1; });
    print_arr(a, SIZE);
    shellSort(a, SIZE);
    print_arr(a, SIZE);

    return 0;
}
