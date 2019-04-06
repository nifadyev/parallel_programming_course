// Copyright 2019 Novikova Alisa
#include <iostream>
#include <ctime>
#include <random>

int ShowMatrix(double* A, int N) {
        std::cout << std::endl;
        for (int i = 0; i < N*N; i += N) {
                for (int j = 0; j< N; j++)
                        std::cout << A[i + j] << " ";
                std::cout << std::endl;
        }
        std::cout << std::endl;
        return 1;
}
double* GetRandMatrix(int N) {
        double *result = new double[N*N];
        for (int i = 0; i< N*N; i++)
                result[i] = (std::rand() % 10000) / 1000.0f;
        return result;
}
void MultiplayMatrix(double* A, double* B, double* C, int blockSize, int n) {
        for (int i = 0; i < blockSize; ++i)
                for (int j = 0; j < blockSize; ++j)
                        for (int k = 0; k < blockSize; ++k) {
                                C[i * n + j] += A[i *n + k] * B[k *n + j];
                        }
}
void Cannon(double *A, double *B, double* C, int n, int q) {
        int blockSize = n / q;
        for (int i = 0; i < q; ++i) {
                for (int j = 0; j < q; ++j) {
                        for (int k = 0; k < q; ++k) {
                                MultiplayMatrix(&A[(i*n + (j + i + k) % q)*blockSize],
                                                &B[(((i + j + k) % q)*n + j)*blockSize],
                                                &C[(i*n + j)*blockSize], blockSize, n);
                        }
                }
        }
}
int CheckMatrixForEqual(double *A, double *C, int N) {
        for (int i = 0; i < N*N; i++)
                if (std::abs(A[i] - C[i]) > 0.000001)
                        return 0;
        return 1;
}
int main(int argc, char** argv) {
        srand((unsigned int)time(0));
        double *A, *B, *C1, *C2;
        int N = 1000, q = 2;
        if (argc == 3) {
                N = atoi(argv[1]);
                q = atoi(argv[2]);
        }
        A = GetRandMatrix(N);
        B = GetRandMatrix(N);
        // ShowMatrix(A, N);
        // ShowMatrix(B, N);
        C1 = new double[N*N];
        C2 = new double[N*N];
        for (int i = 0; i < N*N; i++) {
                C1[i] = 0.0;
                C2[i] = 0.0;
        }
        clock_t start = clock();
        Cannon(A, B, C1, N, q);
        MultiplayMatrix(A, B, C2, N, N);
        clock_t end = clock();
        // ShowMatrix(C1, N);
        // ShowMatrix(C2, N);
        if (CheckMatrixForEqual(C1, C2, N) == 1)
                std::cout << "Matrices are equal" << std::endl;
        else
                std::cout << "Matrices are not equal" << std::endl;
        std::cout << "Time = " << end - start << std::endl;
        return 0;
}
