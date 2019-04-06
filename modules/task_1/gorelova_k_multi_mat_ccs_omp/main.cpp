//  Copyright 2019 Gorelova Ksenia

#include <iostream>
#include <vector>
#include <random>
#include <cstring>
#include <ctime>

using std::vector;
struct ccsMatrix {
    int N;
    int NZ;
    double* val;
    int* row;
    int* ColIndex;
};

void initMatrix(ccsMatrix * mtx, int nz, int n) {
    mtx->NZ = nz;
    mtx->N = n;
    mtx->val = new double[nz];
    mtx->row = new int[nz];
    mtx->ColIndex = new int[n + 1];
}

void freeMatrix(ccsMatrix * mtx) {
    delete[] mtx->val;
    delete[] mtx->row;
    delete[] mtx->ColIndex;
}

void transMatrix(ccsMatrix * A, ccsMatrix * AT) {
    int sum = 0;
    int index, colIn;
    double v;
    memset(AT->ColIndex, 0, (AT->N + 1) * sizeof(int));
    for (int i = 0; i < A->NZ; i++)
        AT->ColIndex[A->row[i] + 1]++;
    for (int i = 1; i <= A->N; i++) {
        int tmp = AT->ColIndex[i];
        AT->ColIndex[i] = sum;
        sum += tmp;
    }
    for (int i = 0; i < A->N; i++) {
        int row = i;
        for (int j = A->ColIndex[i]; j < A->ColIndex[i + 1]; j++) {
            v = A->val[j];
            colIn = A->row[j];
            index = AT->ColIndex[colIn + 1];
            AT->val[index] = v;
            AT->row[index] = row;
            AT->ColIndex[colIn + 1]++;
        }
    }
}

void generateMatrix(ccsMatrix * mtx, int n, int cntInCol) {
    int rows;
    bool exist;
    int nz = cntInCol * n;
    initMatrix(mtx, nz, n);
    for (int col = 0; col < n; col++) {
        for (int row = 0; row < cntInCol; row++) {
            do {
                rows = std::rand() % n;
                exist = false;
                for (int i = 0; i < row; i++)
                    if (rows == mtx->row[col*cntInCol + i])
                        exist = true;
                break;
            } while (exist == true);
            mtx->row[col*cntInCol + row] = rows;
        }
        for (int row = 0; row < cntInCol - 1; row++)
            for (int val = 0; val < cntInCol - 1; val++)
                if (mtx->row[col*cntInCol + val] >
                    mtx->row[col*cntInCol + val + 1]) {
                    int tmp = mtx->row[col*cntInCol + val];
                    mtx->row[col*cntInCol + val] =
                        mtx->row[col*cntInCol + val + 1];
                    mtx->row[col*cntInCol + val + 1] = tmp;
                }
    }
    for (int val = 0; val < nz; val++)
        mtx->val[val] = (std::rand() % 10000) / 1000.0f;
    int c = 0;
    for (int col = 0; col <= n; col++) {
        mtx->ColIndex[col] = c;
        c += cntInCol;
    }
}

void printMatrix(int n, ccsMatrix *mtx) {
    int i;
    int k = mtx[0].NZ;
    std::cout << "\n Value: ";
    for (i = 0; i < k; i++)
        std::cout << " "<< mtx[0].val[i];
    std::cout << "\n Row: ";
    for (i = 0; i < k; i++)
        std::cout << " "<< mtx[0].row[i];
    std::cout << "\n ColIndex: ";
    for (i = 0; i < n + 1; i++)
        std::cout << " " << mtx[0].ColIndex[i];
    std::cout << std::endl;
}

int multMatrix(ccsMatrix * A, ccsMatrix * B, ccsMatrix * C, int size) {
    int n = A->N;
    int NZ = 0, indexA, k, end;
    vector<int> colIndex, row;
    colIndex.push_back(0);
    int *tmp = new int[size];
    for (int i = 0; i < n; i++) {
        memset(tmp, -1, size * sizeof(int));
        for (int j = A->ColIndex[i]; j < A->ColIndex[i + 1]; j++)
            tmp[A->row[j]] = j;
        for (int j = 0; j < size; j++) {
            indexA = -1;
            k = B->ColIndex[j];
            end = B->ColIndex[j + 1];
            while (indexA == -1 && k < end) {
                indexA = tmp[B->row[k]];
                k++;
            }
            if (indexA != -1) {
                row.push_back(j);
                NZ++;
            }
        }
        colIndex.push_back(NZ);
    }
    initMatrix(C, NZ, n);
    for (int j = 0; j < NZ; ++j)
        C->row[j] = row[j];
    for (int i = 0; i <= n; ++i)
        C->ColIndex[i] = colIndex[i];
    delete[] tmp;
    return 0;
}

int numericMult(int size, const ccsMatrix * A, const ccsMatrix * B, ccsMatrix * C) {
    int n = A->N;
    int index = 0, indexA, startC, finishC, colC;
    int *tmp = new int[size];
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        startC = C->ColIndex[i];
        finishC = C->ColIndex[i + 1];
        if (finishC > startC) {
            memset(tmp, -1, sizeof(int)* size);
            for (int j = A->ColIndex[i]; j < A->ColIndex[i + 1]; ++j)
                tmp[A->row[j]] = j;
            for (int j = startC; j < finishC; ++j, sum = 0) {
                colC = C->row[j];
                for (int k = B->ColIndex[colC]; k < B->ColIndex[colC + 1]; ++k) {
                    indexA = tmp[B->row[k]];
                    if (indexA != -1)
                        sum += A->val[indexA] * B->val[k];
                }
                C->val[index] = sum;
                index++;
            }
        }
    }
    delete[] tmp;
    return 0;
}

int main(int argc, char** argv) {
    int N = (argc != 1) ? atoi(argv[1]) : 5;
    int cntInCol = (argc != 1) ? atoi(argv[2]) : 2;
    ccsMatrix A, B, AT, C;
    generateMatrix(&A, N, cntInCol);
    std::cout << "Matrix A: ";
    printMatrix(N, &A);
    initMatrix(&AT, A.NZ, A.N);
    transMatrix(&A, &AT);
    freeMatrix(&A);
    generateMatrix(&B, N, cntInCol);
    std::cout << "Matrix B: ";
    printMatrix(N, &B);
    initMatrix(&C, cntInCol, N);
    multMatrix(&AT, &B, &C, N);
    numericMult(N, &AT, &B, &C);
    std::cout << "Matrix C: ";
    printMatrix(N, &C);
    std::cout << A.N;
    freeMatrix(&AT);
    freeMatrix(&B);
    freeMatrix(&C);
    return 0;
}
