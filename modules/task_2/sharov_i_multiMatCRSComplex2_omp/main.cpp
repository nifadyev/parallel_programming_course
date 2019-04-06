// Copyright 2019 Sharov Ivan
#include <omp.h>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <random>
#include <cstring>
#include <cmath>

const double zerocol = 0.000001;

int random(int lower_bound, int upper_bound) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(lower_bound, upper_bound);
    return uni(rng);
}

class TComplex {
    int Re;
    int Im;

 public:
    explicit TComplex(int a = 0, int b = 0);
    ~TComplex();
    void Print();
    void Input();

    const TComplex& operator= (TComplex);
    TComplex operator +(const TComplex &c);
    TComplex operator *(const TComplex &c);
    TComplex operator -(const TComplex &c);
    TComplex operator /(const TComplex &c);
    double modul() {
        double tmp;
        tmp = sqrt(static_cast<double>((Re*Re) + (Im*Im)));
        return static_cast<double>(tmp);
    }
};

TComplex::TComplex(int a, int b) {
    Re = a;
    Im = b;
}

TComplex::~TComplex(void) {}

const TComplex& TComplex::operator= (TComplex CopyBuffer) {
    Re = CopyBuffer.Re;
    Im = CopyBuffer.Im;
    return *this;
}


void TComplex::Print() {
    if ((Re == 0) && (Im == 0))
        std::cout << "0";
    if ((Re == 0) && (Im != 0))
        printf(" %di ", Im);
    if ((Re != 0) && (Im == 0))
        printf(" %d ", Re);
    if ((Re != 0) && (Im != 0) && (Im > 0))
        printf(" %d+%di ", Re, Im);
    if ((Re != 0) && (Im != 0) && (Im < 0))
        printf(" %d%di ", Re, Im);
}

void TComplex::Input() {
    Re = random(1, 5);
    Im = random(0, 5);
}

TComplex TComplex:: operator + (const TComplex &X) {
    TComplex res;
    res.Re = Re + X.Re;
    res.Im = Im + X.Im;
    return res;
}

TComplex TComplex::operator -(const TComplex &C) {
    TComplex res;
    res.Re = Re - C.Re;
    res.Im = Im - C.Im;
    return res;
}

TComplex TComplex::operator *(const TComplex &D) {
    TComplex res;
    res.Re = Re*D.Re - Im*D.Im;
    res.Im = Re*D.Im + D.Re*Im;
    return res;
}

TComplex TComplex::operator /(const TComplex &R) {
    TComplex res;
    if ((R.Re == 0) && (R.Im == 0)) {
      std::cout << "you can`t divide by 0";
    } else {
        res.Re = (Re*R.Re + Im*R.Im) / (R.Re*R.Re + R.Im*R.Im);
        res.Im = (Im*R.Re - R.Im*Re) / (R.Re*R.Re + R.Im*R.Im);
    }
    return res;
}

struct Matrix {
    int N;
    int NZ;
    TComplex *Value;
    int*Col;
    int*RowIndex;

    Matrix(int n, int nz) : N(n), NZ(nz) {
        Value = new TComplex[NZ];
        Col = new int[NZ];
        RowIndex = new int[N + 1];
    }

    Matrix(Matrix &mat) : N(mat.N), NZ(mat.NZ) {
        Value = new TComplex[mat.NZ];
        for (int i = 0; i < mat.NZ; i++)
            Value[i] = mat.Value[i];
        Col = new int[mat.NZ];
        for (int i = 0; i < mat.NZ; i++)
            Col[i] = mat.Col[i];
        RowIndex = new int[mat.N + 1];
        for (int i = 0; i < mat.N + 1; i++)
            RowIndex[i] = mat.RowIndex[i];
    }
};

void Print(int N, const Matrix *mat) {
    int k = 0;
    int col;

    for (int i = 0; i < N; i++) {
        col = 0;
        for (int j = 0; j < N; j++) {
            if (((mat  -> Col[k]) == j) && (col < (mat -> RowIndex[i + 1] - mat -> RowIndex[i]))) {
                (mat -> Value[k]).Print();
                k++;
                col++;
            } else {
                printf("0 ");
            }
            if (j == (N - 1)) printf("\n");
        }
    }
    printf("\n");
}

void GetMatrix(int N, int countinrow, const Matrix *mat) {
    int tmp;
    bool flag;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < countinrow; j++) {
            do {
                mat -> Col[i*countinrow + j] = random(0, 32767) % N;
                flag = false;
                for (int k = 0; k < j; k++)
                    if (mat -> Col[i*countinrow + j] == mat -> Col[i * countinrow + k]) flag = true;
            } while (flag == true);
        }

        for (int j = 0; j < (countinrow - 1); j++)
            for (int k = 0; k < (countinrow - 1); k++)
                if (mat -> Col[i*countinrow + k] > mat -> Col[i * countinrow + k + 1]) {
                    tmp = mat -> Col[i*countinrow + k];
                    mat -> Col[i*countinrow + k] = mat -> Col[i * countinrow + k + 1];
                    mat -> Col[i * countinrow + k + 1] = tmp;
                }
    }

    TComplex Z;
    for (int i = 0; i < (countinrow*N); i++) {
        Z.Input();
        mat -> Value[i] = Z;
    }

    int startIndex = 0;
    for (int i = 0; i < (N + 1); i++) {
        mat -> RowIndex[i] = startIndex;
        startIndex = startIndex + countinrow;
    }
}

void Transposing(Matrix *B, const Matrix *BT) {
    for (int i = 0; i < B -> NZ; i++) {
        BT -> Value[i] = TComplex(0, 0);
    }

    for (int i = 0; i < B -> NZ; i++) {
        BT -> Col[i] = 0;
    }

    for (int i = 0; i < (B -> N + 1); i++) {
        BT -> RowIndex[i] = 0;
    }

    for (int i = 0; i< B -> NZ; i++)
        BT -> RowIndex[B -> Col[i] + 1]++;

    int startIndex = 0;
    int tmp;
    for (int i = 1; i < ((B -> N) + 1); i++) {
        tmp = BT -> RowIndex[i];
        BT -> RowIndex[i] = startIndex;
        startIndex = startIndex + tmp;
    }

    TComplex z1;
    int j1, j2, colm, index1, index2;
    for (int i = 0; i < B -> N; i++) {
        j1 = B -> RowIndex[i];
        j2 = B -> RowIndex[i + 1];
        colm = i;
        for (int j = j1; j < j2; j++) {
            z1 = B -> Value[j];
            index1 = B -> Col[j];
            index2 = BT -> RowIndex[index1 + 1];
            BT -> Value[index2] = z1;
            BT -> Col[index2] = colm;
            (BT -> RowIndex[index1 + 1])++;
        }
    }

    if (BT -> N < 15) {
        printf("\n matrix Bt \n");
        Print(BT -> N, BT);
    }
}

void multiplication(const Matrix *A, const Matrix *BT) {
    double t1, t2;
    double modul;
    t1 = omp_get_wtime();

    std::vector<TComplex> val;
    std::vector<int> columns;
    std::vector<int> indexrow;

    int  i = 0, j = 0, index1 = 0, index2 = 0, index3 = 0, index4 = 0;

    int *tmp = new int[A->N];
    int NZ = 0;
    indexrow.push_back(0);
    for (i = 0; i < A->N; i++) {
        memset(tmp, -1, A->N * sizeof(int));
        index1 = A->RowIndex[i];
        index2 = A->RowIndex[i + 1];
        for (j = index1; j < index2; j++) {
            int col = A->Col[j];
            tmp[col] = j;
        }

        for (j = 0; j < A->N; j++) {
             TComplex summa = TComplex(0, 0);
             index3 = BT->RowIndex[j];
             index4 = BT->RowIndex[j + 1];
             for (int k = index3; k < index4; k++) {
                int bcol = BT->Col[k];
                int aind = tmp[bcol];
                if (aind != -1)
                    summa = summa + A->Value[aind] * BT->Value[k];
             }
             modul = (summa.modul());
             if (modul > zerocol) {
                columns.push_back(j);
                val.push_back(summa);
                NZ++;
             }
        }
        indexrow.push_back(NZ);
    }
    delete[] tmp;
    Matrix C(A->N, NZ);

    for (i = 0; i < NZ; i++) {
        C.Col[i] = columns[i];
        C.Value[i] = val[i];
    }
    for (i = 0; i < A->N + 1; i++)
        C.RowIndex[i] = indexrow[i];

    t2 = omp_get_wtime();


    if (C.N < 50) {
        printf("Result matrix C: \n");
        Print(C.N, &C);
        printf(" \n");
    }

    printf("Runtime:  %f\n", t2 - t1);
}


void multiplicationOMP(const Matrix *A, const Matrix *BT) {
    double t1, t2;
    double modul;

    t1 = omp_get_wtime();
    std::vector<TComplex>* val = new std::vector<TComplex>[A->N];
    std::vector<int>* columns = new std::vector<int>[A->N];
    int* indexrow = new int[A->N + 1];
    memset(indexrow, 0, sizeof(int) * A->N);
    int  i = 0, j = 0, index1 = 0, index2 = 0, index3 = 0, index4 = 0;
    int chunk = 5;

#pragma omp parallel
{
    int *tmp = new int[A->N];
#pragma omp for private(j, index1, index2, index3, index4)schedule(static, chunk)
    for (i = 0; i < A->N; i++) {
        memset(tmp, -1, A->N * sizeof(int));
        index1 = A->RowIndex[i];
        index2 = A->RowIndex[i + 1];
        for (j = index1; j < index2; j++) {
             int col = A->Col[j];
             tmp[col] = j;
        }

        for (j = 0; j < A->N; j++) {
            TComplex summa = TComplex(0, 0);
            index3 = BT->RowIndex[j];
            index4 = BT->RowIndex[j + 1];
            for (int k = index3; k < index4; k++) {
                int bcol = BT->Col[k];
                int aind = tmp[bcol];
                if (aind != -1)
                    summa = summa + A->Value[aind] * BT->Value[k];
            }
            modul = (summa.modul());
            if (modul > zerocol) {
                columns[i].push_back(j);
                val[i].push_back(summa);
                indexrow[i]++;
            }
        }
    }
    delete[] tmp;
}
    int NZ = 0;
    for (i = 0; i < A->N; i++) {
        int temp = indexrow[i];
        indexrow[i] = NZ;
        NZ += temp;
    }
    indexrow[A->N] = NZ;

    Matrix C(A->N, NZ);

    int count = 0;
    for (i = 0; i < A->N; i++) {
        int size = static_cast<int>(columns[i].size());
        memcpy(&C.Col[count], &columns[i][0], size * sizeof(int));
        memcpy(&C.Value[count], &val[i][0], size * sizeof(TComplex));
        count += size;
    }
    memcpy(C.RowIndex, &indexrow[0], (A->N + 1) * sizeof(int));

    delete[] indexrow;
    delete[] columns;
    delete[] val;

    t2 = omp_get_wtime();

    if (C.N < 50) {
    printf("Result matrix C: \n");
    Print(C.N, &C);
    printf(" \n");
    }

    printf("Runtime:  %f\n", t2 - t1);
}

int main(int argc, char* argv[]) {
    int n, nzInRow;

    if (argc == 3) {
        n = atoi(argv[1]);
        nzInRow = atoi(argv[2]);
    } else {
        n = 5;
        nzInRow = 2;
    }

    if (nzInRow > n) {
        printf("Invalid input count of notnull in string.\n");
        exit(1);
    }

    printf("Zadacha 1. The multiplication of sparse matrices. CRS. Complex type.\n");
    Matrix A(n, nzInRow*n);
    Matrix B(n, nzInRow*n);
    Matrix BT(n, nzInRow*n);
    GetMatrix(n, nzInRow, &A);

    if (A.N < 15) {
        printf("matrix A \n");
        Print(n, &A);
    }

    GetMatrix(n, nzInRow, &B);

    if (B.N < 15) {
        printf("matrix B \n");
        Print(n, &B);
    }

    Transposing(&B, &BT);
    printf("Posledovatelnaya\n");
    multiplication(&A, &BT);
    printf("\n openMP\n");
    multiplicationOMP(&A, &BT);
    return 0;
}
