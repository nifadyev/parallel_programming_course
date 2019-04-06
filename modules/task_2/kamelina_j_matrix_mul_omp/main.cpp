//  Copyright 2019 Kamelina Julia

#include <omp.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <complex>
#include <random>
#include <cstring>
#include <algorithm>

using std::vector;
using std::complex;

struct Matrix {
  complex<double>* val = 0;
  int* rows = 0;
  int* pointer = 0;
  int nz_size;

  explicit Matrix(int N = 0, int nz = 0): nz_size(nz) {
    if (N != 0 && nz_size != 0) {
      val = new complex<double>[nz];
      pointer = new int[N + 1];
      rows = new int[nz];
    }
  }

  ~Matrix() {
    delete[] val;
    delete[] rows;
    delete[] pointer;
  }
};

int createMatrix(int N, complex<double>** mat) {
  std::default_random_engine gen;
  std::uniform_int_distribution<int> dist(0, 10);
  *mat = new complex<double>[N*N];
  int nz_size = 0;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double re = (dist(gen) % 2)*dist(gen);
      double im = re*dist(gen);
      complex<double> val = complex<double>(re, im);
      (*mat)[i*N + j] = val;
      if (re != 0 || im != 0) {
        nz_size++;
      }
    }
  }
  return nz_size;
}

void toCCS(int N, const complex<double>* mat, Matrix* ccs) {
  for (int i = 0; i <= N; i++) {
    ccs->pointer[i] = -1;
  }

  int k = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (mat[j * N + i].real() != 0 || mat[j * N + i].imag() != 0) {
        ccs->val[k] = mat[j * N + i];
        ccs->rows[k] = j;

        if (ccs->pointer[i] < 0) {
          ccs->pointer[i] = k;
        }
        k++;
      }
    }
  }
  ccs->pointer[N] = k;

  for (int i = N; i >= 1; i--) {
    if (ccs->pointer[i - 1] < 0) {
      ccs->pointer[i - 1] = ccs->pointer[i];
    }
  }
}

void Transp(Matrix* a, int col) {
  int nz_size = a->nz_size;
  Matrix a_T(col, nz_size);

  for (int i = 0; i < nz_size; i++) {
    a_T.val[i] = complex<double>(0, 0);
  }
  for (int i = 0; i <= col; i++) {
    a_T.pointer[i] = 0;
  }

  for (int i = 0; i < nz_size; i++) {
    a_T.pointer[a->rows[i] + 1]++;
  }
  for (int i = 1; i < col + 1; i++) {
    a_T.pointer[i] += a_T.pointer[i - 1];
  }

  for (int i = 0; i < col; i++) {
    int j = a->pointer[i];
    int k = a->pointer[i + 1];

    for (; j < k; j++) {
      int pos = a_T.pointer[a->rows[j]];

      while (a_T.val[pos] != complex<double>(0, 0)) {
        pos++;
      }
      a_T.val[pos] = a->val[j];
      a_T.rows[pos] = i;
    }
  }

  std::memcpy(a->val, a_T.val, nz_size * sizeof(complex<double>));
  std::memcpy(a->rows, a_T.rows, nz_size * sizeof(int));
  std::memcpy(a->pointer, a_T.pointer, (col + 1) * sizeof(int));
}

void linMult(const Matrix& a, const Matrix& b, Matrix* c, int N) {
  c->pointer = new int[N + 1];
  for (int i = 0; i <= N; i++) {
    c->pointer[i] = 0;
  }
  vector<int> row;
  vector<complex<double> > vals;
  for (int i = 0; i < N; i++) {    //  cols b
    for (int j = 0; j < N; j++) {  //  cols a_T
      complex<double> sum = complex<double>(0, 0);

      for (int k = b.pointer[i]; k < b.pointer[i + 1]; k++) {
        for (int p = a.pointer[j]; p < a.pointer[j + 1]; p++) {
          if (a.rows[p] == b.rows[k]) {
            sum += a.val[p] * b.val[k];
            break;
          }
        }
      }


      if (sum != complex<double>(0, 0)) {
        vals.push_back(sum);
        row.push_back(j);
        c->pointer[i + 1]++;
      }
    }

    c->pointer[i + 1] += c->pointer[i];
  }

  int size = row.size();
  c->rows = new int[size];
  c->val = new complex<double>[size];
  for (int i = 0; i < size; i++) {
      c->rows[i] = row[i];
      c->val[i] = vals[i];
  }
  // std::copy(row.begin(), row.end(), c->rows);
  // std::copy(vals.begin(), vals.end(), c->val);
  c->nz_size = size;
}

int matrixMult(const complex<double>* a, const complex<double>* b, int N, complex<double>** c) {
  *c = new complex<double>[N*N];
  int nz_size = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      complex<double> sum = complex<double>(0, 0);
      for (int k = 0; k < N; k++) {
        sum += a[i * N + k] * b[k * N + j];
      }
      (*c)[i * N + j] = sum;

      if (sum != complex<double>(0, 0)) {
        nz_size++;
      }
    }
  }
  return nz_size;
}

bool compareCCS(const Matrix* a, const Matrix* b) {
  if (a->nz_size != b->nz_size) {
    return false;
  }
  for (int i = 0; i < a->nz_size; i++) {
    if (a->val[i] != b->val[i] || a->rows[i] != b->rows[i]) {
      return false;
    }
  }
  return true;
}

void ompMult(const Matrix* a, const Matrix* b, Matrix *c, int N) {
  c->pointer = new int[N + 1];
  for (int i = 0; i <= N; i++) {
      c->pointer[i] = 0;
  }
  vector<vector<int> > rows(N);
  vector<vector<complex<double> > > vals(N);

#pragma omp parallel for
  for (int i = 0; i < N; i++) {  //  cols b
    for (int j = 0; j < N; j++) {  //  cols a.T
      complex<double> sum = complex<double>(0, 0);

      for (int k = b->pointer[i]; k < b->pointer[i + 1]; k++) {
        for (int p = a->pointer[j]; p < a->pointer[j + 1]; p++) {
          if (a->rows[p] == b->rows[k]) {
            sum += a->val[p] * b->val[k];
            break;
          }
        }
      }

      if (sum != complex<double>(0, 0)) {
        vals[i].push_back(sum);
        rows[i].push_back(j);
        c->pointer[i + 1]++;
      }
    }
  }

  for (int i = 0; i < N; i++) {
    c->pointer[i + 1] += c->pointer[i];
  }

  c->val = new complex<double>[c->pointer[N]];
  c->rows = new int[c->pointer[N]];
  c->nz_size = c->pointer[N];
  size_t k = 0;
  for (int i = 0; i < N; i++) {
    size_t size = rows[i].size();
    std::memcpy(&(*c).rows[k], &rows[i][0], size * sizeof(int));
    std::memcpy(&(*c).val[k], &vals[i][0], size * sizeof(complex<double>));
    k += size;
  }
}

int main(int argc, char** argv) {
  int N = 5;
  int num_threads = 5;
  if (argc >= 3) {
    N = atoi(argv[1]);
    num_threads = atoi(argv[2]);
  }
  double sereal_time, parallel_time;
  std::cout << "Matrix A" << std::endl;
  complex<double>* a = 0, *b = 0, *c = 0;
  int nz_a = createMatrix(N, &a);
  std::cout << "Matrix B" << std::endl;
  int nz_b = createMatrix(N, &b);
  Matrix ccs_a(N, nz_a);
  Matrix ccs_b(N, nz_b);
  toCCS(N, a, &ccs_a);
  toCCS(N, b, &ccs_b);
  std::cout << "number of nonzero elements: " << ccs_a.pointer[N] << std::endl;

  std::cout << "###### Linear CCS version ######" << std::endl;
  sereal_time = omp_get_wtime();
  Transp(&ccs_a, N);
  Matrix ccs_c;
  linMult(ccs_a, ccs_b, &ccs_c, N);
  sereal_time = omp_get_wtime() - sereal_time;
  std::cout << "Linear time: " << sereal_time << std::endl;

  std::cout << "###### Linear standart version ######" << std::endl;
  int nz_c = matrixMult(a, b, N, &c);
  Matrix std_c(N, nz_c);
  toCCS(N, c, &std_c);
  std::cout << "Compare CCS and standart version " << compareCCS(&ccs_c, &std_c) << std::endl;

  std::cout << "###### Parallel version ######" << std::endl;
  omp_set_num_threads(num_threads);
  Matrix omp_c;
  parallel_time = omp_get_wtime();
  ompMult(&ccs_a, &ccs_b, &omp_c, N);
  parallel_time = omp_get_wtime() - parallel_time;
  std::cout << "Parallel time: " << parallel_time << std::endl;
  std::cout << "Compare parallel and linear version " << compareCCS(&ccs_c, &omp_c) << std::endl;
  std::cout << "Speedup: " << sereal_time/parallel_time << std::endl;
  return 0;
}
