// Copyright 2019 Panova Elena
#include <iostream>
#include <random>
#include "tbb/task_scheduler_init.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/tick_count.h"

const int N = 1024;  // default matrix size
const double EPS = 1e-5;  // default error

static int GRAIN_SIZE = 32;
static int NUM_THREADS = 0;

class Vector {
    int size = 0;
    double *data = 0;

 public:
    Vector() {
        data = 0;
    }
    explicit Vector(int _size) {
        Initialize(_size);
    }
    Vector(const Vector& v) {
        size = v.size;
        data = new double[size];
        for (int i = 0; i < size; i++)
            data[i] = v.data[i];
    }
    ~Vector() {
        Clear();
    }
    void Initialize(int _size) {
        size = _size;
        Clear();
        data = new double[size];
        for (int i = 0; i < size; i++)
            data[i] = 0;
    }
    void Clear() {
        if (data)
            delete[] data;
        data = 0;
    }
    int getSize() const {
        return size;
    }
    double& operator[](int i) {
        return data[i];
    }
    double get(int i) const {
        return data[i];
    }
    Vector& operator=(const Vector& v) {
        if (this != &v) {
            if (v.size != size) {
                Clear();
                data = new double[v.size];
                size = v.size;
            }
            for (int i = 0; i < size; i++)
                data[i] = v.data[i];
        }
        return *this;
    }
    Vector operator+(const Vector& b) const {
        Vector res(size);
        for (int i = 0; i < size; i++)
            res[i] = data[i] + b.data[i];
        return res;
    }
    friend Vector operator*(double a, const Vector& b) {
        Vector res(b.size);
        for (int i = 0; i < b.size; i++)
            res[i] = a * b.data[i];
        return res;
    }
    Vector operator-(const Vector& b) const {
        return (*this) + (-1)*b;
    }
    friend double Dot(const Vector& a, const Vector& b) {
        double res = 0;
        for (int i = 0; i < a.size; i++)
            res = res + a.data[i] * b.data[i];
        return res;
    }
    double Norm() const {
        double max = 0;
        for (int i = 0; i < size; i++)
            if (max < std::abs(data[i]))
                max = std::abs(data[i]);
        return max;
    }
};

class Matrix {
    int size = 0;
    double *data = 0;

 public:
    Matrix() {
        data = 0;
    }
    explicit Matrix(int _size) {
        Initialize(_size);
    }
    Matrix(const Matrix& m) {
        Initialize(m.size);
        for (int i = 0; i < size*size; i++)
            data[i] = m.data[i];
    }
    ~Matrix() {
        Clear();
    }
    void Initialize(int _size) {
        size = _size;
        Clear();
        data = new double[size*size];
        for (int i = 0; i < size*size; i++)
            data[i] = 0;
    }
    void Clear() {
        if (data)
            delete[] data;
        data = 0;
    }
    Matrix& operator=(const Matrix& v) {
        if (this != &v) {
            if (v.size != size) {
                Clear();
                data = new double[v.size*v.size];
                size = v.size;
            }
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    data[i] = v.data[i];
        }
        return *this;
    }
    double& operator()(int i, int j) const {
        return data[i*size + j];
    }

    Vector operator*(const Vector& v) const;
};

class LoopBodyFunctor {
    const Matrix* m;
    const Vector* v;
    Vector* const res;

 public:
    LoopBodyFunctor(const Matrix* _m, const Vector* _v, Vector* _res) :
        m(_m), v(_v), res(_res) {}

    double dot(int i) const {
        double r = 0;
        for (int j = 0; j < v->getSize(); j++)
            r += (*m)(i, j)*v->get(j);
        return r;
    }

    void operator()(const tbb::blocked_range<int>& r) const {
        for (int i = r.begin(); i < r.end(); i++)
            (*res)[i] = dot(i);
    }
};

Vector Matrix::operator*(const Vector& v) const {
    Vector res(size);
    tbb::parallel_for(tbb::blocked_range<int>(0, size, GRAIN_SIZE),
        LoopBodyFunctor(this, &v, &res));
    return res;
}

Vector conjGradientMethod(const Matrix& A, const Vector& b, const Vector& x0,
    double eps, Vector* rn, int* counter) {
    Vector x = x0;
    Vector r = b - A * x0;
    Vector z = r;
    while (r.Norm() / b.Norm() > eps) {
        Vector Az = A * z;
        double alpha = Dot(r, r) / Dot(Az, z);
        x = x + alpha * z;
        r = r - alpha * Az;
        Vector r1 = b - A * x;
        double beta = -Dot(Az, r) / Dot(Az, z);
        z = r + beta * z;
        (*counter)++;
    }
    *rn = r;
    return x;
}

const int NUM_ARGS_MIN = 5;
bool parseArgs(int argc, char* argv[], Matrix* A, Vector* b, double* eps,
    tbb::task_scheduler_init* init) {
    int size = 0;
    if (argc < NUM_ARGS_MIN) {
        NUM_THREADS = tbb::task_scheduler_init::automatic;
        init->initialize(NUM_THREADS);
        size = N;
        *eps = EPS;
    } else {
        NUM_THREADS = atoi(argv[1]);
        init->initialize(NUM_THREADS);
        GRAIN_SIZE = atoi(argv[2]);
        size = atoi(argv[3]);
        *eps = atof(argv[4]);
    }
    A->Initialize(size);
    b->Initialize(size);
    if (argc > NUM_ARGS_MIN) {  // read matrix and vector
        if (argc - NUM_ARGS_MIN != size * size + size) {
            std::cout << "wrong args\n";
            return false;
        }
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                (*A)(i, j) = atof(argv[i*size + j + NUM_ARGS_MIN]);
        for (int i = 0; i < size; i++)
            (*b)[i] = atof(argv[size*size + i + NUM_ARGS_MIN]);
    } else {  // random matrix and vector
        std::uniform_real_distribution<double> undist(0, 1);
        std::default_random_engine re;
        for (int i = 0; i < size; i++)
            for (int j = 0; j <= i; j++)
                (*A)(j, i) = (*A)(i, j) = undist(re);
        for (int i = 0; i < size; i++)
            (*b)[i] = undist(re);
    }
    return true;
}

void printRes(const Vector& res, const Matrix& A, const Vector& b,
    const Vector& r, int c, int argc) {
    std::cout << "Number of threads is " << NUM_THREADS << "\n";
    std::cout << "Grain size is " << GRAIN_SIZE << "\n";
    std::cout << "Matrix size is " << b.getSize() << "\n";
    std::cout << "Relative discrepancy is " <<
        r.Norm() / b.Norm() << "\n";
    std::cout << "Number of iterations is " << c << "\n";
    if (argc > NUM_ARGS_MIN) {
        for (int i = 0; i < res.getSize(); i++)
            std::cout << res.get(i) << " ";
        std::cout << "\n";
    }
}

int main(int argc, char* argv[]) {
    Matrix A;
    Vector b;
    double eps;
    tbb::task_scheduler_init init(tbb::task_scheduler_init::deferred);
    if (parseArgs(argc, argv, &A, &b, &eps, &init)) {
        Vector x0(b.getSize());
        int counter = 0;
        Vector r(b.getSize());
        tbb::tick_count t1 = tbb::tick_count::now();
        Vector res = conjGradientMethod(A, b, x0, eps, &r, &counter);
        tbb::tick_count t2 = tbb::tick_count::now();
        std::cout << "Time of TBB version is " << (t2-t1).seconds() << std::endl;
        printRes(res, A, b, r, counter, argc);
    }
}
