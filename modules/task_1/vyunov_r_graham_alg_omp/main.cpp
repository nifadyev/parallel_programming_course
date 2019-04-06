// Copyright 2019 Roman Vyunov
#include<omp.h>
#include<iostream>
#include<random>
#include<cstring>

struct Point {
    double x;
    double y;
};

struct Stack {
    int *s;
    int size;
    int max_size;
    explicit Stack(int max_size) {
        if (max_size <= 0) {
            std::cerr << "Stack create error: negative size" << std::endl;
        } else {
            this->size = 0;
            this->max_size = max_size;
            this->s = new int[max_size];
        }
    }

    bool push(int val) {
        if (size < max_size) {
            s[size++] = val;
            return true;
        } else {
            std::cerr << "Stack error: stack overflow" << std::endl;
            return false;
        }
    }

    int pop() {
        return s[size--];
    }

    int operator[] (int n) {
        return s[n];
    }
};

inline double rotate(Point a, Point b, Point c) {
    return (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x);
}

Stack graham(const Point *a, const int N) {
    int *p = new int[N];
    for (int i = 0; i < N; i++) {
        p[i] = i;
    }
    for (int i = 1; i < N; i++) {
        if (a[p[i]].x < a[p[0]].x) {
            int tmp = p[i];
            p[i] = p[0];
            p[0] = tmp;
        }
    }
    for (int i = 2; i < N; i++) {
        int j = i;
        while (j > 1 && rotate(a[p[0]], a[p[j - 1]], a[p[j]]) < 0) {
            int tmp = p[j];
            p[j] = p[j - 1];
            p[j - 1] = tmp;
            j--;
        }
    }
    Stack s(N);
    s.push(p[0]);
    s.push(p[1]);
    for (int i = 2; i < N; i++) {
        while (rotate(a[s[s.size - 2]], a[s[s.size - 1]], a[p[i]]) < 0) {
            s.pop();
        }
        s.push(p[i]);
    }
    return s;
}

void generateArrayOfPoints(Point *p, int N, double field) {
    std::random_device rd;
    unsigned int rd_range = rd.max() - rd.min();
    for (int i = 0; i < N; i++) {
        p[i].x = ((static_cast<double>(rd()) - rd.max() / 2) * 2) / rd_range * field;
        p[i].y = ((static_cast<double>(rd()) - rd.max() / 2) * 2) / rd_range * field;
    }
}

int main(int argc, char *argv[]) {
    bool verbose_out = false;
    Point *p = nullptr;
    unsigned int N = 0;
    std::random_device rd;
    if (argc == 1) {
        N = 30;
    } else if (argc > 1) {
        N = atoi(argv[1]);
        if (argc > 2 && std::strcmp(argv[2], "-v") == 0) {
            verbose_out = true;
        }
    } else {
        std::cerr << "Incorrect number of arguments" << std::endl;
        return -1;
    }
    p = new Point[N];
    generateArrayOfPoints(p, N, 10);
    if (verbose_out) {
        std::cout << "Number of points: " << N << std::endl;
        std::cout << "Input array of points:" << std::endl;
        for (unsigned int i = 0; i < N; i++) {
            std::cout << i + 1 << " : " << p[i].x << "\t" << p[i].y << std::endl;
        }
    }
    double start = omp_get_wtime();
    Stack result = graham(p, N);
    double end = omp_get_wtime();
    if (verbose_out) {
        std::cout << "Resuts:" << std::endl;
        for (int i = 0; i < result.size; i++) {
            int tmp = result[i];
            std::cout << i + 1 << " : " << tmp + 1 << " : " << p[tmp].x << "\t" << p[tmp].y << std::endl;
        }
    }
    std::cout << "Execution time: " << end - start << std::endl;
    delete[] p;
    return 0;
}
