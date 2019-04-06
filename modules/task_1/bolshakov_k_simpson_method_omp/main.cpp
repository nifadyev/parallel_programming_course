// Copyright 2019 Bolshakov Konstantin
#include <stdio.h>
#include <math.h>
#include <iostream>

double func1(double x, double y) {
    return x + y;
}

double func2(double x, double y) {
    return x * x + y * y;
}

double func3(double x, double y) {
    return 2;
}

double func4(double x, double y) {
    return sin(x * y) +  exp(cos(y * x));
}

double func_simpson(double f(double, double), double xp, double yp, double x, double y,
            double xn, double yn) {
    return (f(xp, yp) + 4 * f(x, y) + f(xn, yn));
}

double integrate(double f(double, double),
                      double x1, double y1, double hx, double hy, int n) {
    double* x = new double[n + 1];
    double* y = new double[n + 1];
    double sum = 0;
    for (int i = 0; i <= n; i++) {
        x[i] = x1 + i * hx;
        y[i] = y1 + i * hy;
    }
    for (int i = 1; i < n; i+=2) {
        for (int j = 1; j < n; j+=2) {
            sum += func_simpson(f, x[i - 1], y[j - 1], x[i - 1], y[j], x[i - 1], y[j + 1])  +
                    4 * func_simpson(f, x[i], y[j - 1], x[i], y[j], x[i], y[j + 1]) +
                    func_simpson(f, x[i + 1], y[j - 1], x[i + 1], y[j], x[i + 1], y[j + 1]);
        }
    }
    return sum * hx * hy / 9;
}


int main(int argc, char* argv[]) {
    std::cout << integrate(func4, 0, 0, 0.1, 0.1, 30) << std::endl;
    return 0;
}
