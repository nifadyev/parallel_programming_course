// Copyright Ermachenko Boris 2019
// #pragma comment(linker, "/STACK:8000000")
// #pragma comment(linker, "/HEAP:8000000")
#include <math.h>
#include <omp.h>
#include <iostream>
#include <ctime>
#define THREADS 7
const double PI = 3.1415;
int blockX = 0;
int blockY = 0;
double** getNewArr(int width, int height) {
    double** new_arrImage = new double*[height];
    for (int i = 0; i < height; i++) {
        new_arrImage[i] = new double[width];
    }
    return new_arrImage;
}
double** LinierFilterGauss(double** arrImage, double** new_arrImage,
    int width, int height) {
    int n = 3;
    double** w = getNewArr(n, n);
    double sigma = 1.0;
    double r;
    double s = 2.0 * sigma * sigma;
    for (int i = -n / 2; i < n / 2; i++) {
        for (int j = -n / 2; j < n / 2; j++) {
            r = sqrt(i * i + j * j);
            w[i + n / 2][j + n / 2] = (exp(-(r*r) / s)) / (PI * s);
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            new_arrImage[i][j] = arrImage[i][j];
        }
    }
    for (int yi = 0; yi < height ; yi++) {
        for (int xj = 0; xj < width ; xj++) {
            double color = 0;
            double kSum = 0;
            for (int i = -n / 2; i < n / 2; i++) {
                for (int j = -n / 2; j < n / 2; j++) {
                    kSum += w[i + n / 2][j + n / 2];
                    if ((i + yi) >= 0 && (i + yi) < height
                        && (j + xj) >= 0 && (j + xj) < width)
            color +=arrImage[static_cast<int>(i + yi)][static_cast<int>(j + xj)]
                        * w[i + n / 2][j + n / 2];
                }
            }
            if (kSum <= 0) kSum = 1;
            color /= kSum;
            if (color < 0) color = 0;
            if (color > 255) color = 255;
            if ((yi) >= 0 && (yi) < height && (xj) >= 0 && (xj) < width)
            new_arrImage[static_cast<int>(yi)][static_cast<int>(xj)] = color;
        }
    }
    return new_arrImage;
}
void getBlockNum(int throwds) {
    int tmp = static_cast<int>(sqrt(throwds));
    while (throwds % tmp != 0) {
        tmp--;
    }
    blockX = tmp;
    blockY = throwds / tmp;
}
int* arrayLengthBlocksX(int width) {
    int* rezArrX = new int[blockX];
    int rest = width % blockX;
    for (int i = 0; i < blockX; i++) {
        rezArrX[i] = width / blockX;
    }
    if (rest != 0) {
        int k = 0;
        while (rest != 0) {
            rezArrX[k%rest]++;
            rest--;
            k++;
        }
    }
    return rezArrX;
}
int * arrayLengthBlocksY(int height) {
    int* rezArrY = new int[blockY];
    int rest = height % blockY;
    for (int i = 0; i < blockY; i++) {
        rezArrY[i] = height / blockY;
    }
    if (rest != 0) {
        int k = 0;
        while (rest != 0) {
            rezArrY[k%rest]++;
            rest--;
            k++;
        }
    }
    return rezArrY;
}
double** ParallelFilterGauss(double** arrImage,
    double** new_arrImage, int width, int height) {
    int n = 3;
    double** w = getNewArr(n, n);
    double sigma = 1.0;
    double r;
    double s = 2.0 * sigma * sigma;
    for (int i = -n / 2; i < n / 2; i++) {
        for (int j = -n / 2; j < n / 2; j++) {
            r = sqrt(i * i + j * j);
            w[i + n / 2][j + n / 2] = (exp(-(r*r) / s)) / (PI * s);
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            new_arrImage[i][j] = arrImage[i][j];
        }
    }
    getBlockNum(THREADS);
    int* lengthX = arrayLengthBlocksX(width);
    int* lengthY = arrayLengthBlocksY(height);
#pragma omp parallel
    {
        int numberThroat = omp_get_thread_num();
        int startX = 0;
        int startY = 0;
        int row = 0;
        int column = 0;
        int count = 0;
        for (int i = 0; i < blockY; i++) {
            if (i > 0) startY += lengthY[i-1];
            for (int j = 0; j < blockX; j++) {
                if (j > 0) startX += lengthX[j-1];
                count++;
                if (count > numberThroat)break;
                column++;
            }
            if (count > numberThroat)break;
            row++;
            column = 0;
            startX = 0;
        }
        for (int yi = startY; yi < startY + lengthY[row]; yi++) {
            for (int xj = startX; xj < startX + lengthX[column]; xj++) {
                        double color = 0;
                        double kSum = 0;
                        for (int i = -n / 2; i < n / 2; i++) {
                            for (int j = -n / 2; j < n / 2; j++) {
                                kSum += w[i + n / 2][j + n / 2];
                                if ((i + yi) >= 0 && (i + yi) < height
                                    && (j + xj) >= 0 && (j + xj) < width)
            color += arrImage[static_cast<int>(i + yi)][static_cast<int>(j + xj)]
                                    * w[i + n / 2][j + n / 2];
                            }
                        }
                        if (kSum <= 0) kSum = 1;
                        color /= kSum;
                        if (color < 0) color = 0;
                        if (color > 255) color = 255;
                        if ((yi) >= 0 && (yi) < height &&
                            (xj) >= 0 && (xj) < width)
                        new_arrImage[static_cast<int>(yi)][static_cast<int>(xj)]
                            = color;
            }
        }
    }
        return new_arrImage;
}
double** getImage(int width, int height) {
    srand((unsigned int)time(0));
    double** arrImage = getNewArr(width, height);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            arrImage[i][j] = static_cast<int>(std::rand()) % 256;
        }
    }
    return arrImage;
}
void ShowArr(double** arr, int width, int height) {
    std::cout << "\n";
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << "\n";
    }
}
int main() {
    double tStart, tEnd;
    int height = 10;
    int width = 20;
    double** arrImage = getImage(width, height);
    double** new_arrImage_Liner = getNewArr(width, height);
    omp_set_num_threads(THREADS);
    tStart = omp_get_wtime();
    new_arrImage_Liner = LinierFilterGauss(arrImage,
        new_arrImage_Liner, width, height);
    tEnd = omp_get_wtime();
    printf("Time linier :  %.4lf \n", tEnd - tStart);
    double** new_arrImage_Parallel = getNewArr(width, height);
    tStart = omp_get_wtime();
    new_arrImage_Parallel = ParallelFilterGauss(arrImage,
        new_arrImage_Parallel, width, height);
    tEnd = omp_get_wtime();
    printf("Time parallel :  %.4lf \n", tEnd - tStart);
    // ShowArr(arrImage, width, height);
    // ShowArr(new_arrImage_Liner, width, height);
    // ShowArr(new_arrImage_Parallel, width, height);
    return 0;
}
