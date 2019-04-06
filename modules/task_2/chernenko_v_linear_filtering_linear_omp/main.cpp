// Copyright 2019 Chernenko Valery
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
// #include <opencv2\core\core.hpp>
// #include <opencv2\highgui\highgui.hpp>

/*
unsigned char** getArrayInImage(cv::Mat image, int flag = 0) {
    unsigned char** result = new unsigned char*[image.rows];
    for (int i = 0; i < image.rows; i++)
        result[i] = new unsigned char[image.cols];

    if (flag < 0 || flag > 2) flag = 0;
    for (int i = 0; i < image.rows; i++)
        for(int j = 0; j < image.cols; j++)
            result[i][j] = image.data[(i * image.cols + j) * 3 + flag];
    return result;
}
*/

unsigned char** getRandImage(int width, int height) {
    unsigned char** image = new unsigned char*[height];
    for (int i = 0; i < height; i++) {
        image[i] = new unsigned char[width];
        for (int j = 0; j < width; j++)
            image[i][j] = static_cast<unsigned char>(std::rand() % 256);
    }
    return image;
}

int clamp(int x, int min, int max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

double** getGausseFilter(int n) {
    int rad = n / 2;
    double** result = new double*[n];
    for (int i = 0; i < n; i++)
        result[i] = new double[n];
    double coeff1 = 1 / (2 * 3.14159265 * rad * rad);
    for (int i = -rad; i <= rad; i++) {
        for (int j = -rad; j <= rad; j++) {
            result[i+rad][j+rad] = coeff1*exp((-(i*i + j*j))/(2 * rad * rad));
        }
    }

    double coeff2 = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            coeff2 += result[i][j];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i][j] /= coeff2;

    return result;
}

unsigned char getResultPixel(unsigned char** image, int width, int height, int x, int y, double** filter, int n) {
    double result = 0;
    int rad = n / 2;
    for (int i = -rad; i <= rad; i++)
        for (int j = -rad; j <= rad; j++) {
            double r = image[clamp(x + i, 0, height - 1)][clamp(y + j, 0, width - 1)] * filter[i + rad][j + rad];
            result += r;
        }
    return clamp(static_cast<int>(result), 0, 255);
}

unsigned char* getResultRows(unsigned char** image, int width, int height, int x, double** filter, int n) {
    unsigned char* resultRow = new unsigned char[width];
    for (int j = 0; j < width; j++) {
        resultRow[j] = getResultPixel(image, width, height, x, j, filter, n);
    }
    return resultRow;
}

unsigned char** getResultImage(unsigned char** image, int width, int height, double** filter, int n) {
    unsigned char** resultImage = new unsigned char*[height];
    int i = 0;

    #pragma omp parallel num_threads(4) shared(image, width, height, filter, n) private(i)
    #pragma omp for schedule(static)
    for (i = 0; i < height; i++)
        resultImage[i] = getResultRows(image, width, height, i, filter, n);
    return resultImage;
}

/*
void getResultImage(unsigned char** imageR, unsigned char**& resultImageR, unsigned char** imageG, unsigned char**& resultImageG, unsigned char** imageB, unsigned char**& resultImageB, int width, int height, double** filter, int n) {
    resultImageR = new unsigned char*[height];
    resultImageG = new unsigned char*[height];
    resultImageB = new unsigned char*[height];
    int i = 0;

    #pragma omp parallel num_threads(4) shared(imageR, imageG, imageB, resultImageR, resultImageG, resultImageB, width, height, filter, n) private(i) 
    #pragma omp for schedule(static)
    for (i = 0; i < height; i++) {
        resultImageR[i] = getResultRows(imageR, width, height, i, filter, n);
        resultImageG[i] = getResultRows(imageG, width, height, i, filter, n);
        resultImageB[i] = getResultRows(imageB, width, height, i, filter, n);
    }
}
*/

int main(int argv, char** argc) {
    srand((unsigned int)time(0));
    setlocale(LC_ALL, "Russian");

    unsigned char** imageR;
    unsigned char** imageG;
    unsigned char** imageB;
    int width, height;

    int n = 5;
    if (n % 2 == 0) n++;
    double** filter = getGausseFilter(n);

//    char* path = "V:\\Projects\\C++\\Files\\image.jpg";
//    cv::Mat image = cv::imread(path, 1);

    if (/*image.cols == 0 || image.rows == 0*/true) {
        std::cout << "File not found! Generate rand image" << std::endl;
        width = 20;
        height = 20;
        imageR = getRandImage(width, height);
        imageG = getRandImage(width, height);
        imageB = getRandImage(width, height);
        std::cout << "Image generated" << std::endl;
    } else {
//        width = image.cols;
//        height = image.rows;
//        imageR = getArrayInImage(image, 0);
//        imageG = getArrayInImage(image, 1);
//        imageB = getArrayInImage(image, 2);
    }

    std::cout << "Start work" << std::endl;

//    unsigned char** resultImageR;
//    unsigned char** resultImageG;
//    unsigned char** resultImageB;

    double startTime = omp_get_wtime();
    // getResultImage(imageR, resultImageR, imageG, resultImageG, imageB, resultImageB, width, height, filter, n);

    /*unsigned char** resultImageR = */getResultImage(imageR, width, height, filter, n);
    /*unsigned char** resultImageG = */getResultImage(imageG, width, height, filter, n);
    /*unsigned char** resultImageB = */getResultImage(imageB, width, height, filter, n);
    double finishTime = omp_get_wtime();
    std::cout << "Finish work. Elapsed time = " << finishTime - startTime << std::endl;

/*
    cv::Mat resultImage = image.clone();
    std::cout << ResultImageR[0][0];

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            resultImage.data[3 * (i*width + j)] = resultImageR[i][j];
            resultImage.data[3 * (i*width + j) + 1] = resultImageG[i][j];
            resultImage.data[3 * (i*width + j) + 2] = resultImageB[i][j];
        }
    }

    cv::imwrite("V:\\Projects\\C++\\Files\\resultImage.bmp", resultImage);
*/
}
