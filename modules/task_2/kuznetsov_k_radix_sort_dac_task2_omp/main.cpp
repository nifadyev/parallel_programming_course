// Copyright 2019 Kuznetsov Konstantin
// параллельная поразрядная сортировка целых положительных чисел
// сортировка разными потоками со слиянием "Разделяй и властвуй"

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>
#include <vector>
#include <utility>

void print_data(const std::vector<unsigned int>& data, size_t n = 10) {
    if (data.size() < n) return;

    for (size_t i = 0; i < n; ++i)
        std::cout << data[i] << " ";
    std::cout << std::endl;
}

static std::vector<unsigned int> generate_data(size_t size) {
    using value_type = unsigned int;

    static std::uniform_int_distribution<value_type> distribution(
                std::numeric_limits<value_type>::min(),
                std::numeric_limits<value_type>::max());
    static std::default_random_engine generator;

    std::vector<value_type> data(size);
    std::generate(data.begin(), data.end(), []() {
        return distribution(generator); });

    return data;
}

void create_counters(unsigned int* data, unsigned int* counters,
                     size_t size) {
    unsigned char* byte_ptr = reinterpret_cast<unsigned char*>(data);
    unsigned char* data_end = reinterpret_cast<unsigned char*>(data + size);
    for (size_t i = 0; i < sizeof(int) * 256; ++i)
        counters[i] = 0;

    size_t i;
    while (byte_ptr != data_end)
        for (i = 0; i < sizeof(int); ++i)
            counters[256 * i + *byte_ptr++]++;
}

void byte_sort(unsigned int* data, unsigned int* temp,
               unsigned int* counter, size_t byte, size_t size) {
    unsigned char* byte_ptr;
    unsigned int c = 0;

    unsigned int sum = 0;
    for (size_t i = 0; i < 256; ++i) {
        c = counter[i];
        counter[i] = sum;
        sum += c;
    }

    byte_ptr = reinterpret_cast<unsigned char*>(data) + byte;
    for (size_t i = 0; i < size; ++i, byte_ptr += sizeof(int))
        temp[counter[*byte_ptr]++] = data[i];
}

// поразрядная сортировка целых положительных чисел
// сортировка выполняется по байтам
// для сортировки по каждому байту используется сортировка подсчетом
void radix_sort(unsigned int* data, size_t size) {
    unsigned int* temp = new unsigned int[size];
    unsigned int* counters = new unsigned int[sizeof(int) * 256];
    unsigned int* counter;

    create_counters(data, counters, size);
    for (size_t i = 0; i < sizeof(int); ++i) {
        counter = counters + 256 * i;
        byte_sort(data, temp, counter, i, size);
        std::swap(data, temp);
    }

    delete[] temp;
    delete[] counters;
}

void merge(unsigned int* a, unsigned int* b, unsigned int* c,
           size_t n, size_t m) {
    size_t i, j, k;
    i = j = k = 0;

    while (i < n && j < m) {
        if ((a[i]) < (b[j]))
            c[k++] = a[i++];
        else
            c[k++] = b[j++];
    }

    while (i < n) c[k++] = a[i++];
    while (j < m) c[k++] = b[j++];
}

size_t bin_search(unsigned int* data, size_t l, size_t r, unsigned int x) {
    if (l == r) return l;
    if (l + 1 == r) {
        if (x < data[l]) return l;
        else
            return r;
    }

    size_t m = (l + r) / 2;
    if (x < data[m]) {
        r = m;
    } else {
        if (x > data[m]) l = m;
        else
            return m;
    }

    return bin_search(data, l, r, x);
}

void dac_radix_sort(unsigned int* data, size_t size, size_t num_threads) {
    // массив разбивается на части для сортировки на разных потоках
    size_t portion = size / num_threads;
    omp_set_num_threads(static_cast<int>(num_threads));

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < static_cast<int>(num_threads); ++i)
        radix_sort(data + i*static_cast<int>(portion), portion);

    // слияние "Разделяй и властвуй"
    size_t k_max = static_cast<size_t>(std::log2(num_threads));
    for (size_t k = 0; k < k_max; ++k) {
        if (k != 0) portion = portion * 2;
        for (size_t j = 0; j < std::pow(2, k_max - k - 1); ++j) {
            unsigned int* curr = data + portion * j * 2;
            unsigned int* temp = new unsigned int[portion * 2];
            size_t s = portion;
            s /= num_threads;

            std::vector<size_t> l1(num_threads, 0);
            std::vector<size_t> r1(num_threads + 1, s);
            std::vector<size_t> l2(num_threads, 0);
            std::vector<size_t> r2(num_threads + 1, 0);
            r1[num_threads - 1] = portion;
            r2[num_threads - 1] = portion;

            for (size_t i = 0; i < num_threads - 1; ++i) {
                unsigned int x = curr[r1[i]];
                r2[i] = bin_search(curr + portion, 0, portion * 2 - portion, x);

                l1[i+1] = l1[i] + s;
                r1[i+1] = r1[i] + s;
                l2[i+1] = r2[i];
            }

            #pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i <= static_cast<int>(num_threads) - 1; ++i)
                merge(curr + l1[i], curr + portion + l2[i],
                    temp + l1[i] + l2[i], r1[i] - l1[i], r2[i] - l2[i]);

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(portion) * 2; ++i)
                curr[i] = temp[i];

            delete[] temp;
        }
    }
}

bool check_sorting(unsigned int *data1, unsigned int *data2, size_t size) {
    for (size_t i = 0; i < size - 1; ++i) {
        if ((data1[i] > data1[i + 1]) || (data1[i] != data2[i]))
            return false;
    }

    if (data1[size - 1] != data2[size - 1])
        return false;

    return true;
}

int main(int argc, char* argv[]) {
    size_t data_size = 160000;
    size_t num_threads = 4;
    double t1, t2;

    std::vector<unsigned int> data1 = generate_data(data_size);
    std::vector<unsigned int> data2 = data1;

    t1 = omp_get_wtime();
    dac_radix_sort(data1.data(), data_size, 1);
    t2 = omp_get_wtime();
    std::cout << "non-parallel time = " << t2 - t1 << std::endl;

    t1 = omp_get_wtime();
    dac_radix_sort(data2.data(), data_size, num_threads);
    t2 = omp_get_wtime();
    std::cout << "parallel time (" << num_threads << " threads) = " <<
                 t2 - t1 << std::endl;

    if (check_sorting(data1.data(), data2.data(), data_size))
        std::cout << "correct" << std::endl;
    print_data(data1);
    print_data(data2);
    return 0;
}

