// Copyright 2019 Muravev Denis
// Hoar sorting with simple merge
#include <omp.h>
#include <ctime>
#include <iostream>
#include <climits>
#include <cstring>
#define TYPE_MAS double
#define TYPE_SIZE int

TYPE_MAS* Create_mas(TYPE_SIZE _size) {
    if (_size < 1)
        return NULL;
    TYPE_MAS* _mas = new TYPE_MAS[_size];
    TYPE_SIZE left = 0, right = INT_MAX;
    TYPE_SIZE i;
    for (i = 0; i < _size; i++)
        _mas[i] = std::rand() % (left - right + 1) + left;
    return  _mas;
}

void Show_mas(TYPE_MAS* _mas, TYPE_SIZE _size) {
    if (_mas == NULL || _size < 1)
        return;
    TYPE_SIZE i;
    int delta = 7;
    if (_size < 50) {
        std::cout << "Array [" << _size << "]" << std::endl;
        for (i = 0; i < _size; i++) {
            std::cout.width(delta);
            std::cout << _mas[i];
        }
    }
}

int Check_not_decreasing(TYPE_MAS* _mas, TYPE_SIZE _size) {
    TYPE_SIZE i;
    for (i = 1; i < _size; i++)
        if (_mas[i - 1] > _mas[i])
            return 0;
    return 1;
}

void Quick_sort(TYPE_MAS * mas, TYPE_SIZE l, TYPE_SIZE r) {
    TYPE_SIZE i = l, j = r;
    TYPE_MAS e, t;
    e = mas[(r + l) / 2];

    while (i <= j) {
        while (mas[i] < e)
            i++;
        while (mas[j] > e)
            j--;

        if (i <= j) {
            if (i < j) {
                t = mas[i];
                mas[i] = mas[j];
                mas[j] = t;
            }
            i++;
            j--;
        }
    }

    if (j > l)
        Quick_sort(mas, l, j);
    if (r > i)
        Quick_sort(mas, i, r);
}

int* Create_reverse_notation(int _n, int _k) {
    int higher = -1;
    int* _rev = new int[_n];
    _rev[0] = 0;
    for (int i = 1; i < _n; i++) {
        if ((i & (i - 1)) == 0)
            higher++;
        _rev[i] = _rev[i ^ (1 << higher)];
        _rev[i] |= 1 << (_k - higher - 1);
    }
    return _rev;
}

struct BorderList {
    TYPE_SIZE l = 0;
    TYPE_SIZE r = 0;
};

BorderList* Create_border_list(TYPE_SIZE _size, int _num_thr, int _n, int _k) {
    BorderList* _border = new BorderList[_num_thr];
    TYPE_SIZE delta = _size / _num_thr;
    int *rev = Create_reverse_notation(_n, _k);

    int k = 0;
    for (int i = 0; i < _n; i++) {
        if (rev[i] >= _num_thr)
            continue;
        _border[rev[i]].l = k * delta;
        k++;
        _border[rev[i]].r = k * delta - 1;
    }

    if (_n == _num_thr)
        _border[_n - 1].r = _size - 1;
    else
        _border[rev[_n - 2]].r = _size - 1;

    delete[] rev;
    return _border;
}

void Merge(TYPE_MAS* mas, TYPE_MAS* tmp_mas, TYPE_SIZE l1, TYPE_SIZE r1, TYPE_SIZE l2, TYPE_SIZE r2) {
    TYPE_SIZE i = l1, j = l2, k = l1;
    while ((i <= r1) && (j <= r2)) {
        if (mas[i] < mas[j])
            tmp_mas[k++] = mas[i++];
        else
            tmp_mas[k++] = mas[j++];
    }
    if (i > r1) {
        for (; j <= r2; j++)
            tmp_mas[k++] = mas[j];
    } else {
        for (; i <= r1; i++)
            tmp_mas[k++] = mas[i];
    }
    for (i = l1; i <= r2; i++)
        mas[i] = tmp_mas[i];
}

void Task_sort(TYPE_MAS* _mas, TYPE_SIZE _size, TYPE_MAS* _tmp_mas, int _num_thr) {
    int numd_threads;
    int n = 1, k = 0;
    BorderList *border;
    #pragma omp parallel num_threads(_num_thr)
    {
    #pragma omp single
    {
        numd_threads = omp_get_num_threads();
        if (numd_threads != _num_thr) {
            std::cout << "Warning! Maximum number of threads is " << numd_threads;
            std::cout << ". Further calculations will be made with this value." << std::endl;
        }
        while (n < numd_threads) {
            n = n << 1;
            k++;
        }
        border = Create_border_list(_size, numd_threads, n, k);

        std::cout << "Size array: " << _size << "; Numb Threads: " << omp_get_num_threads() << std::endl;
    }
    #pragma omp barrier
    int rank = omp_get_thread_num();
    Quick_sort(_mas, border[rank].l, border[rank].r);
    #pragma omp barrier
    int max_rank = 0;
    int source;
    if (k > 0) {
        max_rank = 1 << (k - 1);
        source = rank ^ max_rank;
        if ((rank < max_rank) && (source < numd_threads)) {
                Merge(_mas, _tmp_mas, border[rank].l, border[rank].r, border[source].l, border[source].r);
                border[rank].r = border[source].r;
        }
    }
    for (int i = k - 2; i >= 0; i--) {
        #pragma omp barrier
        source = rank ^ (1 << i);
        max_rank = max_rank >> 1;
        if (rank < max_rank) {
            Merge(_mas, _tmp_mas, border[rank].l, border[rank].r, border[source].l, border[source].r);
            border[rank].r = border[source].r;
        }
    }
    }
}

int main(int argc, char * argv[]) {
    TYPE_MAS* mas = NULL;
    TYPE_SIZE size = (1 << 20);
    double time_sort = 0;
    int num_thr = 4;

    if (argc == 3 || argc == 5) {
        if (strcmp(argv[1], "-size") == 0)
            size = atoi(argv[2]);
        if (strcmp(argv[1], "-n") == 0)
            num_thr = atoi(argv[2]);
        if (argc == 5) {
            if (strcmp(argv[3], "-size") == 0)
                size = atoi(argv[4]);
            if (strcmp(argv[3], "-n") == 0)
                num_thr = atoi(argv[4]);
        }
    } else {
        std::cout << "Initial data is missing or entered incorrectly" << std::endl;
    }

    std::srand((unsigned)time(NULL));
    mas = Create_mas(size);
    if (mas == NULL) {
        std::cout << "Error! Incorrect input data for array";
        return -1;
    }

    TYPE_MAS* tmp_mas = new TYPE_MAS[size];

    time_sort = omp_get_wtime();
    if (size < 2 * num_thr) {
        Quick_sort(mas, 0, size - 1);
        std::cout << "Size array: " << size << "; Numb Threads: 1" << std::endl;
    } else {
        Task_sort(mas, size, tmp_mas, num_thr);
    }
    time_sort = omp_get_wtime() - time_sort;

    std::cout << "Spend time algorithm (omp parallel version Hoar sorting): " << time_sort << "sec" << std::endl;
    if (Check_not_decreasing(mas, size))
        std::cout << "Array is sorted by not decreasing" << std::endl;
    else
        std::cout << "Array isn't sorted by not decreasing" << std::endl;
    Show_mas(mas, size);

    delete[] mas;
    delete[] tmp_mas;
    return 0;
}
