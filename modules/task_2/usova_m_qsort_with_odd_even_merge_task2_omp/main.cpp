// Copyright 2019 Usova Marina

#include <omp.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <vector>
#include <cmath>

// FUNCTION OF SORTING HOARA

void quicksort(int* mass, int iStart, int iFinish) {
  if (iFinish > iStart) {
    int s = iStart, f = iFinish;
    // as the reference element take the middle of the array
    int middle = mass[(iFinish + iStart) / 2];

    // we carry out division into elements smaller than the reference
    // and larger ones of the reference, and then recursively run
    // the same function for the generated sets
    do {
      while (mass[s] < middle)
        s++;
      while (mass[f] > middle)
        f--;
      if (s <= f) {
        int tmp = mass[s];
        mass[s] = mass[f];
        mass[f] = tmp;
        s++;
        f--;
      }
    } while (s < f);

    if (iStart < f)
      quicksort(mass, iStart, f);
    if (s < iFinish)
      quicksort(mass, s, iFinish);
  }
}

// FUNCTION OF ODD-EVEN MERGE

// function that performs splitting into even and odd elements
void splitUnsplit(int* mass, int* leftMass, int* rightMass, int left, int right, int flag) {
  for (int i = flag; i < left; i += 2) {
    leftMass[i] = mass[i];
  }

  int l_ = flag;
  int r_ = flag;
  int i = flag;

  while ((l_ < left) && (r_ < right)) {
    if (leftMass[l_] <= rightMass[r_]) {
      mass[i] = leftMass[l_];
      l_ += 2;
    } else {
      mass[i] = rightMass[r_];
      r_ += 2;
    }
    i += 2;
  }

  if (l_ == left) {
    for (int j = r_; j < right; j += 2, i += 2) {
      mass[i] = rightMass[j];
    }
  } else {
    for (int j = l_; j < left; j += 2, i += 2) {
      mass[i] = leftMass[j];
    }
  }
}

void compare_exchange(int* x1, int* x2) {
  if (*x1 < *x2) {
    int tmp = *x1;
    *x1 = *x2;
    *x2 = tmp;
  }
}

// function merges ordered arrays
void oddEvenMerge(int* mass, int size, int *array_of_sizes, int threads) {
  int id;
  int sortLevel = threads;
  int* tmp1 = new int[size];
  int* tmp2 = new int[threads + 1];

  while ((sortLevel != 1)) {
#pragma omp parallel for private(id)
    for (id = 0; id < sortLevel; id++) {
      int t = id % 2;
      splitUnsplit(mass + array_of_sizes[id - t],
        tmp1 + array_of_sizes[id],
        mass + array_of_sizes[id + (1 - t)],
        array_of_sizes[id + (1 - t)] - array_of_sizes[id - t],
        array_of_sizes[id + 1 + (1 - t)] - array_of_sizes[id + (1 - t)], t);
    }
#pragma omp parallel for private(id)
    for (id = 0; id < sortLevel; id++) {
      if (id % 2 != 0) {
        for (int j = 1; j < (array_of_sizes[id + 1] - array_of_sizes[id - 1] + 1) / 2; j++) {
          compare_exchange(&(mass + array_of_sizes[id - 1])[2 * j],
            &(mass + array_of_sizes[id - 1])[2 * j - 1]);
        }
      }
    }
    sortLevel = sortLevel / 2;

    tmp2 = new int[sortLevel];
    tmp2[0] = array_of_sizes[0];
    for (int i = 1; i < sortLevel; i++) {
      tmp2[i] = array_of_sizes[i * 2];
    }
    array_of_sizes = new int[sortLevel];
    for (int i = 0; i < sortLevel; i++) {
      array_of_sizes[i] = tmp2[i];
    }
    array_of_sizes[sortLevel] = size;
  }

  delete[] tmp1;
  delete[] tmp2;
}

void parallelQuicksort(int* mass, int threads, int size) {
  int id;
  int *array_of_sizes = new int[threads + 1];

  for (int i = 0; i < threads; i++) {
    array_of_sizes[i] = i * size / threads + (i * size / threads) % 2;
  }
  array_of_sizes[threads] = size;

  omp_set_num_threads(threads);

#pragma omp parallel for private(id)
  for (id = 0; id < threads; id++) {
    quicksort(mass, array_of_sizes[id], array_of_sizes[id + 1] - 1);
  }
  if (threads > 1) {
    oddEvenMerge(mass, size, array_of_sizes, threads);
  }

  delete[] array_of_sizes;
}

// SUPPORTING FUNCTIONS

void create_array(int* array, int* copy, int size) {
  srand((unsigned int)time(NULL));
  const int max_elem = 1000;
  for (int i = 0; i < size; ++i) {
    array[i] = std::rand() % max_elem;
    copy[i] = array[i];
  }
}

void print_array(int* array, int size) {
  for (int i = 0; i < size; ++i)
    std::cout << array[i] << " ";
  std::cout << std::endl;
}

bool check(int* A, int* B, int size) {
  for (int i = 0; i < size; ++i)
    if (std::abs(A[i] - B[i]) != 0) {
      return false;
    }
  return true;
}

/*int compare(const int *i, const int *j) {
  return *i - *j;
}*/

// MAIN FUNCTION

int main(int argc, char** argv) {
  int size, threads = 4;

  srand((unsigned int)time(NULL));

  if (argc > 1) {
    size = atoi(argv[1]);
  } else {
    // size = 500 + std::rand() % 1000;
    // size = 12;
    size = 1000000;
  }

  std::cout << "Array size: n = " << size << std::endl;

  int* data = new int[size];
  int* copy = new int[size];

  // int data[12] = { 3, 7, 10, 1, 10, 17, 4, 8, 5, 9, 11, 10 };
  // int copy[12] = { 3, 7, 10, 1, 10, 17, 4, 8, 5, 9, 11, 10 };

  std::cout << "Start massive generating..." << std::endl;

  create_array(data, copy, size);

  std::cout << "massive is generated!" << std::endl << std::endl;

  // print_array(data, size);

  std::cout << "Start sequential sort..." << std::endl;

  auto timeWork_ = omp_get_wtime();
  quicksort(copy, 0, size - 1);
  timeWork_ = omp_get_wtime() - timeWork_;

  // print_array(copy, size);

  std::cout << "sorting done!" << std::endl << std::endl;
  std::cout << "Start parallel sort..." << std::endl;

  auto timeWork = omp_get_wtime();
  parallelQuicksort(data, threads, size);
  timeWork = omp_get_wtime() - timeWork;

  // print_array(data, size);

  std::cout << "sorting done!" << std::endl << std::endl;

  std::cout.precision(16);
  if (check(copy, data, size))
    std::cout << "- good sorting: results of sequential sort and parallel sort coincide"

    << std::endl << "- time of sequential sorting: " << timeWork_
    << std::endl << "- time of parallel sorting: " << timeWork
    << std::endl << "- acceleration of parallel sorting: " << timeWork_ / timeWork << std::endl;
  else
    std::cout << " error of sorting!!!" << std::endl;

  delete[]data;
  delete[]copy;

  return 0;
}
