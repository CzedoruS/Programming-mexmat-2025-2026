#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "init.h"
#include "solver.h"
#include "print.h"

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " n m k [filename]\n";
        return 1;
    }

    int n = std::atoi(argv[1]);
    int m = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);

    double* A = new double[n * n];
    double* b = new double[n];
    double* x = new double[n];

    bool ok = false;

    if (k == 0) {
        if (argc != 5) {
            std::cerr << "Filename required when k = 0\n";
            delete[] A; delete[] b; delete[] x;
            return 1;
        }
        const char* filename = argv[4];
        ok = read_matrix_from_file(n, A, filename);
        if (!ok) {
            std::cerr << "Cannot open or read file: " << filename << "\n";
            delete[] A; delete[] b; delete[] x;
            return 1;
        }
    } else if (k >= 1 && k <= 4) {
        ok = init_matrix_formula(k, n, A);
        if (!ok) {
            std::cerr << "Error initializing matrix with formula " << k << "\n";
            delete[] A; delete[] b; delete[] x;
            return 1;
        }
    } else {
        std::cerr << "Invalid formula k: " << k << "\n";
        delete[] A; delete[] b; delete[] x;
        return 1;
    }

    compute_b(n, A, b);

    // СОХРАНЯЕМ исходный вектор b перед решением для корректной невязки
    double* b_original = new double[n];
    for (int i = 0; i < n; ++i) {
        b_original[i] = b[i];
    }

    std::cout << "Matrix A:\n";
    print_matrix(n, n, A, m);
    std::cout << "Vector b:\n";
    print_vector(b, n, m);

    clock_t start = clock();
    ok = jordan_solve(n, A, b, x);  // A и b портятся - невязка вычисляется некорректно
    clock_t end = clock();

    if (!ok) {
        std::cerr << "Solver failed\n";
    } else {
        std::cout << "Solution x:\n";
        print_vector(x, n, m);

        // ПЕРЕИНИЦИАЛИЗИРУЕМ матрицу A для вычисления невязки
        if (k == 0) {
            const char* filename = argv[4];
            read_matrix_from_file(n, A, filename);
        } else {
            init_matrix_formula(k, n, A);
        }
        // Вектор b_original уже сохранен

        // Вычисление нормы погрешности
        double error_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            double exact_x;
            if (i % 2 == 0) {
                exact_x = 1.0;
            } else {
                exact_x = 0.0;
            }
            double difference = x[i] - exact_x;
            error_norm += difference * difference;
        }
        error_norm = sqrt(error_norm);
        std::cout << "Error norm ||x_exact - x|| = " << error_norm << "\n";        

        // Вычисляем невязку с восстановленной A и сохраненным b_original
        double res_norm = compute_residual_norm(n, A, x, b_original);
        std::cout << "Residual norm ||Ax-b||/||b|| = " << res_norm << "\n";
    }

    std::cout << "Time elapsed: " << double(end - start) / CLOCKS_PER_SEC << " sec\n";

    delete[] A;
    delete[] b;
    delete[] x;
    delete[] b_original;  
    return 0;

    //Память: n^2 (из матрицы A) + n (полученный вектор b) + n (вектор решения x) + n (сохраненный исходный вектор b)
    // = n^2 + 3n = n^2 + O(n) - удовлетворяет требованию
    //Время: на одну итерацию k: O(n) (нормировка) + O(n^2) (исключение строк) = O(n^2)
    //итого n x O(n^2) = O(n^3) - удовлетворяет требованиям
}