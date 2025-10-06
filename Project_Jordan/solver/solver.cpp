#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include "solver.h"


bool jordan_solve(int n, double* A, double* b, double* x) {


    // Метод Жордана (работаем напрямую с переданными массивами)
    for (int k = 0; k < n; ++k) {
        if (std::fabs(A[k*n + k]) < std::numeric_limits<double>::epsilon()) {
            return false;
        }
        
        // Нормировка k-й строки
        double pivot = A[k*n + k];
        for (int j = k; j < n; ++j) {
            A[k*n + j] /= pivot;
        }
        b[k] /= pivot;
        
        // Исключение из всех строк кроме k-й
        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            double factor = A[i*n + k];
            for (int j = k; j < n; ++j) {
                A[i*n + j] -= factor * A[k*n + j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    // Копируем решение (b теперь содержит решение, потому что слева единичная матрица)
    for (int i = 0; i < n; ++i) {
        x[i] = b[i];
    }
    
    return true;
}