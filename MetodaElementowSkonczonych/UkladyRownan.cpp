#include "UkladyRownan.h"



double* ElimGauss(double** A, double* b, int n)
{
    double** AB = new double* [n];

    for (int i = 0; i < n; i++)
    {
        AB[i] = new double[n + 1];
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            AB[i][j] = A[i][j];
        }
        AB[i][n] = b[i];
    }

    for (int k = 0; k < n; k++) {

        int rmax = k;
        double vmax = std::abs(AB[k][k]);

        // SprawdŸ, czy element diagonalny nie jest bliski zeru
        if (std::abs(vmax) < 1e-10) {
            // Wymieñ wiersze, jeœli element diagonalny jest bliski zeru
            for (int i = k + 1; i < n; i++) {
                if (std::abs(AB[i][k]) > vmax) {
                    rmax = i;
                    vmax = std::abs(AB[i][k]);
                }
            }

            if (std::abs(AB[rmax][k]) < 1e-10) {
                // Element diagonalny jest nadal bliski zeru, program nie mo¿e kontynuowaæ
                std::cerr << "Eliminacja Gaussa: Macierz osobliwa." << std::endl;
                return nullptr;
            }

            // Wymieñ wiersze
            for (int j = k; j < n + 1; j++) {
                std::swap(AB[k][j], AB[rmax][j]);
            }
        }

        for (int i = k + 1; i < n; i++) {
            double factor = AB[i][k] / AB[k][k];
            for (int j = k + 1; j < n + 1; j++) {
                AB[i][j] -= factor * AB[k][j];
            }
            AB[i][k] = 0;
        }
    }

    double* x = new double[n];
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += AB[i][j] * x[j];
        }

        // SprawdŸ, czy element diagonalny nie jest bliski zeru
        if (std::abs(AB[i][i]) < 1e-10) {
            // Element diagonalny jest bliski zeru, program nie mo¿e kontynuowaæ
            std::cerr << "Eliminacja Gaussa: Macierz osobliwa." << std::endl;
            delete[] x;
            for (int i = 0; i < n; i++) {
                delete[] AB[i];
            }
            delete[] AB;
            return nullptr;
        }

        x[i] = (AB[i][n] - sum) / AB[i][i];
    }

    for (int i = 0; i < n; i++) {
        delete[] AB[i];
    }
    delete[] AB;

    return x;
}

void RozkladLU(double** A, int n)
{
    double tmp = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i <= j)
            {
                tmp = 0;
                for (int k = 0; k < i; k++)
                    tmp += A[i][k] * A[k][j];

                A[i][j] = A[i][j] - tmp;
            }
            else
            {
                tmp = 0;
                for (int k = 0; k < j; k++)
                    tmp += A[i][k] * A[k][j];
                A[i][j] = (1 / A[j][j]) * (A[i][j] - tmp);
            }
        }
    }

}

double* URRozkladLU(double** A, double* b, int n)
{
    RozkladLU(A, n);

    double* x = new double[n];

    double tmp = 0;
    //Obliczanie y

    x[0] = b[0];
    for (int i = 1; i < n; i++)
    {
        tmp = b[i];
        for (int j = 0; j < i; j++)
        {
            tmp -= A[i][j] * x[j];
        }
        x[i] = tmp;
    }

    //Obliczanie x

    x[n - 1] /= A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        tmp = x[i];
        for (int j = n - 1; j > i; j--)
        {
            tmp -= A[i][j] * x[j];
        }
        x[i] = tmp / A[i][i];
    }

    return x;
}