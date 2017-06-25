//
// Created by aksol on 25.06.2017.
//

#ifndef DIPLOMSOLDATOVEIGEN_DEBUG_H
#define DIPLOMSOLDATOVEIGEN_DEBUG_H


#include <iomanip>
#include <fstream>

#include "constants/constants.h"

typedef Eigen::SparseMatrix<double> SpMat;


ofstream result_matrix("result_matrix.dat", std::ofstream::out);
ofstream result_m("result_m.dat", std::ofstream::out);
ofstream triangle_temp("triangle_temp.dat", std::ofstream::out);

void ShowRunParameters()
{
    cout << "Number of dots = " << DOTS_NUMBER << endl;
    cout << "Number of triangles = " << TRIANGLES_NUMBER << endl;
    cout << "Space step h = " << STEP_X << endl;
    cout << "Time step tau = " << TAU << endl;
    cout << endl;
}


void ShowMatrixN(double_t matrix[DOTS_NUMBER][DOTS_NUMBER])
{
    result_matrix << endl;
    result_matrix << "---------------------------------------------------" << endl;
    result_matrix << endl;
    for (auto i = 0; i < DOTS_NUMBER; i++) {
        for (auto j = 0; j < DOTS_NUMBER; j++) {
            result_matrix << setw(6) << setprecision(4) << matrix[i][j] << " ";
        }
        result_matrix << endl;
    }
}


void ShowTemperature()
{
    /**
     * Debug вывод распределения температуры на треугольнике во времени
     * для проверки полученного решения
     */
    int_fast32_t dots_shift = DOTS_NUMBER - 3;
    uint_fast32_t level = 3;

    triangle_temp << "Number of dots = " << DOTS_NUMBER << endl;
    triangle_temp << "Number of triangles = " << TRIANGLES_NUMBER << endl;
    triangle_temp << "Space step h = " << STEP_X << endl;
    triangle_temp << "Time step tau = " << TAU << endl;
    triangle_temp << "Thermal_Conductivity = " << Thermal_Conductivity << endl;
    triangle_temp << "flight time = " << flight_time << endl;
    triangle_temp << "heat flow = " << heat_flow << endl;

    triangle_temp << setw(10) << setprecision(6)<< Temperature[DOTS_NUMBER - 1] << endl;

    for (int_fast32_t i = DOTS_NUMBER - 2; level < (TRIANGLE_BASE / STEP_X + 3); i-=level) {
        for (int j = 0; j < level - 1; j++) {
            triangle_temp << setw(10) << setprecision(6)<< Temperature[dots_shift + j] << " ";
        }

        dots_shift -= level;
        level += 1;
        triangle_temp << endl;
    }
    triangle_temp << endl;
}
#endif //DIPLOMSOLDATOVEIGEN_DEBUG_H
