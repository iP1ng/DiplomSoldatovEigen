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
ofstream triangle_temp("triangle_temp.dat", std::ofstream::out | std::fstream::app);

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

    uint_fast32_t dots_num = floor(TRIANGLE_BASE / STEP_X) + 1;
    uint_fast32_t dots_num_next = dots_num;
    uint_fast32_t dots_shift = 0;

    triangle_temp << "Number of dots = " << DOTS_NUMBER << endl;
    triangle_temp << "Number of triangles = " << TRIANGLES_NUMBER << endl;
    triangle_temp << "Space step h = " << STEP_X << endl;
    triangle_temp << "Time step tau = " << TAU << endl;
    triangle_temp << "Thermal_Conductivity = " << Thermal_Conductivity << endl;

    dots_num = floor(TRIANGLE_BASE / STEP_X) + 1;
    dots_num_next = dots_num;
    dots_shift = 0;
    triangle_temp << endl;
    triangle_temp << "Flight time = " << flight_time << endl;
    triangle_temp << "Heat flow = " << heat_flow << endl;
    triangle_temp << "Average temperature = " << average_temp << endl;
    for (auto strs = 0; strs < dots_num; strs++)
        triangle_temp << "**********";
    triangle_temp << endl;
    for (auto d = 0; d < DOTS_NUMBER; d++) {
        if (d == dots_num) {
            dots_num_next = dots_num_next - 1;
            dots_num = dots_num + dots_num_next;
            dots_shift ++;
            triangle_temp << endl;
            for (auto st = 0; st < dots_shift; st++)
                triangle_temp << setw(6) << "       ";
        }
        triangle_temp << setw(6) << setprecision(4)<< Temperature[d] << " ";
    }


}
#endif //DIPLOMSOLDATOVEIGEN_DEBUG_H
