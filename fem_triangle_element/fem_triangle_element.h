//
// Created by aksol on 24.06.2017.
//

#ifndef DIPLOMSOLDATOVEIGEN_FEM_TRIANGLE_ELEMENT_H
#define DIPLOMSOLDATOVEIGEN_FEM_TRIANGLE_ELEMENT_H


// Для sqrt, pow, fabs
#include <math.h>

// Подключаем int фиксированного размера
#include <cstdint>

// Подключаем узлы points
#include "../structures/structures.h"

// Подключаем константы
// Coeff_Mass_c
// Coeff_rho
// Thermal_Conductivity
// EPS_T
// RATIO_Y_TO_X
// PI
#include "../constants/constants.h"


/**
 * Треугольный конечный элемент из трех узлов
 * Методы вычисления коэффициентов функций формы
 * Методы вычисления элементных матриц C, K
 * Метод вычисления элементной правой части F
 *
 * Граничные условия задаются в этом классе
 */
//TODO убрать зависимости от констант, убрать зависимости от основной области
class fem_triangle_element
{
public:
    // узлы элемента
    points first_point;
    points second_point;
    points third_point;

    // число узлов элемента
    static const uint8_t DIMENSION = 3;

    // Конструктор
    fem_triangle_element() = default;
    fem_triangle_element(points _first_point, points _second_point, points _third_point) :
            first_point(_first_point), second_point(_second_point), third_point(_third_point) {}

    // Вычисление коэффициентов функций формы
    void coef_a(double_t a[DIMENSION]);
    void coef_b(double_t b[DIMENSION]);
    void coef_c(double_t c[DIMENSION]);

    // Вычисление длин сторон
    double_t GetFirstSideLength();
    double_t GetSecondSideLength();
    double_t GetThirdSideLength();

    // Вычисление площади треугольника по формуле Герона
    double_t GetSquareTriangleArea();

    // Вычисление определителя матрицы A (для функций формы)
    double_t GetMatrixADeterminant();

    // Вычисление скалярной величины R
    // Сегерлинд, формула (10.23), стр. 189
    double_t GetR();

    // Вспомогательная матрица, используемая при граничных условиях 1 и 3 рода
    // Сегерлинд, формулы (10.31), (10.32), (10.33) стр. 192 -193
    // side задает сторону элемента, на которую действует граничное условие
    // f - ij , s - jk, t - ki
    // TODO сейчас реализовано только для стороны f
    void GetRForDirichle (double_t Rd[DIMENSION][DIMENSION], char side);

    // Матрица K
    // Сегерлин, формула (10.22), стр. 189
    // TODO реализовано только граничное условие второго рода
    void Matrix_K(double_t K[DIMENSION][DIMENSION], char element_side);

    // Матрица С
    // Сегерлинд, формула (11.17), стр. 205
    void Matrix_C(double_t C[DIMENSION][DIMENSION]);

    // Столбец F
    // Сегерлинд, формула (10.34), стр. 193
    // side задает сторону элемента, на которую действует граничное условие
    // f - ij , s - jk, t - ki
    // TODO сейчас реализовано только для стороны t
    // TODO сейчас реализовано только граничное условие второго рода
    // TODO сейчас захардкожена связь с основной областью (большим треугольником)
    void Column_F(double_t F[DIMENSION], double_t q, char element_side);

    // Вычисление коэффициентов функций формы
    // TODO набросок метода, нужна нормальная реализация
    void GetN(double_t * N, double_t * Phi);

};


#endif //DIPLOMSOLDATOVEIGEN_FEM_TRIANGLE_ELEMENT_H
