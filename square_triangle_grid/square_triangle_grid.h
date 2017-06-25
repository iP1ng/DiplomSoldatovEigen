//
// Created by aksol on 24.06.2017.
//

#ifndef DIPLOMSOLDATOVEIGEN_SQUARE_TRIANGLE_GRID_H
#define DIPLOMSOLDATOVEIGEN_SQUARE_TRIANGLE_GRID_H


// Подключаем быстрые int фиксированного размера
#include <cstdint>

// Подключаем векторы
#include <vector>
using std::vector;

// Подключаем double_t тип для оптимальной работы на любой архитектуре компьютера.
// также подключаем fabs
// также подключаем exp
#include <math.h>

// Подключаем константы
// STEP_X отвечает за число узлов в сетке
#include "../constants/constants.h"

// Подключаем узлы points
#include "../structures/structures.h"

// Подключаем конечные элементы, из которых будем строить сетку на заданной области
#include "../fem_triangle_element/fem_triangle_element.h"


// Debug output
#include <iostream>
using namespace std;


class SquareTriangleGrid {

private:
     // first_point_x_coordinate - координата x точки треугольника, находящейся в нуле оси координат
    points first_point_x_coordinate;

public:
    vector<fem_triangle_element *> triangles_array;

    /*
     * Конструктор класса IsoscelesTriangleGrid
     */
    SquareTriangleGrid(points a)
    {
        first_point_x_coordinate = a;
    }

    /*
     * Метод построения сетки
     * step_x - шаг по оси x
     * щаг по оси y определяется уравнением, описывающим гипотенузу
     * debug = true для отладочных сообщений
     * На выходе возвращает число узлов
     */
    uint_fast32_t GetGreed(double_t step_x, points *cordinates, bool debug);

    void GetGreedDebugMessages(fem_triangle_element *triangle, uint_fast32_t triangle_num);
    /*
     * Функция, вычисляющая значение координаты y по заданной координате x на прямой AB треугольника
     * или на прямой, сдвинутой относительно AB на step
     * (структура треугольника описана в constants.h)
     */
    double_t LineFunction_ab(double_t x, double_t step);

    /*
     * Функция, вычисляющая значение координаты y на прямой BC треугольника
     * (структура треугольника описана в constants.h)
     */
    double_t LineFunction_bc();
};


#endif //DIPLOMSOLDATOVEIGEN_SQUARE_TRIANGLE_GRID_H
