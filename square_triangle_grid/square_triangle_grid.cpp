//
// Created by aksol on 24.06.2017.
//

#include "square_triangle_grid.h"


double_t SquareTriangleGrid::LineFunction_ab(double_t x, double_t step)
{
    double_t y = RATIO_Y_TO_X * (x - step);
    return y;
}

double_t SquareTriangleGrid::LineFunction_bc()
{
    return TRIANGLE_BASE;
}

uint_fast32_t SquareTriangleGrid::GetGreed(double_t step_x, points *cordinates, bool debug)
{

    // Координаты узлов сетки
    double_t x = first_point_x_coordinate.x;
    double_t y = 0;

    /*
     * t1 - разница между координатой y правой точки маленького треугольника и правой границей
     * основного треугольника. Если разница меньше eps, считаем что мы дошли до границы.
     * t2 - то же самое, для перевернутых маленьких треугольников.
     * Если выполнены условия для t1 и t2, значит этот уровень триангулирован и надо идти вверх по оси y
     * на следующий уровень.
     */
    double_t eps = EPS_T;
    double_t t1 = 0;
    double_t t2 = 0;
    points first_point, second_point, third_point;


    /*
     * Счетчики
     * i - счетчик прямых, параллельных левой грани треугольника (идет через два шага, одинаково на всех уровнях)
     * обнуляется в конце уровня
     * j - счетчик треугольников в массиве треугольников
     * ii - счетчик сдвига по x при переходе на ряд выше ( счетчик иксов на уровне)
     * l - счетчик рядов по высоте (он идет до середины треугольника)
     * k - нумератор узлов
     * n - число узлов на уровне (на сколько узлов делится горизонтальная сторона треугольника)
     * length - длина уровня
     */
    uint_fast32_t i = 0;
    uint_fast32_t j = 0;
    uint_fast32_t ii = 0;
    uint_fast32_t l = 0;
    uint_fast32_t k = 0;
    uint_fast32_t n = 0;
    double_t length = LineFunction_bc() + STEP_X;
    uint_fast32_t  coord_num = 0;

    /* Переход с ряда на ряд (вверх), пока не дойдем до середины треугольника */
    while (y < TRIANGLE_HEIGHT) {
        length -= STEP_X;
        n += (length / step_x) + 1;
        /* Цикл по горизонтальному ряду */
        for (;;)
        {
            x = ii * step_x;
            y = LineFunction_ab(x, i * step_x);
            first_point = { x, y, k };

            cordinates[coord_num].x = x;
            cordinates[coord_num].y = y;
            coord_num ++;

            x = (ii + 1) * step_x;
            y = LineFunction_ab(x, (i + 1) * step_x);
            second_point = { x, y, k + 1 };

            t1 = abs(x - LineFunction_bc());

            x = (ii + 1) * step_x;
            y = LineFunction_ab(x, i * step_x);
            third_point = { x, y, k + ((length / step_x) + 1) };

            triangles_array.push_back(new fem_triangle_element(first_point, second_point, third_point));

            if (debug == true)
            {
                GetGreedDebugMessages(triangles_array.at(j), j);
            }

            t2 = abs(x - LineFunction_bc());

            // Условие выхода на следующий уровень
            if ((t2 < EPS_T) && (t1 < EPS_T)) {
                ii = l + 1;
                i = 0;
                j++;
                k = n;
                cordinates[coord_num].x = second_point.x;
                cordinates[coord_num].y = second_point.y;
                coord_num ++;
                if (debug == true)
                {
                    cout << "BREAK CONDITION." << endl;
                    cout << "n = " << n << endl;
                    cout << "k = " << k << endl;
                }
                break;
            }

            first_point.x = triangles_array.back()->third_point.x;
            first_point.y = triangles_array.back()->third_point.y;
            first_point.point_num = triangles_array.back()->third_point.point_num;
            second_point.x = triangles_array.back()->second_point.x;
            second_point.y = triangles_array.back()->second_point.y;
            second_point.point_num = triangles_array.back()->second_point.point_num;

            x = (ii + 2) * step_x;
            y = LineFunction_ab(x, (i + 1) * step_x);

            third_point = { x, y, k + 1 + ((length / step_x) + 1) };
            triangles_array.push_back(new fem_triangle_element(first_point, second_point, third_point));

            if (debug == true)
            {
                GetGreedDebugMessages(triangles_array.at(j+1), j + 1);
            }

            i = i + 1;
            ii = ii + 1;
            j = j + 2;
            k++;
        }
        l++;
        if (debug == true)
        {
            cout << "New level l = " << l << endl;
        }
    }
    cordinates[coord_num].x = cordinates[coord_num - 1].x;
    cordinates[coord_num].y =  LineFunction_ab(cordinates[coord_num].x, i * step_x);

    return triangles_array.back()->third_point.point_num;
}

void SquareTriangleGrid::GetGreedDebugMessages(fem_triangle_element *triangle, uint_fast32_t triangle_num)
{
    cout << "First point of triange " << triangle_num << ": " << endl;
    cout << "x: " << triangle->first_point.x << endl;
    cout << "y: " << triangle->first_point.y << endl;
    cout << "k: " << triangle->first_point.point_num << endl;
    cout << "Second point of triange " << triangle_num << ": " << endl;
    cout << "x: " << triangle->second_point.x << endl;
    cout << "y: " << triangle->second_point.y << endl;
    cout << "k: " << triangle->second_point.point_num << endl;
    cout << "Third point of triange " << triangle_num << ": " << endl;
    cout << "x: " << triangle->third_point.x << endl;
    cout << "y: " << triangle->third_point.y << endl;
    cout << "k: " << triangle->third_point.point_num << endl;
}
