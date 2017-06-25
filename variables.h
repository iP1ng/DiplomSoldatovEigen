//
// Created by aksol on 25.06.2017.
//

#ifndef DIPLOMSOLDATOVEIGEN_VARIABLES_H
#define DIPLOMSOLDATOVEIGEN_VARIABLES_H

/**
 * Тепловой поток, заданный на границе
 */
double_t heat_flow = 0;


/**
 * Вектор F правых частей для элемента
 */
double_t F_elem[DIMENSION];
/**
 * Вектор F правой части для всех элементов
 */
double_t *F = new double[DOTS_NUMBER];
/**
 * Вектор F правой части для всех элементов на предыдущем временном слое
 */
double_t *F_old = new double_t[DOTS_NUMBER];
/**
 * 2/tau [C_elem] - [K_elem]
 */
double_t **R = new double_t *[DOTS_NUMBER];


/**
 * Итоговое распределение температуры во всем треугольнике
 */
//double_t *Temperature = new double_t[DOTS_NUMBER];


/**
 * Матрица коэффициентов элемента K_elem
 */
double_t K_elem[DIMENSION][DIMENSION];
/**
 * Матрица коэффициентов элемента C_elem
 */
double_t C_elem[DIMENSION][DIMENSION];


/**
 * Массив номеров строк в итоговом векторе правых частей Result,
 * к которому нужно будет прибавить полученные правые части F_elem на элементе
 */
uint_fast32_t index[DIMENSION];
/**
 * Средняя температура на наконечнике 
 */
double_t average_temp = 0;
/* 
 * Время, которое уже успела пролететь ракета 
 */
double_t flight_time = - TAU;


void InitDimensionN()
{
    for (auto i = 0; i < DOTS_NUMBER; i++) {
        F[i] = 0;
        F_old[i] = F[i];
        R[i] = new double_t[DOTS_NUMBER];
        for (auto j = 0; j < DOTS_NUMBER; j++) {
            R[i][j] = 0;
        }
    }
}

void InitDimension3()
{
    for (auto i = 0; i < DIMENSION; i++) {
        F_elem[i] = 0;
        for (auto j = 0; j < DIMENSION; j++) {
            K_elem[i][j] = 0;
            C_elem[i][j] = 0;
        }
    }
}

#endif //DIPLOMSOLDATOVEIGEN_VARIABLES_H
