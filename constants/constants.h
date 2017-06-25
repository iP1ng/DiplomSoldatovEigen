//
// Created by aksol on 24.06.2017.
//

#ifndef DIPLOMSOLDATOVEIGEN_CONSTANTS_H
#define DIPLOMSOLDATOVEIGEN_CONSTANTS_H

// TODO RATIO_Y_TO_X расчитывать автоматически, возможно убрать вообще

/**
 * Сетка
 */
// Шаг по пространству
const double_t STEP_X = 0.01;
// Шаг по времени
const double_t TAU = 0.001;


/**
 * Треугольник
 *
 *         B
 *        /|
 *       / | TRIANGLE_HEIGHT
 *      /__|
 *     A    C
 *  TRIANGLE_BASE
 */
// Координаты треугольника
const double_t A_x = 0;
// Высота треугольника
const double_t TRIANGLE_HEIGHT = 2;
// Основание треугольника
const double_t TRIANGLE_BASE = 1;
// Коэффициент вытянутости треугольника по высоте (в работе 2:1)
// Задает шаг по y и форму прямоугольних треугольников (и углы)
const double_t RATIO_Y_TO_X = 2.0;
// Точность пересечения границы прямоугольного треугольника
// (используется для контроля выхода за границу области при триангуляции)
const double_t EPS_T = STEP_X * 0.001;
// Число сторон треугольника
const uint_fast32_t DIMENSION = 3;
// Число точек в прямоугольном трегольнике (натуральный ряд)
const uint_fast32_t DOTS_NUMBER = ((TRIANGLE_BASE / STEP_X) + 1) * ((TRIANGLE_BASE / STEP_X) + 2) / 2;
// Число треугольников после триангуляции
const uint_fast32_t TRIANGLES_NUMBER = (TRIANGLE_HEIGHT / (RATIO_Y_TO_X * STEP_X)) * (TRIANGLE_HEIGHT / (RATIO_Y_TO_X * STEP_X));


/**
 * Уравнение теплопроводности
 */
// Коэффициент теплопроводности (Вт/(м*К))
// (в работе это титан)
//  Взят с http://ispu.ru/files/u2/SP._bez_nomera_-_Spravochn._materialy_dlya_resheniya_zadach_po_kursu_Teplomassoobmen..pdf
const double_t Thermal_Conductivity = 15.1;
// Удельная массовая теплоемкость (Дж/кг*К)
const double_t Coeff_Mass_c = 532.0;
// Плотность титана (кг/м^3)
const double_t Coeff_rho = 4540.0;


/**
 * Начальные и граничные условия
 */
//коэффициент конвективного теплообмена на ребре (Вт/м^2)
const double_t Coeff_H = 0.0;
// Коэффициент сопротивления формы конуса 2:1 (острием к потоку)
// (безразмерная величина, определяющая реакцию среды на движение в ней тела)
// Взят с https://ru.wikipedia.org/wiki/%D0%9A%D0%BE%D1%8D%D1%84%D1%84%D0%B8%D1%86%D0%B8%D0%B5%D0%BD%D1%82_%D1%81%D0%BE%D0%BF%D1%80%D0%BE%D1%82%D0%B8%D0%B2%D0%BB%D0%B5%D0%BD%D0%B8%D1%8F_%D1%84%D0%BE%D1%80%D0%BC%D1%8B
const double_t Coeff_S = 0.5;
// Коэффициент плотности воздуха (м^-1)
const double_t Coeff_K = 0.00013;
// Начальная плотность воздуха (кг/м^3)
const double_t RHO_0 = 1.2;
// Начальное распределение температур (К)
const double_t INITIAL_TEMPERATURE = 273.0;
// Температура среды (К)
const double_t OUTER_TEMPERATURE = 273.0;
// Ускорение ракеты (м/с^2)
const double_t ACCELERATION = 90;
// Температура плавления (К)
// (в задаче температура плавления титана)
// Взято с https://ru.wikipedia.org/wiki/%D0%A2%D0%B8%D1%82%D0%B0%D0%BD_(%D1%8D%D0%BB%D0%B5%D0%BC%D0%B5%D0%BD%D1%82)
// Ракета летит до того момента, пока средняя температура на наконечнике ниже температуры плавления (по постановке)
const double_t MELTING_TEMPERATURE = 1941;


/**
 * Константы
 */
// Число Пи
const double_t PI = 3.1415926535897932384626433832795;

#endif //DIPLOMSOLDATOVEIGEN_CONSTANTS_H
