//
// Created by aksol on 24.06.2017.
//

#ifndef DIPLOMSOLDATOVEIGEN_STRUCTURES_H
#define DIPLOMSOLDATOVEIGEN_STRUCTURES_H

/**
 * Структура, описывающая координаты и номер узла, а также радиус-вектор узла.
 */
struct points {
    double_t x;
    double_t y;
    uint_fast32_t point_num;
    double_t rad_vector() { return x; }
};

#endif //DIPLOMSOLDATOVEIGEN_STRUCTURES_H
