//
// Created by aksol on 25.06.2017.
//

#ifndef DIPLOMSOLDATOVEIGEN_AIR_PARAMETERS_H
#define DIPLOMSOLDATOVEIGEN_AIR_PARAMETERS_H

// Плотность воздуха на высоте at^2/2 над Землей
double_t func_calculate_rho(double_t t) {
    return (RHO_0 * exp(- Coeff_K * (ACCELERATION * t * t / 2)));
}

// Тепловой поток на ребро наконечника ракеты на времени полета t
double_t func_calculate_q(double_t t) {
    double_t rho = func_calculate_rho(t);
    double_t V = ACCELERATION * t;

    return (( Coeff_S * rho * V * V * V) /
            (100 * Thermal_Conductivity * PI * TRIANGLE_BASE * sqrt(TRIANGLE_BASE * TRIANGLE_BASE + TRIANGLE_HEIGHT * TRIANGLE_HEIGHT)));
}

#endif //DIPLOMSOLDATOVEIGEN_AIR_PARAMETERS_H
