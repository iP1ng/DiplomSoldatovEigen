//
// Created by aksol on 24.06.2017.
//

#include "fem_triangle_element.h"

void fem_triangle_element::coef_a(double_t a[DIMENSION])
{
    a[0] = second_point.x * third_point.y - third_point.x * second_point.y;
    a[1] = third_point.x * first_point.y - first_point.x * third_point.y;
    a[2] = first_point.x * second_point.y - second_point.x * first_point.y;
}


void fem_triangle_element::coef_b(double_t b[DIMENSION])
{
    b[0] = second_point.y - third_point.y;
    b[1] = third_point.y - first_point.y;
    b[2] = first_point.y - second_point.y;
}


void fem_triangle_element::coef_c(double_t c[DIMENSION])
{
    c[0] = third_point.x - second_point.x;
    c[1] = first_point.x - third_point.x;
    c[2] = second_point.x - first_point.x;
}


double_t fem_triangle_element::GetFirstSideLength()
{
    return sqrt(pow(first_point.x - second_point.x, 2) + pow(first_point.y - second_point.y, 2));
};


double_t fem_triangle_element::GetSecondSideLength()
{
    return sqrt(pow(second_point.x - third_point.x, 2) + pow(second_point.y - third_point.y, 2));
};


double_t fem_triangle_element::GetThirdSideLength()
{
    return sqrt(pow(third_point.x - first_point.x, 2) + pow(third_point.y - first_point.y, 2));
};


double_t fem_triangle_element::GetSquareTriangleArea()
{
    double_t first_side_length = GetFirstSideLength();
    double_t second_side_length = GetSecondSideLength();
    double_t third_side_length = GetThirdSideLength();

    double_t p = (first_side_length + second_side_length + third_side_length) / 2;
    return sqrt(p * (p - first_side_length) * (p - second_side_length) * (p - third_side_length));
}


double_t fem_triangle_element::GetMatrixADeterminant()
{
    return 0.5 * (second_point.x * third_point.y
                  - third_point.x * second_point.y
                  - first_point.x * third_point.y
                  + first_point.x * second_point.y
                  + third_point.x * first_point.y
                  - second_point.x * first_point.y);
}


double_t fem_triangle_element::GetR() {
    points p[DIMENSION];

    p[0].x = first_point.x;
    p[0].y = first_point.y;
    p[1].x = second_point.x;
    p[1].y = second_point.y;
    p[2].x = third_point.x;
    p[2].y = third_point.y;

    double_t R = 1 / 12.0
                 * ((2 * p[0].rad_vector() + p[1].rad_vector() + p[2].rad_vector()) * p[0].rad_vector()
                    + (p[0].rad_vector() + 2 * p[1].rad_vector() + p[2].rad_vector()) * p[1].rad_vector()
                    + (p[0].rad_vector() + p[1].rad_vector() + 2 * p[2].rad_vector()) * p[2].rad_vector());
    return R;
}


void fem_triangle_element::GetRForDirichle (double_t Rd[DIMENSION][DIMENSION], char element_side) {
    points p[DIMENSION];

    p[0].x = first_point.x;
    p[0].y = first_point.y;
    p[1].x = second_point.x;
    p[1].y = second_point.y;

    double_t coeff = 0;

    switch(element_side) {
        case 'f':
            coeff = 2.0 * PI * Coeff_H * GetFirstSideLength() / 12.0;

            for (auto i = 0; i < DIMENSION; i++)
                for (auto j = 0; j < DIMENSION; j++)
                    Rd[i][j] = 0;


            Rd[0][0] = coeff * (3 * p[0].rad_vector() + p[1].rad_vector());
            Rd[0][1] = coeff * (p[0].rad_vector() + p[1].rad_vector());
            Rd[1][0] = coeff * (p[0].rad_vector() + p[1].rad_vector());
            Rd[1][1] = coeff * (p[0].rad_vector() + 3 * p[1].rad_vector());
            break;
    }

}


void fem_triangle_element::Matrix_K(double_t K[DIMENSION][DIMENSION], char element_side)
{
    double_t a[DIMENSION];
    double_t b[DIMENSION];
    double_t c[DIMENSION];
    coef_a(a);
    coef_b(b);
    coef_c(c);

    double_t Rd[DIMENSION][DIMENSION];

    GetRForDirichle(Rd, element_side);

    double_t R = GetR();

    for (int i = 0; i < DIMENSION; i++)
        for (int j = 0; j < DIMENSION; j++)
            K[i][j] = 0;

        for (int i = 0; i < DIMENSION; i++) {
            for (int j = 0; j < DIMENSION; j++) {
                K[i][j] = ( 2.0 * PI * R * Thermal_Conductivity / (4.0 * GetMatrixADeterminant()))
                          * (b[i] * b[j] + c[i] * c[j]);
            }
        }

    /**
     * Вычисляем K так, только если задан кожффициент теплообмена Coeff_H
     * т.е. задано граничное условие первого или третьего рода */
    /* TODO доделать для всех сторон и граничных условий */
    /*if ((fabs(first_point.y) < EPS_T) && (fabs(second_point.y) < EPS_T) && (Coeff_H > EPS_T)) {
        for (int i = 0; i < DIMENSION; i++) {
            for (int j = 0; j < DIMENSION; j++) {
                K[i][j] = (Thermal_Conductivity  * 2.0 * PI * R / (4.0 * GetMatrixADeterminant()))
                          * (b[i] * b[j] + c[i] * c[j])
                          + Rd[i][j];
            }
        }
    }*/
}


void fem_triangle_element::Matrix_C(double_t C[DIMENSION][DIMENSION])
{
    points p[DIMENSION];
    p[0].x = first_point.x;
    p[0].y = first_point.y;
    p[1].x = second_point.x;
    p[1].y = second_point.y;
    p[2].x = third_point.x;
    p[2].y = third_point.y;

    double_t D = 2.0 * PI * GetMatrixADeterminant() *  Coeff_Mass_c * Coeff_rho  /  180.0;
    double_t R[DIMENSION];

    for (int i = 0; i < DIMENSION; i++) {
        R[i] = p[i].rad_vector();
    }

    for (int i = 0; i < DIMENSION; i++) {
        for (int j = 0; j < DIMENSION; j++) {
            C[i][j] = 0;
        }
    }

    C[0][0] = D * (12 * R[0] * R[0]
                   + 2 * R[1] * R[1]
                   + 2 * R[2] * R[2]
                   + 6 * R[0] * R[1]
                   + 6 * R[0] * R[2]
                   + 2 * R[1] * R[2]);
    C[0][1] = D * (3 * R[0] * R[0]
                   + 3 * R[1] * R[1]
                   + R[2] * R[2]
                   + 4 * R[0] * R[1]
                   + 2 * R[0] * R[2]
                   + 2 * R[1] * R[2]);
    C[1][0] = C[0][1];
    C[0][2] = D * (3 * R[0] * R[0]
                   + R[1] * R[1]
                   + 3 * R[2] * R[2]
                   + 2 * R[0] * R[1]
                   + 4 * R[0] * R[2]
                   + 2 * R[1] * R[2]);
    C[2][0] = C[0][2];
    C[1][1] = D * (2 * R[0] * R[0]
                   + 12 * R[1] * R[1]
                   + 2 * R[2] * R[2]
                   + 6 * R[0] * R[1]
                   + 2 * R[0] * R[2]
                   + 6 * R[1] * R[2]);
    C[1][2] = D * (R[0] * R[0]
                   + 3 * R[1] * R[1]
                   + 3 * R[2] * R[2]
                   + 2 * R[0] * R[1]
                   + 2 * R[0] * R[2]
                   + 4 * R[1] * R[2]);
    C[2][1] = C[1][2];
    C[2][2] = D * (2 * R[0] * R[0]
                   + 2 * R[1] * R[1]
                   + 12 * R[2] * R[2]
                   + 2 * R[0] * R[1]
                   + 6 * R[0] * R[2]
                   + 6 * R[1] * R[2]);
}


void fem_triangle_element::Column_F(double_t F[DIMENSION], double_t q, char element_side)
{
    double_t k = 0;

    F[0] = 0;
    F[1] = 0;
    F[2] = 0;

    switch(element_side)
    {
        case 't':
            /*k = 2.0 * PI * GetSecondSideLength() * q / 6.0;*/


            /* Задаем граничное условие второго рода на стороне t */
            /*if fabs(- RATIO_Y_TO_X * (second_point.x - TRIANGLE_BASE) - second_point.y) < EPS_T
                && fabs(- RATIO_Y_TO_X * (third_point.x - TRIANGLE_BASE) - third_point.y) < EPS_T  )
                    (third_point.y == TRIANGLE_HEIGHT) {
                F[1] = k * (2.0 * second_point.rad_vector() + third_point.rad_vector());
                //F[2] = k * (second_point.rad_vector() + 2.0 * third_point.rad_vector());
            }
            break;*/

            double_t k1 = 2.0 * PI * GetFirstSideLength() * q / 6.0;
            double_t k2 = 2.0 * PI * GetSecondSideLength() * q / 6.0;
            double_t k3 = 2.0 * PI * GetThirdSideLength() * q / 6.0;

            if (fabs(- RATIO_Y_TO_X * (second_point.x - TRIANGLE_BASE) - second_point.y) < EPS_T
                && fabs(- RATIO_Y_TO_X * (third_point.x - TRIANGLE_BASE) - third_point.y) < EPS_T  )
            {
                F[1] += k2 * (2.0 * second_point.rad_vector() + third_point.rad_vector());
                F[2] += k2 * (second_point.rad_vector() + 2.0 * third_point.rad_vector());
            }

            /*if (first_point.y < EPS_T * STEP_X )
            {
                F[0] += k1 * (2.0 * first_point.rad_vector() + second_point.rad_vector());
                F[1] += k1 * (first_point.rad_vector() + 2.0 * second_point.rad_vector());
            }

            if (first_point.x < EPS_T * STEP_X)
            {
                F[0] += k3 * (2.0 * first_point.rad_vector() + third_point.rad_vector());
                F[2] += k3 * (first_point.rad_vector() + 2.0 * third_point.rad_vector());
            }*/
    }

    /* TODO пример реализации граничного условия первого или третьего рода на стороне f */
    /*if (Coeff_H > EPS_T) {
        double_t k1 = GetFirstSideLength() * OUTER_TEMPERATURE * 2.0 * PI * Coeff_H / 6;
        if (fabs(first_point.y) < EPS_T && fabs(second_point.y) < EPS_T ) {
            F[0] = - k1 * (2.0 * first_point.rad_vector() + second_point.rad_vector());
            F[1] = - k1 * (first_point.rad_vector() + 2.0 * second_point.rad_vector());
        }
    }*/

}


void fem_triangle_element::GetN(double_t * N, double_t * Phi) {
    double_t a[3];
    double_t b[3];
    double_t c[3];
    coef_a(a);
    coef_b(b);
    coef_c(c);

    for (auto i = 0; i < DIMENSION; i++)
        N[i] = 0;

    for (auto i = 0; i < DIMENSION; i++) {
        N[0] += 1 / (2 * GetMatrixADeterminant()) * (a[i]) * Phi[i];
        N[1] += 1 / (2 * GetMatrixADeterminant()) * (b[i]) * Phi[i]; // при rho
        N[2] += 1 / (2 * GetMatrixADeterminant()) * (c[i]) * Phi[i]; // при z
    }
}