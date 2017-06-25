#include <iostream>

#include <Eigen/Sparse>
#include <vector>

#include "square_triangle_grid/square_triangle_grid.h"

#include "Rocket/air_parameters.h"

#include "variables.h"

#include "debug.h"


typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


int main(int argc, char** argv)
{
    // производим триангуляцию прямоугольного треугольника
    SquareTriangleGrid triangle_grid((points){A_x});
    //Массив координат узлов треугольников, полученных после триангуляции
    points *coordinates = new points[DOTS_NUMBER];
    triangle_grid.GetGreed(STEP_X, coordinates, false);


    SpMat Result_M(DOTS_NUMBER, DOTS_NUMBER);
    std::vector<T> tripletList;
    tripletList.reserve(DOTS_NUMBER);
    Eigen::VectorXd b(DOTS_NUMBER);

    ShowRunParameters();

    // инициализируем переменные из variables.h
    InitDimension3();
    InitDimensionN();


    /************************************************************************************/
    /*
     * Идем по всем элементам
     */
    /************************************************************************************/


    for (int element = 0; element < TRIANGLES_NUMBER; element++) {
        /*
         * Вычисляем матрицы коэффициентов K_elem, C_elem и вектор правых частей F_elem для элемента
         * Предыдущие значения обнуляются в используемых функциях
         */

        triangle_grid.triangles_array[element]->Matrix_K(K_elem, 'f');
        triangle_grid.triangles_array[element]->Matrix_C(C_elem);


        /**
         * Запоминаем индексы глобальной матрицы и глобального вектора,
         * с которыми потом будем объединять значения, полученные на элементе
         */
        index[0] = triangle_grid.triangles_array[element]->first_point.point_num;
        index[1] = triangle_grid.triangles_array[element]->second_point.point_num;
        index[2] = triangle_grid.triangles_array[element]->third_point.point_num;

        /* Накапливаем коэффициенты */
        for (int i = 0; i < DIMENSION; i++) {
            /* по формуле 11.21 */
            //F[index[i]] += (F_elem[i] + F_old[index[i]]);
            for (int j = 0; j < DIMENSION; j++) {
                /* левая часть */
                //Result_matrix[index[i]][index[j]] += ((C_elem[i][j] * 2.0 / TAU) + K_elem[i][j]);
                tripletList.push_back(T(index[i],index[j], ((C_elem[i][j] * 2.0 / TAU) + K_elem[i][j])));
                /* правая часть */
                R[index[i]][index[j]] += ((C_elem[i][j] * 2.0 / TAU) - K_elem[i][j]);

            }
        }
    } // Конец цикла по элементам

    Result_M.setFromTriplets(tripletList.begin(), tripletList.end());
    Result_M.makeCompressed();



/*    result_matrix << endl;
    result_matrix << "---------------------------------------------------" << endl;
    result_matrix << endl;
    for (auto i = 0; i < DOTS_NUMBER; i++) {
        for (auto j = 0; j < DOTS_NUMBER; j++) {
            result_matrix <<  Result_matrix[i][j] << " ";
        }
        result_matrix << endl;
    }

    result_m << endl;
    result_m << "---------------------------------------------------" << endl;
    result_m << endl;
    result_m << Result_M;
    cout << "------------" << endl;*/

    /************************************************************************************/
    /*
     * Ищем решения по временным слоям
     * Пока не начнется плавление наконечника
     */
    /************************************************************************************/
    while (Temperature[0] < MELTING_TEMPERATURE) {

        average_temp = 0;
        flight_time += TAU;

        for(auto i = 0; i < DOTS_NUMBER; i++)
            average_temp +=  Temperature[i] / DOTS_NUMBER;

        /**
         * Вычисляем тепловой поток на границе в зависимости от времени полета
         */
        //heat_flow = - func_calculate_q(flight_time);
        /*if (heat_flow >= -100) {
            flight_time += 0.05;
            heat_flow = - func_calculate_q(flight_time);
        }*/
        heat_flow = - 100;
        

        /**
         * Обнуляем матрицы коэффициентов и вектор правой части элементов
         * перед заходом на новый временной слой
         */
        for (auto i = 0; i < DIMENSION; i++) {
            F_elem[i] = 0;
        }

        /**
         * Обнуление матрицы жесткости и вектора правой части итоговой системы
         * перед заходом на новый временной слой
         * Также запоминаем значение вектора правой части с предыдущего временного слоя
         */
        for (auto i = 0; i < DOTS_NUMBER; i++) {
            F_old[i] = F[i];
            F[i] = 0;
            //Result_matrix[i][DOTS_NUMBER] = 0;
            b[i] = 0;
        }

        for (auto element = 0; element < TRIANGLES_NUMBER; element++) {
            triangle_grid.triangles_array[element]->Column_F(F_elem, heat_flow, 't');

            index[0] = triangle_grid.triangles_array[element]->first_point.point_num;
            index[1] = triangle_grid.triangles_array[element]->second_point.point_num;
            index[2] = triangle_grid.triangles_array[element]->third_point.point_num;

            for (auto i = 0; i < DIMENSION; i++) {
                /* по формуле 11.21 */
                F[index[i]] += (F_elem[i] + F_old[index[i]]);
                //cout << index[i] << " " << F[index[i]] << endl;
            }
        } // конец цикла по элементам
        
        /*
         * Заполняыем  итоговый вектор правой части уравнения 11.23
         * Он состоит из поэлементно накопленных значений R + F
         */
        for (auto i = 0; i < DOTS_NUMBER; i++) {
            for (auto j = 0; j <DOTS_NUMBER; j++)
                //Result_matrix[i][DOTS_NUMBER] += (R[i][j] * Temperature[j]);
                b[i] += (R[i][j] * Temperature[j]);
        }
        for (auto i = 0; i < DOTS_NUMBER; i++) {
            //Result_matrix[i][DOTS_NUMBER] -= F[i];
            b[i] -= F[i];
            //cout << Result_matrix[i][DOTS_NUMBER] << endl;
        }
        //for (auto i = 0; i < DOTS_NUMBER; i++) {
        //    b[i] = Result_matrix[i][DOTS_NUMBER];
        //}

        /**
         * Вычисляем температуру на новом слое по времени
         */
        //equation.Solve(Temperature, Result_matrix);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        solver.compute(Result_M);
        //Eigen::SimplicialCholesky<SpMat> chol(Result_M);  // performs a Cholesky factorization of A
        //Eigen::VectorXd x = chol.solve(b);
        Eigen::VectorXd x = solver.solve(b);

        for (auto i = 0; i < DOTS_NUMBER; i++) {
            Temperature[i] = x[i];
            cout << Temperature[i] << endl;
        }
        cout <<" ===== " << endl;
        ShowTemperature();


    } // Конец цикла по времени

    result_matrix.close();
    result_m.close();
    triangle_temp.close();
    return 0;
}