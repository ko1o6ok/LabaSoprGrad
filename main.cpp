#include <iostream>
#include <cmath>
#include <vector>
# define M_PI           3.14159265358979323846
#include <chrono>
#include <queue>
using namespace std::chrono;
using namespace std;



// Границы прямоугольной области
const double a = 0.0;
const double b = 1.0;
const double c = 0.0;
const double d = 1.0;



// Вариант - 1
// Функция f - плотность источников/стоков
// Основная задача
double f_main(double x,double y){
    double temp = sin(M_PI * x * y);
    return temp * temp;
}
// Основные гранусловия
double mu1_main(double y){
    return sin(M_PI * y);
}
double mu2_main(double y){
    return sin(M_PI * y);
}
double mu3_main(double x){
    return x-x*x;
}
double mu4_main(double x){

    return x-x*x;
}
// Функция u для тестовой задачи
// u* = exp(sin^2(pi * x * y))
double u_test(double x, double y){
    return exp(f_main(x,y));
}
// Функция f для тестовой задачи
double f_test(double x, double y){
    return 0.5 * u_test(x,y) * M_PI * M_PI * (x * x + y * y) * (-1-4 * cos(2 * M_PI * x * y) + cos(4 * M_PI * x * y));
}

// Тестовые гранусловия
double mu1_test(double y){
    return 1;
}
double mu2_test(double y){
    double temp = sin(M_PI * y);
    return exp(temp * temp);
}
double mu3_test(double x){
    return 1;
}
double mu4_test(double x){
    double temp = sin(M_PI * x);
    return exp(temp * temp);
}

// Совокупность граничных функций
// x - абсцисса точки
// y - ордината точки
// test - является ли задача тестовой
double gran_func(double x, double y, int i,int j,int n,int m, bool test = true){
    if(test){
        if(i==0)
            return mu1_test(y);
        if(i==n)
            return mu2_test(y);
        if(j==0)
            return mu3_test(x);
        if(j==m)
            return mu4_test(x);
    }

    if(i==0)
        return mu1_main(y);
    if(i == n)
        return mu2_main(y);
    if(j == 0)
        return mu3_main(x);

    if(j == m)
        return mu4_main(x);
}


int main() {
    system("chcp 65001");

    // Применение метода сопряжённых градиентов к решению
    // Задачи Дирихле для уравнения Пуассона

    int n = 290; // Число разбиений по оси x
    int m = 290; // Число разбиений по оси y


    double h = (b-a)/n; // Шаг по оси x
    double k = (d-c)/m; // Шаг по оси y

    //double laplacian[n+1][m+1]; // Кэширование лапласиана

    double v[n+1][m+1]; // Текущее состояние
    // Взять вместо них векторы, чтобы ДИНАМИЧЕСКИ формировалась память
    // Тогда оно разрешит жрать больше
    double h_s[n+1][m+1]; // Текущий вектор сопряжённого направления
    double r_s[n+1][m+1]; // Текущая невязка

    //double real_sol[n+1][m+1]; // Точное решение




    string str = ""; // Сообщение при окончании работы

    int S = 0; // Количество итераций метода
    int Nmax = 5000; // Максимальное число итераций

    double eps = 0.000005; // Требуемая точность для основной/тестовой задачи

    bool test = true; // Является ли задача основной/тестовой

    double epsMax = 0; // Точность на текущей итерации
    double epsCur = 0; // Вспомогательная точность

    double v_old; // Старое значение в узле
    double v_new; // Новое значение в узле

    double h2 =  ((n/(b-a)) * (n / (b-a))); // Вспомогательная величина 1/h^2
    double k2 =  ((m/(d-c)) * (m/(d-c))); // Вспомогательная величина 1/k^2

    double a2 = - 2 * (h2 + k2); // Вспомогательная величина -2(1/h^2+1/k^2)

    vector<double> x(n+1); // Абсциссы узлов
    vector<double> y(m+1); // Ординаты узлов

    x[0] = a;
    #pragma omp parallel for default(none) shared(n,x,h,a)
    for (int i = 1; i <= n; ++i)
        x[i] = a + i*h;
        //x[i] = x[i-1] + h;

    y[0] = c;
    #pragma omp parallel for default(none) shared(m,y,k,c)
    for (int j = 1; j <= m; ++j)
        y[j] = c + j*k;
        //y[j] = y[j-1] + k;
    // Заполнение
    for(int j = 0;j<=m;j++){
        #pragma omp parallel for default(none) shared(n,v,j,h_s,r_s,x,y,m,test)
        for(int i = 0;i<=n;i++){
            v[i][j] = 0;
            h_s[i][j] = 0;
            r_s[i][j] = 0;
            //real_sol[i][j] = u_test(x[i],y[j]);
            //laplacian[i][j] = f_test(x[i],y[j]);
        }
    }

    // Учёт гран условий
    #pragma omp parallel for default(none) shared(m,n,v,x,y,test)
    for(int j = 0;j<=m;j++)
        for(int i = 0;i<=n;i++)
            if((i==0)||(i==n)||(j==0)||(j==m))
                v[i][j] = gran_func(x[i],y[j],i,j,n,m,test);

    double alpha ; // alpha(s)
    double beta = 0.0; // beta(s)
    //int system_size = (n-1)*(m-1); // Размерность системы

    // Начинаем итерации метода сопр. градиентов

    double current_norm; // Текущая норма невязки
    // Вспомогательные суммы
    double s_up = 0.0;
    double s_down = 0.0;

    double prev_eps;

//    cout << "Посчитанное решение"<<" на "<<S<<"-той итерации"<<": \n";
//    for(int j = 0;j<=m;j++){
//        for(int i = 0;i<=n;i++){
//            cout << v[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl << "Реальное решение: \n";
//    for(int j = 0;j<=m;j++){
//        for(int i = 0;i<=n;i++){
//            cout << u_test(x[i],y[j]) << " ";
//        }
//        cout << endl;
//    }


    cout << "\n-----------------------------------------"<<endl;

    auto start = high_resolution_clock::now(); // Для замера времени

    bool alert_enabled = false; // beta слишком маленькая - > - > вычисл. погрешность
    while (true){

        if(S>= Nmax){
            cout << "Текущая точность "<<epsMax<<endl;
            cout << "Норма невязки = "<<current_norm<<endl;
//            cout << "Посчитанное решение"<<" на "<<S<<"-той итерации"<<": \n";
//            for(int j = 0;j<=m;j++){
//                for(int i = 0;i<=n;i++){
//                    cout <<v[i][j] << " ";
//                }
//                cout << endl;
//            }
//            cout << endl << "Реальное решение: \n";
//            for(int j = 0;j<=m;j++){
//                for(int i = 0;i<=n;i++){
//                    cout << u_test(x[i],y[j]) << " ";
//                }
//                cout << endl;
//            }
            cout << "\n-----------------------------------------"<<endl;

            str = "Выход по числу итераций!";
            break;
        }

        epsMax = 0;
        current_norm = -10;
        // Проведём итерацию методом сопряжённых градиентов.
        // Начинаем с подсчёта текущей невязки
        // r_s = A x_s-F


        #pragma omp parallel for default(none) shared(m,n,r_s,a2,v,h2,k2,x,y) reduction(max:current_norm)
        for(int i = 1;i<n;i++){
            for(int j = 1;j<m;j++){
                // Проходимся по вектору v
                // Разреженность помогает!

                r_s[i][j] = -a2 * v[i][j]-h2*(v[i+1][j]+v[i-1][j])-k2*(v[i][j+1]+v[i][j-1])- f_test(x[i],y[j]);

                current_norm = fmax(current_norm, fabs(r_s[i][j]));
            }
        }


        if(S>=1){

                // Вспомогательные суммы
                s_up = 0.0;
                s_down = 0.0;
                #pragma omp parallel for default(none) shared(m,n,a2,h_s,r_s,h2,k2,s_up,s_down)
                for(int j = 1;j<m;j++){
                    double temp ; // Вспомогательная переменная = A h_s
                    for(int i = 1;i<n;i++){
                        temp = -a2 * h_s[i][j]-h2*(h_s[i+1][j]+h_s[i-1][j])-k2*(h_s[i][j+1]+h_s[i][j-1]);
                        s_up += temp * r_s[i][j];
                        s_down += temp * h_s[i][j];
                    }
                }

                beta = s_up/s_down;
//                if((abs(beta) < 1e-10)&&(S!=0)&&(abs(beta)>0)){
//
//                    cout << "Alert! Производится перезапуск метода" << endl;
//                    beta = 0;
//                    for(int j = 0;j<=m;j++){
//                        //#pragma omp parallel for default(none) shared(n,v,j,h_s,r_s,x,y,m,test)
//                        for(int i = 0;i<=n;i++){
//                            h_s[i][j] = 0;
//                            //real_sol[i][j] = u_test(x[i],y[j]);
//                            //laplacian[i][j] = f_test(x[i],y[j]);
//                        }
//                    }
//                }



        }

            // Пересчитываем значение h_s
        #pragma omp parallel for default(none) shared(m,n,h_s,beta,r_s)
        for(int j = 1;j<m;j++){
            for(int i = 1;i<n;i++){
                h_s[i][j] = beta * h_s[i][j] - r_s[i][j];
            }
        }

        s_up = 0.0;
        s_down = 0.0;
        #pragma omp parallel for default(none) shared(m,n,h_s,r_s,a2,h2,k2) reduction(+:s_up,s_down)
        for(int j = 1;j<m;j++)
            for(int i = 1;i<n;i++){
                s_up += h_s[i][j] * r_s[i][j];
                s_down += (-a2 * h_s[i][j]-h2*(h_s[i+1][j]+h_s[i-1][j])-k2*(h_s[i][j+1]+h_s[i][j-1])) * h_s[i][j];
            }
        alpha = -s_up/s_down;

        double local_epsMax = epsMax; // Thread-local variable

#pragma omp parallel for default(none) shared(m, n, v, h_s, alpha, test, x,y) reduction(max:local_epsMax)
        for(int j = 1; j < m; j++) {
            for(int i = 1; i < n; i++) {
                double v_old_local = v[i][j];
                double v_new_local = v_old_local + alpha * h_s[i][j];

                v[i][j] = v_new_local;

                if(test)
                    local_epsMax = fmax(local_epsMax, fabs(u_test(x[i],y[j]) - v_new_local));
                else
                    local_epsMax = fmax(local_epsMax, fabs(v_old_local - v_new_local));
            }
        }

// Update shared epsMax after all calculations are done
        epsMax = local_epsMax;
        //cout << "beta на "<< S<<"-й итерации = "<<beta << " Текущая точность "<<epsMax<< " Норма невязки = "<<current_norm<<endl;

        // Из-за вычислительной погрешности метод приходится перезапускать


        ++S;

        if(S==1 || S==2){

            cout << "Текущая точность "<<epsMax<<endl;
        }

        if((epsMax < eps)||((epsMax>prev_eps)&&(S>10))){
            cout << "Норма невязки = "<<current_norm<<endl;
//            cout << "Посчитанное решение"<<" на "<<S<<"-той итерации"<<": \n";
//            for(int j = 0;j<=m;j++){
//                for(int i = 0;i<=n;i++){
//                    cout << v[i][j] << " ";
//                }
//                cout << endl;
//            }
//            cout << endl << "Реальное решение: \n";
//            for(int j = 0;j<=m;j++){
//                for(int i = 0;i<=n;i++){
//                    cout << u_test(x[i],y[j]) << " ";
//                }
//                cout << endl;
//            }
            cout << "\n-----------------------------------------"<<endl;
            str = "Выход по точности!";
            break;
        }
        prev_eps = epsMax;
        if(S==1||S==2){
            cout << "Норма невязки = "<<current_norm<<endl;
//            cout << "Посчитанное решение"<<" на "<<S<<"-той итерации"<<": \n";
//            for(int j = 0;j<=m;j++){
//                for(int i = 0;i<=n;i++){
//                    cout << v[i][j] << " ";
//                }
//                cout << endl;
//            }
//            cout << endl << "Реальное решение: \n";
//            for(int j = 0;j<=m;j++){
//                for(int i = 0;i<=n;i++){
//                    cout << u_test(x[i],y[j]) << " ";
//                }
//                cout << endl;
//            }
            cout << "\n-----------------------------------------"<<endl;
        }

    }
    auto stop = high_resolution_clock::now(); // Для замера времени
    auto duration = duration_cast<milliseconds>(stop - start);

    cout << "Проведено "<<S<<" итераций"<<endl;
    cout << "Достигнута погрешность eps = " << epsMax <<endl;
    cout << "Достигнута норма невязки "<<current_norm<<endl;
    cout << str << endl;
    cout << "Время просчёта "<<duration.count() << " мс"<<endl;
    return 0;
}
