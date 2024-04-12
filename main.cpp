//#include "boost/multiprecision/cpp_dec_float.hpp"
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
double gran_func_test(double x, double y, int i,int j,int n,int m){
    if(i==0)
        return mu1_test(y);
    if(i==n)
        return mu2_test(y);
    if(j==0)
        return mu3_test(x);
    if(j==m)
        return mu4_test(x);
}
double gran_func_main(double x, double y, int i,int j,int n,int m) {
    if(i==0)
        return mu1_main(y);
    if(i == n)
        return mu2_main(y);
    if(j == 0)
        return mu3_main(x);

    if(j == m)
        return mu4_main(x);
}
//double gran_func(double x, double y, int i,int j,int n,int m, bool test = true){
//    if(test){
//        if(i==0)
//            return mu1_test(y);
//        if(i==n)
//            return mu2_test(y);
//        if(j==0)
//            return mu3_test(x);
//        if(j==m)
//            return mu4_test(x);
//    }
//
//    if(i==0)
//        return mu1_main(y);
//    if(i == n)
//        return mu2_main(y);
//    if(j == 0)
//        return mu3_main(x);
//
//    if(j == m)
//        return mu4_main(x);
//}

// Функция для решения тестовой задачи
// n - Число разбиений по оси x
// m - Число разбиений по оси y
// Nmax - Максимальное количество итераций
// eps - требуемая погрешность
void solve_test(int n, int m,int Nmax,double eps){
    double h = (b-a)/n; // Шаг по оси x
    double k = (d-c)/m; // Шаг по оси y

    vector<vector<double>> v,h_s,r_s;

    string str = ""; // Сообщение при окончании работы

    int S = 0; // Количество итераций метода


    double epsMax = 0; // Погрешность на текущей итерации

    double h2 =  ((n/(b-a)) * (n / (b-a))); // Вспомогательная величина 1/h^2
    double k2 =  ((m/(d-c)) * (m/(d-c))); // Вспомогательная величина 1/k^2

    double a2 = - 2 * (h2 + k2); // Вспомогательная величина -2(1/h^2+1/k^2)

    vector<double> x(n+1); // Абсциссы узлов
    vector<double> y(m+1); // Ординаты узлов

    x[0] = a;
#pragma omp parallel for default(none) shared(n,x,h,a)
    for (int i = 1; i <= n; ++i)
        x[i] = a + i*h;

    y[0] = c;
#pragma omp parallel for default(none) shared(m,y,k,c)
    for (int j = 1; j <= m; ++j)
        y[j] = c + j*k;
    // Заполнение
    vector<double> vec;
    vec.reserve(m+1);

#pragma omp parallel for default(none) shared(m,vec)
    for(int t = 0;t<=m;t++){
        vec.push_back(0);
    }

// Нельзя параллелить из-за особенностей работы памяти
    for (int i = 0; i <= n; ++i) {
        v.push_back(vec);
        h_s.push_back(vec);
        r_s.push_back(vec);
    }

    // Учёт гран условий
#pragma omp parallel for default(none) shared(m,n,v,x,y)
    for(int j = 0;j<=m;j++)
        for(int i = 0;i<=n;i++)
            if((i==0)||(i==n)||(j==0)||(j==m))
                v[i][j] = gran_func_test(x[i],y[j],i,j,n,m);

    double alpha ; // alpha(s)
    double beta = 0.0; // beta(s)
    //int system_size = (n-1)*(m-1); // Размерность системы

    // Начинаем итерации метода сопр. градиентов

    double current_norm; // Текущая норма невязки
    // Вспомогательные суммы
    double s_up;
    double s_down;

    double prev_eps;

    cout << "\n-----------------------------------------"<<endl;

    auto start = high_resolution_clock::now(); // Для замера времени

    while (true){

        if(S>= Nmax){
            cout << "Текущая погрешность "<<epsMax<<endl;
            cout << "Норма невязки = "<<current_norm<<endl;
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

        double local_epsMax = epsMax;

#pragma omp parallel for default(none) shared(m, n, v, h_s, alpha, x,y) reduction(max:local_epsMax)
        for(int j = 1; j < m; j++) {
            for(int i = 1; i < n; i++) {
                double v_old_local = v[i][j];
                double v_new_local = v_old_local + alpha * h_s[i][j];
                v[i][j] = v_new_local;
                local_epsMax = fmax(local_epsMax, fabs(u_test(x[i],y[j]) - v_new_local));
            }
        }

        epsMax = local_epsMax;
        if(S%200 == 0)
            cout << "beta на "<< S<<"-й итерации = "<<beta << " Текущая погрешность "<<epsMax<< " Норма невязки = "<<current_norm<<endl;
        ++S;

        if(S==1 || S==2){

            cout << "Текущая погрешность "<<epsMax<<endl;
        }
        if((epsMax < eps)||((S>7001)&&(epsMax>prev_eps))){
            cout << "Норма невязки = "<<current_norm<<endl;
            cout << "\n-----------------------------------------"<<endl;
            str = "Выход по точности!";
            break;
        }
        prev_eps = epsMax;
        if(S==1||S==2){
            cout << "Норма невязки = "<<current_norm<<endl;
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
}
// Функция для решения основной задачи
// n - Число разбиений по оси x
// m - Число разбиений по оси y
// Nmax - Максимальное количество итераций
// eps - требуемая погрешность
void solve_main(int n, int m,int Nmax,double eps){
    double h = (b-a)/n; // Шаг по оси x
    double k = (d-c)/m; // Шаг по оси y

    vector<vector<double>> v,h_s,r_s;

    string str = ""; // Сообщение при окончании работы

    int S = 0; // Количество итераций метода


    double epsMax = 0; // Точность на текущей итерации

    double h2 =  ((n/(b-a)) * (n / (b-a))); // Вспомогательная величина 1/h^2
    double k2 =  ((m/(d-c)) * (m/(d-c))); // Вспомогательная величина 1/k^2

    double a2 = - 2 * (h2 + k2); // Вспомогательная величина -2(1/h^2+1/k^2)

    vector<double> x(n+1); // Абсциссы узлов
    vector<double> y(m+1); // Ординаты узлов

    x[0] = a;
#pragma omp parallel for default(none) shared(n,x,h,a)
    for (int i = 1; i <= n; ++i)
        x[i] = a + i*h;

    y[0] = c;
#pragma omp parallel for default(none) shared(m,y,k,c)
    for (int j = 1; j <= m; ++j)
        y[j] = c + j*k;
    // Заполнение
    vector<double> vec;
    vec.reserve(m+1);

#pragma omp parallel for default(none) shared(m,vec)
    for(int t = 0;t<=m;t++){
        vec.push_back(0);
    }

// Нельзя параллелить из-за особенностей работы операции вставки в вектор
    for (int i = 0; i <= n; ++i) {
        v.push_back(vec);
        h_s.push_back(vec);
        r_s.push_back(vec);
    }

    // Учёт гран условий
#pragma omp parallel for default(none) shared(m,n,v,x,y)
    for(int j = 0;j<=m;j++)
        for(int i = 0;i<=n;i++)
            if((i==0)||(i==n)||(j==0)||(j==m))
                v[i][j] = gran_func_main(x[i],y[j],i,j,n,m);

    double alpha ; // alpha(s)
    double beta = 0.0; // beta(s)
    //int system_size = (n-1)*(m-1); // Размерность системы

    // Начинаем итерации метода сопр. градиентов

    double current_norm; // Текущая норма невязки
    // Вспомогательные суммы
    double s_up;
    double s_down;

    double prev_eps;

    cout << "\n-----------------------------------------"<<endl;

    auto start = high_resolution_clock::now(); // Для замера времени

    while (true){

        if(S>= Nmax){
            cout << "Текущая точность "<<epsMax<<endl;
            cout << "Норма невязки = "<<current_norm<<endl;
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
                r_s[i][j] = -a2 * v[i][j]-h2*(v[i+1][j]+v[i-1][j])-k2*(v[i][j+1]+v[i][j-1])- f_main(x[i],y[j]);
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

        double local_epsMax = epsMax;

#pragma omp parallel for default(none) shared(m, n, v, h_s, alpha, x,y) reduction(max:local_epsMax)
        for(int j = 1; j < m; j++) {
            for(int i = 1; i < n; i++) {
                double v_old_local = v[i][j];
                double v_new_local = v_old_local + alpha * h_s[i][j];
                v[i][j] = v_new_local;
                local_epsMax = fmax(local_epsMax, fabs(v_old_local - v_new_local));
            }
        }

        epsMax = local_epsMax;
        if(S%200 == 0)
            cout << "beta на "<< S<<"-й итерации = "<<beta << " Текущая точность "<<epsMax<< " Норма невязки = "<<current_norm<<endl;
        ++S;

        if(S==1 || S==2){

            cout << "Текущая точность "<<epsMax<<endl;
        }
        if((epsMax < eps)||((S>7001)&&(epsMax>prev_eps))){
            cout << "Норма невязки = "<<current_norm<<endl;
            cout << "\n-----------------------------------------"<<endl;
            str = "Выход по точности!";
            break;
        }
        prev_eps = epsMax;
        if(S==1||S==2){
            cout << "Норма невязки = "<<current_norm<<endl;
            cout << "\n-----------------------------------------"<<endl;
        }

    }
    auto stop = high_resolution_clock::now(); // Для замера времени
    auto duration = duration_cast<milliseconds>(stop - start);

    cout << "Проведено "<<S<<" итераций"<<endl;
    cout << "Достигнута точность eps = " << epsMax <<endl;
    cout << "Достигнута норма невязки "<<current_norm<<endl;
    cout << str << endl;
    cout << "Время просчёта "<<duration.count() << " мс"<<endl;
}
// Одна итерация метода
void iterate_conjugate_gradients()
int main() {
    system("chcp 65001");
    cout << "<<МЕТОД СОПРЯЖЁННЫХ ГРАДИЕНТОВ>>\nЩербаков Павел, гр. 3821Б1ПМоп2, Вар. 1"<<endl;
    // Применение метода сопряжённых градиентов к решению
    // Задачи Дирихле для уравнения Пуассона

    int n = 15; // Число разбиений по оси x 1500
    int m = 15; // Число разбиений по оси y 1501

    int Nmax = 20; // Максимальное число итераций 20 000

    double eps = 0.0000005; // Требуемая точность/погрешность для основной/тестовой задачи

    bool test; // Является ли задача основной/тестовой

    cout << "Выберите задачу: 0 - тестовая, 1 - основная"<<endl;
    cin >> test;

    if(!test)
        solve_test(n,m,Nmax,eps);
    else
        solve_main(n,m,Nmax,eps);
    return 0;
}