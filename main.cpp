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
void solve_main(int n, int m,int Nmax,double eps,int Nmax2,double eps2){
    double h = (b-a)/n; // Шаг по оси x
    double k = (d-c)/m; // Шаг по оси y

    // Для контрольной сетки
    int n_ = 2*n;
    int m_ = 2*m;
    double h_ = (b-a)/n_; // Шаг по оси x
    double k_ = (d-c)/m_; // Шаг по оси y
    // ------------------------------
    vector<vector<double>> v,h_s,r_s;
    //
    vector<vector<double>> v_,h_s_,r_s_;

    int Nmax_ = Nmax2;
    double eps_ = eps2;

    string str = ""; // Сообщение при окончании работы
    string str_ = ""; // Сообщение при окончании работы

    int S = 0; // Количество итераций метода

    int S_ = 0;

    double precision; // Точность решения основной задачи

    double epsMax = 0; // Точность метода на текущей итерации
    double epsMax_ = 0;

    double h2 =  ((n/(b-a)) * (n / (b-a))); // Вспомогательная величина 1/h^2
    double k2 =  ((m/(d-c)) * (m/(d-c))); // Вспомогательная величина 1/k^2

    double h2_ = ((n_/(b-a)) * (n_ / (b-a)));
    double k2_ =  ((m_/(d-c)) * (m_/(d-c)));

    double a2 = - 2 * (h2 + k2); // Вспомогательная величина -2(1/h^2+1/k^2)

    double a2_  = - 2 * (h2_ + k2_);

    vector<double> x(n+1); // Абсциссы узлов
    vector<double> y(m+1); // Ординаты узлов

    vector<double> x_(n_+1); // Абсциссы узлов
    vector<double> y_(m_+1); // Ординаты узлов

    x[0] = a;

#pragma omp parallel for default(none) shared(n,x,h,a)
    for (int i = 1; i <= n; ++i){
        x[i] = a + i*h;
    }
    x_[0] = a;
#pragma omp parallel for default(none) shared(n_,x_,h_,a)
    for (int i = 1; i <= n_; ++i){
        x_[i] = a + i*h_;
    }

    y[0] = c;
#pragma omp parallel for default(none) shared(m,y,k,c)
    for (int j = 1; j <= m; ++j)
        y[j] = c + j*k;

    y_[0] = c;
#pragma omp parallel for default(none) shared(m_,y_,k_,c)
    for (int j = 1; j <= m_; ++j)
        y_[j] = c + j*k_;
    // Заполнение
    vector<double> vec;
    vec.reserve(m+1);

#pragma omp parallel for default(none) shared(m,vec)
    for(int t = 0;t<=m;t++){
        vec.push_back(0);
    }

    vector<double> vec_;
    vec_.reserve(m_+1);

#pragma omp parallel for default(none) shared(m_,vec_)
    for(int t = 0;t<=m_;t++){
        vec_.push_back(0);
    }

// Нельзя параллелить из-за особенностей работы операции вставки в вектор
    for (int i = 0; i <= n; ++i) {
        v.push_back(vec);
        h_s.push_back(vec);
        r_s.push_back(vec);
    }

    for (int i = 0; i <= n_; ++i) {
        v_.push_back(vec_);
        h_s_.push_back(vec_);
        r_s_.push_back(vec_);
    }

    // Учёт гран условий
#pragma omp parallel for default(none) shared(m,n,v,x,y)
    for(int j = 0;j<=m;j++)
        for(int i = 0;i<=n;i++)
            if((i==0)||(i==n)||(j==0)||(j==m))
                v[i][j] = gran_func_main(x[i],y[j],i,j,n,m);

#pragma omp parallel for default(none) shared(m_,n_,v_,x_,y_)
    for(int j = 0;j<=m_;j++)
        for(int i = 0;i<=n_;i++)
            if((i==0)||(i==n_)||(j==0)||(j==m_))
                v_[i][j] = gran_func_main(x_[i],y_[j],i,j,n_,m_);

    double alpha ; // alpha(s)
    double beta = 0.0; // beta(s)

    double alpha_ ; // alpha(s)
    double beta_ = 0.0; // beta(s)

    //int system_size = (n-1)*(m-1); // Размерность системы

    // Начинаем итерации метода сопр. градиентов

    double current_norm; // Текущая норма невязки
    // Вспомогательные суммы
    double s_up;
    double s_down;

    //double prev_eps;

    double current_norm_; // Текущая норма невязки
    // Вспомогательные суммы
    double s_up_;
    double s_down_;

    //double prev_precision;

    cout << "\n-----------------------------------------"<<endl;

    auto start = high_resolution_clock::now(); // Для замера времени

    while (true) {


        // Стандартная сетка
        // ------------------------------------------------------
        epsMax = 0;
        current_norm = -10;
        // Проведём итерацию методом сопряжённых градиентов.
        // Начинаем с подсчёта текущей невязки
        // r_s = A x_s-F
#pragma omp parallel for default(none) shared(m, n, r_s, a2, v, h2, k2, x, y) reduction(max:current_norm)
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < m; j++) {
                // Проходимся по вектору v
                r_s[i][j] = -a2 * v[i][j] - h2 * (v[i + 1][j] + v[i - 1][j]) - k2 * (v[i][j + 1] + v[i][j - 1]) -
                            f_main(x[i], y[j]);
                current_norm = fmax(current_norm, fabs(r_s[i][j]));
            }
        }
        if (S >= 1) {
            // Вспомогательные суммы
            s_up = 0.0;
            s_down = 0.0;
#pragma omp parallel for default(none) shared(m, n, a2, h_s, r_s, h2, k2, s_up, s_down)
            for (int j = 1; j < m; j++) {
                double temp; // Вспомогательная переменная = A h_s
                for (int i = 1; i < n; i++) {
                    temp = -a2 * h_s[i][j] - h2 * (h_s[i + 1][j] + h_s[i - 1][j]) -
                           k2 * (h_s[i][j + 1] + h_s[i][j - 1]);
                    s_up += temp * r_s[i][j];
                    s_down += temp * h_s[i][j];
                }
            }
            beta = s_up / s_down;
        }
        // Пересчитываем значение h_s
#pragma omp parallel for default(none) shared(m, n, h_s, beta, r_s)
        for (int j = 1; j < m; j++) {
            for (int i = 1; i < n; i++) {
                h_s[i][j] = beta * h_s[i][j] - r_s[i][j];
            }
        }
        s_up = 0.0;
        s_down = 0.0;
#pragma omp parallel for default(none) shared(m, n, h_s, r_s, a2, h2, k2) reduction(+:s_up, s_down)
        for (int j = 1; j < m; j++)
            for (int i = 1; i < n; i++) {
                s_up += h_s[i][j] * r_s[i][j];
                s_down += (-a2 * h_s[i][j] - h2 * (h_s[i + 1][j] + h_s[i - 1][j]) -
                           k2 * (h_s[i][j + 1] + h_s[i][j - 1])) * h_s[i][j];
            }
        alpha = -s_up / s_down;

        double local_epsMax = epsMax;

#pragma omp parallel for default(none) shared(m, n, v, h_s, alpha, x, y) reduction(max:local_epsMax)
        for (int j = 1; j < m; j++) {
            for (int i = 1; i < n; i++) {
                double v_old_local = v[i][j];
                double v_new_local = v_old_local + alpha * h_s[i][j];
                v[i][j] = v_new_local;
                local_epsMax = fmax(local_epsMax, fabs(v_old_local - v_new_local));
            }
        }
        epsMax = local_epsMax;

        ++S;
        if (S >= Nmax) {
            str = "Выход по числу итераций!";
            break;
        }
        if(epsMax < eps){
            str = "Выход по точности метода!";
            break;
        }
    }
    while(true) {

        // Задача на контрольной сетке
        epsMax_ = 0;
        current_norm_ = -10;
        // Проведём итерацию методом сопряжённых градиентов.
        // Начинаем с подсчёта текущей невязки
        // r_s = A x_s-F
#pragma omp parallel for default(none) shared(m_, n_, r_s_, a2_, v_, h2_, k2_, x_, y_) reduction(max:current_norm_)
        for (int i = 1; i < n_; i++) {
            for (int j = 1; j < m_; j++) {
                // Проходимся по вектору v
                r_s_[i][j] =
                        -a2_ * v_[i][j] - h2_ * (v_[i + 1][j] + v_[i - 1][j]) - k2_ * (v_[i][j + 1] + v_[i][j - 1]) -
                        f_main(x_[i], y_[j]);
                current_norm_ = fmax(current_norm_, fabs(r_s_[i][j]));
            }
        }
        if (S_ >= 1) {
            // Вспомогательные суммы
            s_up_ = 0.0;
            s_down_ = 0.0;
#pragma omp parallel for default(none) shared(m_, n_, a2_, h_s_, r_s_, h2_, k2_, s_up_, s_down_)
            for (int j = 1; j < m_; j++) {
                double temp; // Вспомогательная переменная = A h_s
                for (int i = 1; i < n_; i++) {
                    temp = -a2_ * h_s_[i][j] - h2_ * (h_s_[i + 1][j] + h_s_[i - 1][j]) -
                           k2_ * (h_s_[i][j + 1] + h_s_[i][j - 1]);
                    s_up_ += temp * r_s_[i][j];
                    s_down_ += temp * h_s_[i][j];
                }
            }
            beta_ = s_up_ / s_down_;
        }
        // Пересчитываем значение h_s
#pragma omp parallel for default(none) shared(m_, n_, h_s_, beta_, r_s_)
        for (int j = 1; j < m_; j++) {
            for (int i = 1; i < n_; i++) {
                h_s_[i][j] = beta_ * h_s_[i][j] - r_s_[i][j];
            }
        }
        s_up_ = 0.0;
        s_down_ = 0.0;
#pragma omp parallel for default(none) shared(m_, n_, h_s_, r_s_, a2_, h2_, k2_) reduction(+:s_up_, s_down_)
        for (int j = 1; j < m_; j++)
            for (int i = 1; i < n_; i++) {
                s_up_ += h_s_[i][j] * r_s_[i][j];
                s_down_ += (-a2_ * h_s_[i][j] - h2_ * (h_s_[i + 1][j] + h_s_[i - 1][j]) -
                            k2_ * (h_s_[i][j + 1] + h_s_[i][j - 1])) * h_s_[i][j];
            }
        alpha_ = -s_up_ / s_down_;

        double local_epsMax_ = epsMax_;

#pragma omp parallel for default(none) shared(m_, n_, v_, h_s_, alpha_, x_, y_) reduction(max:local_epsMax_)
        for (int j = 1; j < m_; j++) {
            for (int i = 1; i < n_; i++) {
                double v_old_local = v_[i][j];
                double v_new_local = v_old_local + alpha_ * h_s_[i][j];
                v_[i][j] = v_new_local;
                local_epsMax_ = fmax(local_epsMax_, fabs(v_old_local - v_new_local));
            }
        }

        epsMax_ = local_epsMax_;
        ++S_;

        if (S_ >= Nmax_) {
            str_ = "Выход по числу итераций!";
            break;
        }
        if(epsMax_ < eps_){
            str_ = "Выход по точности метода!";
            break;
        }
    }


    // Подсчёт точности решения основной задачи
    precision = 0;
    // Ищем максимум модуля разности по общим узлам
    #pragma parallel for default(none) shared(n, m, precision, v, v_) reduction(max:precision)
    for (int i = 1; i < n; i++)
        for (int j = 1; j < m; j++) {
            precision = fmax(precision, abs(v[i][j] - v_[2 * i][2 * j]));
        }





    auto stop = high_resolution_clock::now(); // Для замера времени
    auto duration = duration_cast<milliseconds>(stop - start);

    cout << "На основной сетке проведено "<<S<<" итераций"<<endl;
    cout << "На контрольной сетке проведено "<<S_<<" итераций"<<endl;
    cout << "На основной сетке достигнута точность метода eps = " << epsMax <<endl;
    cout << "На контрольной сетке достигнута точность метода eps = " << epsMax_ <<endl;
    cout << "Точность решения задачи = "<<precision<<endl;
    cout << "Достигнута норма невязки "<<current_norm<<endl;
    cout << "Основная сетка: "<<str << endl;
    cout << "Контрольная сетка: "<<str_ << endl;
    cout << "Время просчёта "<<duration.count() << " мс"<<endl;
}

int main() {
    system("chcp 65001");
    cout << "<<МЕТОД СОПРЯЖЁННЫХ ГРАДИЕНТОВ>>\nЩербаков Павел, гр. 3821Б1ПМоп2, Вар. 1"<<endl;
    // Применение метода сопряжённых градиентов к решению
    // Задачи Дирихле для уравнения Пуассона

    int n = 300; // Число разбиений по оси x 1500
    int m = 300; // Число разбиений по оси y 1501

    int Nmax = 4000; // Максимальное число итераций 20 000

    double eps = 0.0000005; // Требуемая погрешность для тестовой задачи

    double eps_met =  0.00000000001; // Целевая точность метода для основной задачи на основной сетке
    double eps_met2 = 0.0000000001; // Целевая точность метода для основной задачи на контрольной сетке

    int Nmax2 = 8000; // Максимальное число итераций на контрольной сетке

    bool test; // Является ли задача основной/тестовой

    cout << "Выберите задачу: 0 - тестовая, 1 - основная"<<endl;
    //cin >> test;
    test = 1;
    if(!test)
        solve_test(n,m,Nmax,eps);
    else
        solve_main(n,m,Nmax,eps_met,Nmax2,eps_met2);
    return 0;
}