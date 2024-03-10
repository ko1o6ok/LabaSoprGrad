#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

// Скалярное произведение
double scalar_product(vector<double> v1, vector<double> v2){
    double S = 0.0; // Сумма
    for (int i = 0; i < v1.size(); ++i)
        S += v1[i] * v2[i];
    return S;
}
// Чебышевская норма
double norm_inf(const vector<double>& v){
    double mx = -100;
    for (auto& t:v)
        if(abs(t)>mx)
            mx = abs(t);

    return mx;
}
// Умножение матрицы на вектор
// A * v
vector<double> mat_mult(const vector<vector<double>>& matrix, const vector<double>& vec){
    vector<double> res;
    res.reserve(matrix.size());
    for (auto& v:matrix) {
        res.push_back(scalar_product(v,vec));
    }
    return res;
}
// Линейная комбинация векторов v1 и v2:
// v = a * v1 + b * v2
vector<double> lin_comb(vector<double> v1,vector<double> v2,double a,double b){
    vector<double> res;

    res.reserve(v1.size());
    for (int i = 0; i < v1.size(); ++i)
            res.push_back(a*v1[i]+b*v2[i]);

    return res;
}
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
    return mu1_main(y);
}
double mu3_main(double x){
    return x-x*x;
}
double mu4_main(double x){

    return mu3_main(x);
}
// Функция u для тестовой задачи
// u* = exp(sin^2(pi * x * y))
double u_test(double x, double y){
    return exp(f_main(x,y));
}
// Функция f для тестовой задачи
double f_test(double x, double y){
    return -0.5 * u_test(x,y) * M_PI * M_PI * (x * x + y * y) * (-1-4 * cos(2 * M_PI * x * y) + cos(4 * M_PI * x * y));
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

// Массив для прохода узлов в правильном порядке
// n - число разбиений по оси x
// m - число разбиений по оси y
vector<pair<int,int>> uzly_region(int n, int m){
    vector<pair<int,int>> res;
    for (int j = 0; j <= m; ++j)
        for (int i = 0; i <= n; ++i)
            res.emplace_back(i,j);
    return res;
}

// Проверка того, что узел находится на границе
// i - номер узла по x
// j - номер узла по y
// n - число разбиений по оси x
// m - число разбиений по оси y
bool on_border(int i, int j, int n, int m){
    return (i == 0) || (i == n) || ( j == 0) || (j == m) ;
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
void fill_zeros(vector<double> v){
    for(auto& t:v)
        t = 0.0;
}
void print_vector(const vector<double>& v){
    for (auto& t:v) {
        cout << t<<" ";
    }
    cout << endl;
}
// Максимальный модуль разности
double accuracy(vector<double> v1, vector<double> v2){
    double diffMax = -100;
    double diffCur;
    for (int i = 0; i < v1.size(); ++i) {
        diffCur = abs(v1[i]-v2[i]);
        if(diffCur > diffMax)
            diffMax = diffCur;
    }
    return diffMax;
}
// Функция для решения системы Ax = b
// Методом сопряжённых градиентов
// Nmax - макс. число итераций
// eps - целевая точность
//vector<double> conjugate_grad(vector<vector<double>> A,vector<double> b, int Nmax,double eps){
//    int n = b.size();// Размерность системы
//    // Метод найдёт точное решение не более, чем за n итераций
//
//    vector<double> x(n); // Текущее решение
//    fill_zeros(x);
//    double alpha; // Текущий коэффициент
//    double beta = 0.0; // Текущий коэффициент
//
//    int S = 0; // Номер текущей итерации
//
//    vector<double> r; // Текущая невязка r_s
//    vector<double> h(n);// Текущий вектор h_s
//    fill_zeros(h);
//    cout << "Численное решение на шаге 0 есть: \n";print_vector(x);
//
//    double r_norm; // Текущая норма невязки
//
//    double epsMax = 0; // Точность на текущей итерации
//    double epsCur = 0; // Вспомогательная точность
//
//    vector<double> x_old; // Предыдущее решение, с которым сравниваем
//    double acc; // Текущая точность
//    while(true) {
//        x_old = x;
//        // Итерация метода по формулам
//
//        r = lin_comb(mat_mult(A,x),b,1.0,-1.0); // Текущая невязка r_s
//        r_norm = norm_inf(r);
//        cout<<"Норма невязки = "<<r_norm<<endl;
//        if(S>0)
//            beta = scalar_product(mat_mult(A,h),r)/ scalar_product(mat_mult(A,h),h);
//        h = lin_comb(h,r,beta,-1.0); // Текущий вектор h_s
//        alpha = -scalar_product(r,h)/ scalar_product(mat_mult(A,h),h);
//        x = lin_comb(x,h,1.0,alpha);
//
//        S++;
//        cout << "Численное решение на шаге "<<S<<" есть: \n";
//
//        print_vector(x);
//        acc = accuracy(x,x_old);
//        cout << "Текущая точность: "<<acc<<endl;
//        if(S>= min(Nmax,n)){
//            cout<<"Выход по числу итераций!"<<endl;
//            return x;
//        }
//
//        if(acc<eps){
//            cout << "Выход по точности!"<<endl;
//            return x;
//        }
//    }
//}
int main() {
    system("chcp 65001");

//    vector<vector<double>> matrix = {{10,1},{1,5}};
//    vector<double> F = {12,11};
//    conjugate_grad(matrix,F,10,0.01);




    // Применение метода сопряжённых градиентов к решению
    // Задачи Дирихле для уравнения Пуассона

    int n = 40; // Число разбиений по оси x
    int m = 40; // Число разбиений по оси y


    double h = (b-a)/n; // Шаг по оси x
    double k = (d-c)/m; // Шаг по оси y



    string str = ""; // Сообщение при окончании работы

    int S = 0; // Количество итераций метода
    int Nmax = 1000; // Максимальное число итераций

    double eps = 0.001; // Требуемая точность для основной/тестовой задачи

    bool test = true; // Является ли задача основной/тестовой

    double epsMax = 0; // Точность на текущей итерации
    double epsCur = 0; // Вспомогательная точность

//    double v_old; // Старое значение в узле
//    double v_new; // Новое значение в узле

    double h2 =  ((n/(b-a)) * (n / (b-a))); // Вспомогательная величина 1/h^2
    double k2 =  ((m/(d-c)) * (m/(d-c))); // Вспомогательная величина 1/k^2

    double a2 = - 2 * (h2 + k2); // Вспомогательная величина -2(1/h^2+1/k^2)

    vector<double> x(n+1); // Абсциссы узлов
    vector<double> y(m+1); // Ординаты узлов

    x[0] = a;
    for (int i = 1; i <= n; ++i)
        x[i] = x[i-1] + h;

    y[0] = c;
    for (int j = 1; j <= m; ++j)
        y[j] = y[j-1] + k;


    // Сгенерируем узлы области
    auto uzly  = uzly_region(n, m);

//    vector<double> inserted(m+1);
//    for(auto& t:inserted)
//        t = 0.0;
//    vector<vector<double>> v(n+1); // Численное решение на s-той итерации
//    for(auto& vect:v)
//        vect = inserted;
//
//    // Учёт гран условий
//    for(auto& [i,j]:uzly)
//        if(on_border(i,j,n,m))
//            v[i][j] = gran_func(x[i],y[j],test);

    // К сожалению, матрицу здесь придётся хранить
    vector<vector<double>> A; // Матрица для системы -Ax = -F
    vector<double> F; // Правая часть - F
    // Чтобы матрица была положительно определена, нам нужно решать -Ax = -F

    // Заполним матрицу
    // -A

    for (auto& [i,j]:uzly)
        if(!on_border(i,j,n,m)){
            // Делаем двойной проход
            vector<double> add;
            double add_F = 0.0;

            for (auto& [i_,j_]:uzly)
            {
                if((i==i_)&&(j==j_)){
                    if(!on_border(i_,j_,n,m))
                        add.push_back(-a2);
                }
                else if((i==i_)&& (abs(j-j_)==1)){
                    if(!on_border(i_,j_,n,m))
                        add.push_back(-k2);
                    if(on_border(i_,j_,n,m))
                        add_F -= gran_func(x[i_],y[j_],i_,j_,n,m)*k2; // Подставляем гранусловия
                }

                else if((j==j_)&& (abs(i-i_)==1)){
                    if(!on_border(i_,j_,n,m))
                        add.push_back(-h2);
                    if(on_border(i_,j_,n,m))
                        add_F -= gran_func(x[i_],y[j_],i_,j_,n,m)*h2; // Подставляем гранусловия
                }

                else
                    if(!on_border(i_,j_,n,m))
                        add.push_back(0.0);
            }
            double fir = x[i];
            double sec = y[j];
            add_F += f_test(x[i],y[j]);

            add_F = -add_F; // Нам нужно -F
            F.push_back(add_F);
            A.push_back(add);
        }

//    cout<<"Получена матрица"<<endl;
//    for(auto& vec:A)
//        print_vector(vec);
//    cout << "Правая часть: ";
//    print_vector(F);

    vector<double> x_s; // Текущий вектор x(s)
    vector<double> h_s; // Текущий вектор h(s)
    x_s.reserve(A.size());
    h_s.reserve(A.size());
    for(auto t:A){
        x_s.push_back(0.0);
        h_s.push_back(0.0);
    }

    // Для простоты берём нулевой начальный вектор

    vector<double> r_s; // Текущая невязка
    double alpha ; // alpha(s)
    double beta = 0.0; // beta(s)

    vector<double> x_s_new; // Новое значение

    vector<double> real_sol; // Реальное решение

    // Начинаем итерации метода сопр. градиентов
    int counter ; // Счётчик для прохода
    double current_norm; // Текущая норма невязки
    while (true){

        if(S>= min(Nmax,(int)A.size())){
            cout << "Норма невязки = "<<current_norm<<endl;
            cout << "Посчитанное решение: ";
            for(auto& tmp:x_s)
                cout << tmp <<" ";
            cout << endl << "Реальное решение: ";
            for(auto& temp:real_sol)
                cout << temp <<" ";
            cout << "\n-----------------------------------------"<<endl;
            str = "Выход по числу итераций!";
            break;
        }
        real_sol.clear();
        epsMax = 0;
        // Проведём итерацию методом сопряжённых градиентов.
        // Начинаем с подсчёта текущей невязки
        // r_s = A x_s-F
        r_s = lin_comb(mat_mult(A,x_s),F,1.0,-1.0); // Текущая невязка
        current_norm = norm_inf(r_s); // Текущая норма невязки
        if(S>=1)
            beta = scalar_product(mat_mult(A,h_s),r_s)/scalar_product(mat_mult(A,h_s),h_s); // Так считаем коэффициент

        h_s = lin_comb(h_s,r_s,beta,-1.0);

        alpha = -scalar_product(r_s,h_s)/ scalar_product(mat_mult(A,h_s),h_s);

        x_s_new = lin_comb(x_s,h_s,1.0,alpha); // Получаем новый вектор

        // Теперь смотрим его погрешность/точность
        //
        counter = 0;
        for (auto& [i,j]:uzly)
            if(!on_border(i,j,n,m)){

//                v_old = v[i][j];
//                // Здесь будет метод сопряжённых градиентов
//                v_new = 0;


                // В случае тестовой задачи смотрим погрешность в норме бесконечность
                if(test){
                    //epsCur = abs(u_test(x[i],y[j])-x_s_new[counter]);
                    real_sol.push_back(u_test(x[i],y[j]));
                }

                else
                    // В случае основной задачи смотрим точность на s-той итерации
                    epsCur = abs(x_s[counter]-x_s_new[counter]);

                // Используем норму бесконечность

//                if(epsCur > epsMax)
//                    epsMax = epsCur;

//                v[i][j] = v_new;
            }
        x_s = x_s_new;
        ++S;
        epsMax = accuracy(x_s,real_sol);
        if(S==1 || S==2){

            cout << "Текущая точность "<<epsMax<<endl;
        }

        if(epsMax < eps){
            str = "Выход по точности!";
            break;
        }
        if(S==1||S==2){
            cout << "Норма невязки = "<<current_norm<<endl;
            cout << "Посчитанное решение: ";
            for(auto& tmp:x_s)
                cout << tmp <<" ";
            cout << endl << "Реальное решение: ";
            for(auto& temp:real_sol)
                cout << temp <<" ";
            cout << "\n-----------------------------------------"<<endl;
        }

    }


    cout << "Проведено "<<S<<" итераций"<<endl;
    cout << "Достигнута точность eps = " << epsMax <<endl;
    cout << str << endl;

    return 0;
}
