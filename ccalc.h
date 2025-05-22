/**
 *
 * Классическая конечно-разностная схема
 *
**/
#ifndef CCALC_H
#define CCALC_H

#include <QObject>
#include <complex>
#include <vector>

#define Nlines 4
#define NU 10
#define MaxThreadCount 10

using namespace std;

typedef float T;
/// размеры реальной системы
//const int Nx = 40;
//const T Lx = 0.204;//[nm]  => dx = 0.0051
//const int Ny = 1000;//
//const T Ly = 3200.0;//[nm] если dx -> dy => Ny = 627'451 !!! => конечно-разностная схема не вариант !
/// размеры моделируемой системы: Nx*Ny = 40.000 точек

const T Ly2Lx = 2.0;
const int Np = 100;
const T L = 2.0;

//const T Ly2Lx = 2.0;
//const int Np = 100;
//const T L = 1.0;

const int Nx = Np*Ly2Lx;// = 400
const int Ny = Np/Ly2Lx;// = 100
const T Lx = L*Ly2Lx;
const T Ly = L/Ly2Lx;
const T dx = Lx/(Nx-1);
const T dy = Ly/(Ny-1);
const T p2 = dx*dx/dy/dy;
const T dx2 = dx*dx;
const int N = 20;// число Фурье-компонент в Фурье-образе
/// константы
const T kappa = 58.0;//58.0;
const T sigma = 1.42;
const complex<T> im = {0.0f, 1.0f};
const complex<T> czero = {0.0f, 0.0f};
///
/// \brief структура модели
///
struct SModel
{
    int counter = 0;
    const int counter_max = 5000;
    T dt = .01f;
    T t;
    T B = .0f;
    T J = .0f;
    complex<T> px = {0.0f, 0.5f*dx*B};
    complex<T> py = {0.0f, 0.5f*dy*B};
    T C[5] = { dt/kappa/kappa/dx/dx, dt/kappa/kappa/dy/dy,
                    .5f*B*dt/kappa/dx, .5f*B*dt/kappa/dy, .25f*dt*B*B };
    const string Cstr[5] = {"dt/kappa^2/dx^2", "dt/kappa^2/dy^2",
                           "i*.5*dt*B/kappa/dx", "i*.5*dt*B/kappa/dy", ".25*dt*B^2"};

    std::complex<T>** psi = nullptr;
    std::complex<T>** npsi = nullptr;
    std::complex<T>** Ux = nullptr;
    std::complex<T>** Uy = nullptr;
    std::complex<T>** Vx = nullptr;
    std::complex<T>** Vy = nullptr;
    T** divJs = nullptr;
    T** phi = nullptr;
    T Uinstant;
    T Uavg, Sum4Uavg, dUavg;
    T container[NU];
    const float ParamCr = 0.5e-7;
    int countU;
    /// вспомогательные величины
    T* x = nullptr;
    T* y = nullptr;
    T* ky = nullptr;        /**< набор значений волнового вектора */
    T** Fphi = nullptr;     /**< Фурье-образ потенциала */
    T** FdivJs = nullptr;   /**< Фурье-образ дивергенции сверхпровод. тока */
    T** Cos = nullptr;      /**< матрица Фурье-мод (собств.функций) по которым выполняется разложение */
    T delta_psi;

    bool boolElectricField = false;
};

class CCalc : public QObject
{
    Q_OBJECT
public:
    explicit CCalc(QObject *parent = nullptr);
    ~CCalc();
    void init();
//    int getCounter()
//        { return counter; }
    T getTime()
        { return m.t; }
    void setBoundCond(bool periodic_bound_cond)
        { this->periodic_bound_cond = periodic_bound_cond; }
    void setBoolConinue(bool boolContinue)
        { this->boolContinue = boolContinue; }
    void updateMagneticField(T B);
    void updateInjectionСurrent(T J);
    void updateTimeStep(T dt);

    bool isTerminated;
    bool isFinished;
    int divider = 100;
    SModel m;
    float* data;

signals:
    void dataReady();
    void resetChart();
    void updateChart();
    void updateScatter(int, float);
    void nextValue(double);
    void processFinished();

    void NumerSolReady(bool);
    void AnalytSolReady(bool);

public slots:
    void process(); /**< без link variables */
    void process1();/**< c link variables однопоточные вычисления */
    void process2();/**< c link variables многопоточные вычисления */

    void getNumerSol(bool);
    void getAnalytSol(bool);

private:
    QList<double> listJ;
    bool boolListJ = true;
    bool periodic_bound_cond = false;
    bool boolContinue = false;

    void getListJ();
    void get_parallel_psi();
    void get_parallel_divJs();
    void get_parallel_FdivJs();
    void get_parallel_Fphi();
    void get_parallel_Potential();

    void get_divJs();
    void get_FdivJs();
    void get_Fphi();
    void get_Potential();

    T U(float x, float y);
    T V(float x, float y);
    T U2(float x, float y);
    T V2(float x, float y);
};

#endif // CCALC_H
