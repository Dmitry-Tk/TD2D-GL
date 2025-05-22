#include "ccalc.h"
#include "unistd.h"
#include <iostream>
#include <chrono>
#include <random>
#include <execution>
#include <QThreadPool>
#include <QDateTime>
#include <fstream>
#include <iomanip>//setw(..)

using namespace std;

/// генератор случайных чиселл
random_device rd;
mt19937 engine( rd() );// каждый раз случайная
//mt19937 engine( 5 ); // фиксированная последовательность
uniform_real_distribution<float> real_distribution(0.0f, 1.0f);
auto Rand = bind(real_distribution, engine);
/// для временных измерений
QDateTime timer;
quint64 start_time;
///
/// \brief класс СTask1 для параллельного вычисления параметра порядка
///
class СTask1 : public QRunnable
{
public:
    СTask1(int k)
    {
        this->k = k;
    }
    ~СTask1()
    {
    }
    static void setModel(SModel* model)
    {
        m = model;
        Sum[0] = czero.real();
        Sum[1] = czero.imag();
    }
    static int Counter;
    static float Sum[2];

private:
    void run() override
    {
        /// граничные условия для параметра порядка:
        ///   без электродов: d(U_mu*Psi)/dmu |_{\vec{n}} = 0
        ///   с электродами, но без учета потенциала: conj(U_mu)*d(U_mu*Psi)/dmu |_{\vec{n}} = i*Jx/conj(Psi)|_{\vec{n}}
        if (k == 0)
        {   /// y = -.5*Ly : нижняя граница
            for (int n = 0; n < Nx; n++)
            {
                m->npsi[n][k] = m->Uy[n][k+1]/m->Uy[n][k]*m->psi[n][k+1];
                m->delta_psi += abs(m->npsi[n][k] - m->psi[n][k])*abs(m->npsi[n][k] - m->psi[n][k]);
            }
        }
        else if (k == Ny-1)
        {   /// y = +.5*Ly : верхняя граница
            for (int n = 0; n < Nx; n++)
            {
                m->npsi[n][k] = m->Uy[n][k-1]/m->Uy[n][k]*m->psi[n][k-1];
                m->delta_psi += abs(m->npsi[n][k] - m->psi[n][k])*abs(m->npsi[n][k] - m->psi[n][k]);
            }
        }
        else
        {   /// -.5*Ly < y < .5*Ly
            /// x = -Lx/2 : левая граница
            /// старый вариант: учёт электрического поля через инжекцию тока вдоль оси x
            if (m->boolElectricField)
            {
                /// I
//                m->npsi[0][k] = m->Ux[1][k]/m->Ux[0][k]*m->psi[1][k] - im*m->J0/sigma/conj(m->psi[0][k])*dx;
                /// II
                m->npsi[0][k] = m->npsi[Nx-1][k] = czero;
            }
            else
                m->npsi[0][k] = m->Ux[1][k]/m->Ux[0][k]*m->psi[1][k];
            m->delta_psi += abs(m->npsi[0][k] - m->psi[0][k])*abs(m->npsi[0][k] - m->psi[0][k]);
            complex<T> dpsi;
            for (int n = 1; n < Nx-1; n++)
            {
                float mod2 = m->psi[n][k].real()*m->psi[n][k].real() + m->psi[n][k].imag()*m->psi[n][k].imag();
                dpsi = m->C[0]*(conj(m->Ux[n][k])*(m->Ux[n+1][k]*m->psi[n+1][k] + m->Ux[n-1][k]*m->psi[n-1][k]) - m->psi[n][k] - m->psi[n][k])
                        + m->C[1]*(conj(m->Uy[n][k])*(m->Uy[n][k+1]*m->psi[n][k+1] + m->Uy[n][k-1]*m->psi[n][k-1]) - m->psi[n][k] - m->psi[n][k])
                        + m->dt*(1.0f - mod2)*m->psi[n][k];
//                        - im*m->dt*kappa*m->phi[n][k]*m->psi[n][k];
                m->npsi[n][k] = (m->psi[n][k] + dpsi)*exp(-im*m->phi[n][k]*m->dt);
                m->delta_psi += abs(m->npsi[n][k] - m->psi[n][k])*abs(m->npsi[n][k] - m->psi[n][k]);
            }
            /// x = +Lx/2 : правая граница
            /// старый вариант: учёт электрического поля через инжекцию тока вдоль оси x
            if (m->boolElectricField)
            {
//                m->npsi[Nx-1][k] = m->Ux[Nx-2][k]/m->Ux[Nx-1][k]*m->psi[Nx-2][k] + im*m->J0/sigma/conj(m->psi[Nx-1][k])*dx;
            }
            else
                m->npsi[Nx-1][k] = m->Ux[Nx-2][k]/m->Ux[Nx-1][k]*m->psi[Nx-2][k];
            m->delta_psi += abs(m->npsi[Nx-1][k] - m->psi[Nx-1][k])*abs(m->npsi[Nx-1][k] - m->psi[Nx-1][k]);
//            if (m->boolElectricField)
//                for (int k = 0; k < Ny; k++)
//                {
//                    Sum[0] += m->npsi[Nx-1][k].real();
//                    Sum[1] += m->npsi[Nx-1][k].imag();
//                }
        }
//        cout << "  " << counter;
//        if (Counter == 0)
        {
            /// контроль
//            cout << "\nThreadPool has been  finished... k :" << k //<< endl
//                 << "\tCounter: " << Counter << endl;
        }
        Counter--;
    }
    int k;/**< индекс строки обрабатываемой отдельным потоком */
    static SModel* m;
};
SModel* СTask1::m = nullptr;
int СTask1::Counter = 0;
float СTask1::Sum[2] = {0.0f, 0.0f};
///
/// \brief класс СTask2 для параллельного вычисления div{Js}*dx^2
///
class СTask2 : public QRunnable
{
public:
    СTask2(int k)
    {
        this->k = k;
    }
    ~СTask2()
    {
    }
    static void setModel(SModel* model)
    {
        m = model;
    }
    static int Counter;
private:
    void run() override
    {
        complex<T> dJsdx, dJsdy;
        if (k == 0)
        {   /// y = -.5*Ly: нижняя граница
            for (int n = 1; n < Nx-1; n++)
            {
                dJsdx = conj(m->Ux[n][k]*m->psi[n][k])*
                        (m->Ux[n+1][k]*m->psi[n+1][k] + m->Ux[n-1][k]*m->psi[n-1][k]);
                dJsdy = conj(m->Uy[n][k]*m->psi[n][k])*
                    (-m->Uy[n][k+3]*m->psi[n][k+3] + 4.f*m->Uy[n][k+2]*m->psi[n][k+2]- 5.f*m->Uy[n][k+1]*m->psi[n][k+1]);
                m->divJs[n][k] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
            }
        }
        else if (k == Ny-1)
        {   /// y = +.5*Ly: верхняя граница
            for (int n = 1; n < Nx-1; n++)
            {
                dJsdx = conj(m->Ux[n][k]*m->psi[n][k])*
                        (m->Ux[n+1][k]*m->psi[n+1][k] + m->Ux[n-1][k]*m->psi[n-1][k]);
                dJsdy = conj(m->Uy[n][k])*conj(m->psi[n][k])*
                    (-m->Uy[n][k-3]*m->psi[n][k-3] + 4.f*m->Uy[n][k-2]*m->psi[n][k-2]- 5.f*m->Uy[n][k-1]*m->psi[n][k-1]);
                m->divJs[n][k] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
            }
        }
        else
        {
            /// x = -.5*Lx: левая граница
            dJsdx = conj(m->Ux[0][k]*m->psi[0][k])*
                (-m->Ux[3][k]*m->psi[3][k] + 4.f*m->Ux[2][k]*m->psi[2][k]- 5.f*m->Ux[1][k]*m->psi[1][k]);
            dJsdy = conj(m->Uy[0][k]*m->psi[0][k])*
                    (m->Uy[0][k+1]*m->psi[0][k+1] + m->Uy[0][k-1]*m->psi[0][k-1]);
            m->divJs[0][k] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
            /// -.5*Lx < x < .5*Lx
            for (int n = 1; n < Nx-1; n++)
            {
                dJsdx = conj(m->Ux[n][k]*m->psi[n][k])*
                        (m->Ux[n+1][k]*m->psi[n+1][k] + m->Ux[n-1][k]*m->psi[n-1][k]);
                dJsdy = conj(m->Uy[n][k])*conj(m->psi[n][k])*
                        (m->Uy[n][k+1]*m->psi[n][k+1] + m->Uy[n][k-1]*m->psi[n][k-1]);
                m->divJs[n][k] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
            }
            /// x = +.5*Lx: правая граница
            dJsdx = conj(m->Ux[Nx-1][k]*m->psi[Nx-1][k])*
                (-m->Ux[Nx-4][k]*m->psi[Nx-4][k] + 4.f*m->Ux[Nx-3][k]*m->psi[Nx-3][k] - 5.f*m->Ux[Nx-2][k]*m->psi[Nx-2][k]);
            dJsdy = conj(m->Uy[Nx-1][k]*m->psi[Nx-1][k])*
                    (m->Uy[Nx-1][k+1]*m->psi[Nx-1][k+1] + m->Uy[Nx-1][k-1]*m->psi[Nx-1][k-1]);
            m->divJs[Nx-1][k] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
        }
        Counter--;
    }
    int k;/**< индекс строки обрабатываемой отдельным потоком */
    static SModel* m;
};
SModel* СTask2::m = nullptr;
int СTask2::Counter = 0;
///
/// \brief класс СTask3 для параллельного вычисления Фурье-коэф-тов дивергенции сверхпровод. тока
///
class СTask3 : public QRunnable
{
public:
    СTask3(int n)
    {
        this->n = n;
    }
    ~СTask3()
    {
    }
    static void setModel(SModel* model)
    {
        m = model;
    }
    static int Counter;
private:
    void run() override
    {
        for (int i = 0; i < Nx; i++)
        {
            m->FdivJs[i][n] = .0f;
            for (int j = 1; j < Ny; j++)
                m->FdivJs[i][n] += (m->divJs[i][j]*m->Cos[n][j] + m->divJs[i][j-1]*m->Cos[n][j-1]);
            m->FdivJs[i][n] /= (Ny-1);
        }
        Counter--;
    }
    int n;/**< индекс моды, для которой вычисляются Фурье-коэф-ты  div{Js} для всех x-координат */
    static SModel* m;
};
SModel* СTask3::m = nullptr;
int СTask3::Counter = 0;
///
/// \brief класс СTask4 для параллельного вычисления Фурье-коэф-тов потенциала
///
class СTask4 : public QRunnable
{
public:
    СTask4(int n)
    {
        this->n = n;
    }
    ~СTask4()
    {
    }
    static void setModel(SModel* model)
    {
        m = model;
    }
    static int Counter;
private:
    void run() override
    {
        static const float c1 = 9.0/8.0,  c2 = 3.0/8.0;
        static const float c3 = 1.0/24.0, c4 = 11.0/12.0;
        static const float c5 = 1.0/3.0,  c6 = 5.0/24.0;
        T f[Nx];
        /// Фурье-коэффициент для нулевой моды  (ky[n=0] = 0)
        /// решаем f'(x,ky[0]) = div{Js}(x,ky[0])*dx^2 = D(x), где f(x)=phi'(x,ky[0])
        f[0] = m->J/sigma;
        f[Nx-1] = -m->J/sigma;// ??? согласно формуле д.б. "-", согласно физич. соображениям - "+"
        /// 1. учитываем оба граничных условия
//        for (int i = 2; i < Nx-1; i=i+2)
//        {
//            f[i] = f[i-2] + 2*m->FdivJs[i-1][0]/dx;// + O(dx^2)
//            f[Nx-1-i] = f[Nx+1-i] - 2*m->FdivJs[Nx-i][0]/dx;// + O(dx^2)
//        }
        /// 2. учитываем только левое граничное условие
        f[1] = f[0] + (c1*m->FdivJs[1][0] + c2*m->FdivJs[0][0])*dx;// + O(dx^2)
        for (int i = 2; i < Nx-1; i++)
            f[i] = f[i-1] + (c3*m->FdivJs[i+1][0] + c4*m->FdivJs[i][0] + c3*m->FdivJs[i-1][0])/dx;// + O(dx^2)
        f[Nx-1] = f[Nx-2] + (c5*m->FdivJs[Nx-1][0] + c6*m->FdivJs[Nx-2][0] - c3*m->FdivJs[Nx-3][0])/dx;// + O(dx^2)
        /// решаем phi'(x,ky[0]) = f(x,ky[0])
        m->Fphi[0][0] = .0f;
        m->Fphi[1][0] = m->Fphi[0][0] + (c1*f[1] + c2*f[0])*dx;// + O(dx^2)
        for (int i = 2; i < Nx-1; i++)
            m->Fphi[i][0] = m->Fphi[i-1][0] + (c3*f[i+1] + c4*f[i] + c3*f[i-1])*dx;// + O(dx^2)
        m->Fphi[Nx-1][0] = m->Fphi[Nx-2][0] + (c5*f[Nx-1] + c6*f[Nx-2] - c3*f[Nx-3])*dx;// + O(dx^2)
        /// Фурье-коэфф-ты для остальных мод: ky[n] > 0 -> решение СЛАУ: A*Fphi = FdivJs*dx^2 методом прогонки
        /// Fphi, FdivJs*dx^2 - вектора Фурье-коэфф-тов потенциала и дивергенции сверхпр. тока
        /// A - тридиагональная матрица; a - главная диагональ; b(=1) и c(=1) - над и поддиагонали соот-но
        for (int n = 1; n < N; n++)
        {
            T D[Nx], a[Nx], b[Nx-1], a_const = -(2.0f + m->ky[n]*m->ky[n]*dx2);
            for (int i = 0; i < Nx; i++)
            {
                D[i] = m->FdivJs[i][n];
                a[i] = a_const;
                if (i < Nx-1)
                    b[i] = 1.f;
            }
            a[0] += 1.0f;   a[Nx-1] += 1.0f;
            /// первый прогон - исключение нижней поддиагонали с учётом, что  c[i] = b[i] = 1
            /// и приведение главного диагонального элемента к 1.0
            /// делением текущей строки на a[i] для устранения больших значений
            for (int i = 1; i < Nx; i++)
            {
                a[i] = a[i-1] - 1.f/a[i];
                if (i < Nx-1)
                    b[i] = a[i-1]/a[i];
                D[i] = b[i]*D[i] - D[i-1]/a[i];
                a[i] = 1.f;
            }
            /// второй прогон - вычисление корней СЛАУ c учётом, что a[i] = 1
            m->Fphi[Nx-1][n] = D[Nx-1];
            for (int i = Nx-2; i > -1; i--)
                m->Fphi[i][n] = D[i] - b[i]*m->Fphi[i+1][n];
        }
        Counter--;
    }
    int n;/**< индекс Фурье-коэф-тов потенциала вычисляемых для всех x-координат */
    static SModel* m;
};
SModel* СTask4::m = nullptr;
int СTask4::Counter = 0;
///
/// \brief класс СTask5 для параллельного вычисления потенциала с помощью его Фурье-разложения
///
class СTask5 : public QRunnable
{
public:
    СTask5(int n)
    {
        this->n = n;
    }
    ~СTask5()
    {
    }
    static void setModel(SModel* model)
    {
        m = model;
    }
    static int Counter;
private:
    void run() override
    {
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
                m->phi[i][j] = (n > 0) ? m->phi[i][j] + m->Fphi[i][n]*m->Cos[n][j] : m->Fphi[i][0];
        Counter--;
    }
    int n;/**< индекс моды добавляемой отдельным потоком в потенциал во всех точках */
    static SModel* m;
};
SModel* СTask5::m = nullptr;
int СTask5::Counter = 0;

CCalc::CCalc(QObject* /*parent*/)
{
    m.psi = new complex<T>*[Nx];
    m.npsi = new complex<T>*[Nx];
    m.Ux = new complex<T>*[Nx];
    m.Uy = new complex<T>*[Nx];
    m.Vx = new complex<T>*[Nx];
    m.Vy = new complex<T>*[Nx];
    m.phi = new T*[Nx];
    m.divJs = new T*[Nx];
    m.Fphi = new T*[Nx];
    m.FdivJs = new T*[Nx];
    m.x = new T[Nx];
    m.y = new T[Ny];
    for (int i = 0; i < Nx; i++)
    {
        m.psi[i] = new complex<T>[Ny];
        m.npsi[i] = new complex<T>[Ny];
        m.Ux[i] = new complex<T>[Ny];
        m.Uy[i] = new complex<T>[Ny];
        m.Vx[i] = new complex<T>[Ny];
        m.Vy[i] = new complex<T>[Ny];
        m.phi[i] = new T[Ny];
        m.divJs[i] = new T[Ny];
        m.Fphi[i] = new T[N];
        m.FdivJs[i] = new T[N];
    }
    m.ky = new T[N];
    m.Cos = new T*[N];
    for (int n = 0; n < N; n++)
        m.Cos[n] = new T[Ny];
    data = new float[Nlines];
}

CCalc::~CCalc()
{
    for (int i = 0; i < Nx; i++)
    {
        delete [] m.psi[i];
        delete [] m.npsi[i];
        delete [] m.Ux[i];
        delete [] m.Uy[i];
        delete [] m.Vx[i];
        delete [] m.Vy[i];
        delete [] m.phi[i];
        delete [] m.divJs[i];
        delete [] m.Fphi[i];
        delete [] m.FdivJs[i];
    }
    for (int n = 0; n < N; n++)
        delete [] m.Cos[n];
    delete [] m.psi;
    delete [] m.npsi;
    delete [] m.Ux;
    delete [] m.Uy;
    delete [] m.Vx;
    delete [] m.Vy;
    delete [] m.phi;
    delete [] m.Fphi;
    delete [] m.divJs;
    delete [] m.FdivJs;
    delete [] m.Cos;
    delete [] m.ky;
    delete [] m.x;
    delete [] m.y;
    delete data;
}

void CCalc::init()
{
    cout << "Размеры моделируемой системы - Nx*Ny = " << Nx*Ny << " точек" << endl <<endl;
    float sum = 0.0f;
    complex<T> sum_psi = czero;
    /// случайные данные для Psi
    float Ax, Ay;//, pAx = 0, pAy = 0;
    for (int i = 0; i < Nx; i++)
    {
        m.x[i] = (i > 0) ? m.x[i-1] + dx : -.5*Lx;
        for (int j = 0; j < Ny; j++)
        {
            m.y[j] = (j > 0) ? m.y[j-1] + dy : -.5*Ly;
            float r = Rand(), alpha = 2*M_PI*Rand();
            sum += r + alpha;
            m.psi[i][j] = complex<T>( r*cos(alpha), r*sin(alpha) );
            sum_psi += m.psi[i][j];
            Ax =-0.5f*m.B*m.y[j];
            Ay = 0.5f*m.B*m.x[i];
            /// I.
            /// конечные формулы для link variables, как в статье
            /// проблема -> все вихри постепенно исчезают
//            Ux[i][j] = exp(-im*kappa*dx*Ax);
//            Uy[i][j] = exp(-im*kappa*dy*Ay);
            /// II.
            /// конечные приближенные формулы для link variables в соответствии с определением,
//            Ux[i][j] = (i > 0) ? Ux[i-1][j]*exp(-im*kappa*dx*.5f*(Ax+pAx)) : complex<T>(1.0,0.0);
//            Uy[i][j] = (j > 0) ? Uy[i][j-1]*exp(-im*kappa*dy*.5f*(Ay+pAy)) : complex<T>(1.0,0.0);
            /// III.
            /// точные аналитические выражения согласно определению link variables
            m.Ux[i][j] = exp(-im*kappa*Ax*m.x[i]);
            m.Uy[i][j] = exp(-im*kappa*Ay*m.y[j]);
            /// !!!
//            m.Vx[i][j] = exp(-im*kappa*Ax*m.x[i]);
//            m.Vy[i][j] = exp(-im*kappa*Ay*m.y[j]);
            m.phi[i][j] = m.divJs[i][j] = .0f;
//            pAy = Ay;
        }
//        pAx = Ax;
    }
    for (int n = 0; n < N; n++)
    {
        m.ky[n] = n*M_PI/Ly;
        for (int i = 0; i < Nx; i++)
            m.FdivJs[i][n] = m.Fphi[i][n] = .0f;
        for (int j = 0; j < Ny; j++)
            m.Cos[n][j] = (n > 0) ? cos(m.ky[n]*(m.y[j] - .5f*Ly)) : .5f;
    }
    for (int n = 0; n < Nlines; n++)
        data[n] = 0.0f;
    for (int i = 0; i < NU; i++)
        m.container[i] = 0.0f;

    cout << "Sum of all Rand: " << sum << endl
         << "Sum of all Psi: " << sum_psi << endl;
    /// отображаем параметры модели
    cout << "dt = " << m.dt << endl
         << "dx = " << dx << "\tdy = " << dy << endl
         << "B  = " << m.B << endl
         << "px = " << m.px << "\tpy = " << m.py << endl;
    for (ulong i = 0; i < sizeof(m.Cstr)/sizeof(m.Cstr[0]); i++)
    {
         cout << "C[" << i << "] = " << m.Cstr[i]
                << " = "  << m.C[i] << endl;
    }

    emit resetChart();
}

void CCalc::updateMagneticField(T B)
{
    m.B = B;
    m.C[2] = 0.5f*B*m.dt/kappa/dx;
    m.C[3] = 0.5f*B*m.dt/kappa/dy;
    m.C[4] = 0.25f*m.dt*B*B;
    m.px = {0.0f, 0.5f*dx*B};
    m.py = {0.0f, 0.5f*dy*B};
    if (m.boolElectricField)
    {
        m.countU = 0;
        m.Sum4Uavg = .0f;
        for (int i = 0; i < NU; i++)
            m.container[i] = 1.0f;
    }
}

void CCalc::updateInjectionСurrent(T J)
{
    m.J = J;
    if (fabs(m.J) < 0.001f)
        m.boolElectricField = false;
    else
    {
        m.boolElectricField = true;
        m.countU = 0;
        m.Sum4Uavg = .0f;
        for (int i = 0; i < NU; i++)
            m.container[i] = 1.0f;
    }
    cout << "J : " << m.J << "\telectric field: "
         << ((m.boolElectricField) ? " true " : " false ") << endl;
}

void CCalc::updateTimeStep(T dt)
{
    m.dt = dt;
    m.C[0] = dt/kappa/kappa/dx/dx;
    m.C[1] = dt/kappa/kappa/dy/dy;
    m.C[2] = 0.5f*m.B*dt/kappa/dx;
    m.C[3] = 0.5f*m.B*dt/kappa/dy;
    m.C[4] = 0.25f*dt*m.B*m.B;
}
void CCalc::getListJ()
{
    listJ.clear();
    if (boolListJ)
    {
        /// подготовка списка заданий моделирования -> заполняем список значений параметра
        for (int i = 1; i < 10; i++)
            listJ.append(0.001f*i);
        for (int i = 1; i < 10; i++)
            listJ.append(0.01f*i);
        for (int i = 1; i < 2; i++)
            listJ.append(0.1f*i);
//        float J = 1.0f;
//        while (J <= 2.0f)
//        {
//            listJ.append(J);
//            J += 0.25f;
//        }
    }
    /// контроль
    cout << "listJ = (";
    for (QList<double>::Iterator it = listJ.begin(); it!=listJ.end(); ++it)
        cout << *it << ((it != listJ.end()-1) ? ", " : "");
    cout << ")" << endl;
}
///
/// \brief CCalc::process - расчёт в простейшем наивном варианте
///
void CCalc::process()
{
    start_time = timer.currentMSecsSinceEpoch();
    cout << "start the process without the link variables" << endl;
    isTerminated = isFinished = false;
    if (!boolContinue)
    {
        m.delta_psi = 1.0;// если используется для прерывания расчетов
        m.counter = 0;
        m.t = 0.0;
    }
    ///
    int np, nm, kp, km, bc;
    ///
    bc = (periodic_bound_cond == true) ? 0 : 1;// периодические г.у.: bc=0 иначе bc=1
    while (!isTerminated && !isFinished && !isnan(m.delta_psi) && !isinf(m.delta_psi))// && delta_psi > 1e-6)
    {
        m.delta_psi = 0.0;
        /// обновляем данные внутри массива
//        if (periodic_bound_cond == false)
//            x += dx;
        for (int n = bc; n < Nx-bc; n++)
        {
            if (periodic_bound_cond)
            {
                np = (n+1 == Nx) ? 0 : n+1;
                nm = (n-1 == -1) ? Nx-1 : n-1;
            }
            else
            {
                np = n+1; nm = n-1;
            }
//            if (periodic_bound_cond == false)
//                y += dy;
            float Ey = 0.0f;
            for (int k = bc; k < Ny-bc; k++)
            {
                if (periodic_bound_cond)
                {
                    kp = (k+1 == Ny) ? 0 : k+1;
                    km = (k-1 == -1) ? Ny-1 : k-1;
                }
                else
                {
                    kp = k+1; km = k-1;
                }
                Ey = (m.counter < 1000) ? Ey + 0.00005 : Ey;
                m.phi[n][k] = -Ey*m.y[k];
                float mod2 = m.psi[n][k].real()*m.psi[n][k].real() + m.psi[n][k].imag()*m.psi[n][k].imag();
                m.npsi[n][k] = m.psi[n][k] + m.C[0]*(m.psi[np][k] - m.psi[n][k] - m.psi[n][k] + m.psi[nm][k])
                        + m.C[1]*(m.psi[n][kp] - m.psi[n][k] - m.psi[n][k] + m.psi[n][km])
                        + m.C[2]*im*m.y[k]*(m.psi[np][k] - m.psi[nm][k]) - m.C[3]*im*m.x[n]*(m.psi[n][kp] - m.psi[n][km])
                        - m.C[4]*(m.x[n]*m.x[n] + m.y[k]*m.y[k])*m.psi[n][k]
                        + m.dt*(1.0f - mod2)*m.psi[n][k]
                        - im*m.dt*kappa*m.phi[n][k]*m.psi[n][k];

                m.delta_psi += abs(m.npsi[n][k] - m.psi[n][k])*abs(m.npsi[n][k] - m.psi[n][k]);
            }
        }
        if (periodic_bound_cond == false)
        {
        /// граничные условия
        /// I. согласно алгоритму
        for (int k = bc; k < Ny-bc; k++)
        {
            // (i)
            m.npsi[0][k] = m.psi[1][k]/(1.0f + m.px*m.y[k]);        // bound: x = -Lx/2
            m.npsi[Nx-1][k] = m.psi[Nx-2][k]/(1.0f + m.px*m.y[k]);  // bound: x = +Lx/2
            // (ii)
//            npsi[0][k] = (1.0f - px*y)*psi[1][k];        // bound: x = -Lx/2
//            npsi[Nx-1][k] = (1.0f - px*y)*psi[Nx-2][k];  // bound: x = +Lx/2
        }
        for (int n = bc; n < Nx-bc; n++)
        {
            // (i)
            m.npsi[n][0] = m.psi[n][1]/(1.0f - m.py*m.x[n]);        // bound: y = -Ly/2
            m.npsi[n][Ny-1] = m.psi[n][Ny-2]/(1.0f - m.py*m.x[n]);  // bound: y = +Ly/2
            // (ii)
//            npsi[n][0] = (1.0f + py*x)*psi[n][1];        // bound: y = -Ly/2
//            npsi[n][Ny-1] = (1.0f + py*x)*psi[n][Ny-2];  // bound: y = +Ly/2
        }
        m.npsi[0][0] = 0.5f*(m.npsi[0][1] + m.npsi[1][0]);
        m.npsi[0][Ny-1] = 0.5f*(m.npsi[0][Ny-2] + m.npsi[1][Ny-1]);
        m.npsi[Nx-1][0] = 0.5f*(m.npsi[Nx-2][0] + m.npsi[Nx-1][1]);
        m.npsi[Nx-1][Ny-1] = 0.5f*(m.npsi[Nx-1][Ny-2] + m.npsi[Nx-2][Ny-1]);
        /// II. градиент на границе нулевой
//        for (int n = 0; n < Nx; n++)
//        {
//            npsi[n][0] = psi[n][1];       // dΨ/dy|_{y = 0} = 0
//            npsi[n][Ny-1] = psi[n][Ny-2]; // dΨ/dy|_{y = Ly} = 0
//        }
//        for (int k = 0; k < Ny; k++)
//        {
//            npsi[0][k] = psi[1][k];       // dΨ/dx|_{x = 0} = 0
//            npsi[Nx-1][k] = psi[Nx-2][k]; // dΨ/dx|_{x = Lx} = 0
//        }
        }
        m.t += m.dt;
//        cout << "iteration: " << m.counter << "\tdelta_psi: " << delta_psi << endl;
        if (m.counter % divider == 0)
        {
            cout << "iteration: " << m.counter << "\tdelta_psi: " << m.delta_psi << endl;
//            for (int n = 0; n < Nx; n++)
//            {
//                cout << "\n";
//                for (int k= 0; k < Ny; k++)
//                    cout << psi[n][k] << "\t";
//            }
//            cin.get();

            emit dataReady();
        }
        /// перенос данных psi <- npsi
        for (int n = 0; n < Nx; n++)
            for (int k = 0; k < Ny; k++)
                m.psi[n][k] = m.npsi[n][k];
        m.counter++;
        if (m.counter == m.counter_max)
            isTerminated = true;
    }
    cout << "iteration: " << m.counter << "\tdelta_psi: " << m.delta_psi << endl;
    cout << "time of calculation: "
         << (timer.currentMSecsSinceEpoch() - start_time)/1000 << "[s]" << endl;
    ///
    isTerminated = true;
    if (isFinished)
        emit processFinished();
}
///
/// используются link variables
/// iterations: 10040 -> time of calculation: 236[s] =~42[it/sec] (1 thread)
void CCalc::process1()
{
    start_time = timer.currentMSecsSinceEpoch();
    cout << "start the process with the link variables" << endl;
    isTerminated = isFinished = false;
    if (!boolContinue)
    {
        m.delta_psi = 1.0;// если используется для прерывания расчетов
        m.counter = 0;
        m.t = 0.0;
    }
    ///
    while (!isTerminated && !isFinished && !isnan(m.delta_psi) && !isinf(m.delta_psi))// && delta_psi > 1e-6)
    {
        m.delta_psi = 0.0;
        /// обновляем данные внутри массива последовательно
        for (int n = 1; n < Nx-1; n++)
        {
            for (int k = 1; k < Ny-1; k++)
            {
                float mod2 = m.psi[n][k].real()*m.psi[n][k].real() + m.psi[n][k].imag()*m.psi[n][k].imag();
                m.npsi[n][k] = m.psi[n][k]
                        + m.C[0]*(conj(m.Ux[n][k])*(m.Ux[n+1][k]*m.psi[n+1][k] + m.Ux[n-1][k]*m.psi[n-1][k]) - m.psi[n][k] - m.psi[n][k])
                        + m.C[1]*(conj(m.Uy[n][k])*(m.Uy[n][k+1]*m.psi[n][k+1] + m.Uy[n][k-1]*m.psi[n][k-1]) - m.psi[n][k] - m.psi[n][k])
                        + m.dt*(1.0f - mod2)*m.psi[n][k]
                        - im*m.dt*kappa*m.phi[n][k]*m.psi[n][k];

                m.delta_psi += abs(m.npsi[n][k] - m.psi[n][k])*abs(m.npsi[n][k] - m.psi[n][k]);
            }
        }
        /// граничные условия d(Psi*U_mu)/dmu |_{граница} = 0
        for (int k = 1; k < Ny-1; k++)//вдоль оси y
        {
            m.npsi[0][k] = m.psi[1][k]*m.Ux[1][k]/m.Ux[0][k];               // bound: x = -Lx/2
            m.npsi[Nx-1][k] = m.psi[Nx-2][k]*m.Ux[Nx-2][k]/m.Ux[Nx-1][k];   // bound: x = +Lx/2
            m.delta_psi += abs(m.npsi[0][k] - m.psi[0][k])*abs(m.npsi[0][k] - m.psi[0][k])
                   + abs(m.npsi[Nx-1][k] - m.psi[Nx-1][k])*abs(m.npsi[Nx-1][k] - m.psi[Nx-1][k]);
        }
        for (int n = 1; n < Nx-1; n++)//вдоль оси x
        {
            m.npsi[n][0] = m.psi[n][1]*m.Uy[n][1]/m.Uy[n][0];               // bound: y = -Ly/2
            m.npsi[n][Ny-1] = m.psi[n][Ny-2]*m.Uy[n][Ny-2]/m.Uy[n][Ny-1];   // bound: y = +Ly/2
            m.delta_psi += abs(m.npsi[n][0] - m.psi[n][0])*abs(m.npsi[n][0] - m.psi[n][0])
                   + abs(m.npsi[n][Ny-1] - m.psi[n][Ny-1])*abs(m.npsi[n][Ny-1] - m.psi[n][Ny-1]);
        }
        /// углы
        m.npsi[0][0] = 0.5f*(m.npsi[0][1] + m.npsi[1][0]);
        m.npsi[0][Ny-1] = 0.5f*(m.npsi[0][Ny-2] + m.npsi[1][Ny-1]);
        m.npsi[Nx-1][0] = 0.5f*(m.npsi[Nx-2][0] + m.npsi[Nx-1][1]);
        m.npsi[Nx-1][Ny-1] = 0.5f*(m.npsi[Nx-1][Ny-2] + m.npsi[Nx-2][Ny-1]);
        /// учёт электрического поля - инжекция тока
        if (m.boolElectricField)
        {
            for (int n = 0; n < Nx; n++)
            {
                m.npsi[n][0] -= im*m.J*dy*conj(m.psi[n][0])/m.Uy[n][0];// bound: y = -Ly/2
//                npsi[n][0] *= exp(im*kappa*(-0.5f*Ly)*mp.J0/npsi[n][0]/conj(npsi[n][0])*dy);        // bound: y = -Ly/2
                m.npsi[n][Ny-1] += im*m.J*dy*conj(m.psi[n][Ny-1])/m.Uy[n][Ny-1];// bound: y = +Ly/2
//                npsi[n][Ny-1] *= exp(im*kappa*(0.5f*Ly)*mp.J0/npsi[n][Ny-1]/conj(npsi[n][Ny-1])*dy);// bound: y = +Ly/2
            }
        }
        /// перенос данных: psi <- npsi
        for (int n = 0; n < Nx; n++)
            for (int k = 0; k < Ny; k++)
                m.psi[n][k] = m.npsi[n][k];
        ///
        if (m.counter % divider == 0)
        {
            cout << "iteration: " << m.counter << "\tdelta_psi: " << m.delta_psi << endl;
            emit dataReady();
        }
        m.t += m.dt;
        m.counter++;
    }
    cout << "iteration: " << m.counter << "\tdelta_psi: " << m.delta_psi << endl;
    cout << "time of calculation: "
         << (timer.currentMSecsSinceEpoch() - start_time)/1000 << "[s]" << endl;
    ///
    isTerminated = true;
    if (isFinished)
        emit processFinished();
}
///
/// \brief CCalc::process2
/// iteration: 10077 -> time of calculation: 42[s] =~240[it/sec] (8 threads)
/// iteration: 10117 -> time of calculation: 28[s] =~361[it/sec] (16 threads)
/// iteration: 29969 -> time of calculation: 85[s] =~352[it/sec] (16 threads)
/// iteration: 30096 -> time of calculation: 87[s] =~345[it/sec] (16 threads)
/// iteration: 30125 -> time of calculation: 86[s] =~350[it/sec] (16 threads)
/// ускорение: ~5.7(8); ~8.3(16)
/// в однопотоковом варианте этот процесс выдаёт:
/// iteration: 5021 -> time of calculation: 114[s] ~44[it/sec] (1 threads)
void CCalc::process2()
{
    static QThreadPool* ThreadPool = QThreadPool::globalInstance();
    start_time = timer.currentSecsSinceEpoch();
    cout << "Start the parallel process of calculations with the link variables" << endl
         << "on the maxThreadCount(): " << ThreadPool->maxThreadCount()
         << " threads" << endl;
    isTerminated = isFinished = false;
    if (!boolContinue)
    {
        m.delta_psi = 1.0;// если используется для прерывания расчетов
        m.t = 0.0;
        m.counter = 0;
    }
    m.countU = 0;
    m.Sum4Uavg = .0f;
    СTask1::setModel( &m );
    СTask2::setModel( &m );
    СTask3::setModel( &m );
    СTask4::setModel( &m );
    СTask5::setModel( &m );
//    quint64 msec, min, hour;
    quint64 sec, psec = timer.currentSecsSinceEpoch() - start_time;
    int container_counter = 0;
    float pUavg = .0f;
    /// подготовка файла
    std::string FileName = std::string( "J-U_B=0.5.dat" );// c расширением файла
    std::ofstream outFile( FileName, std::ofstream::out );
    /// очистка файла
    outFile.precision(5);
    outFile.setf( std::ios::fixed, std::ios::floatfield );
    outFile.close();
    /// получаем список заданий
    getListJ();
    QList<double>::Iterator itJ = listJ.end()-1;
//    QList<double>::Iterator itJ = listJ.begin();
    while (!isTerminated && !isFinished)
    {
        /// однопоточный вариант вычислений
        /// 2. вычисление div{Js}
        //get_divJs();
        /// 3. вычисление фурье-образа div{Js} методом трапеций
        //get_FdivJs();
        /// 4. решение уравнения Пуассона для фурье-образа потенциала методом прогонки
        //get_Fphi();
        /// 5. восстановление потенциала из его фурье-образа
        //get_Potential();
        /// многопоточный вариант вычислений
        /// 1. многопоточное вычисление параметра порядка
        m.delta_psi = 0.0;
        get_parallel_psi();
//        cout << "(" << СTask1::Sum[0] << ", " << СTask1::Sum[1] << ")" << endl;
        /// 2. вычисление div{Js}
        get_parallel_divJs();
        /// 3. вычисление фурье-образа div{Js} методом трапеций
        get_parallel_FdivJs();
        /// 4. решение уравнения Пуассона для фурье-образа потенциала методом прогонки
        get_parallel_Fphi();
        /// 5. восстановление потенциала из его фурье-образа
        get_parallel_Potential();
        /// вычисляем мгновенное, усредненное по времени напряжение и его производную по времени
        m.Uinstant = .0f;
        for (int j = 0; j < Ny; j++)
            m.Uinstant += m.phi[Nx-1][j] - m.phi[0][j];
        m.Uinstant /= (Ny-1);
        m.countU++;
        m.Sum4Uavg += m.Uinstant;
        pUavg = m.Uavg;
        m.Uavg = m.Sum4Uavg/m.countU;
        m.dUavg = fabs(m.Uavg - pUavg);
        /// перенос данных psi <- npsi
        for (int n = 0; n < Nx; n++)
            for (int k = 0; k < Ny; k++)
                m.psi[n][k] = m.npsi[n][k];
        ///
        m.t += m.dt;
        m.counter++;
        /// вычисление статистических данных для напряжения
        if (m.countU % 10000 == 0)
        {
            /// обновляем массив маркеров отслеживаемого параметра
            emit updateScatter(m.counter, m.dUavg);
            /// I.
//            m.container[count1] = m.Uavg;
//            count1 = (count1 == NU) ? 0 : count1+1;/// циклический счётчик с шагом 10000
//            /// вычисляем среднее и дисперсию напряжения по выборке хранимой в контейнере
//            /// на основании последних NU значений измеренных с шагом 10.000
//            T SumUavg = 0.0f, Uavg, DUavg = 0.0f;
//            for (size_t i = 0; i < NU; i++)
//                SumUavg += m.container[i];
//            Uavg = SumUavg / NU;
//            for (size_t i = 0; i < NU; i++)
//                DUavg += (m.container[i] - Uavg)*(m.container[i] - Uavg);
//            DUavg /= NU;
//            cout << "iteration: " << m.counter << "\tdispersion of Uavg:\t" << DUavg << endl;
//            /// проверка условия окончания моделирования для текущего значения параметра
//            if (DUavg < 1e-6)
//            {
//                /// проверка условия перехода к следующему значению списка
//                if (listJ.size() > 0)// если имеется список параметров
//                {
//                    /// дозапись текущего результата моделирования в файл
//                    oFile.open( FileName, std::ios::app );
//                    oFile << *itJ << "\t" << m.Uavg << "\n";
//                    oFile.close();
//                    /// переход к следующему значению списка
//                    if (itJ != listJ.begin())
//                        ++itJ;
//                    /// проверка условия достижения конца списка
//                    if (itJ == listJ.end())
//                        isFinished = true;
//                    else
//                    {
//                        emit resetChart();
//                        emit nextValue(*itJ);
//                    }
//                }
//                else// списка нет -> моделирование одиночное
//                    isFinished = true;
//            }
            /// II.
            /// заполняем контейнер анализируемого параметра,
            /// используя циклический счётчик с шагом 10.000
            m.container[container_counter] = m.dUavg;
            container_counter = (container_counter < NU-1) ? container_counter+1 : 0;
            /// проверка условия заполненности контейнера анализируемого параметра
            if (m.countU > 100000)// m.countU > NU*10.000 = 100.000
            {
                /// определение условия окончания моделирования для текущего значения параметра
                bool key = true;
                for (int i = 0; i < NU; i++)
                    if (m.container[i] > m.ParamCr)
                    {
                        key = false;
                        break;
                    }
                if (key)
                {
                    /// контроль
                    for (int i = 0; i < NU; i++)
                        cout << m.container[i] << ((i < NU-1) ? "\t":"\n");
                    /// проверка условия перехода к следующему значению списка
                    if (listJ.size() > 0)// если имеется список параметров
                    {
                        /// дозапись текущего результата моделирования в файл
                        outFile.open( FileName, std::ios::app );
                        outFile << m.J << "\t" << m.Uavg << "\n";
                        outFile.close();
                        /// проверка условия перехода к следующему элементу списка
                        if (itJ != listJ.end())
                        {
                            container_counter = 0;
                            emit resetChart();
                            emit nextValue(*itJ);
                        }
                        /// переход к следующему значению списка &
                        /// проверка условия достижения конца списка
                        /// и окончания моделирования
                        if (itJ == listJ.begin())
                            isFinished = false;
                        --itJ;
//                        ++itJ;
//                        if (itJ == listJ.end())
//                            isFinished = false;
                    }
                    else// списка нет -> моделирование одиночное
                        isFinished = true;
                }
            }
        }
        /// вывод данных в терминал и на графики
        if (m.counter % divider == 0)
        {
            sec = timer.currentSecsSinceEpoch() - start_time;
            if (psec != sec)
            {
                cout << "iteration: " << m.counter << "\tdelta_psi: " << m.delta_psi; //<< endl;
                cout << "\ttime: ";
    //            hour = sec/3600;
    //            min = sec/60 - hour*3600;
    //            sec = sec - hour*3600 - min*60;
    //            if (hour > 0)
    //                 cout << hour << "[h] ";
    //            if (min > 0 || hour > 0)
    //                cout << min << "[m] " << endl;
                cout << sec << "[s]" << endl;
                data[0] = m.delta_psi;
                data[1] = m.Uinstant;
                data[2] = m.Uavg;
                data[3] = m.dUavg;// = |Uavg(t+dt)-Uavg(t)|
                psec = sec;
                emit updateChart();
            }
            emit dataReady();
        }
//        cin.get();
    }
    cout << "iteration: " << m.counter << "\tdelta_psi: " << m.delta_psi << endl;
    cout << "rate of calculation: "
         << m.counter/(timer.currentSecsSinceEpoch() - start_time) << "[iter/s]" << endl;
    ///
    isTerminated = true;
    if (isFinished)
        emit processFinished();
}
///
/// многопоточный вариант вычислений Psi
void CCalc::get_parallel_psi()
{
    static QThreadPool* ThreadPool = QThreadPool::globalInstance();
    ThreadPool->setMaxThreadCount(MaxThreadCount);
    СTask1::Counter = N;
    for (int j = 0; j < Ny;)
    {
        if (ThreadPool->activeThreadCount() < ThreadPool->maxThreadCount())
        {
            /// каждый поток вычисляет параметр порядка с учётом г.у.
            /// для всех x-координат и заданной y-координаты: Psi[0..Nx-1][j]
            СTask1* task = new СTask1(j++);
            ThreadPool->start(task);
        }
    }
    /// ожидание завершения вычисления для всего массива -> НЕ РАБОТАЕТ !
//        while (СTask::CounterActiveTask > 0 && !isTerminated)
    /// альтернативный вариант -> тоже НЕ РАБОТАЕТ !
//        ThreadPool->waitForDone(-1);
    /// граничные условия в углах -> среднее по двум ближайшим точкам
//    m.npsi[0][0] = .5f*(m.npsi[0][1] + m.npsi[1][0]);
//    m.npsi[0][Ny-1] = .5f*(m.npsi[0][Ny-2] + m.npsi[1][Ny-1]);
//    m.npsi[Nx-1][0] = .5f*(m.npsi[Nx-2][0] + m.npsi[Nx-1][1]);
//    m.npsi[Nx-1][Ny-1] = .5f*(m.npsi[Nx-1][Ny-2] + m.npsi[Nx-2][Ny-1]);
}
///
/// многопоточный вариант вычислений div{Js}
void CCalc::get_parallel_divJs()
{
    static QThreadPool* ThreadPool = QThreadPool::globalInstance();
    ThreadPool->setMaxThreadCount(MaxThreadCount);
    СTask2::Counter = Ny;
    for (int j = 0; j < Ny;)
    {
        if (ThreadPool->activeThreadCount() < ThreadPool->maxThreadCount())
        {
            /// каждый поток вычисляет div{Js}
            /// для всех x-координат и заданной y-координаты: div{Js}[0..Nx-1][j]
            СTask2* task = new СTask2(j++);
            ThreadPool->start(task);
        }
    }
    /// ожидание завершения
    /// ...
    /// г.у. для углов - арифмет. среднее двух ближайших соседей
    m.divJs[0][0] = .5f*(m.divJs[0][1] + m.divJs[1][0]);
    m.divJs[0][Ny-1] = .5f*(m.divJs[0][Ny-2] + m.divJs[1][Ny-1]);
    m.divJs[Nx-1][0] = .5f*(m.divJs[Nx-2][0] + m.divJs[Nx-1][1]);
    m.divJs[Nx-1][Ny-1] = .5f*(m.divJs[Nx-1][Ny-2] + m.divJs[Nx-2][Ny-1]);
}
///
/// однопоточный вариант вычислений div{Js}
void CCalc::get_divJs()
{
    complex<T> dJsdx, dJsdy;
    for (int k = 1; k < Ny-1; k++)
        for (int n = 1; n < Nx-1; n++)
        {
            dJsdx = conj(m.Ux[n][k])*conj(m.psi[n][k])*
                    (m.Ux[n+1][k]*m.psi[n+1][k] + m.Ux[n-1][k]*m.psi[n-1][k]);
            dJsdy = conj(m.Uy[n][k])*conj(m.psi[n][k])*
                    (m.Uy[n][k+1]*m.psi[n][k+1] + m.Uy[n][k-1]*m.psi[n][k-1]);
            m.divJs[n][k] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
        }
    for (int n = 1; n < Nx-1; n++)
    {
        /// y = -.5*Ly:
        dJsdx = conj(m.Ux[n][0])*conj(m.psi[n][0])*
                (m.Ux[n+1][0]*m.psi[n+1][0] + m.Ux[n-1][0]*m.psi[n-1][0]);
        dJsdy = conj(m.Uy[n][0])*conj(m.psi[n][0])*
            (-m.Uy[n][3]*m.psi[n][3] + 4.f*m.Uy[n][2]*m.psi[n][2]- 5.f*m.Uy[n][1]*m.psi[n][1]);
        m.divJs[n][0] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
        /// y = +.5*Ly:
        dJsdx = conj(m.Ux[n][Ny-1])*conj(m.psi[n][Ny-1])*
                (m.Ux[n+1][Ny-1]*m.psi[n+1][Ny-1] + m.Ux[n-1][Ny-1]*m.psi[n-1][Ny-1]);
        dJsdy = conj(m.Uy[n][Ny-1])*conj(m.psi[n][Ny-1])*
            (-m.Uy[n][Ny-4]*m.psi[n][Ny-4] + 4.f*m.Uy[n][Ny-3]*m.psi[n][Ny-3]- 5.f*m.Uy[n][Ny-2]*m.psi[n][Ny-2]);
        m.divJs[n][Ny-1] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
    }
    for (int k = 1; k < Ny-1; k++)
    {
        /// x = -.5*Lx:
        dJsdx = conj(m.Ux[0][k])*conj(m.psi[0][k])*
            (-m.Ux[3][k]*m.psi[3][k] + 4.f*m.Ux[2][k]*m.psi[2][k]- 5.f*m.Ux[1][k]*m.psi[1][k]);
        dJsdy = conj(m.Uy[0][k])*conj(m.psi[0][k])*
                (m.Uy[0][k+1]*m.psi[0][k+1] + m.Uy[0][k-1]*m.psi[0][k-1]);
        m.divJs[0][k] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
        /// x = +.5*Lx:
        dJsdx = conj(m.Ux[Nx-1][k])*conj(m.psi[Nx-1][k])*
            (-m.Ux[Nx-4][k]*m.psi[Nx-4][k] + 4.f*m.Ux[Nx-3][k]*m.psi[Nx-3][k]- 5.f*m.Ux[Nx-2][k]*m.psi[Nx-2][k]);
        dJsdy = conj(m.Uy[Nx-1][k])*conj(m.psi[Nx-1][k])*
                (m.Uy[Nx-1][k+1]*m.psi[Nx-1][k+1] + m.Uy[Nx-1][k-1]*m.psi[Nx-1][k-1]);
        m.divJs[0][k] = (dJsdx.imag() + p2*dJsdy.imag())/kappa;
    }
    /// г.у. для углов - арифмет. среднее двух ближайших соседей
    m.divJs[0][0] = .5f*(m.divJs[0][1] + m.divJs[1][0]);
    m.divJs[0][Ny-1] = .5f*(m.divJs[0][Ny-2] + m.divJs[1][Ny-1]);
    m.divJs[Nx-1][0] = .5f*(m.divJs[Nx-2][0] + m.divJs[Nx-1][1]);
    m.divJs[Nx-1][Ny-1] = .5f*(m.divJs[Nx-1][Ny-2] + m.divJs[Nx-2][Ny-1]);
}

void CCalc::get_parallel_FdivJs()
{
    static QThreadPool* ThreadPool = QThreadPool::globalInstance();
    ThreadPool->setMaxThreadCount(MaxThreadCount);
    СTask3::Counter = N;
    for (int n = 0; n < N;)
    {
        if (ThreadPool->activeThreadCount() < ThreadPool->maxThreadCount())
        {
            /// каждый поток вычисляет Фурье-коэф-ты дивергенции сверхпровод. тока
            /// для всех x-координат n-ой моды: FdivJs[0..Nx-1][n]
            СTask3* task = new СTask3(n++);
            ThreadPool->start(task);
        }
    }
    /// ожидание завершения
    /// ...
}

void CCalc::get_FdivJs()
{
    for (int n = 0; n < N; n++)
        for (int i = 0; i < Nx; i++)
        {
            m.FdivJs[i][n] = .0f;
            for (int j = 1; j < Ny; j++)
                m.FdivJs[i][n] += (m.divJs[i][j]*m.Cos[n][j] + m.divJs[i][j-1]*m.Cos[n][j-1]);
            m.FdivJs[i][n] /= (Ny-1);
        }
}

void CCalc::get_parallel_Fphi()
{
    static QThreadPool* ThreadPool = QThreadPool::globalInstance();
    ThreadPool->setMaxThreadCount(MaxThreadCount);
    СTask4::Counter = N;
    for (int n = 0; n < N;)
    {
        if (ThreadPool->activeThreadCount() < ThreadPool->maxThreadCount())
        {
            /// каждый поток вычисляет Фурье-коэф-ты потенциала
            /// для всех x-координат n-ой моды: Fphi[0..Nx-1][n]
            СTask4* task = new СTask4(n++);
            ThreadPool->start(task);
        }
    }
    /// ожидание завершения
    /// ...
}
///
/// \brief вычисляем Фурье-коэффициенты потенциала
/// \return
///
void CCalc::get_Fphi()
{
    static const float c1 = 9.0/8.0,  c2 = 3.0/8.0;
    static const float c3 = 1.0/24.0, c4 = 11.0/12.0;
    static const float c5 = 1.0/3.0,  c6 = 5.0/24.0;
    T f[Nx];
    /// Фурье-коэффициент для нулевой моды  (ky[n=0] = 0)
    /// решаем f'(x,ky[0]) = div{Js}(x,ky[0])*dx^2 = D(x), где f(x)=phi'(x,ky[0])
    f[0] = m.J/sigma; f[Nx-1] =-m.J/sigma;
    /// 1. учитываем оба граничных условия
    for (int i = 2; i < Nx-1; i=i+2)
    {
        f[i] = f[i-2] + 2*m.FdivJs[i-1][0]/dx;// + O(dx^2)
        f[Nx-1-i] = f[Nx+1-i] + 2*m.FdivJs[Nx-i][0]/dx;// + O(dx^2)
    }
    /// 2. учитываем левое граничное условие
//        f[1] = f[0] + (c1*m->FdivJs[1][0] + c2*m->FdivJs[0][0])*dx;// + O(dx^2)
//        for (int i = 2; i < Nx-1; i++)
//            f[i] = f[i-1] + (c3*m->FdivJs[i+1][0] + c4*m->FdivJs[i][0] + c3*m->FdivJs[i-1][0])/dx;// + O(dx^2)
//        f[Nx-1] = f[Nx-2] + (c5*m->FdivJs[Nx-1][0] + c6*m->FdivJs[Nx-2][0] - c3*m->FdivJs[Nx-3][0])/dx;// + O(dx^2)
    /// решаем phi'(x,ky[0]) = f(x)
    m.Fphi[0][0] = .0f;
    m.Fphi[1][0] = m.Fphi[0][0] + (c1*f[1] + c2*f[0])*dx;// + O(dx^2)
    for (int i = 2; i < Nx-1; i++)
        m.Fphi[i][0] = m.Fphi[i-1][0] + (c3*f[i+1] + c4*f[i] + c3*f[i-1])*dx;// + O(dx^2)
    m.Fphi[Nx-1][0] = m.Fphi[Nx-2][0] + (c5*f[Nx-1] + c6*f[Nx-2] - c3*f[Nx-3])*dx;// + O(dx^2)
    /// Фурье-коэфф-ты для остальных мод: ky[n] > 0 -> решение СЛАУ: A*Fphi = FdivJs методом прогонки
    /// Fphi, FdivJs - вектора Фурье-коэфф-тов потенциала и дивергенции сверхпр. тока
    /// A - тридиагональная матрица; a - главная диагональ; b(=1) и c(=1) - над и поддиагонали соот-но
    T D[Nx], a[Nx], b[Nx-1], a_const;
    for (int n = 1; n < N; n++)
    {
        a_const = -(2.0f + m.ky[n]*m.ky[n]*dx2);
        for (int i = 0; i < Nx; i++)
        {
            D[i] = m.FdivJs[i][n]*dx2;
            a[i] = a_const;
            if (i < Nx-1)
                b[i] = 1.f;
        }
        a[0] += 1.0f;   a[Nx-1] += 1.0f;
        /// первый прогон - исключение нижней поддиагонали с учётом, что  c[i] = b[i] = 1
        /// и дополнено делением текущей строки на a[i] для устранения больших значений
        for (int i = 1; i < Nx; i++)
        {                                       /// исходные формулы:
            a[i] = a[i-1] - 1.f/a[i];           /// a[i] = a[i-1]*a[i] - c[i-1]*b[i-1]
            if (i < Nx-1)
                b[i] = a[i-1]/a[i];             /// b[i] = a[i-1]*b[i]
            D[i] = (a[i-1]*D[i] - D[i-1])/a[i]; /// D[i] = a[i-1]*D[i] - c[i-1]*D[i-1]
            a[i] = 1.f;
        }
        /// второй прогон - вычисление корней СЛАУ с учётом, что a[i] = 1
        m.Fphi[Nx-1][n] = D[Nx-1];              /// sol[Nx-1] = D[Nx-1] / a[Nx-1];
        for (int i = Nx-2; i > -1; i--)
            m.Fphi[i][n] = D[i] - b[i]*m.Fphi[i+1][n];/// sol[i] = (D[i] - b[i]*sol[i+1]) / a[i]
    }
}

void CCalc::get_parallel_Potential()
{
    static QThreadPool* ThreadPool = QThreadPool::globalInstance();
    ThreadPool->setMaxThreadCount(MaxThreadCount);
    СTask5::Counter = N;
    for (int n = 0; n < N;)
    {
        if (ThreadPool->activeThreadCount() < ThreadPool->maxThreadCount())
        {
            /// каждый поток добавляет вклад n-ой моды во все точки phi[0..Nx-1][0..Ny-1]
            СTask5* task = new СTask5(n++);
            ThreadPool->start(task);
        }
    }
    /// ожидание завершения
    /// ...
}

void CCalc::get_Potential()
{
    for (int n = 0; n < N; n++)
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
                m.phi[i][j] = (n > 0) ? m.phi[i][j] + m.Fphi[i][n]*m.Cos[n][j] : m.Fphi[i][0];
}
///
/// ниже функции для тестирования Фурье-метода применительно к уравнению Пуассона
///
/// d2U/dx^2 + d2U/dy^2 = V(x,y)
/// dU/dy|{y=±.5*Ly} = 0;    => с.ф.: cos(k[n]*(y+.5*Ly)); k[n] = n*Pi/Ly
/// dU/dx|{x=±.5*Lx} = ∓J/sigma;
/// Подобрать правильные U & V функции для тестов !!!
T CCalc::U(float x, float y)
{
    return (m.J/Lx*x*x*x/3 + 0.0) + y*(y*y/3 - Ly*Ly/4);
}
T CCalc::V(float x, float y)
{
    return  2.f*m.J*x/Lx + 2.f*y;
}

T CCalc::U2(float x, float y)
{
    float res = .0f;
    int i = round((x + .5*Lx)/dx);
    for (int n = 0; n < N; n++)
    {
        if (m.Fphi[i][n] == .0f)
            break;
        res += m.Fphi[i][n]*((n > 0) ? cos(m.ky[n]*(y - .5f*Ly)) : .5f);
    }

    return res;
}
T CCalc::V2(float x, float y)
{
    int i = round((x + .5*Lx)/dx);
    T res;
    if (0 < i && i < Nx-1)
        res = (m.Fphi[i+1][0] - 2*m.Fphi[i][0] + m.Fphi[i-1][0])/dx2;
    else
        if (i == 0)
            res = (-m.Fphi[3][0] + 4*m.Fphi[2][0] - 5*m.Fphi[1][0] + 2*m.Fphi[0][0])/dx2;
        else// i == Nx-1
            res = (-m.Fphi[Nx-4][0] + 4*m.Fphi[Nx-3][0] - 5*m.Fphi[Nx-2][0] + 2*m.Fphi[Nx-1][0])/dx2;
    for (int n = 1; n < N; n++)
    {
        if (m.Fphi[i][n] == .0f)
            break;
        /// V(x,y) += d^2U/dy^2
        res += -m.ky[n]*m.ky[n]*m.Fphi[i][n]
                *cos(m.ky[n]*(y - .5f*Ly));
        /// V(x,y) += d^2U/dx^2
        if (0 < i && i < Nx-1)
            res += (m.Fphi[i+1][n] - 2*m.Fphi[i][n] + m.Fphi[i-1][n])/dx2
                    *cos(m.ky[n]*(y - .5f*Ly));
        else
            if (i == 0)
                res += (-m.Fphi[3][n] + 4*m.Fphi[2][n] - 5*m.Fphi[1][n] + 2*m.Fphi[0][n])/dx2
                        *cos(m.ky[n]*(y - .5f*Ly));
            else// i == Nx-1
                res += (-m.Fphi[Nx-4][n] + 4*m.Fphi[Nx-3][n] - 5*m.Fphi[Nx-2][n] + 2*m.Fphi[Nx-1][n])/dx2
                        *cos(m.ky[n]*(y - .5f*Ly));
    }

    return res;
}

void print_matrix(string name, float** matrix, int rows, int cols, int Nmax)
{
    cout << "mode | " << name << endl;
    for (int n = 0; n < rows; n++)
    {
        cout << n << ":  ";
        for (int i = 0; i < cols; i++)
            cout << matrix[i][n] << ((i < cols-1) ? "   " : "\n");
    }
    cout << "\t\t...\n";
    for (int n = Nmax-rows; n < Nmax; n++)
    {
        cout << n << ":  ";
        for (int i = 0; i < cols; i++)
            cout << matrix[i][n] << ((i < cols-1) ? "   " : "\n");
    }
    cout << "\n\n";
}
///
/// \brief тест метода Фурье для уравнения Пуассона
///
void CCalc::getNumerSol(bool key)
{
    if (key)
    {
        /// многопоточный вариант вычислений
        /// 2. вычисление div{Js} (правой части уравнения Пуассона): divJs[i][j]
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                m.divJs[i][j] = V(m.x[i], m.y[j])*dx2;
        m.J = .0f;// sigma*dU/dx_{x=±.5*Lx}
        m.J = m.J*Lx*sigma;
        /// 3. вычисление Фурье-коэффициентов div{Js} x-координат методом трапеций
        get_FdivJs();
        print_matrix("Fdiv{Js}", m.FdivJs, 5, 5, N);
        /// 4. решение уравнения Пуассона для Фурье-образа потенциала методом прогонки
        get_Fphi();
        print_matrix("Fphi", m.Fphi, 5, 5, N);
        /// 5. восстановление потенциала из его Фурье-образа
        get_Potential();
        print_matrix("phi", m.phi, 5, 5, Ny);
        ///
        emit dataReady();
        emit AnalytSolReady(false);
    }
}
///
/// \brief аналитическое решение уравнения Пуассона
///
void CCalc::getAnalytSol(bool key)
{
    if (key)
    {
        /// генерируем спектр потенциала для U2 & V2
//        for (int i = 0; i < Nx; i++)
//            for (int n = 0; n < N; n++)
////                m.Fphi[i][n] = (n < N) ? 1.f/fabs(((fabs(m.x[i]) < 1e-8) ? 1.f : m.x[i])) : .0f;
//                m.Fphi[i][n] = (n < 5) ? 5.f - 1.f*n : .0f;
        /// аналитическое решение
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                m.phi[i][j] = U(m.x[i], m.y[j]);
        emit dataReady();
        emit NumerSolReady(false);
    }
}
