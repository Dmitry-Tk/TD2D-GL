#include "csurfacegraph.h"
#include <QtDataVisualization/QValue3DAxis>
#include <QtDataVisualization/Q3DTheme>
#include <QtCore/qmath.h>
#include <QtGui/QImage>
#include <QGraphicsTextItem>
#include <iostream>

using namespace std;
using namespace QtDataVisualization;

static int countPNG = 1;

CSurfaceGraph::CSurfaceGraph(Q3DSurface *surface, QLabel *info) : m_graph(surface), m_info(info)
{
    m_graph->setAxisX(new QValue3DAxis);
    m_graph->setAxisY(new QValue3DAxis);
    m_graph->setAxisZ(new QValue3DAxis);

    // подписи к осям
    m_graph->axisX()->setTitle("X");
    m_graph->axisX()->setTitleVisible(true);
    m_graph->axisY()->setTitle("Z");
    m_graph->axisY()->setTitleVisible(true);
    m_graph->axisZ()->setTitle("Y");
    m_graph->axisZ()->setTitleVisible(true);

    m_surfaceDataProxy = new QSurfaceDataProxy();
    m_surface3DSeries = new QSurface3DSeries(m_surfaceDataProxy);

    /// отображаем по умолчанию |Psi|^2
    what_data = _abs2;
    /// устанавливаем ортогональную проекцию
    toggleProj();
    /// задвигаем модель в отдельный поток
    calc.moveToThread( &CalcThread );
    if (boolLinkVars)
        QObject::connect( &CalcThread, &QThread::started, &calc, &CCalc::process1 );
    else
        QObject::connect( &CalcThread, &QThread::started, &calc, &CCalc::process );
    QObject::connect( &calc, &CCalc::processFinished, this, &CSurfaceGraph::stopThread );
    QObject::connect( &calc, &CCalc::dataReady, this,  &CSurfaceGraph::drawSurf );
    /// добавляем метку
    label = new QCustom3DLabel(this);
    label->setText("counter");
    label->setRotationAxisAndAngle(QVector3D(0.01f, 0.f, 0.f), -60.0);
    label->setPosition(QVector3D(0.f, 1.f, .6f));
    label->setScaling(QVector3D(1.25f, 1.25f, 1.25f));
    label->setBorderEnabled( false );
    m_graph->addCustomItem(label);

    calc.init();
}

CSurfaceGraph::~CSurfaceGraph()
{
    delete m_graph;
}

void CSurfaceGraph::initSurface()
{
//    m_surface3DSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
//    m_surface3DSeries->setDrawMode(QSurface3DSeries::DrawWireframe);
    m_surface3DSeries->setDrawMode(QSurface3DSeries::DrawSurface);

    m_graph->axisX()->setLabelFormat("%.2f");
    m_graph->axisY()->setLabelFormat("%.2f");

    m_graph->axisX()->setRange(-.5*Lx, .5*Lx);
    m_graph->axisZ()->setRange(-.5*Ly,  0.75*Ly);

    m_graph->axisX()->setLabelAutoRotation(90);
    m_graph->axisY()->setLabelAutoRotation(90);
//    m_graph->axisZ()->setLabelAutoRotation(30);

    // поверхностное плоское затенение (surface flat shading)
    m_surface3DSeries->setFlatShadingEnabled(true);
    m_surface3DSeries->setFlatShadingEnabled(false);

    m_graph->addSeries(m_surface3DSeries);

    // Reset range sliders
//    m_rangeMinX = sampleMin;
//    m_stepX = (sampleMax - sampleMin) / float(sampleCountX - 1);
//    m_axisMinSliderX->setMaximum(sampleCountX - 2);
//    m_axisMinSliderX->setValue(0);
//    m_axisMaxSliderX->setMaximum(sampleCountX - 1);
//    m_axisMaxSliderX->setValue(sampleCountX - 1);
}

void CSurfaceGraph::changeView(int view)
{
    m_graph->scene()->activeCamera()->setCameraPreset(Q3DCamera::CameraPreset(view));
}

void CSurfaceGraph::setBlackToYellowGradient()
{
    QLinearGradient gr;
    gr.setColorAt(0.0, Qt::black);
    gr.setColorAt(0.33, Qt::blue);
    gr.setColorAt(0.67, Qt::red);
    gr.setColorAt(1.0, Qt::yellow);

    m_graph->seriesList().at(0)->setBaseGradient(gr);
    m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
}

void CSurfaceGraph::setGreenToRedGradient()
{
    QLinearGradient gr;
    gr.setColorAt(0.0, Qt::darkGreen);
    gr.setColorAt(0.5, Qt::yellow);
    gr.setColorAt(0.8, Qt::red);
    gr.setColorAt(1.0, Qt::darkRed);

    m_graph->seriesList().at(0)->setBaseGradient(gr);
    m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
}

void CSurfaceGraph::start_stopThread()
{
    if (!CalcThread.isRunning())
    {
        if (!boolContinue)
        {
            countPNG = 1;
            calc.init();
        }
        CalcThread.start();
    }
    else
    {
        calc.isTerminated = true;
        CalcThread.quit();
        //CalcThread.wait();
    }
}

void CSurfaceGraph::stopThread()
{
    CalcThread.quit();
    //CalcThread.wait();
    emit threadStopped();
}

void CSurfaceGraph::drawSurf()
{
    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve(Ny);
    float z;
    for (int j = 0; j < Ny; j++)
    {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(Nx);
        int index = 0;
        for (int i = 0; i < Nx; i++)
        {
            float re = calc.m.psi[i][j].real();
            float im = calc.m.psi[i][j].imag();
            switch (what_data)
            {
            case _re:
                z = re;
                break;
            case _im:
                z = im;
                break;
            case _arg:// фаза параметра порядка
                z = atan(im/re) + M_PI;
                break;
            case _divJs:// дивергенция сверхпроводящего тока
                z = calc.m.divJs[i][j]/dx2;
                break;
            case _phi:// распределение потенциала
                z = calc.m.phi[i][j];
                break;
            case _kappa:// распределение kappa
                z = 0.0;
                break;
            default:  // плотность параметра порядка
                z = re*re + im*im;
            }
            if (i + j == 0)
            {
                minZ = maxZ = z;
            }
            if (maxZ < z && !isnan(z) && !isinf(z))
                maxZ = z;
            if (minZ > z && !isnan(z) && !isinf(z))
                minZ = z;

            (*newRow)[index++].setPosition(QVector3D(calc.m.x[i], z, calc.m.y[j]));
        }
        *dataArray << newRow;
    }
    if (autoZoom)
        m_graph->axisY()->setRange(minZ, maxZ);

    m_surfaceDataProxy->resetArray(dataArray);

    /// отображаем текущее значение счётчика и стат. показателей
    QString str1 =
//            QString("time = ") + QString::number(calc.getTime()) +
//            QString("counter = ") + QString::number(calc.m.counter) +
            QString("U = ") + QString::number(calc.m.Uinstant) +
            QString("\t\tavg{ U } = ") + QString::number(calc.m.Uavg) +
            QString("\t\tΔU = ") + QString::number(fabs((calc.m.Uavg - calc.m.Uinstant)/calc.m.Uavg));
    m_info->setText(str1);
    ///
    QString str2 = QString("counter = ") + QString::number(calc.m.counter) +
                    QString("\tJ = ") + QString::number(calc.m.J) +
                    QString("\tB = ") + QString::number(calc.m.B);
    QCustom3DLabel* nlabel = new QCustom3DLabel(nullptr);
    nlabel->setScaling( label->scaling() );
    nlabel->setPosition( label->position() );
    nlabel->setRotation( label->rotation() );
    nlabel->setText( str2 );
    m_graph->releaseCustomItem( label );
    label = nlabel;
    m_graph->addCustomItem( label );
    /// скриншоты
//    savePng();
}

void CSurfaceGraph::setAutoZoom(bool autoZoom)
{
    this->autoZoom = autoZoom;
}

void CSurfaceGraph::setMaxZ(double maxZ)
{
    if (!autoZoom)
    {
        this->minZ = 0.0;
        this->maxZ = maxZ;
        m_graph->axisY()->setRange(this->minZ, this->maxZ);
    }
}

void CSurfaceGraph::setDivider(int divider)
{
    calc.divider = divider;
}

void CSurfaceGraph::setMagneticField(double B)
{
    calc.updateMagneticField(static_cast<float>(B));
}

void CSurfaceGraph::setElectricField(double J)
{
    calc.updateInjectionСurrent(static_cast<float>(J));
}

void CSurfaceGraph::setTimeStep(double dt)
{
    calc.updateTimeStep(static_cast<float>(dt));
}

void CSurfaceGraph::setBoundCond(bool periodic_bound_cond)
{
    calc.setBoundCond(periodic_bound_cond);
}

void CSurfaceGraph::setBoolContinue(bool boolContinue)
{
    this->boolContinue = boolContinue;
    calc.setBoolConinue(boolContinue);
}

void CSurfaceGraph::setBoolLinkVars(bool boolLinkVars)
{
    QObject::disconnect(&CalcThread, 0, &calc, 0);
    if (boolLinkVars)
        QObject::connect( &CalcThread, &QThread::started, &calc, &CCalc::process2 );
    else
        QObject::connect( &CalcThread, &QThread::started, &calc, &CCalc::process );
}

void CSurfaceGraph::toggleProj()
{
    static bool proj = true;
    QAbstract3DGraph* surf = static_cast<QAbstract3DGraph*>(m_graph);
    surf->setOrthoProjection(proj);
    proj = !proj;
}

void CSurfaceGraph::toggleDataType(int index)
{
    switch (index)
    {
    case 0:
        what_data = type_data::_re;
        break;
    case 1:
        what_data = type_data::_im;
        break;
    case 3:
        what_data = type_data::_arg;
        break;
    case 4:
        what_data = type_data::_divJs;
        break;
    case 5:
        what_data = type_data::_phi;
        break;
    case 6:
        what_data = type_data::_kappa;
        break;
    default:
        what_data = type_data::_abs2;
    }
    if (!CalcThread.isRunning())
        drawSurf();
}

void CSurfaceGraph::savePng()
{
    if (CalcThread.isRunning())
    {
        QString filename = tr("/home/dmitry/Work/projects/_Qt/TD2D-GL/png/")
                + QString::number(countPNG) + tr(".png");
        QImage image = m_graph->renderToImage(0, QSize(900, 450));
        image.save(filename);
        countPNG++;
    }
//    QString filename = tr("/home/dmitry/Work/projects/_Qt/TD2D-GL/png/")
//            + QString::number(countPNG) + tr(".png");
//    QWidget wgt = static_cast<QWidget>(this->parent());
//    QPixmap p = QPixmap::grabWindow(QRect(0, 0, 600, 500));
//    p.save(filename);
}
