#ifndef CSURFACEGRAPH_H
#define CSURFACEGRAPH_H

#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QHeightMapSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtWidgets/QSlider>
#include <QThread>
#include <QLabel>
#include <QCustom3DLabel>
#include "ccalc.h"

using namespace QtDataVisualization;

enum type_data {_re, _im, _abs2, _arg, _divJs, _phi, _kappa};

class CSurfaceGraph : public QObject
{
    Q_OBJECT
public:
    explicit CSurfaceGraph(Q3DSurface *surface, QLabel *info);
    ~CSurfaceGraph();

    void initSurface();
    void savePng();
    //! [0]
    void toggleModeNone() { m_graph->setSelectionMode(QAbstract3DGraph::SelectionNone); }
    void toggleModeItem() { m_graph->setSelectionMode(QAbstract3DGraph::SelectionItem); }
    void toggleModeSliceRow() { m_graph->setSelectionMode(QAbstract3DGraph::SelectionItemAndRow
                                                          | QAbstract3DGraph::SelectionSlice); }
    void toggleModeSliceColumn() { m_graph->setSelectionMode(QAbstract3DGraph::SelectionItemAndColumn
                                                             | QAbstract3DGraph::SelectionSlice); }
    //! [0]
    void toggleProj();
    void toggleDataType(int index);

    void setBlackToYellowGradient();
    void setGreenToRedGradient();

    void setAxisMinSliderX(QSlider *slider) { m_axisMinSliderX = slider; }
    void setAxisMaxSliderX(QSlider *slider) { m_axisMaxSliderX = slider; }
    void setAxisMinSliderZ(QSlider *slider) { m_axisMinSliderZ = slider; }
    void setAxisMaxSliderZ(QSlider *slider) { m_axisMaxSliderZ = slider; }

    void adjustXMin(int min);
    void adjustXMax(int max);
    void adjustZMin(int min);
    void adjustZMax(int max);

    void setCalc(CCalc* calc, QThread* thread);
    CCalc calc;

signals:
    void threadStopped();
    void data4ChartReady();

public Q_SLOTS:
    void changeView(int view);
    void start_stopThread();
    void stopThread();
    void drawSurf();
    void setAutoZoom(bool autoZoom);
    void setMaxZ(double maxZ);
    void setDivider(int divider);
    void setMagneticField(double B);
    void setElectricField(double J);
    void setTimeStep(double dt);
    void setBoundCond(bool periodic_bound_cond);
    void setBoolContinue(bool boolContinue);
    void setBoolLinkVars(bool boolLinkVars);

private:
    QCustom3DLabel* label;
    Q3DSurface *m_graph;
    QLabel *m_info;
//    QHeightMapSurfaceDataProxy *m_heightMapProxy;
    QSurfaceDataProxy *m_surfaceDataProxy;
//    QSurface3DSeries *m_heightMapSeries;
    QSurface3DSeries *m_surface3DSeries;

    QSlider *m_axisMinSliderX;
    QSlider *m_axisMaxSliderX;
    QSlider *m_axisMinSliderZ;
    QSlider *m_axisMaxSliderZ;
    float m_rangeMinX;
    float m_stepX;
    int m_heightMapWidth;
    int m_heightMapHeight;

    void setAxisXRange(float min, float max);
    void setAxisZRange(float min, float max);

    bool autoZoom = false;
    bool boolContinue = false;
    bool boolLinkVars = false;
    double minZ, maxZ;
    type_data what_data = type_data::_abs2;
    QThread CalcThread;
};

#endif // CSURFACEGRAPH_H
