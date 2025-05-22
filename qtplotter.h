#ifndef QTPLOTTER_H
#define QTPLOTTER_H

#include <QObject>
#include <QWidget>
#include <QChartView>
#include <QChart>
#include <QLineSeries>
#include <QScatterSeries>

#define Nlines 4+1
#define NU 10

const float ParamCr = 0.5e-7;

using namespace QtCharts;

class QtPlotter : public QChartView
{
    Q_OBJECT

public:
    QtPlotter(QWidget *parent = nullptr);
    ~QtPlotter()
    {};

    QChart chart;
    QLineSeries line[Nlines];
    QScatterSeries scatter;

    void CreateAxesAndSeries();
    void set_pointers2data(int* counter, float* data)
    {
        this->counter = counter;
        this->data = data;
    }

public slots:
    void onUpdateChart();
    void onUpdateScatter(int, float);
    void onResetChart();
    void setAutoZoom(bool autoZoom);

private:
    bool autoZoom = true;
    float min, max;
    float *data;
    int *counter, counter_start;
    bool trigger;
    int count1;
};

#endif // QTPLOTTER_H
