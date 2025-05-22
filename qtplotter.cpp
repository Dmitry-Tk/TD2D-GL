#include "qtplotter.h"
#include <qdebug.h>

static const QColor colorID[Nlines] = { Qt::green, Qt::red, Qt::blue, Qt::magenta, Qt::black};
//Qt::magenta, Qt::cyan, Qt::yellow,
//Qt::darkGreen, Qt::darkRed, Qt::darkBlue,
//Qt::darkMagenta, Qt::darkCyan, Qt::darkYellow

QtPlotter::QtPlotter(QWidget* parent) : QChartView(parent)
{
    setChart( &chart );
    setRubberBand( QChartView::RectangleRubberBand );
//    chart.setTheme( QChart::ChartTheme::ChartThemeDark );
    CreateAxesAndSeries();
    chart.legend()->setFont(QFont("Times New Roman", 18));
    QList<QAbstractAxis*> axes = chart.axes();
    axes[0]->setTitleFont(QFont("Times New Roman", 20));
    axes[1]->setTitleFont(QFont("Times New Roman", 20));
    min = 1.0/0.0;
    max =-1.0/0.0;
}

void QtPlotter::CreateAxesAndSeries()
{
//    QString l = u8"\u003C", r = u8"\u003D", text = l + "U" + r;
    for (int j = 0; j < Nlines; j++)
    {
        chart.addSeries( &line[j] );
        line[j].setPen( QPen( colorID[j], (j < Nlines-1) ? 3.0 : 1.5) );
        switch (j)
        {
            case 0:
                line[j].setName("|Ψ(t+Δt) - Ψ(t)|^2");
                break;
            case 1:
                line[j].setName("U(t)");
                break;
            case 2:
                line[j].setName("avg{U(t)}");
//                line[j].setName(text);
                break;
            case 3:
                line[j].setName("|avg{U(t+Δt)} - avg{U(t)}|");
                break;
            case 4:
                line[j].setName("crit.param.");
                break;
        }
    }
    line[Nlines-1] << QPointF(ParamCr, 0) << QPointF(ParamCr, 0);
    for (int j = 0; j < 10; j++)
        scatter << QPointF(0, 0);
    scatter.setPen( QPen( Qt::darkMagenta, 3.0 ));
    scatter.setBrush( QBrush( Qt::darkMagenta ) );
    scatter.setBorderColor( Qt::darkMagenta );
    scatter.setMarkerSize( 8.0 );
    scatter.setName("container avg{U(t)}");
    chart.addSeries( &scatter );

    chart.createDefaultAxes();
    QList<QAbstractAxis*> axes = chart.axes();
    axes[0]->setTitleText("counter");
    axes[1]->setTitleText("data");
    /// тест
//    float t;
//    float* data = new float[Nlines];
//    set_pointers2data(data, &t);
//    for (int i = 0; i < 100; i++)
//    {
//        t = static_cast<float>(i);
//        for (int j = 0; j < Nlines; j++)
//            data[j] = sin(2*M_PI*t/100 + j*M_PI/Nlines);
//        updateChart();
//        axes[0]->setMax(t);
//    }
//    axes[1]->setRange(-1.1, 1.1);
}

void QtPlotter::onResetChart()
{
    min = 1.0/0.0;
    max =-1.0/0.0;
    trigger = true;
    for (int j = 0; j < Nlines; j++)
        line[j].clear();
    line[Nlines-1] << QPointF(0.f, ParamCr) << QPointF(0.f, ParamCr);
    for (int j = 0; j < 10; j++)
        scatter.replace(j, QPointF(0, 0));
}

void QtPlotter::onUpdateChart()
{
    /// фиксация старта нового графика в рамках одной и той же симуляции
    if (trigger)
    {
        counter_start = *counter;
        trigger = false;
        count1 = 0;
    }
    /// обновление границ графика
    if (autoZoom)
    {
        QList<QAbstractAxis*> axes = chart.axes();
        axes[0]->setRange(0, *counter - counter_start);
        axes[1]->setRange(min, max);
    }
    /// обновление данных
    for (int j = 0; j < Nlines-1; j++)
    {
        line[j] << QPointF((float)(*counter - counter_start), data[j]);
        if (j == Nlines-1)// масштабирование по dUavg
        {
            if (min > data[j])
                min = data[j];
            if (max < data[j])
                max = data[j];
        }
    }
    line[Nlines-1].replace(1, QPointF((float)(*counter - counter_start), ParamCr));
}

void QtPlotter::onUpdateScatter(int counter, float value)
{
    count1++;
    if (count1 == NU)
        count1 = 0;
    scatter.replace(count1, QPointF((float)(counter - counter_start), value));
}

void QtPlotter::setAutoZoom(bool autoZoom)
{
    this->autoZoom = autoZoom;
}


