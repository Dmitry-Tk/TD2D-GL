/**
 *  Начало проекта:  2025-01-09
 *  Решение время-зависимого 2D уравнения Гинзбурга-Ландау в прямогоульной области
 *  с использованием множества потоков/нитей
 *  Вычисление ВАХ моделируемой системы
 *  Численное решение двумерного нестационарного уравнения Гинзубурга-Ландау для прямоугольной области сверхпроводника 2-го рода методом конечных разностей. Расчёт статистически усредненной ВАХ системы при различных величинах магнитного поля.
**/
#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QSlider>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMessageBox>
#include <QtGui/QScreen>
#include "csurfacegraph.h"
#include "qtplotter.h"

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    QTabWidget *tabWidget = new QTabWidget();
    /// 1.
    Q3DSurface *graph = new Q3DSurface();
    QWidget *container = QWidget::createWindowContainer(graph);
    /// 2.
    QtPlotter plotter;

    tabWidget->addTab(container, QString("3d-model"));
    tabWidget->addTab(&plotter, QString("|Ψ(t+Δt)-Ψ(t)|^2, U, <U>"));

    container->setWindowTitle("Graph");
    if (!graph->hasContext())
    {
        QMessageBox msgBox;
        msgBox.setText("Couldn't initialize the OpenGL context.");
        msgBox.exec();
        return -1;
    }

    QSize screenSize = graph->screen()->size();
    container->setMinimumSize(QSize(screenSize.width() / 2, screenSize.height() / 1.6));
    container->setMaximumSize(screenSize);
    container->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    container->setFocusPolicy(Qt::StrongFocus);

    QWidget *widget = new QWidget();
    widget->setWindowTitle(QStringLiteral("Time-dependent Ginzburg-Landau 2D eq"));

    /// инфо-панель - текущий шаг моделирования
    QLabel *info = new QLabel("Info");
    QFont font = info->font();
    font.setPointSize(18);
    info->setFont(font);
    info->setAlignment(Qt::AlignHCenter);
    info->setAutoFillBackground(true);
    QPalette pal = info->palette();
    pal.setColor(QPalette::Window, QColor(Qt::white));
    info->setPalette(pal);

    QHBoxLayout *hLayout = new QHBoxLayout(widget);
    QVBoxLayout *vLayoutL = new QVBoxLayout();

    vLayoutL->addWidget(tabWidget, 1);
    vLayoutL->addWidget(info);
    hLayout->addLayout(vLayoutL, 1);
    QVBoxLayout *vLayoutR = new QVBoxLayout();
    vLayoutR->setAlignment(Qt::AlignTop);
    hLayout->addLayout(vLayoutR);

    QGroupBox *projectionGroupBox = new QGroupBox(QStringLiteral("Projection"));
    QRadioButton *orthoProj = new QRadioButton(widget);
    orthoProj->setText(QStringLiteral("ortho"));
    orthoProj->setChecked(false);
    QRadioButton *perspecProj = new QRadioButton(widget);
    perspecProj->setText(QStringLiteral("perspec"));
    perspecProj->setChecked(false);
    QHBoxLayout *hboxproj= new QHBoxLayout;
    hboxproj->addWidget(orthoProj);
    hboxproj->addWidget(perspecProj);
    projectionGroupBox->setLayout(hboxproj);

    QGroupBox *zoomGroupBox = new QGroupBox(QStringLiteral("Zoom"));
    QCheckBox *chboxAutoZoom = new QCheckBox(QStringLiteral("autoscaling Z-axis"), widget);
    QDoubleSpinBox *spboxZoom = new QDoubleSpinBox(widget);
    spboxZoom->setLocale(QLocale("C"));
    spboxZoom->setRange(0.1, 10.1);
    spboxZoom->setDecimals(1);
    spboxZoom->setSingleStep(0.1);
    spboxZoom->setValue(1.1);
    spboxZoom->setToolTip(QString("max Z"));
    spboxZoom->setAlignment(Qt::AlignHCenter);
    QCheckBox *chboxAutoZooChart = new QCheckBox(QStringLiteral("autoscaling Chart"), widget);

    QVBoxLayout *zoomVBox = new QVBoxLayout;
    zoomVBox->addWidget(chboxAutoZoom);
    zoomVBox->addWidget(spboxZoom);
    zoomVBox->addWidget(chboxAutoZooChart);
    zoomGroupBox->setLayout(zoomVBox);

    QSpinBox *spboxDivider= new QSpinBox(widget);
    spboxDivider->setRange(1, 1000);
    spboxDivider->setValue(100);
    spboxDivider->setPrefix(QString("every "));
    spboxDivider->setSuffix(QString(" steps"));
    spboxDivider->setAlignment(Qt::AlignHCenter);

    QComboBox *listDataType = new QComboBox(widget);
    listDataType->addItems({"re{ Ψ }", "im{ Ψ }", "|Ψ|^2", "arg{ Ψ }",
            "div{ Js }", "φ", "κ"});

    QGroupBox *grboxModelParams = new QGroupBox(QStringLiteral("Model parameters"));
    QDoubleSpinBox *spboxdt = new QDoubleSpinBox(widget);
    spboxdt->setLocale(QLocale("C"));
    spboxdt->setDecimals(4);
    spboxdt->setRange(0.0001, 0.5);
    spboxdt->setSingleStep(0.025);
    spboxdt->setValue(0.05);
    spboxdt->setToolTip(QString("time step"));
    spboxdt->setAlignment(Qt::AlignHCenter);
    QDoubleSpinBox *spboxB = new QDoubleSpinBox(widget);
    spboxB->setLocale(QLocale("C"));
    spboxB->setDecimals(3);
    spboxB->setRange(0.0, 50.0);
    spboxB->setSingleStep(0.25);
    spboxB->setToolTip(QString("magnetic field"));
    spboxB->setAlignment(Qt::AlignHCenter);
    QDoubleSpinBox *spboxJ0 = new QDoubleSpinBox(widget);
    spboxJ0->setLocale(QLocale("C"));
    spboxJ0->setDecimals(3);
    spboxJ0->setRange(-50.0, 50.0);
    spboxJ0->setSingleStep(0.1);
    spboxJ0->setToolTip(QString("injection current"));
    spboxJ0->setAlignment(Qt::AlignHCenter);
    QCheckBox *chboxBoundCond = new QCheckBox(QStringLiteral("periodic bound cond-s"), widget);
    chboxBoundCond->setChecked(false);
    QCheckBox *chboxContinue = new QCheckBox(QStringLiteral("continue"), widget);
    chboxContinue->setChecked(false);
    QCheckBox *chboxLinkVars = new QCheckBox(QStringLiteral("link variables"), widget);
    chboxLinkVars->setChecked(false);
    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(spboxdt);
    vbox->addWidget(spboxB);
    vbox->addWidget(spboxJ0);
    vbox->addWidget(chboxBoundCond);
    vbox->addWidget(chboxContinue);
    vbox->addWidget(chboxLinkVars);
    grboxModelParams->setLayout(vbox);

    QGroupBox *selectionGroupBox = new QGroupBox(QStringLiteral("Selection Mode"));
    QRadioButton *modeNoneRB = new QRadioButton(widget);
    modeNoneRB->setText(QStringLiteral("No selection"));
    modeNoneRB->setChecked(false);
    QRadioButton *modeItemRB = new QRadioButton(widget);
    modeItemRB->setText(QStringLiteral("Item"));
    modeItemRB->setChecked(false);
    QRadioButton *modeSliceRowRB = new QRadioButton(widget);
    modeSliceRowRB->setText(QStringLiteral("Row Slice"));
    modeSliceRowRB->setChecked(false);
    QRadioButton *modeSliceColumnRB = new QRadioButton(widget);
    modeSliceColumnRB->setText(QStringLiteral("Column Slice"));
    modeSliceColumnRB->setChecked(false);

    QPushButton *start_stop = new QPushButton(widget);
    start_stop->setText("start");
    start_stop->setCheckable(true);

    QPushButton *theor_sol = new QPushButton(widget);
    theor_sol->setText("theor.");
    theor_sol->setCheckable(true);

    QPushButton *numer_sol = new QPushButton(widget);
    numer_sol->setText("numer.");
    numer_sol->setCheckable(true);

    hLayout = new QHBoxLayout;
    hLayout->addWidget(theor_sol);
    hLayout->addWidget(numer_sol);

    QVBoxLayout *selectionVBox = new QVBoxLayout;
    selectionVBox->addWidget(modeNoneRB);
    selectionVBox->addWidget(modeItemRB);
    selectionVBox->addWidget(modeSliceRowRB);
    selectionVBox->addWidget(modeSliceColumnRB);
    selectionGroupBox->setLayout(selectionVBox);

    QSlider *axisMinSliderX = new QSlider(Qt::Horizontal, widget);
    axisMinSliderX->setMinimum(0);
    axisMinSliderX->setTickInterval(1);
    axisMinSliderX->setEnabled(true);
    QSlider *axisMaxSliderX = new QSlider(Qt::Horizontal, widget);
    axisMaxSliderX->setMinimum(1);
    axisMaxSliderX->setTickInterval(1);
    axisMaxSliderX->setEnabled(true);

    QComboBox *viewList = new QComboBox(widget);
    viewList->addItems({"None",
            "FrontLow", "Front", "FrontHigh",
            "LeftLow", "Left", "LeftHigh",
            "RightLow", "Right", "RightHigh",
            "BehindLow", "Behind", "BehindHigh",
            "IsometricLeft", "~(1, 1, 1)", "IsometricRight", "XY plane",
            "DirectlyAbove", "DirectlyAboveCW", "XZ plane",
            "FrontBelow", "LeftBelow", "RightBelow", "BehindBelow", "DirectlyBelow"});

//    "IsometricLeftHigh" = "~(1, 1, 1)"
//    "IsometricRightHigh" = "XY plane"
//    "DirectlyAboveCCW" = "XZ plane"
    QGroupBox *colorGroupBox = new QGroupBox(QStringLiteral("Custom gradient"));

    QLinearGradient grBtoY(0, 0, 1, 100);
    grBtoY.setColorAt(1.0, Qt::black);
    grBtoY.setColorAt(0.67, Qt::blue);
    grBtoY.setColorAt(0.33, Qt::red);
    grBtoY.setColorAt(0.0, Qt::yellow);
    QPixmap pm(24, 100);
    QPainter pmp(&pm);
    pmp.setBrush(QBrush(grBtoY));
    pmp.setPen(Qt::NoPen);
    pmp.drawRect(0, 0, 24, 100);
    QPushButton *gradientBtoYPB = new QPushButton(widget);
    gradientBtoYPB->setIcon(QIcon(pm));
    gradientBtoYPB->setIconSize(QSize(24, 100));

    QLinearGradient grGtoR(0, 0, 1, 100);
    grGtoR.setColorAt(1.0, Qt::darkGreen);
    grGtoR.setColorAt(0.5, Qt::yellow);
    grGtoR.setColorAt(0.2, Qt::red);
    grGtoR.setColorAt(0.0, Qt::darkRed);
    pmp.setBrush(QBrush(grGtoR));
    pmp.drawRect(0, 0, 24, 100);
    QPushButton *gradientGtoRPB = new QPushButton(widget);
    gradientGtoRPB->setIcon(QIcon(pm));
    gradientGtoRPB->setIconSize(QSize(24, 100));

    QHBoxLayout *colorHBox = new QHBoxLayout;
    colorHBox->addWidget(gradientBtoYPB);
    colorHBox->addWidget(gradientGtoRPB);
    colorGroupBox->setLayout(colorHBox);

    vLayoutR->addWidget(start_stop);
    vLayoutR->addLayout(hLayout);
    vLayoutR->addWidget(projectionGroupBox);
    vLayoutR->addWidget(zoomGroupBox);
    vLayoutR->addWidget(new QLabel(QStringLiteral("frames are displayed after")));
    vLayoutR->addWidget(spboxDivider);
    vLayoutR->addWidget(new QLabel(QStringLiteral("type of data")));
    vLayoutR->addWidget(listDataType);
    vLayoutR->addWidget(grboxModelParams);
    vLayoutR->addWidget(selectionGroupBox);
    vLayoutR->addWidget(new QLabel(QStringLiteral("View")));
    vLayoutR->addWidget(viewList);
    vLayoutR->addWidget(colorGroupBox);
    vLayoutR->addWidget(new QLabel(QStringLiteral("Column range")));
    vLayoutR->addWidget(axisMinSliderX);
    vLayoutR->addWidget(axisMaxSliderX);

    widget->show();

    CSurfaceGraph *modifier = new CSurfaceGraph(graph, info);

    QObject::connect(modeNoneRB, &QRadioButton::toggled,
                     modifier, &CSurfaceGraph::toggleModeNone);
    QObject::connect(modeItemRB,  &QRadioButton::toggled,
                     modifier, &CSurfaceGraph::toggleModeItem);
    QObject::connect(modeSliceRowRB,  &QRadioButton::toggled,
                     modifier, &CSurfaceGraph::toggleModeSliceRow);
    QObject::connect(modeSliceColumnRB,  &QRadioButton::toggled,
                     modifier, &CSurfaceGraph::toggleModeSliceColumn);
    QObject::connect(viewList, qOverload<int>(&QComboBox::currentIndexChanged),
                     modifier, &CSurfaceGraph::changeView);
    QObject::connect(gradientBtoYPB, &QPushButton::pressed,
                     modifier, &CSurfaceGraph::setBlackToYellowGradient);
    QObject::connect(gradientGtoRPB, &QPushButton::pressed,
                     modifier, &CSurfaceGraph::setGreenToRedGradient);
    /// добавлено
    /// управление настройками отображения графика
    QObject::connect(start_stop, &QPushButton::pressed,
                     modifier, &CSurfaceGraph::start_stopThread);
    QObject::connect(modifier, &CSurfaceGraph::threadStopped,
                      start_stop, &QPushButton::toggle);
    QObject::connect(orthoProj, &QRadioButton::pressed,
                     modifier, &CSurfaceGraph::toggleProj);
    QObject::connect(perspecProj, &QRadioButton::pressed,
                     modifier, &CSurfaceGraph::toggleProj);
    QObject::connect(chboxAutoZoom, &QCheckBox::clicked,
                     spboxZoom, &QDoubleSpinBox::setDisabled);
    QObject::connect(chboxAutoZoom, &QCheckBox::clicked,
                     modifier, &CSurfaceGraph::setAutoZoom);
    QObject::connect(spboxZoom, qOverload<double>(&QDoubleSpinBox::valueChanged),
                     modifier, &CSurfaceGraph::setMaxZ);
    /// управление частотой отображения фрейма в процессе расчётов
    QObject::connect(spboxDivider, qOverload<int>(&QSpinBox::valueChanged),
                     modifier, &CSurfaceGraph::setDivider);
    /// управление типом отображаемых данных
    QObject::connect(listDataType, qOverload<int>(&QComboBox::currentIndexChanged),
                     modifier, &CSurfaceGraph::toggleDataType);
    /// управление параметрами модели
    QObject::connect(spboxdt, qOverload<double>(&QDoubleSpinBox::valueChanged),
                     modifier, &CSurfaceGraph::setTimeStep);
    QObject::connect(spboxB, qOverload<double>(&QDoubleSpinBox::valueChanged),
                     modifier, &CSurfaceGraph::setMagneticField);
    QObject::connect(spboxJ0, qOverload<double>(&QDoubleSpinBox::valueChanged),
                     modifier, &CSurfaceGraph::setElectricField);
    QObject::connect(chboxBoundCond, &QCheckBox::clicked,
                     modifier, &CSurfaceGraph::setBoundCond);
    QObject::connect(chboxContinue, &QCheckBox::clicked,
                     modifier, &CSurfaceGraph::setBoolContinue);
    QObject::connect(chboxLinkVars, &QCheckBox::clicked,
                     modifier, &CSurfaceGraph::setBoolLinkVars);
    QObject::connect(chboxLinkVars, &QCheckBox::clicked,
                     modifier, &CSurfaceGraph::setBoolLinkVars);
    QObject::connect(theor_sol, &QPushButton::clicked,
                     &modifier->calc, &CCalc::getAnalytSol, Qt::DirectConnection);
    QObject::connect(numer_sol, &QPushButton::clicked,
                     &modifier->calc, &CCalc::getNumerSol, Qt::DirectConnection);
    QObject::connect(&modifier->calc, &CCalc::AnalytSolReady,
                     theor_sol, &QPushButton::setChecked, Qt::DirectConnection);
    QObject::connect(&modifier->calc, &CCalc::NumerSolReady,
                     numer_sol, &QPushButton::setChecked, Qt::DirectConnection);

    modifier->setAxisMinSliderX(axisMinSliderX);
    modifier->setAxisMaxSliderX(axisMaxSliderX);
    /// настройки графика
    modifier->initSurface();
    orthoProj->setChecked(true);
    modeNoneRB->setChecked(true);
    viewList->setCurrentIndex(16);// задаём положение камеры
    gradientBtoYPB->click();
    /// настройки масштабирования
    modifier->setAutoZoom(chboxAutoZoom->isChecked());
    spboxZoom->setDisabled(chboxAutoZoom->isChecked());
    if (!chboxAutoZoom->isChecked())
        modifier->setMaxZ(spboxZoom->value());
    /// количество шагов на отображение одного фрейма
    modifier->setDivider(spboxDivider->value());
    /// магнитное поле
    spboxB->setValue(0.50);
    /// ток инжекции
    spboxJ0->setValue(0.0);
    /// активируем по умолчанию link vars
    chboxLinkVars->click();
    /// отображаем по умолчанию abs{Psi}^2
    listDataType->setCurrentIndex(2);
    /// настройка плоттера, отображающего статистические данные моделирования
    /// передаем указатель на массив со стат. данными моделирования
    plotter.set_pointers2data(&modifier->calc.m.counter, modifier->calc.data);
    /// обновление плоттера при получении новых стат. данных
    QObject::connect(&(modifier->calc), &CCalc::updateChart,
                     &plotter, &QtPlotter::onUpdateChart);
    QObject::connect(&(modifier->calc), &CCalc::updateScatter,
                     &plotter, &QtPlotter::onUpdateScatter);
    /// очистка плоттера и сброс отображаемых данных
    QObject::connect(&(modifier->calc), &CCalc::resetChart,
                     &plotter, &QtPlotter::onResetChart);
    /// связываем сигнал готовности к новому этапу моделирования
    /// при следующем значении управляющего параметра
    QObject::connect(&(modifier->calc), &CCalc::nextValue,
                     spboxJ0,  &QDoubleSpinBox::setValue);
    /// автомасштабирование плоттера
    QObject::connect(chboxAutoZooChart, &QCheckBox::clicked,
                     &plotter, &QtPlotter::setAutoZoom);
    /// активируем автомасштабирование на плоттере
    chboxAutoZooChart->click();

    return app.exec();
}
/*
to do:
    [x] правильная реализация граничного условия с учётом link vars
    [x] ускорение расчётов за счет распараллеливания
    [x] оптимизация кода
    [x] добавить метку со счётчиком шагов к графику

    [ ] добавить опцию: kappa-константа / двумерная модуляция kappa
    [ ] добавить опцию: список задач / одиночное моделирование
    [ ] сохранение/загрузка результатов моделирования в xml
    [ ] добавить опцию: загрузка равновесной конфигурации / вычисление равновесной конфигурации
    [?] возможность скрытия / открытия панели управления
    [?] разместить управление и настройки на отдельных табвиджетах
    [?] возможность указания критического значения параметра остановки моделирования

    [?] метод сверхрелаксации для уравнения Пуассона
    [?] инжекция тока -> правильная реализация граничного условия
*/


