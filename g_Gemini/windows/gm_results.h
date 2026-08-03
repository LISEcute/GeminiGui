#ifndef RESULTS_H
#define RESULTS_H

#include <QDialog>
#include <QString>
#include <QStringList>
#include <QTabWidget>
#include <QUrl>
#include <QVector>

#include "gm_angular_distribution.h"

struct YieldPlotPoint
{
    int z = 0;
    double residualEvents = 0.0;
    double imfEvents = 0.0;
};

struct YieldPlotData
{
    QVector<YieldPlotPoint> points;
};

struct AngularResultTab
{
    QString label;
    QString html;
    std::map<std::pair<int, int>, AngularDistEntry> entries;
    double sigmaTotal = 0.0;
    int nEvents = 0;
    double lowLimitPercent = 0.0;
    double highLimitPercent = 100.0;
    QString title;
    AngularDistEntry neutronEntry;
    AngularDistEntry protonEntry;
    AngularDistEntry alphaEntry;
    AngularDistEntry gammaEntry;
    double compoundExcitationMeV = 0.0;
    int compoundA = 1;
    int compoundZ = 0;
    double recoilBetaCN = 0.0;
    int mdir = 0;
};

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class Result_Widget : public QDialog
{
    Q_OBJECT
public:
    explicit Result_Widget(QString result_string, QWidget *parent = 0);
    explicit Result_Widget(QString result_string,
                           const YieldPlotData &yieldPlot,
                           QWidget *parent = 0);
    explicit Result_Widget(QString result_string,
                           const AngularResultTab &residualAngular,
                           const AngularResultTab &imfAngular,
                           QWidget *parent = 0);
    explicit Result_Widget(QString result_string,
                           const AngularResultTab &residualAngular,
                           const AngularResultTab &imfAngular,
                           const YieldPlotData &yieldPlot,
                           QWidget *parent = 0);

signals:

public slots:
    void save_clicked();
    void print_clicked();
    void link_clicked(QUrl url);
private:
    void init(QString result_string,
              const AngularResultTab *residualAngular,
              const AngularResultTab *imfAngular,
              const YieldPlotData *yieldPlot);
    QString currentHtml() const;
    void openYieldPlotWindow();

    QString result;
    QTabWidget *tabs = nullptr;
    QStringList tabHtml;
    YieldPlotData yieldPlotData;
};

#endif // RESULTS_H
