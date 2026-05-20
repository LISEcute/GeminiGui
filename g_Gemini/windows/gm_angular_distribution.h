#ifndef GM_ANGULAR_DISTRIBUTION_H
#define GM_ANGULAR_DISTRIBUTION_H

#include <QDialog>
#include <QString>
#include <QUrl>

#include <map>
#include <utility>
#include <vector>

struct AngularDistEntry
{
    int z = 0;
    int n = 0;
    std::vector<float> kineticEnergy;
    std::vector<float> thetaDeg;
    std::vector<float> vz;
    std::vector<float> vxy;
};

void addAngularSample(AngularDistEntry &entry,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy);

void addAngularSample(std::map<std::pair<int, int>, AngularDistEntry> &entries,
                      int z,
                      int n,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy);

QString buildAngularDistributionHtmlPACEStyle(
    const std::map<std::pair<int, int>, AngularDistEntry> &entries,
    double sigmaTotalMb,
    int nCascades,
    double lowLimitPercent,
    double highLimitPercent,
    double compoundExcitationMeV,
    int compoundA,
    double recoilBetaCN,
    const QString &title,
    int mdir = 0,
    int inputMode = 1,
    const AngularDistEntry &neutronEntry = AngularDistEntry(),
    const AngularDistEntry &protonEntry = AngularDistEntry(),
    const AngularDistEntry &alphaEntry = AngularDistEntry());

class AngularDistributionWidget : public QDialog
{
    Q_OBJECT

public:
    explicit AngularDistributionWidget(
        const QString &htmlContent,
        const std::map<std::pair<int, int>, AngularDistEntry> &entries = std::map<std::pair<int, int>, AngularDistEntry>(),
        double sigmaTotal = 0.0,
        int nEvents = 0,
        double lowLimitPercent = 0.0,
        double highLimitPercent = 100.0,
        const QString &title = QString(),
        const AngularDistEntry &neutronEntry = AngularDistEntry(),
        const AngularDistEntry &protonEntry = AngularDistEntry(),
        const AngularDistEntry &alphaEntry = AngularDistEntry(),
        double compoundExcitationMeV = 0.0,
        int compoundA = 1,
        double recoilBetaCN = 0.0,
        int mdir = 0,
        QWidget *parent = nullptr);

private slots:
    void save_clicked();
    void print_clicked();
    void link_clicked(const QUrl &url);

private:
    void openPlotWindow(bool plotAllTables, int tableIndex, int plotKind);

    QString html;
    std::map<std::pair<int, int>, AngularDistEntry> entriesForPlots;
    double sigmaTotalForPlots = 0.0;
    int nEventsForPlots = 0;
    double lowLimitForPlots = 0.0;
    double highLimitForPlots = 100.0;
    QString plotTitle;
    AngularDistEntry neutronEntryForPlots;
    AngularDistEntry protonEntryForPlots;
    AngularDistEntry alphaEntryForPlots;
    double compoundExcitationForPlots = 0.0;
    int compoundAForPlots = 1;
    double recoilBetaForPlots = 0.0;
    int mdirForPlots = 0;
};

#endif // GM_ANGULAR_DISTRIBUTION_H
