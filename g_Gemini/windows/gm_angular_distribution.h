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
    std::vector<float> cmEnergy;
    std::vector<float> weight;
};

void addAngularSample(AngularDistEntry &entry,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy);

void addAngularSample(AngularDistEntry &entry,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy,
                      float cmEnergy);

void addAngularSample(AngularDistEntry &entry,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy,
                      float cmEnergy,
                      float weight);

void addAngularSample(std::map<std::pair<int, int>, AngularDistEntry> &entries,
                      int z,
                      int n,
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
                      float vxy,
                      float cmEnergy);

void addAngularSample(std::map<std::pair<int, int>, AngularDistEntry> &entries,
                      int z,
                      int n,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy,
                      float cmEnergy,
                      float weight);

QString buildEmittedParticleCMSpectraHtmlGemini(
    const AngularDistEntry &neutronEntry,
    const AngularDistEntry &protonEntry,
    const AngularDistEntry &alphaEntry,
    const AngularDistEntry &gammaEntry);

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
    double yieldSigmaTotalMb = -1.0,
    const AngularDistEntry &neutronEntry = AngularDistEntry(),
    const AngularDistEntry &protonEntry = AngularDistEntry(),
    const AngularDistEntry &alphaEntry = AngularDistEntry(),
    const AngularDistEntry &gammaEntry = AngularDistEntry());

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
        const AngularDistEntry &gammaEntry = AngularDistEntry(),
        double compoundExcitationMeV = 0.0,
        int compoundA = 1,
        int compoundZ = 0,
        double recoilBetaCN = 0.0,
        int mdir = 0,
        QWidget *parent = nullptr);

private slots:
    void link_clicked(const QUrl &url);

private:
    void openPlotWindow(bool plotAllTables, int tableIndex, int plotKind);
    void openCMSpectraPlotWindow();

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
    AngularDistEntry gammaEntryForPlots;

    double compoundExcitationForPlots = 0.0;
    int compoundAForPlots = 1;
    int compoundZForPlots = 0;
    double recoilBetaForPlots = 0.0;
    int mdirForPlots = 0;
};

#endif // GM_ANGULAR_DISTRIBUTION_H
