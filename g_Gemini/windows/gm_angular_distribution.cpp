// gm_angular_distribution.cpp
#include "gm_angular_distribution.h"

#include <QAction>
#include <QBoxLayout>
#include <QDesktopServices>
#include <QFile>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QIcon>
#include <QLabel>
#include <QPainter>
#include <QPen>
#include <QPixmap>
#include <QPolygonF>
#include <QPrintDialog>
#include <QPrinter>
#include <QPushButton>
#include <QScrollArea>
#include <QTextBrowser>
#include <QTextDocument>
#include <QTextStream>
#include <QToolBar>
#include <QVBoxLayout>
#include <QWidget>

#include "g_Gemini/source/CNucleus.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace {

constexpr int kMaxResiduesToPrint = 15;
constexpr double kDegToRad = 0.01745;
constexpr double kTwoPiPACE = 6.2832;
constexpr double kAMUMeV = 931.49432;

struct PaceSelectedResidue
{
    int z = 0;
    int n = 0;
    int a = 0;
    int count = 0;
    double percent = 0.0;
};

struct PaceAngularAccum
{
    double ELOW = 0.0;
    double DELE = 1.0;
    double DELANG = 1.0;

    std::array<double, 38> ANGLE1{};
    std::array<double, 38> ANGLE2{};

    std::array<double, 5> EWR{};
    std::array<double, 5> EWR1{};
    std::array<double, 5> EWR2{};

    std::array<std::array<int, 38>, 33> NR{};
    std::array<std::array<int, 38>, 5> NRW{};
};

QString fmt0(double x)
{
    return QString::number(x, 'f', 0);
}

QString fmt1(double x)
{
    return QString::number(x, 'f', 1);
}

void setParticleEnergyInterval(PaceAngularAccum &acc)
{
    acc.DELE = 2.0;
}

QString nucleusLabelFromZNPlain(int z, int n)
{
    CNucleus nuc(z, z + n);
    QString s = nuc.getGName();
    s.replace("<sup>", "");
    s.replace("</sup>", "");
    s.replace("<font color=\"darkBlue\">", "");
    s.replace("</font>", "");
    s.replace("<font size=+1>", "");
    s.replace("</font>", "");
    return s.simplified();
}

int isum(const PaceAngularAccum &acc, int energyBin, int l1, int l2)
{
    int sum = 0;
    for (int l = l1; l <= l2; ++l)
        sum += acc.NR[energyBin][l];
    return sum;
}

void initPaceAccum(PaceAngularAccum &acc,
                   int compoundA,
                   double compoundExcitationMeV,
                   double recoilBetaCN,
                   int mdir)
{
    const double A = std::max(compoundA, 1);
    const double EREC = 460.0 * A * recoilBetaCN * recoilBetaCN;

    const int ILOW = int((1.0 - 0.4 * compoundExcitationMeV / A) * EREC + 0.001);
    acc.ELOW = std::max(double(ILOW), 0.0);

    const int IDELE = int((EREC - acc.ELOW) / 16.0 + 0.5);
    acc.DELE = std::max(double(IDELE), 1.0);

    acc.DELANG = 1.0 + 4.0 * double(mdir);
    if (acc.DELANG <= 0.0 || acc.DELANG > 5.0) acc.DELANG = 5.0;

    acc.EWR[1] = acc.ELOW;
    acc.EWR[2] = int(EREC);
    if (acc.EWR[2] < acc.EWR[1] + 1.0) acc.EWR[2] = acc.EWR[1] + 1.0;
    acc.EWR[3] = 2.0 * acc.EWR[2] - acc.ELOW;

    for (int m = 1; m <= 36; ++m)
    {
        const double am = double(m);
        acc.ANGLE1[m] = (am - 1.0) * acc.DELANG;
        acc.ANGLE2[m] = am * acc.DELANG;
    }
    acc.ANGLE1[37] = acc.ANGLE2[36];
    acc.ANGLE2[37] = acc.ANGLE2[36];

    acc.EWR1[1] = 0.0;
    acc.EWR2[1] = acc.EWR[1];
    acc.EWR1[2] = acc.EWR[1];
    acc.EWR2[2] = acc.EWR[2];
    acc.EWR1[3] = acc.EWR[2];
    acc.EWR2[3] = acc.EWR[3];
    acc.EWR1[4] = acc.EWR[3];
}

void addToPaceAccum(PaceAngularAccum &acc, double energyLabMeV, double angleLabDeg)
{
    int K = 1;
    if (energyLabMeV >= acc.ELOW)
        K = std::min(31, int((energyLabMeV - acc.ELOW) / acc.DELE) + 2);

    int L = std::min(37, int(angleLabDeg / acc.DELANG) + 1);
    if (L < 1) L = 1;

    acc.NR[K][L]++;
    acc.NR[32][L]++;

    int KW = 4;
    if (energyLabMeV <= acc.EWR[1]) KW = 1;
    else if (energyLabMeV <= acc.EWR[2]) KW = 2;
    else if (energyLabMeV <= acc.EWR[3]) KW = 3;

    acc.NRW[KW][L]++;
}

std::vector<PaceSelectedResidue> buildSelectedResiduesPACE(
    const std::map<std::pair<int, int>, AngularDistEntry> &entries,
    int nCascades,
    double lowLimitPercent,
    double highLimitPercent)
{
    std::vector<PaceSelectedResidue> allResidues;
    allResidues.reserve(entries.size());

    for (const auto &it : entries)
    {
        const AngularDistEntry &e = it.second;
        const int cnt = int(e.kineticEnergy.size());
        if (cnt <= 0) continue;

        PaceSelectedResidue r;
        r.z = e.z;
        r.n = e.n;
        r.a = e.z + e.n;
        r.count = cnt;
        r.percent = 100.0 * double(cnt) / double(std::max(1, nCascades));
        allResidues.push_back(r);
    }

    // PACE-like yield-table order: descending A, then descending Z
    std::sort(allResidues.begin(), allResidues.end(),
              [](const PaceSelectedResidue &lhs, const PaceSelectedResidue &rhs)
              {
                  if (lhs.a != rhs.a) return lhs.a > rhs.a;
                  return lhs.z > rhs.z;
              });

    std::vector<PaceSelectedResidue> selected;
    selected.reserve(kMaxResiduesToPrint);

    for (const auto &r : allResidues)
    {
        if (r.percent < lowLimitPercent) continue;
        if (r.percent > highLimitPercent) continue;
        if ((int)selected.size() >= kMaxResiduesToPrint) break;
        selected.push_back(r);
    }

    return selected;
}

QString buildVelocityLinePACE(const PaceAngularAccum &acc,
                              int particleA,
                              int inputMode,
                              bool isAll,
                              const QString &velocityLabel)
{
    if (inputMode != 1 || isAll) return QString();

    const double globSum = isum(acc, 32, 1, 37);
    if (globSum <= 2.0) return QString();

    const double A = std::max(particleA, 1);
    const double fff = 2.0 / (kAMUMeV * A);

    double s_vz = 0.0;

    for (int k = 1; k <= 37; ++k)
    {
        const double angleDeg = (k < 37) ? (acc.ANGLE2[k] + acc.ANGLE1[k]) / 2.0 : acc.ANGLE2[k];
        const double angle = angleDeg * kDegToRad;
        const double angcos = std::cos(angle);

        for (int K = 1; K <= 31; ++K)
        {
            if (acc.NR[K][k] <= 0) continue;

            double ener = 0.0;
            if (K == 1) ener = acc.ELOW;
            else if (K < 31) ener = (double(K) - 1.5) * acc.DELE + acc.ELOW;
            else ener = 30.0 * acc.DELE + acc.ELOW;

            const double v = std::sqrt(fff * ener);
            const double vz = v * angcos;
            s_vz += vz * acc.NR[K][k];
        }
    }

    const double vzm = s_vz / globSum;

    s_vz = 0.0;
    double s_vxy = 0.0;

    for (int k = 1; k <= 37; ++k)
    {
        const double angleDeg = (k < 37) ? (acc.ANGLE2[k] + acc.ANGLE1[k]) / 2.0 : acc.ANGLE2[k];
        const double angle = angleDeg * kDegToRad;
        const double angsin = std::sin(angle);
        const double angcos = std::cos(angle);

        for (int K = 1; K <= 31; ++K)
        {
            if (acc.NR[K][k] <= 0) continue;

            double ener = 0.0;
            if (K == 1) ener = acc.ELOW;
            else if (K < 31) ener = (double(K) - 1.5) * acc.DELE + acc.ELOW;
            else ener = 30.0 * acc.DELE + acc.ELOW;

            const double v = std::sqrt(fff * ener);
            const double vz = v * angcos;
            const double vxy = v * angsin;

            s_vz += (vz - vzm) * (vz - vzm) * acc.NR[K][k];
            s_vxy += vxy * vxy * acc.NR[K][k];
        }
    }

    const double dvz = std::sqrt(s_vz / globSum);
    const double dvxy = std::sqrt(s_vxy / globSum);

    QString line;
    line += "<p><em>" + velocityLabel.toHtmlEscaped() + "</em> V<sub>z</sub> = ";
    line += QString::number(vzm, 'e', 2);
    line += " (Rms = ";
    line += QString::number(dvz, 'e', 2);
    line += ") V<sub>xy</sub> = ";
    line += QString::number(dvxy, 'e', 2);
    line += "</p><p>&nbsp;</p>";
    return line;
}

void appendAngleHeader18(QString &html, const PaceAngularAccum &acc, int l1, int l2, bool includeAboveLabel)
{
    html += "<table border=\"1\" cellspacing=\"0\">";
    html += "<tr><th>Energy Range</th><th colspan=\"18\">Angular range (deg)</th></tr>";
    html += "<tr><th>(MeV)</th>";

    for (int k = l1; k <= l2; ++k)
        html += "<td>" + fmt0(acc.ANGLE1[k]) + "</td>";

    if (includeAboveLabel)
        html += "<td>Above</td>";

    html += "</tr><tr><td></td>";

    for (int k = l1; k <= l2; ++k)
        html += "<td>" + fmt0(acc.ANGLE2[k]) + "</td>";

    if (includeAboveLabel)
        html += "<td>" + fmt0(acc.ANGLE2[36]) + "</td>";

    html += "</tr>";
}

void appendMainPaceTable(QString &html,
                         const PaceAngularAccum &acc,
                         int particleA,
                         double sigmaTotalMb,
                         int nCascades,
                         int inputMode,
                         bool isAll,
                         int plotIndex,
                         const QString &velocityLabel)
{
    html += buildVelocityLinePACE(acc, particleA, inputMode, isAll, velocityLabel);

    appendAngleHeader18(html, acc, 1, 18, false);

    double E2 = 0.0;
    if (acc.ELOW >= 0.01)
    {
        E2 = acc.ELOW;
        html += "<tr><td>Below " + fmt0(E2) + "</td>";
        for (int k = 1; k <= 18; ++k)
            html += (acc.NR[1][k] > 0) ? "<td align=\"center\">" + QString::number(acc.NR[1][k]) + "</td>" : "<td></td>";
        html += "</tr>";
    }

    for (int K = 2; K <= 30; ++K)
    {
        const double E1 = (double(K) - 2.0) * acc.DELE + acc.ELOW;
        E2 = E1 + acc.DELE;

        if (isum(acc, K, 1, 18) > 0)
        {
            html += "<tr><td>" + fmt1(E1) + " - " + fmt1(E2) + "</td>";
            for (int k = 1; k <= 18; ++k)
                html += (acc.NR[K][k] > 0) ? "<td align=\"center\">" + QString::number(acc.NR[K][k]) + "</td>" : "<td></td>";
            html += "</tr>";
        }
    }

    html += "<tr><td>Above " + fmt0(E2) + "</td>";
    for (int k = 1; k <= 18; ++k)
        html += (acc.NR[31][k] > 0) ? "<td align=\"center\">" + QString::number(acc.NR[31][k]) + "</td>" : "<td></td>";
    html += "</tr>";

    std::array<double, 37> DSIG{};
    const double FAC = sigmaTotalMb / (kTwoPiPACE * double(std::max(1, nCascades)));
    for (int i = 1; i <= 36; ++i)
    {
        const double TET = (double(i) - 0.5) * kDegToRad;
        DSIG[i] = FAC * double(acc.NR[32][i]) / (std::sin(TET) * kDegToRad);
    }

    html += "<tr><td>Total</td>";
    for (int k = 1; k <= 18; ++k)
        html += (acc.NR[32][k] > 0) ? "<td align=\"center\">" + QString::number(acc.NR[32][k]) + "</td>" : "<td></td>";
    html += "</tr>";

    html += "<tr><td>d&sigma;/d&Omega;</td>";
    for (int k = 1; k <= 18; ++k)
        html += (DSIG[k] > 0.0) ? "<td>" + QString::number(DSIG[k], 'g', 2) + "</td>" : "<td>0.00</td>";
    html += "</tr>";

    for (int i = 1; i <= 4; ++i)
    {
        html += (i != 4)
                    ? "<tr><td>" + fmt0(acc.EWR1[i]) + " - " + fmt0(acc.EWR2[i]) + "</td>"
                    : "<tr><td>Above " + fmt0(acc.EWR1[4]) + "</td>";

        for (int k = 1; k <= 18; ++k)
            html += (acc.NRW[i][k] > 0) ? "<td align=\"center\">" + QString::number(acc.NRW[i][k]) + "</td>" : "<td></td>";
        html += "</tr>";
    }

    html += "</table>";
    html += "<p>";
    html += "<a href=\"gemini://plot_table/" + QString::number(plotIndex) + "\">Plot E vs &theta;</a> &nbsp; ";
    html += "<a href=\"gemini://plot_ntheta/" + QString::number(plotIndex) + "\">Plot N vs &theta;</a> &nbsp; ";
    html += "<a href=\"gemini://plot_ne/" + QString::number(plotIndex) + "\">Plot N vs Energy</a> &nbsp; ";
    html += "<a href=\"gemini://plot_dcsdtheta/" + QString::number(plotIndex) + "\">Plot d&sigma;/d&theta; vs &theta;</a>";
    html += "</p>";
    html += "<p>&nbsp;</p>";
}

void appendSecondPaceTableIfNeeded(QString &html,
                                   const PaceAngularAccum &acc,
                                   int plotIndex)
{
    int IS = 0;
    for (int IJ = 1; IJ <= 32; ++IJ)
        for (int LQ = 19; LQ <= 37; ++LQ)
            IS += acc.NR[IJ][LQ];

    if (IS == 0) return;

    appendAngleHeader18(html, acc, 19, 36, true);

    double E2 = 0.0;
    if (acc.ELOW >= 0.01)
    {
        E2 = acc.ELOW;
        html += "<tr><td>Below " + fmt0(E2) + "</td>";
        for (int k = 19; k <= 36; ++k)
            html += (acc.NR[1][k] > 0) ? "<td align=\"center\">" + QString::number(acc.NR[1][k]) + "</td>" : "<td></td>";
        html += "<td></td></tr>";
    }

    for (int K = 2; K <= 30; ++K)
    {
        const double E1 = (double(K) - 2.0) * acc.DELE + acc.ELOW;
        E2 = E1 + acc.DELE;

        if (isum(acc, K, 19, 37) > 0)
        {
            html += "<tr><td>" + fmt1(E1) + " - " + fmt1(E2) + "</td>";
            for (int k = 19; k <= 36; ++k)
                html += (acc.NR[K][k] > 0) ? "<td align=\"center\">" + QString::number(acc.NR[K][k]) + "</td>" : "<td></td>";
            html += "<td></td></tr>";
        }
    }

    html += "<tr><td>Above " + fmt0(E2) + "</td>";
    for (int k = 19; k <= 36; ++k)
        html += (acc.NR[31][k] > 0) ? "<td align=\"center\">" + QString::number(acc.NR[31][k]) + "</td>" : "<td></td>";
    html += "<td></td></tr>";

    html += "<tr><th>Total</th>";
    for (int k = 19; k <= 36; ++k)
        html += (acc.NR[32][k] > 0) ? "<th>" + QString::number(acc.NR[32][k]) + "</th>" : "<th></th>";
    html += "<th></th></tr>";

    html += "<tr><td>d&sigma;/d&Omega;</td>";
    for (int k = 19; k <= 36; ++k)
    {
        const double TET = (double(k) - 0.5) * kDegToRad;
        const double DS = (std::sin(TET) > 0.0) ? double(acc.NR[32][k]) / (std::sin(TET) * kDegToRad) : 0.0;
        html += (DS > 0.0) ? "<td align=\"center\">" + QString::number(DS, 'g', 2) + "</td>" : "<td></td>";
    }
    html += "<td></td></tr>";

    for (int i = 1; i <= 3; ++i)
    {
        html += "<tr><td>" + fmt0(acc.EWR1[i]) + " - " + fmt0(acc.EWR2[i]) + "</td>";
        for (int k = 19; k <= 36; ++k)
            html += (acc.NRW[i][k] > 0) ? "<td align=\"center\">" + QString::number(acc.NRW[i][k]) + "</td>" : "<td></td>";
        html += "<td></td></tr>";
    }

    html += "<tr><td>Above " + fmt0(acc.EWR1[4]) + "</td>";
    for (int k = 19; k <= 36; ++k)
        html += (acc.NRW[4][k] > 0) ? "<td>" + QString::number(acc.NRW[4][k]) + "</td>" : "<td></td>";
    html += "<td></td></tr></table>";
    html += "<p>";
    html += "<a href=\"gemini://plot_table/" + QString::number(plotIndex) + "\">Plot E vs &theta;</a> &nbsp; ";
    html += "<a href=\"gemini://plot_ntheta/" + QString::number(plotIndex) + "\">Plot N vs &theta;</a> &nbsp; ";
    html += "<a href=\"gemini://plot_ne/" + QString::number(plotIndex) + "\">Plot N vs E</a> &nbsp; ";
    html += "<a href=\"gemini://plot_dcsdtheta/" + QString::number(plotIndex) + "\">Plot d&sigma;/d&theta; vs &theta;</a>";
    html += "</p>";
    html += "<p>&nbsp;</p>";
}

QString buildHeaderForResiduePACE(int displayIndex, int z, int n)
{
    const QString nuc = nucleusLabelFromZNPlain(z, n);
    return QString("<br><h3>%1. Energy and angular distribution of residual nucleus Z = "
                   "<span style=\"color:blue\">%2</span> and N = <span style=\"color:blue\">%3</span> "
                   "(<span style=\"color:blue\">%4</span>)</h3>")
        .arg(displayIndex)
        .arg(z)
        .arg(n)
        .arg(nuc);
}

QString buildHeaderAllPACE(int displayIndex)
{
    return QString("<br><h3>%1. Energy and angular distribution of ALL residual nuclei</h3>")
        .arg(displayIndex);
}

QString buildHeaderForParticlePACE(int displayIndex, const QString &particleLabel)
{
    return QString("<br><h3>%1. Energy and angular distribution of emitted %2</h3>")
        .arg(displayIndex)
        .arg(particleLabel.toHtmlEscaped());
}

bool hasAngularSamples(const AngularDistEntry &entry)
{
    return !entry.kineticEnergy.empty() && !entry.thetaDeg.empty();
}

bool hasAngularSamples(const std::map<std::pair<int, int>, AngularDistEntry> &entries)
{
    for (const auto &it : entries)
    {
        if (hasAngularSamples(it.second)) return true;
    }
    return false;
}

void appendEntryToAccum(PaceAngularAccum &acc, const AngularDistEntry &entry)
{
    const int n = std::min(entry.kineticEnergy.size(), entry.thetaDeg.size());
    for (int i = 0; i < n; ++i)
        addToPaceAccum(acc, entry.kineticEnergy[i], entry.thetaDeg[i]);
}

struct PlotTableEntry
{
    QString label;
    AngularDistEntry entry;
    PaceAngularAccum acc;
};

std::vector<PlotTableEntry> buildPlotTableListPACE(
    const std::map<std::pair<int, int>, AngularDistEntry> &entries,
    const AngularDistEntry &neutronEntry,
    const AngularDistEntry &protonEntry,
    const AngularDistEntry &alphaEntry,
    const AngularDistEntry &gammaEntry,
    int nEvents,
    double lowLimitPercent,
    double highLimitPercent,
    double compoundExcitationMeV,
    int compoundA,
    double recoilBetaCN,
    int mdir)
{
    const std::vector<PaceSelectedResidue> selected =
        buildSelectedResiduesPACE(entries, nEvents, lowLimitPercent, highLimitPercent);

    PaceAngularAccum allAcc;
    initPaceAccum(allAcc, compoundA, compoundExcitationMeV, recoilBetaCN, mdir);

    std::vector<std::pair<PaceSelectedResidue, PaceAngularAccum>> selectedAcc;
    selectedAcc.reserve(selected.size());
    for (const auto &r : selected)
    {
        PaceAngularAccum acc;
        initPaceAccum(acc, compoundA, compoundExcitationMeV, recoilBetaCN, mdir);
        selectedAcc.push_back({r, acc});
    }

    for (const auto &it : entries)
    {
        const AngularDistEntry &e = it.second;
        const int m = std::min(e.kineticEnergy.size(), e.thetaDeg.size());

        int selectedIndex = -1;
        for (int i = 0; i < (int)selected.size(); ++i)
        {
            if (selected[i].z == e.z && selected[i].n == e.n)
            {
                selectedIndex = i;
                break;
            }
        }

        for (int i = 0; i < m; ++i)
        {
            addToPaceAccum(allAcc, e.kineticEnergy[i], e.thetaDeg[i]);
            if (selectedIndex >= 0)
                addToPaceAccum(selectedAcc[selectedIndex].second, e.kineticEnergy[i], e.thetaDeg[i]);
        }
    }

    std::vector<PlotTableEntry> tables;
    tables.reserve(selected.size() + 5);

    int idx = 1;
    for (const auto &item : selectedAcc)
    {
        const PaceSelectedResidue &r = item.first;
        auto it = entries.find({r.z, r.n});
        if (it == entries.end()) continue;

        PlotTableEntry t;
        t.label = QString("%1. Energy and angular distribution of residual nucleus Z = %2 and N = %3 (%4)")
                      .arg(idx)
                      .arg(r.z)
                      .arg(r.n)
                      .arg(nucleusLabelFromZNPlain(r.z, r.n));
        t.entry = it->second;
        t.acc = item.second;
        tables.push_back(t);
        idx++;
    }

    AngularDistEntry all;
    all.z = -1;
    all.n = -1;
    for (const auto &it : entries)
    {
        const AngularDistEntry &e = it.second;
        all.kineticEnergy.insert(all.kineticEnergy.end(), e.kineticEnergy.begin(), e.kineticEnergy.end());
        all.thetaDeg.insert(all.thetaDeg.end(), e.thetaDeg.begin(), e.thetaDeg.end());
        all.vz.insert(all.vz.end(), e.vz.begin(), e.vz.end());
        all.vxy.insert(all.vxy.end(), e.vxy.begin(), e.vxy.end());
        all.cmEnergy.insert(all.cmEnergy.end(), e.cmEnergy.begin(), e.cmEnergy.end());
    }

    PlotTableEntry tAll;
    tAll.label = QString("%1. Energy and angular distribution of ALL residual nuclei").arg(idx);
    tAll.entry = all;
    tAll.acc = allAcc;
    tables.push_back(tAll);
    idx++;

    auto appendParticlePlot = [&](const QString &particleLabel, const AngularDistEntry &entry)
    {
        if (!hasAngularSamples(entry)) return;

        PlotTableEntry table;
        table.label = QString("%1. Energy and angular distribution of emitted %2")
                          .arg(idx)
                          .arg(particleLabel);
        table.entry = entry;
        initPaceAccum(table.acc, compoundA, compoundExcitationMeV, recoilBetaCN, mdir);
        setParticleEnergyInterval(table.acc);
        appendEntryToAccum(table.acc, entry);
        tables.push_back(table);
        idx++;
    };

    appendParticlePlot("neutrons", neutronEntry);
    appendParticlePlot("protons", protonEntry);
    appendParticlePlot("alpha particles", alphaEntry);
    appendParticlePlot("gamma particles", gammaEntry);

    return tables;
}

enum PlotKind
{
    PlotEnergyVsTheta = 0,
    PlotCountsVsTheta = 1,
    PlotCountsVsEnergy = 2,
    PlotCrossSectionVsTheta = 3
};

class ScatterPlotWidget : public QWidget
{
public:
    ScatterPlotWidget(const AngularDistEntry &entry,
                      const PaceAngularAccum &acc,
                      QWidget *parent = nullptr)
        : QWidget(parent), m_entry(entry), m_acc(acc)
    {
        setMinimumSize(980, 560);
        resize(980, 560);
    }

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter p(this);
        p.setRenderHint(QPainter::Antialiasing, true);
        p.fillRect(rect(), Qt::white);

        const int left = 90;
        const int right = 120;
        const int top = 35;
        const int bottom = 80;

        QRect plotRect(left, top, width() - left - right, height() - top - bottom);

        p.setPen(QPen(Qt::black, 1));
        p.drawRect(plotRect);

        if (m_entry.thetaDeg.empty() || m_entry.kineticEnergy.empty())
        {
            p.drawText(plotRect, Qt::AlignCenter, "No data");
            return;
        }

        const int n = int(std::min(m_entry.thetaDeg.size(), m_entry.kineticEnergy.size()));

        double rawMaxTheta = 0.0;
        for (int i = 0; i < n; ++i)
            rawMaxTheta = std::max(rawMaxTheta, double(m_entry.thetaDeg[i]));

        double maxTheta = 18.0;
        if (rawMaxTheta > 18.0) maxTheta = 36.0;
        if (rawMaxTheta > 36.0) maxTheta = std::ceil(rawMaxTheta);

        const double minTheta = 0.0;
        const int thetaBins = int(std::ceil(maxTheta / m_acc.DELANG));
        const int energyBins = 31;

        std::vector<std::vector<int>> counts(
            energyBins,
            std::vector<int>(thetaBins, 0)
            );

        int maxCount = 0;
        int highestPopulatedEnergyBin = 0;

        for (int i = 0; i < n; ++i)
        {
            double theta = m_entry.thetaDeg[i];
            double energy = m_entry.kineticEnergy[i];

            if (theta < minTheta || theta > maxTheta) continue;

            int thetaBin = int(theta / m_acc.DELANG);
            int energyBin = 0;
            if (energy >= m_acc.ELOW)
                energyBin = std::min(30, int((energy - m_acc.ELOW) / m_acc.DELE) + 1);

            if (thetaBin >= thetaBins) thetaBin = thetaBins - 1;
            if (energyBin >= energyBins) energyBin = energyBins - 1;

            if (thetaBin < 0 || energyBin < 0) continue;

            counts[energyBin][thetaBin]++;
            maxCount = std::max(maxCount, counts[energyBin][thetaBin]);
            highestPopulatedEnergyBin = std::max(highestPopulatedEnergyBin, energyBin);
        }

        if (maxCount <= 0)
        {
            p.drawText(plotRect, Qt::AlignCenter, "No histogram data");
            return;
        }

        const double cellW = double(plotRect.width()) / thetaBins;
        const int visibleEnergyBins = std::max(1, highestPopulatedEnergyBin + 1);
        const double cellH = double(plotRect.height()) / visibleEnergyBins;

        // Draw 2D histogram cells
        for (int e = 0; e < visibleEnergyBins; ++e)
        {
            for (int t = 0; t < thetaBins; ++t)
            {
                int count = counts[e][t];

                int x = plotRect.left() + int(t * cellW);
                int y = plotRect.bottom() - int((e + 1) * cellH);

                QRect cell(x, y, int(cellW) + 1, int(cellH) + 1);

                if (count == 0)
                {
                    p.fillRect(cell, QColor(250, 250, 250));
                }
                else
                {
                    double frac = double(count) / double(maxCount);

                    int red = int(255 * frac);
                    int green = int(230 * (1.0 - frac));
                    int blue = 255 - int(220 * frac);

                    p.fillRect(cell, QColor(red, green, blue));
                }

                p.setPen(QPen(QColor(220, 220, 220), 1));
                p.drawRect(cell);

                if (count > 0)
                {
                    p.setPen(Qt::black);
                    p.drawText(cell, Qt::AlignCenter, QString::number(count));
                }
            }
        }

        // Axes
        p.setPen(QPen(Qt::black, 2));
        p.drawLine(plotRect.bottomLeft(), plotRect.bottomRight());
        p.drawLine(plotRect.bottomLeft(), plotRect.topLeft());

        QFont axisFont = p.font();
        axisFont.setPointSize(10);
        p.setFont(axisFont);

        // X-axis theta labels
        int thetaTickStep = 1;
        if (maxTheta > 18.0) thetaTickStep = 2;
        if (maxTheta > 36.0) thetaTickStep = 5;

        for (int i = 0; i <= int(maxTheta); i += thetaTickStep)
        {
            int x = plotRect.left() + int((double(i) / maxTheta) * plotRect.width());
            p.drawLine(x, plotRect.bottom(), x, plotRect.bottom() + 6);
            p.drawText(x - 15,
                       plotRect.bottom() + 24,
                       30,
                       18,
                       Qt::AlignCenter,
                       QString::number(i));
        }

        // Y-axis energy labels: same boundaries as the table, clipped at the highest populated bin.
        for (int e = 0; e <= visibleEnergyBins; ++e)
        {
            const double frac = double(e) / double(visibleEnergyBins);
            double energyValue = 0.0;
            if (e == 0) energyValue = 0.0;
            else if (e == 1) energyValue = m_acc.ELOW;
            else energyValue = m_acc.ELOW + double(e - 1) * m_acc.DELE;

            const int y = plotRect.bottom() - int(frac * plotRect.height());

            if (e == 0 || e == 1 || e == visibleEnergyBins || (e - 1) % 5 == 0)
            {
                p.drawLine(plotRect.left() - 6, y, plotRect.left(), y);
                p.drawText(5,
                           y - 10,
                           left - 15,
                           20,
                           Qt::AlignRight | Qt::AlignVCenter,
                           QString::number(energyValue, 'f', 1));
            }
        }

        // Axis titles
        p.drawText(plotRect.left(),
                   height() - 35,
                   plotRect.width(),
                   24,
                   Qt::AlignCenter,
                   "θ range (deg)");

        p.save();
        p.translate(30, plotRect.top() + plotRect.height() / 2);
        p.rotate(-90);
        p.drawText(QRect(-plotRect.height() / 2, -20, plotRect.height(), 20),
                   Qt::AlignCenter,
                   "Energy range (MeV)");
        p.restore();

        // Legend
        const int legendX = plotRect.right() + 25;
        const int legendY = plotRect.top() + 20;
        const int legendW = 22;
        const int legendH = 180;

        for (int i = 0; i < legendH; ++i)
        {
            double frac = 1.0 - double(i) / double(legendH - 1);

            int red = int(255 * frac);
            int green = int(230 * (1.0 - frac));
            int blue = 255 - int(220 * frac);

            p.setPen(QColor(red, green, blue));
            p.drawLine(legendX, legendY + i, legendX + legendW, legendY + i);
        }

        p.setPen(Qt::black);
        p.drawRect(legendX, legendY, legendW, legendH);

        p.drawText(legendX + 30,
                   legendY,
                   70,
                   20,
                   Qt::AlignLeft,
                   QString::number(maxCount));

        p.drawText(legendX + 30,
                   legendY + legendH - 15,
                   70,
                   20,
                   Qt::AlignLeft,
                   "0");

        p.drawText(legendX - 10,
                   legendY + legendH + 15,
                   90,
                   40,
                   Qt::AlignCenter,
                   "Counts");
    }

private:
    AngularDistEntry m_entry;
    PaceAngularAccum m_acc;
};


struct GaussianOverlayStats
{
    double amplitude = 0.0;        // A in y = A exp(-0.5 ((x - mean) / sigma)^2)
    double mean = 0.0;             // Gaussian center
    double sigma = 0.0;            // Gaussian width
    double chiSquare = 0.0;        // Chi-square
    double reducedChiSquare = 0.0; // Chi-square / degrees of freedom
    double area = 0.0;             // A * sigma * sqrt(2pi)
    bool valid = false;
};

double gaussianValue(double amplitude, double mean, double sigma, double x)
{
    if (sigma <= 1.0e-12) return 0.0;

    const double z = (x - mean) / sigma;
    return amplitude * std::exp(-0.5 * z * z);
}

void computeGaussianQualityAndArea(GaussianOverlayStats &stats,
                                   const std::vector<double> &xCenters,
                                   const std::vector<double> &yValues)
{
    stats.chiSquare = 0.0;
    stats.reducedChiSquare = 0.0;
    stats.area = 0.0;

    if (!stats.valid || stats.sigma <= 1.0e-12) return;

    int usedPoints = 0;

    for (int i = 0; i < int(yValues.size()); ++i)
    {
        const double observed = std::max(0.0, yValues[i]);
        const double expected = gaussianValue(stats.amplitude,
                                              stats.mean,
                                              stats.sigma,
                                              xCenters[i]);

        const double denom = std::max(1.0, expected);
        const double diff = observed - expected;

        stats.chiSquare += (diff * diff) / denom;
        usedPoints++;
    }

    const int dof = std::max(1, usedPoints - 3);
    stats.reducedChiSquare = stats.chiSquare / double(dof);

    const double sqrtTwoPi = std::sqrt(2.0 * 3.14159265358979323846);
    stats.area = stats.amplitude * stats.sigma * sqrtTwoPi;
}

GaussianOverlayStats computeGaussianOverlayStats(const std::vector<double> &xCenters,
                                                 const std::vector<double> &yValues)
{
    GaussianOverlayStats stats;

    if (xCenters.size() != yValues.size() || xCenters.size() < 3)
        return stats;

    double sumW = 0.0;
    double sumX = 0.0;
    double maxY = 0.0;

    for (int i = 0; i < int(yValues.size()); ++i)
    {
        const double w = std::max(0.0, yValues[i]);
        if (w <= 0.0) continue;

        sumW += w;
        sumX += w * xCenters[i];
        maxY = std::max(maxY, w);
    }

    if (sumW <= 0.0 || maxY <= 0.0)
        return stats;

    const double xMin = xCenters.front();
    const double xMax = xCenters.back();
    const double xRange = std::max(1.0e-9, xMax - xMin);

    double A = maxY;
    double mu = sumX / sumW;

    double variance = 0.0;
    for (int i = 0; i < int(yValues.size()); ++i)
    {
        const double w = std::max(0.0, yValues[i]);
        if (w <= 0.0) continue;

        const double dx = xCenters[i] - mu;
        variance += w * dx * dx;
    }

    double sigma = std::sqrt(std::max(variance / sumW, 1.0e-12));
    sigma = std::max(sigma, xRange / 1000.0);

    auto sseFor = [&](double testA, double testMu, double testSigma)
    {
        if (testA <= 0.0 || testSigma <= 1.0e-12)
            return 1.0e300;

        double sse = 0.0;

        for (int i = 0; i < int(yValues.size()); ++i)
        {
            const double yFit = gaussianValue(testA, testMu, testSigma, xCenters[i]);
            const double r = yValues[i] - yFit;
            sse += r * r;
        }

        return sse;
    };

    double lambda = 1.0e-3;
    double bestSse = sseFor(A, mu, sigma);

    for (int iter = 0; iter < 100; ++iter)
    {
        double jtj[3][3] = {};
        double jtr[3] = {};

        for (int i = 0; i < int(yValues.size()); ++i)
        {
            const double x = xCenters[i];
            const double y = yValues[i];

            const double z = (x - mu) / sigma;
            const double e = std::exp(-0.5 * z * z);
            const double yFit = A * e;
            const double r = y - yFit;

            const double dA = e;
            const double dMu = A * e * (z / sigma);
            const double dSigma = A * e * ((z * z) / sigma);

            const double J[3] = {dA, dMu, dSigma};

            for (int a = 0; a < 3; ++a)
            {
                jtr[a] += J[a] * r;
                for (int b = 0; b < 3; ++b)
                    jtj[a][b] += J[a] * J[b];
            }
        }

        for (int d = 0; d < 3; ++d)
            jtj[d][d] += lambda;

        double M[3][4] =
            {
                {jtj[0][0], jtj[0][1], jtj[0][2], jtr[0]},
                {jtj[1][0], jtj[1][1], jtj[1][2], jtr[1]},
                {jtj[2][0], jtj[2][1], jtj[2][2], jtr[2]}
            };

        bool singular = false;

        for (int col = 0; col < 3; ++col)
        {
            int pivot = col;
            for (int row = col + 1; row < 3; ++row)
            {
                if (std::fabs(M[row][col]) > std::fabs(M[pivot][col]))
                    pivot = row;
            }

            if (std::fabs(M[pivot][col]) < 1.0e-18)
            {
                singular = true;
                break;
            }

            if (pivot != col)
            {
                for (int k = col; k < 4; ++k)
                    std::swap(M[col][k], M[pivot][k]);
            }

            const double div = M[col][col];
            for (int k = col; k < 4; ++k)
                M[col][k] /= div;

            for (int row = 0; row < 3; ++row)
            {
                if (row == col) continue;

                const double factor = M[row][col];
                for (int k = col; k < 4; ++k)
                    M[row][k] -= factor * M[col][k];
            }
        }

        if (singular)
            break;

        const double dA = M[0][3];
        const double dMu = M[1][3];
        const double dSigma = M[2][3];

        double newA = A + dA;
        double newMu = mu + dMu;
        double newSigma = sigma + dSigma;

        newA = std::max(0.0, newA);
        newMu = std::max(xMin - xRange, std::min(xMax + xRange, newMu));
        newSigma = std::max(xRange / 1000.0, std::min(xRange * 10.0, newSigma));

        const double newSse = sseFor(newA, newMu, newSigma);

        if (newSse < bestSse)
        {
            A = newA;
            mu = newMu;
            sigma = newSigma;
            bestSse = newSse;
            lambda *= 0.3;

            if (std::fabs(dA) < 1.0e-9 &&
                std::fabs(dMu) < 1.0e-9 &&
                std::fabs(dSigma) < 1.0e-9)
            {
                break;
            }
        }
        else
        {
            lambda *= 10.0;
        }

        if (lambda > 1.0e12)
            break;
    }

    stats.amplitude = A;
    stats.mean = mu;
    stats.sigma = sigma;
    stats.valid = std::isfinite(stats.amplitude) &&
                  std::isfinite(stats.mean) &&
                  std::isfinite(stats.sigma) &&
                  stats.amplitude > 0.0 &&
                  stats.sigma > 1.0e-12;

    computeGaussianQualityAndArea(stats, xCenters, yValues);

    return stats;
}

QString formatGaussianNumber(double value)
{
    if (!std::isfinite(value))
        return "nan";

    if (std::fabs(value) >= 1000.0 ||
        (std::fabs(value) > 0.0 && std::fabs(value) < 0.01))
    {
        return QString::number(value, 'e', 2);
    }

    return QString::number(value, 'g', 4);
}

class OneDPlotWidget : public QWidget
{
public:
    OneDPlotWidget(const AngularDistEntry &entry,
                   const PaceAngularAccum &acc,
                   int plotKind,
                   double sigmaTotal,
                   int nEvents,
                   QWidget *parent = nullptr)
        : QWidget(parent),
        m_entry(entry),
        m_acc(acc),
        m_plotKind(plotKind),
        m_sigmaTotal(sigmaTotal),
        m_nEvents(std::max(1, nEvents))
    {
        setMinimumSize(980, 560);
        resize(980, 560);
    }

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter p(this);
        p.setRenderHint(QPainter::Antialiasing, true);
        p.fillRect(rect(), Qt::white);

        const int left = 90;
        const int right = 60;
        const int top = 45;
        const int bottom = 85;

        QRect plotRect(left, top, width() - left - right, height() - top - bottom);

        p.setPen(QPen(Qt::black, 1));
        p.drawRect(plotRect);

        const int n = int(std::min(m_entry.thetaDeg.size(), m_entry.kineticEnergy.size()));

        if (n <= 0)
        {
            p.drawText(plotRect, Qt::AlignCenter, "No data");
            return;
        }

        std::vector<double> values;
        std::vector<double> xCenters;

        QString title;
        QString xTitle;
        QString yTitle;

        double minX = 0.0;
        double maxX = 1.0;
        int bins = 1;

        if (m_plotKind == PlotCountsVsTheta || m_plotKind == PlotCrossSectionVsTheta)
        {
            double rawMaxTheta = 0.0;

            for (int i = 0; i < n; ++i)
                rawMaxTheta = std::max(rawMaxTheta, double(m_entry.thetaDeg[i]));

            maxX = 18.0;
            if (rawMaxTheta > 18.0) maxX = 36.0;
            if (rawMaxTheta > 36.0) maxX = std::ceil(rawMaxTheta);

            minX = 0.0;
            bins = std::max(1, int(maxX));

            const double binWidth = (maxX - minX) / double(bins);

            std::vector<int> counts(bins, 0);

            for (int i = 0; i < n; ++i)
            {
                const double theta = m_entry.thetaDeg[i];

                if (theta < minX || theta > maxX)
                    continue;

                int bin = int((theta - minX) / (maxX - minX) * bins);

                if (bin >= bins) bin = bins - 1;
                if (bin >= 0) counts[bin]++;
            }

            values.resize(bins, 0.0);
            xCenters.resize(bins, 0.0);

            for (int i = 0; i < bins; ++i)
                xCenters[i] = minX + (double(i) + 0.5) * binWidth;

            for (int i = 0; i < bins; ++i)
            {
                if (m_plotKind == PlotCountsVsTheta)
                {
                    values[i] = double(counts[i]);
                }
                else
                {
                    values[i] = m_sigmaTotal * double(counts[i]) /
                                (double(m_nEvents) * binWidth);
                }
            }

            xTitle = "θ (deg)";

            if (m_plotKind == PlotCountsVsTheta)
            {
                title = "N = f(θ)";
                yTitle = "Counts N";
            }
            else
            {
                title = "dσ/dθ = f(θ)";
                yTitle = "dσ/dθ (mb/deg)";
            }
        }
        else if (m_plotKind == PlotCountsVsEnergy)
        {
            minX = 0.0;
            maxX = m_acc.ELOW + 29.0 * m_acc.DELE;
            bins = 31;

            values.assign(bins, 0.0);
            xCenters.assign(bins, 0.0);

            for (int i = 0; i < bins; ++i)
            {
                if (i == 0)
                    xCenters[i] = (m_acc.ELOW > 0.0) ? 0.5 * m_acc.ELOW : 0.0;
                else if (i < bins - 1)
                    xCenters[i] = m_acc.ELOW + (double(i) - 0.5) * m_acc.DELE;
                else
                    xCenters[i] = m_acc.ELOW + 29.5 * m_acc.DELE;
            }

            for (int i = 0; i < n; ++i)
            {
                const double energy = m_entry.kineticEnergy[i];

                int bin = 0;

                if (energy >= m_acc.ELOW)
                    bin = std::min(30, int((energy - m_acc.ELOW) / m_acc.DELE) + 1);

                if (bin >= 0 && bin < bins)
                    values[bin] += 1.0;
            }

            title = "N = f(E)";
            xTitle = "Energy (MeV)";
            yTitle = "Counts N";
        }

        double maxY = 0.0;

        for (double v : values)
            maxY = std::max(maxY, v);

        const bool showGaussian =
            (m_plotKind == PlotCountsVsTheta || m_plotKind == PlotCountsVsEnergy);

        GaussianOverlayStats gaussianStats;

        if (showGaussian)
            gaussianStats = computeGaussianOverlayStats(xCenters, values);

        if (gaussianStats.valid)
            maxY = std::max(maxY, gaussianStats.amplitude * 1.15);

        if (maxY <= 0.0)
        {
            p.drawText(plotRect, Qt::AlignCenter, "No histogram data");
            return;
        }

        maxY *= 1.10;

        auto mapY = [&](double yValue)
        {
            return plotRect.bottom() - int((yValue / maxY) * plotRect.height());
        };

        QFont titleFont = p.font();
        titleFont.setPointSize(12);
        titleFont.setBold(true);
        p.setFont(titleFont);
        p.setPen(Qt::black);

        p.drawText(0,
                   8,
                   width(),
                   24,
                   Qt::AlignCenter,
                   title);

        QFont axisFont = p.font();
        axisFont.setPointSize(9);
        axisFont.setBold(false);
        p.setFont(axisFont);

        p.setPen(QPen(Qt::black, 2));
        p.drawLine(plotRect.bottomLeft(), plotRect.bottomRight());
        p.drawLine(plotRect.bottomLeft(), plotRect.topLeft());

        const int yTicks = 5;

        for (int i = 0; i <= yTicks; ++i)
        {
            const double frac = double(i) / double(yTicks);
            const double yValue = frac * maxY;
            const int y = plotRect.bottom() - int(frac * plotRect.height());

            p.setPen(QPen(QColor(220, 220, 220), 1));
            p.drawLine(plotRect.left(), y, plotRect.right(), y);

            p.setPen(Qt::black);
            p.drawLine(plotRect.left() - 6, y, plotRect.left(), y);

            QString yText;

            if (m_plotKind == PlotCountsVsTheta || m_plotKind == PlotCountsVsEnergy)
                yText = QString::number(int(std::round(yValue)));
            else
                yText = QString::number(yValue, 'g', 3);

            p.drawText(5,
                       y - 10,
                       left - 15,
                       20,
                       Qt::AlignRight | Qt::AlignVCenter,
                       yText);
        }

        int xTickCount = 6;

        if (m_plotKind == PlotCountsVsTheta || m_plotKind == PlotCrossSectionVsTheta)
            xTickCount = int(maxX / 3.0);

        xTickCount = std::max(3, xTickCount);

        for (int i = 0; i <= xTickCount; ++i)
        {
            const double frac = double(i) / double(xTickCount);
            const double xValue = minX + frac * (maxX - minX);
            const int x = plotRect.left() + int(frac * plotRect.width());

            p.setPen(Qt::black);
            p.drawLine(x, plotRect.bottom(), x, plotRect.bottom() + 6);

            p.drawText(x - 30,
                       plotRect.bottom() + 10,
                       60,
                       20,
                       Qt::AlignCenter,
                       QString::number(xValue, 'g', 3));
        }

        const double barW = double(plotRect.width()) / double(bins);

        for (int i = 0; i < bins; ++i)
        {
            const int x = plotRect.left() + int(i * barW);
            const int y = mapY(values[i]);
            const int h = plotRect.bottom() - y;

            QRect bar(x + 1,
                      y,
                      std::max(1, int(barW) - 2),
                      std::max(0, h));

            if (values[i] > 0.0)
                p.fillRect(bar, QColor(160, 190, 255));
            else
                p.fillRect(bar, QColor(248, 248, 248));

            p.setPen(QPen(Qt::black, 1));
            p.drawRect(bar);
        }

        auto drawOneDLegend = [&](const QPen &samplePen, bool includeGaussian)
        {
            p.save();

            const int legendW = 210;
            const int legendH = includeGaussian ? 140 : 36;
            const int legendX = plotRect.right() - legendW - 12;
            const int legendY = plotRect.top() + 12;

            QRect legendRect(legendX, legendY, legendW, legendH);

            p.fillRect(legendRect, QColor(255, 255, 255, 225));
            p.setPen(QPen(QColor(120, 120, 120), 1));
            p.drawRect(legendRect);

            QFont legendFont = p.font();
            legendFont.setPointSize(8);
            legendFont.setBold(false);
            p.setFont(legendFont);
            p.setPen(Qt::black);

            p.drawText(legendX + 10,
                       legendY + 8,
                       legendW - 20,
                       14,
                       Qt::AlignLeft,
                       "N = " + QString::number(n));

            if (!includeGaussian)
            {
                p.restore();
                return;
            }

            p.setPen(samplePen);
            p.drawLine(legendX + 10, legendY + 33, legendX + 42, legendY + 33);
            p.setPen(Qt::black);

            p.drawText(legendX + 50,
                       legendY + 24,
                       legendW - 58,
                       16,
                       Qt::AlignLeft,
                       "Gaussian fit");

            p.drawText(legendX + 10,
                       legendY + 48,
                       legendW - 20,
                       14,
                       Qt::AlignLeft,
                       "Amp = " + formatGaussianNumber(gaussianStats.amplitude));

            p.drawText(legendX + 10,
                       legendY + 63,
                       legendW - 20,
                       14,
                       Qt::AlignLeft,
                       "Mean = " + formatGaussianNumber(gaussianStats.mean));

            p.drawText(legendX + 10,
                       legendY + 78,
                       legendW - 20,
                       14,
                       Qt::AlignLeft,
                       "Rms = " + formatGaussianNumber(gaussianStats.sigma));

            p.drawText(legendX + 10,
                       legendY + 93,
                       legendW - 20,
                       14,
                       Qt::AlignLeft,
                       "Chi^2 = " + formatGaussianNumber(gaussianStats.chiSquare));

            p.drawText(legendX + 10,
                       legendY + 108,
                       legendW - 20,
                       14,
                       Qt::AlignLeft,
                       "Area = " + formatGaussianNumber(gaussianStats.area));

            p.restore();
        };

        QPen gaussianPen(QColor(220, 40, 40, 155), 2, Qt::DotLine);
        gaussianPen.setCapStyle(Qt::RoundCap);

        if (showGaussian && gaussianStats.valid)
        {
            auto mapX = [&](double xValue)
            {
                const double frac = (xValue - minX) / (maxX - minX);
                return plotRect.left() + int(frac * plotRect.width());
            };

            QPolygonF gaussianCurve;

            const int samples = 500;

            for (int i = 0; i <= samples; ++i)
            {
                const double frac = double(i) / double(samples);
                const double xValue = minX + frac * (maxX - minX);

                const double yValue = gaussianValue(gaussianStats.amplitude,
                                                    gaussianStats.mean,
                                                    gaussianStats.sigma,
                                                    xValue);

                gaussianCurve << QPointF(mapX(xValue), mapY(yValue));
            }

            p.setPen(gaussianPen);
            p.drawPolyline(gaussianCurve);
        }

        drawOneDLegend(gaussianPen, showGaussian && gaussianStats.valid);

        p.setPen(Qt::black);

        p.drawText(plotRect.left(),
                   height() - 35,
                   plotRect.width(),
                   24,
                   Qt::AlignCenter,
                   xTitle);

        p.save();
        p.translate(28, plotRect.top() + plotRect.height() / 2);
        p.rotate(-90);
        p.drawText(QRect(-plotRect.height() / 2,
                         -20,
                         plotRect.height(),
                         20),
                   Qt::AlignCenter,
                   yTitle);
        p.restore();
    }

private:
    AngularDistEntry m_entry;
    PaceAngularAccum m_acc;
    int m_plotKind = PlotCountsVsTheta;
    double m_sigmaTotal = 1.0;
    int m_nEvents = 1;
};
}

void addAngularSample(AngularDistEntry &entry,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy)
{
    addAngularSample(entry, kineticEnergy, thetaDeg, vz, vxy, kineticEnergy);
}

void addAngularSample(AngularDistEntry &entry,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy,
                      float cmEnergy)
{
    entry.kineticEnergy.push_back(kineticEnergy);
    entry.thetaDeg.push_back(thetaDeg);
    entry.vz.push_back(vz);
    entry.vxy.push_back(vxy);
    entry.cmEnergy.push_back(cmEnergy);
}

void addAngularSample(std::map<std::pair<int, int>, AngularDistEntry> &entries,
                      int z,
                      int n,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy)
{
    const std::pair<int, int> key(z, n);
    AngularDistEntry &e = entries[key];
    e.z = z;
    e.n = n;
    addAngularSample(e, kineticEnergy, thetaDeg, vz, vxy);
}

void addAngularSample(std::map<std::pair<int, int>, AngularDistEntry> &entries,
                      int z,
                      int n,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy,
                      float cmEnergy)
{
    const std::pair<int, int> key(z, n);
    AngularDistEntry &e = entries[key];
    e.z = z;
    e.n = n;
    addAngularSample(e, kineticEnergy, thetaDeg, vz, vxy, cmEnergy);
}

QString buildEmittedParticleCMSpectraHtmlGemini(
    const AngularDistEntry &neutronEntry,
    const AngularDistEntry &protonEntry,
    const AngularDistEntry &alphaEntry,
    const AngularDistEntry &gammaEntry)
{
    constexpr int kMaxEnergyBins = 50;
    constexpr int kTotalIndex = 51;

    std::array<std::array<int, kTotalIndex + 1>, 5> NSPC{};
    std::array<double, 5> averageEnergy{};

    auto fillParticle = [&](int mode, const AngularDistEntry &entry)
    {
        const std::vector<float> &energies =
            entry.cmEnergy.empty() ? entry.kineticEnergy : entry.cmEnergy;

        for (float energy : energies)
        {
            if (energy < 0.0f) continue;

            int bin = int(energy) + 1;
            if (bin < 1) bin = 1;
            if (bin > kMaxEnergyBins) bin = kMaxEnergyBins;

            NSPC[mode][bin]++;
        }
    };

    fillParticle(1, neutronEntry);
    fillParticle(2, protonEntry);
    fillParticle(3, alphaEntry);
    fillParticle(4, gammaEntry);

    bool hasAnyCounts = false;
    for (int mode = 1; mode <= 4; ++mode)
    {
        for (int bin = 1; bin <= kMaxEnergyBins; ++bin)
        {
            if (NSPC[mode][bin] > 0)
            {
                hasAnyCounts = true;
                break;
            }
        }
    }

    if (!hasAnyCounts) return QString();

    QString html;
    html += "<p>&nbsp;</p>";
    html += "<h2 align=\"center\">C.M. spectra of emitted particles</h2>";
    html += "<table class=\"cm-spectra\" align=\"center\">";
    html += "<tr><th>Ex(MeV)</th><th>Neut</th><th>Prot</th><th>Alpha</th><th>Gamma</th></tr>";

    for (int bin = 1; bin <= kMaxEnergyBins; ++bin)
    {
        const bool printRow =
            NSPC[1][bin] != 0 ||
            NSPC[2][bin] != 0 ||
            NSPC[3][bin] != 0 ||
            NSPC[4][bin] != 0;

        const double particleEnergy = double(bin) - 0.5;
        for (int mode = 1; mode <= 4; ++mode)
        {
            NSPC[mode][kTotalIndex] += NSPC[mode][bin];
            averageEnergy[mode] += particleEnergy * double(NSPC[mode][bin]);
        }

        if (!printRow) continue;

        const int eLow = bin - 1;
        const int eHigh = (bin == kMaxEnergyBins) ? 99 : bin;

        html += "<tr><td align=\"center\">" +
                QString::number(eLow) + " - " + QString::number(eHigh) +
                "</td>";

        for (int mode = 1; mode <= 4; ++mode)
        {
            if (NSPC[mode][bin] > 0)
                html += "<td align=\"center\">" + QString::number(NSPC[mode][bin]) + "</td>";
            else
                html += "<td></td>";
        }

        html += "</tr>";
    }

    for (int mode = 1; mode <= 4; ++mode)
        averageEnergy[mode] /= double(NSPC[mode][kTotalIndex]) + 1.0e-9;

    html += "<tr><th>Total</th>";
    for (int mode = 1; mode <= 4; ++mode)
        html += "<th>" + QString::number(NSPC[mode][kTotalIndex]) + "</th>";
    html += "</tr>";

    html += "<tr><th>Average Energy</th>";
    for (int mode = 1; mode <= 4; ++mode)
        html += "<td align=\"center\">" + QString::number(averageEnergy[mode], 'f', 2) + "</td>";
    html += "</tr>";

    html += "</table>";
    return html;
}

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
    int mdir,
    int inputMode,
    const AngularDistEntry &neutronEntry,
    const AngularDistEntry &protonEntry,
    const AngularDistEntry &alphaEntry,
    const AngularDistEntry &gammaEntry)
{
    QString html;
    html += "<!DOCTYPE html><html><head><meta charset=\"utf-8\">";
    html += "<style>"
            "body{font-family:Sans-Serif;}"
            "table{border-collapse:collapse; margin:10px 0;}"
            "th,td{border:1px solid #909090; padding:4px 6px;}"
            "th{background:#ececec;}"
            "h2{color:#1d4f91;}"
            "h3{margin-top:25px;}"
            ".small{font-size:12px; color:#666;}"
            ".cm-spectra{margin:24px auto; min-width:520px;}"
            ".cm-spectra th{background:#dfeaf7; color:#1d4f91;}"
            ".cm-spectra td,.cm-spectra th{padding:5px 10px; text-align:center;}"
            "</style></head><body>";

    html += "<p class=\"small\">" + title + "</p>";

    if ((!hasAngularSamples(entries)
         && !hasAngularSamples(neutronEntry)
         && !hasAngularSamples(protonEntry)
         && !hasAngularSamples(alphaEntry)
         && !hasAngularSamples(gammaEntry))
        || nCascades <= 0)
    {
        html += "<p>No angular-distribution events were found for output.</p></body></html>";
        return html;
    }

    const std::vector<PaceSelectedResidue> selected =
        buildSelectedResiduesPACE(entries, nCascades, lowLimitPercent, highLimitPercent);

    html += "<p>&nbsp;</p><h2 align=\"center\">Angular distribution results</h2>";
    html += "<p>";
    html += "<a href=\"gemini://plot_all\">Plot All: E vs &theta;</a> &nbsp; ";
    html += "<a href=\"gemini://plot_all_ntheta\">Plot All: N vs &theta;</a> &nbsp; ";
    html += "<a href=\"gemini://plot_all_ne\">Plot All: N vs Energy</a> &nbsp; ";
    html += "<a href=\"gemini://plot_all_dcsdtheta\">Plot All: d&sigma;/d&theta; vs &theta;</a>";
    html += "</p>";

    if (mdir == 0)
        html += "<p>*** Spin alignment perpendicular to recoil axis - standard compound nucleus angular distribution</p>";
    else if (mdir == 1)
        html += "<p>*** Spin alignment perpendicular to reaction plane - angular distribution is around Z axis perpendicular to<br>*** reaction plane</p>";

    PaceAngularAccum allAcc;
    initPaceAccum(allAcc, compoundA, compoundExcitationMeV, recoilBetaCN, mdir);

    std::vector<std::pair<PaceSelectedResidue, PaceAngularAccum>> selectedAcc;
    selectedAcc.reserve(selected.size());
    for (const auto &r : selected)
    {
        PaceAngularAccum acc;
        initPaceAccum(acc, compoundA, compoundExcitationMeV, recoilBetaCN, mdir);
        selectedAcc.push_back({r, acc});
    }

    for (const auto &it : entries)
    {
        const AngularDistEntry &e = it.second;
        const int m = std::min(e.kineticEnergy.size(), e.thetaDeg.size());

        int selectedIndex = -1;
        for (int i = 0; i < (int)selected.size(); ++i)
        {
            if (selected[i].z == e.z && selected[i].n == e.n)
            {
                selectedIndex = i;
                break;
            }
        }

        for (int i = 0; i < m; ++i)
        {
            addToPaceAccum(allAcc, e.kineticEnergy[i], e.thetaDeg[i]);
            if (selectedIndex >= 0)
                addToPaceAccum(selectedAcc[selectedIndex].second, e.kineticEnergy[i], e.thetaDeg[i]);
        }
    }

    int displayIndex = 1;
    int plotIndex = 0;

    for (const auto &item : selectedAcc)
    {
        const PaceSelectedResidue &r = item.first;
        const PaceAngularAccum &acc = item.second;

        html += buildHeaderForResiduePACE(displayIndex, r.z, r.n);
        appendMainPaceTable(html, acc, r.a, sigmaTotalMb, nCascades, inputMode, false, plotIndex,
                            "Residual velocity/c");
        appendSecondPaceTableIfNeeded(html, acc, plotIndex);
        displayIndex++;
        plotIndex++;
    }

    html += buildHeaderAllPACE(displayIndex);
    appendMainPaceTable(html, allAcc, compoundA, sigmaTotalMb, nCascades, inputMode, true, plotIndex,
                        "Residual velocity/c");
    appendSecondPaceTableIfNeeded(html, allAcc, plotIndex);
    displayIndex++;
    plotIndex++;

    auto appendParticleSection = [&](const QString &particleLabel,
                                     const QString &velocityLabel,
                                     const AngularDistEntry &entry,
                                     int particleA,
                                     bool showVelocity)
    {
        if (!hasAngularSamples(entry)) return;

        PaceAngularAccum acc;
        initPaceAccum(acc, compoundA, compoundExcitationMeV, recoilBetaCN, mdir);
        setParticleEnergyInterval(acc);
        appendEntryToAccum(acc, entry);

        html += buildHeaderForParticlePACE(displayIndex, particleLabel);
        appendMainPaceTable(html, acc, particleA, sigmaTotalMb, nCascades,
                            showVelocity ? inputMode : 0,
                            false, plotIndex, velocityLabel);

        appendSecondPaceTableIfNeeded(html, acc, plotIndex);
        displayIndex++;
        plotIndex++;
    };

    appendParticleSection("neutrons", "Neutron velocity/c", neutronEntry, 1, true);
    appendParticleSection("protons", "Proton velocity/c", protonEntry, 1, true);
    appendParticleSection("alpha particles", "Alpha velocity/c", alphaEntry, 4, true);
    appendParticleSection("gamma particles", "", gammaEntry, 1, false);

    html += buildEmittedParticleCMSpectraHtmlGemini(neutronEntry,
                                                    protonEntry,
                                                    alphaEntry,
                                                    gammaEntry);

    html += "</body></html>";
    return html;
}

// ----------------------- Widget -----------------------

AngularDistributionWidget::AngularDistributionWidget(const QString &htmlContent,
                                                     const std::map<std::pair<int, int>, AngularDistEntry> &entries,
                                                     double sigmaTotal,
                                                     int nEvents,
                                                     double lowLimitPercent,
                                                     double highLimitPercent,
                                                     const QString &title,
                                                     const AngularDistEntry &neutronEntry,
                                                     const AngularDistEntry &protonEntry,
                                                     const AngularDistEntry &alphaEntry,
                                                     const AngularDistEntry &gammaEntry,
                                                     double compoundExcitationMeV,
                                                     int compoundA,
                                                     double recoilBetaCN,
                                                     int mdir,
                                                     QWidget *parent)
    : QDialog(parent),
    html(htmlContent),
    entriesForPlots(entries),
    sigmaTotalForPlots(sigmaTotal),
    nEventsForPlots(nEvents),
    lowLimitForPlots(lowLimitPercent),
    highLimitForPlots(highLimitPercent),
    plotTitle(title),
    neutronEntryForPlots(neutronEntry),
    protonEntryForPlots(protonEntry),
    alphaEntryForPlots(alphaEntry),
    gammaEntryForPlots(gammaEntry),
    compoundExcitationForPlots(compoundExcitationMeV),
    compoundAForPlots(compoundA),
    recoilBetaForPlots(recoilBetaCN),
    mdirForPlots(mdir)
{
    setWindowTitle("Gemini: Angular distribution");
    setWindowIcon(QIcon(":/Gemini_logo.png"));
    setModal(false);

    QVBoxLayout *layout = new QVBoxLayout;

    QToolBar *toolbar = new QToolBar;
    toolbar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);

    QAction *saveAction = new QAction(QIcon(":/save29.png"), "&Save", this);
    toolbar->addAction(saveAction);
    connect(saveAction, SIGNAL(triggered()), this, SLOT(save_clicked()));

    QAction *printAction = new QAction(QIcon(":/printer70.png"), "&Print", this);
    toolbar->addAction(printAction);
    connect(printAction, SIGNAL(triggered()), this, SLOT(print_clicked()));

    layout->addWidget(toolbar);

    QTextBrowser *text = new QTextBrowser;
    text->setOpenLinks(false);
    text->setHtml(htmlContent);
    connect(text, SIGNAL(anchorClicked(QUrl)), this, SLOT(link_clicked(QUrl)));
    layout->addWidget(text);

    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    setLayout(layout);
    setMinimumSize(1200, 760);
}

void AngularDistributionWidget::link_clicked(const QUrl &url)
{
    if (url.scheme() != "gemini")
    {
        QDesktopServices::openUrl(url);
        return;
    }

    if (url.host() == "plot_all")
    {
        openPlotWindow(true, -1, PlotEnergyVsTheta);
        return;
    }

    if (url.host() == "plot_all_ntheta")
    {
        openPlotWindow(true, -1, PlotCountsVsTheta);
        return;
    }

    if (url.host() == "plot_all_ne")
    {
        openPlotWindow(true, -1, PlotCountsVsEnergy);
        return;
    }

    if (url.host() == "plot_all_dcsdtheta")
    {
        openPlotWindow(true, -1, PlotCrossSectionVsTheta);
        return;
    }

    if (url.host() == "plot_table" ||
        url.host() == "plot_ntheta" ||
        url.host() == "plot_ne" ||
        url.host() == "plot_dcsdtheta")
    {
        bool ok = false;
        int idx = url.path().mid(1).toInt(&ok);
        if (!ok) return;

        int kind = PlotEnergyVsTheta;
        if (url.host() == "plot_ntheta") kind = PlotCountsVsTheta;
        else if (url.host() == "plot_ne") kind = PlotCountsVsEnergy;
        else if (url.host() == "plot_dcsdtheta") kind = PlotCrossSectionVsTheta;

        openPlotWindow(false, idx, kind);
    }
}

void AngularDistributionWidget::openPlotWindow(bool plotAllTables, int tableIndex, int plotKind)
{
    if ((!hasAngularSamples(entriesForPlots)
         && !hasAngularSamples(neutronEntryForPlots)
         && !hasAngularSamples(protonEntryForPlots)
         && !hasAngularSamples(alphaEntryForPlots)
         && !hasAngularSamples(gammaEntryForPlots))
        || nEventsForPlots <= 0) return;

    const std::vector<PlotTableEntry> tables =
        buildPlotTableListPACE(entriesForPlots,
                               neutronEntryForPlots,
                               protonEntryForPlots,
                               alphaEntryForPlots,
                               gammaEntryForPlots,
                               nEventsForPlots,
                               lowLimitForPlots,
                               highLimitForPlots,
                               compoundExcitationForPlots,
                               compoundAForPlots,
                               recoilBetaForPlots,
                               mdirForPlots);

    QDialog *dlg = new QDialog(this);
    dlg->setAttribute(Qt::WA_DeleteOnClose, true);
    dlg->setWindowIcon(QIcon(":/Gemini_logo.png"));
    QString plotKindLabel = "E vs θ";
    if (plotKind == PlotCountsVsTheta) plotKindLabel = "N vs θ";
    else if (plotKind == PlotCountsVsEnergy) plotKindLabel = "N vs Energy";
    else if (plotKind == PlotCrossSectionVsTheta) plotKindLabel = "dσ/dθ vs θ";

    dlg->setWindowTitle(plotAllTables
                            ? "Gemini: Angular distribution plots (all) - " + plotKindLabel
                            : "Gemini: Angular distribution plot - " + plotKindLabel);
    dlg->setModal(false);

    QVBoxLayout *mainLayout = new QVBoxLayout(dlg);

    QHBoxLayout *topRow = new QHBoxLayout;

    QLabel *topLabel = new QLabel("<h2 style='color:#1d4f91; margin:0;'>Angular distribution plots: "
                                  + plotKindLabel.toHtmlEscaped() + "</h2>"
                                                                    "<p style='margin:2px 0 0 0;'><em>" + plotTitle.toHtmlEscaped() + "</em></p>");
    topLabel->setTextFormat(Qt::RichText);
    topRow->addWidget(topLabel);
    topRow->addStretch();

    QPushButton *savePngButton = new QPushButton("Save as PNG");
    savePngButton->setFixedHeight(26);
    savePngButton->setMaximumWidth(110);
    topRow->addWidget(savePngButton, 0, Qt::AlignTop | Qt::AlignRight);

    mainLayout->addLayout(topRow);

    QScrollArea *scroll = new QScrollArea;
    scroll->setWidgetResizable(true);

    QWidget *container = new QWidget;
    QVBoxLayout *plotsLayout = new QVBoxLayout(container);

    int shownCount = 0;
    for (int i = 0; i < (int)tables.size(); ++i)
    {
        if (!plotAllTables && i != tableIndex) continue;

        QLabel *lbl = new QLabel("<b>" + tables[i].label.toHtmlEscaped() + "</b>");
        lbl->setWordWrap(true);
        lbl->setTextFormat(Qt::RichText);
        plotsLayout->addWidget(lbl);

        QWidget *plot = nullptr;
        if (plotKind == PlotEnergyVsTheta)
            plot = new ScatterPlotWidget(tables[i].entry, tables[i].acc);
        else
            plot = new OneDPlotWidget(tables[i].entry,
                                      tables[i].acc,
                                      plotKind,
                                      sigmaTotalForPlots,
                                      nEventsForPlots);

        plotsLayout->addWidget(plot);

        shownCount++;
    }

    plotsLayout->addStretch();
    scroll->setWidget(container);
    mainLayout->addWidget(scroll);

    QObject::connect(savePngButton, &QPushButton::clicked, dlg, [container, dlg]()
                     {
                         const QString fileName = QFileDialog::getSaveFileName(
                             dlg,
                             QObject::tr("Save Plot as PNG"),
                             QString(),
                             QObject::tr("PNG Image (*.png)")
                             );

                         if (fileName.isEmpty()) return;

                         QPixmap pixmap(container->size());
                         pixmap.fill(Qt::white);
                         container->render(&pixmap);
                         pixmap.save(fileName, "PNG");
                     });

    dlg->setLayout(mainLayout);

    const int dialogWidth = 1080;
    int dialogHeight = 220 + shownCount * 620;
    if (dialogHeight > 920) dialogHeight = 920;
    if (dialogHeight < 700) dialogHeight = 700;

    dlg->resize(dialogWidth, dialogHeight);
    dlg->show();
}

void AngularDistributionWidget::save_clicked()
{
    const QString fileName = QFileDialog::getSaveFileName(this,
                                                          tr("Save File"),
                                                          QString(),
                                                          tr("Gemini angular distribution (*.html)"));
    if (fileName.isEmpty()) return;

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly)) return;

    QTextStream stream(&file);
    stream << html;
    file.close();
}

void AngularDistributionWidget::print_clicked()
{
    QPrinter printer;
    QPrintDialog printDialog(&printer, this);
    printDialog.setWindowTitle(tr("Print Angular Distribution"));

    if (printDialog.exec() != QDialog::Accepted) return;

    printer.setFullPage(true);
    QTextDocument textDoc;
    textDoc.setHtml(html);
    textDoc.print(&printer);
}
