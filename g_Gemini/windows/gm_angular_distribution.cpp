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

QString nucleusLabelFromZNHtml(int z, int n)
{
    CNucleus nuc(z, z + n);
    QString s = nuc.getGName();
    s.replace("<font color=\"darkBlue\">", "");
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

int totalAngularCount(const PaceAngularAccum &acc)
{
    return isum(acc, 32, 1, 37);
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

    html += "<tr class=\"total-row\"><th>Total</th>";
    for (int k = 1; k <= 18; ++k)
        html += (acc.NR[32][k] > 0) ? "<th align=\"center\">" + QString::number(acc.NR[32][k]) + "</th>" : "<th></th>";
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

    html += "<tr class=\"total-row\"><th>Total</th>";
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

QString buildHeaderForResiduePACE(int displayIndex,
                                  int z,
                                  int n,
                                  const QString &fragmentKind = "residual nucleus")
{
    const QString nuc = nucleusLabelFromZNHtml(z, n);
    return QString("<br><h3>%1. Energy and angular distribution of %2 Z = "
                   "<span style=\"color:blue\">%3</span> and N = <span style=\"color:blue\">%4</span> "
                   "(<span style=\"color:blue\">%5</span>)</h3>")
        .arg(displayIndex)
        .arg(fragmentKind.toHtmlEscaped())
        .arg(z)
        .arg(n)
        .arg(nuc);
}

QString buildHeaderAllPACE(int displayIndex, const QString &fragmentKindPlural = "residual nuclei")
{
    return QString("<br><h3>%1. Energy and angular distribution of ALL %2</h3>")
        .arg(displayIndex)
        .arg(fragmentKindPlural.toHtmlEscaped());
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

enum OverlayParticleKind
{
    OverlayNone = 0,
    OverlayNeutron,
    OverlayProton,
    OverlayAlpha
};

struct PlotTableEntry
{
    QString label;
    QString htmlLabel;
    AngularDistEntry entry;
    PaceAngularAccum acc;
    bool allowGaussianEnergyOverlay = true;
    bool allowBoltzmannEnergyFit = false;
    OverlayParticleKind boltzmannKind = OverlayNone;
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
    int mdir,
    bool isImf = false)
{
    const QString fragmentKind = isImf ? "IMF fragment" : "residual nucleus";
    const QString fragmentKindPlural = isImf ? "IMF fragments" : "residual nuclei";

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
        if (totalAngularCount(item.second) <= 0) continue;
        t.label = QString("%1. Energy and angular distribution of %2 Z = %3 and N = %4 (%5)")
                      .arg(idx)
                      .arg(fragmentKind)
                      .arg(r.z)
                      .arg(r.n)
                      .arg(nucleusLabelFromZNPlain(r.z, r.n));
        t.htmlLabel = QString("%1. Energy and angular distribution of %2 Z = %3 and N = %4 (%5)")
                          .arg(idx)
                          .arg(fragmentKind.toHtmlEscaped())
                          .arg(r.z)
                          .arg(r.n)
                          .arg(nucleusLabelFromZNHtml(r.z, r.n));
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
    if (totalAngularCount(allAcc) > 0)
    {
        tAll.label = QString("%1. Energy and angular distribution of ALL %2")
                         .arg(idx)
                         .arg(fragmentKindPlural);
        tAll.entry = all;
        tAll.acc = allAcc;
        tables.push_back(tAll);
        idx++;
    }

    auto appendParticlePlot = [&](const QString &particleLabel,
                                  const AngularDistEntry &entry,
                                  OverlayParticleKind overlayParticleKind)
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
        if (totalAngularCount(table.acc) <= 0) return;
        if (overlayParticleKind != OverlayNone)
        {
            table.allowBoltzmannEnergyFit = true;
            table.boltzmannKind = overlayParticleKind;
        }
        else
        {
            table.allowGaussianEnergyOverlay = false;
        }
        tables.push_back(table);
        idx++;
    };

    appendParticlePlot("neutrons", neutronEntry, OverlayNeutron);
    appendParticlePlot("protons", protonEntry, OverlayProton);
    appendParticlePlot("alpha particles", alphaEntry, OverlayAlpha);
    appendParticlePlot("gamma particles", gammaEntry, OverlayNone);

    return tables;
}

enum PlotKind
{
    PlotEnergyVsTheta = 0,
    PlotCountsVsTheta = 1,
    PlotCountsVsEnergy = 2,
    PlotCrossSectionVsTheta = 3
};

constexpr int kCMSpectraMaxEnergyBins = 50;
constexpr int kCMSpectraTotalIndex = 51;

struct CMSpectraSeries
{
    QString label;
    QColor color;
    std::array<int, kCMSpectraTotalIndex + 1> counts{};
    double averageEnergy = 0.0;
};

void fillCMSpectraSeries(CMSpectraSeries &series, const AngularDistEntry &entry)
{
    const std::vector<float> &energies =
        entry.cmEnergy.empty() ? entry.kineticEnergy : entry.cmEnergy;

    double weightedEnergy = 0.0;

    for (float energy : energies)
    {
        if (energy < 0.0f) continue;

        int bin = int(energy) + 1;
        if (bin < 1) bin = 1;
        if (bin > kCMSpectraMaxEnergyBins) bin = kCMSpectraMaxEnergyBins;

        const double particleEnergy = double(bin) - 0.5;
        series.counts[bin]++;
        series.counts[kCMSpectraTotalIndex]++;
        weightedEnergy += particleEnergy;
    }

    if (series.counts[kCMSpectraTotalIndex] > 0)
        series.averageEnergy = weightedEnergy / double(series.counts[kCMSpectraTotalIndex]);
}

std::vector<CMSpectraSeries> buildCMSpectraSeries(
    const AngularDistEntry &neutronEntry,
    const AngularDistEntry &protonEntry,
    const AngularDistEntry &alphaEntry,
    const AngularDistEntry &gammaEntry)
{
    std::vector<CMSpectraSeries> series;
    series.reserve(4);

    auto append = [&](const QString &label, const QColor &color, const AngularDistEntry &entry)
    {
        CMSpectraSeries s;
        s.label = label;
        s.color = color;
        fillCMSpectraSeries(s, entry);

        series.push_back(s);
    };

    append("Neutrons", QColor(175, 205, 238), neutronEntry);
    append("Protons", QColor(246, 190, 166), protonEntry);
    append("Alpha", QColor(178, 224, 188), alphaEntry);
    append("Gamma", QColor(216, 190, 235), gammaEntry);

    return series;
}

bool hasCMSpectraSamples(const AngularDistEntry &entry)
{
    return !entry.cmEnergy.empty() || !entry.kineticEnergy.empty();
}

bool hasCMSpectraSamples(const AngularDistEntry &neutronEntry,
                         const AngularDistEntry &protonEntry,
                         const AngularDistEntry &alphaEntry,
                         const AngularDistEntry &gammaEntry)
{
    return hasCMSpectraSamples(neutronEntry) ||
           hasCMSpectraSamples(protonEntry) ||
           hasCMSpectraSamples(alphaEntry) ||
           hasCMSpectraSamples(gammaEntry);
}

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
                                                 const std::vector<double> &yValues,
                                                 bool trimToPopulatedSpan = false)
{
    GaussianOverlayStats stats;

    if (xCenters.size() != yValues.size() || xCenters.empty())
        return stats;

    int firstPopulated = -1;
    int lastPopulated = -1;
    int populatedBins = 0;

    for (int i = 0; i < int(yValues.size()); ++i)
    {
        if (yValues[i] <= 0.0) continue;

        if (firstPopulated < 0)
            firstPopulated = i;
        lastPopulated = i;
        populatedBins++;
    }

    if (firstPopulated < 0 || lastPopulated < 0)
        return stats;

    const int fitFirst = trimToPopulatedSpan ? firstPopulated : 0;
    const int fitLast = trimToPopulatedSpan ? lastPopulated : int(yValues.size()) - 1;

    std::vector<double> fitX;
    std::vector<double> fitY;
    fitX.reserve(fitLast - fitFirst + 1);
    fitY.reserve(fitLast - fitFirst + 1);

    for (int i = fitFirst; i <= fitLast; ++i)
    {
        fitX.push_back(xCenters[i]);
        fitY.push_back(std::max(0.0, yValues[i]));
    }

    double sumW = 0.0;
    double sumX = 0.0;
    double maxY = 0.0;

    for (int i = 0; i < int(fitY.size()); ++i)
    {
        const double w = fitY[i];
        if (w <= 0.0) continue;

        sumW += w;
        sumX += w * fitX[i];
        maxY = std::max(maxY, w);
    }

    if (sumW <= 0.0 || maxY <= 0.0)
        return stats;

    const double xMin = fitX.front();
    const double xMax = fitX.back();
    const double xRange = std::max(1.0e-9, xMax - xMin);
    const double fallbackSigma =
        (fitX.size() >= 2)
            ? std::max(1.0e-9, std::fabs(fitX.back() - fitX.front()) / double(fitX.size() - 1))
            : std::max(1.0, xRange);

    double A = maxY;
    double mu = sumX / sumW;

    double variance = 0.0;
    for (int i = 0; i < int(fitY.size()); ++i)
    {
        const double w = fitY[i];
        if (w <= 0.0) continue;

        const double dx = fitX[i] - mu;
        variance += w * dx * dx;
    }

    double sigma = std::sqrt(std::max(variance / sumW, 1.0e-12));
    sigma = std::max(sigma, fallbackSigma);

    if (populatedBins < 3 || fitX.size() < 3)
    {
        stats.amplitude = A;
        stats.mean = mu;
        stats.sigma = sigma;
        stats.valid = std::isfinite(stats.amplitude) &&
                      std::isfinite(stats.mean) &&
                      std::isfinite(stats.sigma) &&
                      stats.amplitude > 0.0 &&
                      stats.sigma > 1.0e-12;
        computeGaussianQualityAndArea(stats, fitX, fitY);
        return stats;
    }

    auto sseFor = [&](double testA, double testMu, double testSigma)
    {
        if (testA <= 0.0 || testSigma <= 1.0e-12)
            return 1.0e300;

        double sse = 0.0;

        for (int i = 0; i < int(fitY.size()); ++i)
        {
            const double yFit = gaussianValue(testA, testMu, testSigma, fitX[i]);
            const double r = fitY[i] - yFit;
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

        for (int i = 0; i < int(fitY.size()); ++i)
        {
            const double x = fitX[i];
            const double y = fitY[i];

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
        newSigma = std::max(fallbackSigma, std::min(std::max(fallbackSigma, xRange * 10.0), newSigma));

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

    computeGaussianQualityAndArea(stats, fitX, fitY);

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

struct BoltzmannOverlayCurve
{
    bool valid = false;
    bool drawBarrier = false;
    bool unavailable = false;
    QString unavailableReason;
    double amplitude = 0.0;
    double temperatureMeV = 1.0;
    double barrierMeV = 0.0;
    bool barrierEstimated = false;
    double chiSquare = 0.0;
    double reducedChiSquare = 0.0;
    double maxY = 0.0;
    QPolygonF points;
};

struct BoltzmannFitResult
{
    bool valid = false;
    bool unavailable = false;
    QString unavailableReason;
    double amplitude = 0.0;
    double temperatureMeV = 0.0;
    double barrierMeV = 0.0;
    double chiSquare = 0.0;
    double reducedChiSquare = 0.0;
};

struct CoulombBarrierInfo
{
    double valueMeV = 0.0;
    bool estimated = false;
};

CoulombBarrierInfo coulombBarrierForParticle(OverlayParticleKind particleKind,
                                             int sourceZ,
                                             int sourceA)
{
    CoulombBarrierInfo barrier;

    if (particleKind == OverlayNeutron)
        return barrier;

    int zParticle = 0;
    int aParticle = 0;
    double fallback = 0.0;
    double minBarrier = 0.0;
    double maxBarrier = 0.0;

    if (particleKind == OverlayProton)
    {
        zParticle = 1;
        aParticle = 1;
        fallback = 5.0;
        minBarrier = 1.0;
        maxBarrier = 30.0;
    }
    else if (particleKind == OverlayAlpha)
    {
        zParticle = 2;
        aParticle = 4;
        fallback = 10.0;
        minBarrier = 2.0;
        maxBarrier = 60.0;
    }
    else
    {
        return barrier;
    }

    barrier.valueMeV = fallback;

    const int zDaughter = sourceZ - zParticle;
    const int aDaughter = sourceA - aParticle;
    if (sourceZ <= 0 || sourceA <= 0 || zDaughter <= 0 || aDaughter <= 0)
        return barrier;

    constexpr double r0 = 1.2;
    const double radiusSum = r0 * (std::cbrt(double(aParticle)) + std::cbrt(double(aDaughter)));
    if (radiusSum <= 0.0)
        return barrier;

    const double estimated =
        1.44 * double(zParticle) * double(zDaughter) / radiusSum;
    if (!std::isfinite(estimated) || estimated <= 0.0)
        return barrier;

    barrier.valueMeV = std::clamp(estimated, minBarrier, maxBarrier);
    barrier.estimated = true;
    return barrier;
}

double boltzmannShape(double energyMeV,
                      double temperatureMeV,
                      double barrierMeV,
                      bool charged)
{
    if (temperatureMeV <= 0.0)
        return 0.0;

    if (!charged)
    {
        if (energyMeV <= 0.0)
            return 0.0;
        return std::sqrt(energyMeV) * std::exp(-energyMeV / temperatureMeV);
    }

    if (energyMeV <= barrierMeV)
        return 0.0;

    const double shiftedEnergy = energyMeV - barrierMeV;
    if (shiftedEnergy <= 0.0 || energyMeV <= 1.0e-12)
        return 0.0;

    return (std::pow(shiftedEnergy, 1.5) / energyMeV) *
           std::exp(-shiftedEnergy / temperatureMeV);
}

BoltzmannFitResult fitBoltzmannToHistogram(const std::vector<double> &xCenters,
                                           const std::vector<double> &yValues,
                                           OverlayParticleKind particleKind,
                                           double barrierMeV)
{
    // Boltzmann fitting is used only as a diagnostic for emitted neutron/proton/alpha spectra.
    // It is not applied to residual nuclei or IMF fragments. In GEMINI, IMFs are
    // complex-fragment/asymmetric binary-decay products, so a simple emitted-particle
    // Boltzmann+Coulomb spectrum is not appropriate for IMF plots. The Coulomb barrier
    // used here is an approximate plotting barrier, not the full GEMINI
    // transmission-coefficient treatment.
    //
    // This is a hard Coulomb-barrier approximation using a simple touching-spheres
    // Coulomb estimate. Sub-barrier tunneling, barrier distributions, and full GEMINI
    // transmission coefficients are not included in this plotting fit.
    BoltzmannFitResult result;

    if (xCenters.size() != yValues.size() || xCenters.empty() || particleKind == OverlayNone)
    {
        result.unavailable = true;
        result.unavailableReason = "Boltzmann fit unavailable";
        return result;
    }

    const bool charged = particleKind != OverlayNeutron;
    const double barrier = charged ? barrierMeV : 0.0;
    result.barrierMeV = barrier;

    std::vector<double> fitX;
    std::vector<double> fitY;

    for (int i = 0; i < int(yValues.size()); ++i)
    {
        const double y = yValues[i];
        const double x = xCenters[i];
        if (y <= 0.0) continue;
        if (charged && x <= barrier) continue;

        fitX.push_back(x);
        fitY.push_back(y);
    }

    if (fitX.size() < 3)
    {
        result.unavailable = true;
        result.unavailableReason = "Boltzmann fit unavailable: too few points";
        return result;
    }

    auto evaluateT = [&](double temperature,
                         double &amplitude,
                         double &objective,
                         double &chiSquare,
                         double &reducedChiSquare)
    {
        double numerator = 0.0;
        double denominator = 0.0;
        std::vector<double> shapes;
        shapes.reserve(fitX.size());

        for (double x : fitX)
        {
            const double shape = boltzmannShape(x, temperature, barrier, charged);
            shapes.push_back(shape);
            numerator += fitY[shapes.size() - 1] * shape;
            denominator += shape * shape;
        }

        if (denominator <= 1.0e-18)
            return false;

        amplitude = numerator / denominator;
        if (!std::isfinite(amplitude) || amplitude <= 0.0)
            return false;

        objective = 0.0;
        chiSquare = 0.0;

        for (int i = 0; i < int(fitY.size()); ++i)
        {
            const double expected = amplitude * shapes[i];
            const double diff = fitY[i] - expected;
            objective += diff * diff;
            chiSquare += (diff * diff) / std::max(1.0, expected);
        }

        reducedChiSquare = chiSquare / double(std::max(1, int(fitY.size()) - 2));
        return std::isfinite(objective) && std::isfinite(chiSquare) && std::isfinite(reducedChiSquare);
    };

    constexpr double tMin = 0.2;
    constexpr double tMax = 10.0;
    constexpr int gridSteps = 800;

    double bestA = 0.0;
    double bestT = 0.0;
    double bestObjective = 1.0e300;
    double bestChi = 0.0;
    double bestReducedChi = 0.0;

    auto tryRange = [&](double rangeMin, double rangeMax, int steps)
    {
        for (int i = 0; i <= steps; ++i)
        {
            const double frac = double(i) / double(steps);
            const double temperature = rangeMin + frac * (rangeMax - rangeMin);

            double amplitude = 0.0;
            double objective = 0.0;
            double chi = 0.0;
            double reducedChi = 0.0;
            if (!evaluateT(temperature, amplitude, objective, chi, reducedChi))
                continue;

            if (objective < bestObjective)
            {
                bestObjective = objective;
                bestA = amplitude;
                bestT = temperature;
                bestChi = chi;
                bestReducedChi = reducedChi;
            }
        }
    };

    tryRange(tMin, tMax, gridSteps);

    if (bestT > 0.0)
    {
        const double coarseStep = (tMax - tMin) / double(gridSteps);
        const double localMin = std::max(tMin, bestT - 8.0 * coarseStep);
        const double localMax = std::min(tMax, bestT + 8.0 * coarseStep);
        tryRange(localMin, localMax, 500);
    }

    result.valid = bestT > 0.0 &&
                   std::isfinite(bestA) &&
                   std::isfinite(bestT) &&
                   bestA > 0.0;
    if (!result.valid)
    {
        result.unavailable = true;
        result.unavailableReason = "Boltzmann fit unavailable";
        return result;
    }

    result.amplitude = bestA;
    result.temperatureMeV = bestT;
    result.chiSquare = bestChi;
    result.reducedChiSquare = bestReducedChi;
    return result;
}

class OneDPlotWidget : public QWidget
{
public:
    OneDPlotWidget(const AngularDistEntry &entry,
                   const PaceAngularAccum &acc,
                   int plotKind,
                   double sigmaTotal,
                   int nEvents,
                   int sourceZ,
                   int sourceA,
                   bool allowGaussianEnergyOverlay,
                   bool allowBoltzmannEnergyFit,
                   OverlayParticleKind boltzmannKind,
                   QWidget *parent = nullptr)
        : QWidget(parent),
        m_entry(entry),
        m_acc(acc),
        m_plotKind(plotKind),
        m_sigmaTotal(sigmaTotal),
        m_nEvents(std::max(1, nEvents)),
        m_sourceZ(sourceZ),
        m_sourceA(sourceA),
        m_allowGaussianEnergyOverlay(allowGaussianEnergyOverlay),
        m_allowBoltzmannEnergyFit(allowBoltzmannEnergyFit),
        m_boltzmannKind(boltzmannKind)
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

        const bool useCmEnergy =
            (m_plotKind == PlotCountsVsEnergy && !m_entry.cmEnergy.empty());
        const std::vector<float> &energySamples =
            useCmEnergy ? m_entry.cmEnergy : m_entry.kineticEnergy;

        const int n = (m_plotKind == PlotCountsVsEnergy)
                          ? int(energySamples.size())
                          : int(std::min(m_entry.thetaDeg.size(), m_entry.kineticEnergy.size()));

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
            const bool useDataEnergyBins = useCmEnergy || m_allowBoltzmannEnergyFit;
            double binWidth = 1.0;

            if (useDataEnergyBins)
            {
                double maxEnergy = 0.0;
                for (float energy : energySamples)
                {
                    if (std::isfinite(double(energy)) && energy > 0.0f)
                        maxEnergy = std::max(maxEnergy, double(energy));
                }

                maxX = (maxEnergy > 0.0) ? std::ceil(maxEnergy * 1.10) : 1.0;
                if (maxX <= 0.0)
                    maxX = 1.0;

                binWidth = (maxX <= 80.0) ? 1.0 : 2.0;
                bins = int(std::ceil(maxX / binWidth));
                bins = std::clamp(bins, 20, 80);
                binWidth = maxX / double(bins);

                values.assign(bins, 0.0);
                xCenters.assign(bins, 0.0);

                for (int i = 0; i < bins; ++i)
                    xCenters[i] = minX + (double(i) + 0.5) * binWidth;
            }
            else
            {
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
            }

            for (int i = 0; i < n; ++i)
            {
                const double energy = energySamples[i];
                if (!std::isfinite(energy) || energy < minX)
                    continue;
                int bin = 0;

                if (useDataEnergyBins)
                    bin = int((energy - minX) / binWidth);
                else if (energy >= m_acc.ELOW)
                    bin = std::min(30, int((energy - m_acc.ELOW) / m_acc.DELE) + 1);

                if (bin >= bins)
                    bin = bins - 1;
                if (bin >= 0 && bin < bins)
                    values[bin] += 1.0;
            }

            title = "N = f(E)";
            xTitle = useCmEnergy ? "C.M. / source energy (MeV)" : "Lab energy (MeV)";
            yTitle = "Counts N";
        }

        double maxY = 0.0;

        for (double v : values)
            maxY = std::max(maxY, v);

        const bool showGaussian =
            (m_plotKind == PlotCountsVsTheta ||
             (m_plotKind == PlotCountsVsEnergy && m_allowGaussianEnergyOverlay));

        GaussianOverlayStats gaussianStats;

        if (showGaussian)
            gaussianStats = computeGaussianOverlayStats(xCenters,
                                                        values,
                                                        m_plotKind == PlotCountsVsEnergy);

        if (gaussianStats.valid)
            maxY = std::max(maxY, gaussianStats.amplitude * 1.15);

        BoltzmannOverlayCurve boltzmannCurve;
        const bool showBoltzmann =
            m_plotKind == PlotCountsVsEnergy &&
            m_allowBoltzmannEnergyFit &&
            m_boltzmannKind != OverlayNone;

        if (showBoltzmann)
        {
            const bool charged = m_boltzmannKind != OverlayNeutron;
            const CoulombBarrierInfo barrier =
                coulombBarrierForParticle(m_boltzmannKind, m_sourceZ, m_sourceA);
            const BoltzmannFitResult fit =
                fitBoltzmannToHistogram(xCenters, values, m_boltzmannKind, barrier.valueMeV);

            boltzmannCurve.drawBarrier = charged;
            boltzmannCurve.barrierMeV = barrier.valueMeV;
            boltzmannCurve.barrierEstimated = barrier.estimated;
            boltzmannCurve.unavailable = fit.unavailable;
            boltzmannCurve.unavailableReason = fit.unavailableReason;

            if (fit.valid)
            {
                const double startX = charged ? std::max(minX, barrier.valueMeV + 1.0e-6) : minX;

                boltzmannCurve.valid = startX < maxX;
                boltzmannCurve.amplitude = fit.amplitude;
                boltzmannCurve.temperatureMeV = fit.temperatureMeV;
                boltzmannCurve.chiSquare = fit.chiSquare;
                boltzmannCurve.reducedChiSquare = fit.reducedChiSquare;

                if (boltzmannCurve.valid)
                {
                    const int samples = 700;
                    for (int i = 0; i <= samples; ++i)
                    {
                        const double frac = double(i) / double(samples);
                        const double xValue = startX + frac * (maxX - startX);
                        const double yValue =
                            fit.amplitude *
                            boltzmannShape(xValue, fit.temperatureMeV, barrier.valueMeV, charged);
                        boltzmannCurve.maxY = std::max(boltzmannCurve.maxY, yValue);
                        boltzmannCurve.points << QPointF(xValue, yValue);
                    }
                }
            }
            else if (fit.unavailable)
            {
                boltzmannCurve.unavailable = true;
                boltzmannCurve.unavailableReason = fit.unavailableReason;
            }
        }

        if (boltzmannCurve.valid)
            maxY = std::max(maxY, boltzmannCurve.maxY * 1.15);

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

        auto mapX = [&](double xValue)
        {
            const double frac = (xValue - minX) / (maxX - minX);
            return plotRect.left() + int(frac * plotRect.width());
        };

        QPen gaussianPen(QColor(220, 40, 40, 155), 2, Qt::DotLine);
        gaussianPen.setCapStyle(Qt::RoundCap);
        QPen boltzmannPen(QColor(210, 75, 30), 2, Qt::SolidLine);
        boltzmannPen.setCapStyle(Qt::RoundCap);
        QPen barrierPen(QColor(80, 80, 80), 2, Qt::DashLine);

        auto drawOneDLegend = [&]()
        {
            const bool includeGaussian = showGaussian && gaussianStats.valid;
            const bool includeBoltzmann = boltzmannCurve.valid;
            const bool includeBoltzmannUnavailable =
                showBoltzmann && boltzmannCurve.unavailable && !boltzmannCurve.unavailableReason.isEmpty();

            p.save();

            const int legendW = (includeBoltzmann || includeBoltzmannUnavailable) ? 270 : 210;
            int legendH = 36;
            if (includeGaussian)
                legendH += 105;
            if (includeBoltzmann)
                legendH += boltzmannCurve.drawBarrier ? 100 : 85;
            else if (includeBoltzmannUnavailable)
                legendH += 35;
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

            if (!includeGaussian && !includeBoltzmann && !includeBoltzmannUnavailable)
            {
                p.restore();
                return;
            }

            int cursorY = legendY + 33;

            if (includeGaussian)
            {
                p.setPen(gaussianPen);
                p.drawLine(legendX + 10, cursorY, legendX + 42, cursorY);
                p.setPen(Qt::black);

                p.drawText(legendX + 10,
                           cursorY + 15,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Amp = " + formatGaussianNumber(gaussianStats.amplitude));

                p.drawText(legendX + 10,
                           cursorY + 30,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Mean = " + formatGaussianNumber(gaussianStats.mean));

                p.drawText(legendX + 10,
                           cursorY + 45,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Sigma = " + formatGaussianNumber(gaussianStats.sigma));

                p.drawText(legendX + 10,
                           cursorY + 60,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Chi^2 = " + formatGaussianNumber(gaussianStats.chiSquare));

                p.drawText(legendX + 10,
                           cursorY + 75,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Area = " + formatGaussianNumber(gaussianStats.area));

                p.drawText(legendX + 50,
                           cursorY - 9,
                           legendW - 58,
                           16,
                           Qt::AlignLeft,
                           "Gaussian fit");

                cursorY += 105;
            }

            if (includeBoltzmann)
            {
                p.setPen(boltzmannPen);
                p.drawLine(legendX + 10, cursorY, legendX + 42, cursorY);
                p.setPen(Qt::black);

                p.drawText(legendX + 50,
                           cursorY - 9,
                           legendW - 58,
                           16,
                           Qt::AlignLeft,
                           "Boltzmann fit");

                p.drawText(legendX + 10,
                           cursorY + 15,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "A = " + formatGaussianNumber(boltzmannCurve.amplitude));

                p.drawText(legendX + 10,
                           cursorY + 30,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "T_fit = " + QString::number(boltzmannCurve.temperatureMeV, 'g', 3) + " MeV");

                int infoY = cursorY + 45;
                if (boltzmannCurve.drawBarrier)
                {
                    p.drawText(legendX + 10,
                               infoY,
                               legendW - 20,
                               14,
                               Qt::AlignLeft,
                               "Coulomb barrier B = " +
                                   QString::number(boltzmannCurve.barrierMeV, 'g', 3) + " MeV (" +
                                   QString(boltzmannCurve.barrierEstimated ? "estimated" : "fallback") + ")");
                    infoY += 15;
                }

                p.drawText(legendX + 10,
                           infoY,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Reduced Chi^2 = " +
                               formatGaussianNumber(boltzmannCurve.reducedChiSquare));
            }
            else if (includeBoltzmannUnavailable)
            {
                p.setPen(boltzmannPen);
                p.drawLine(legendX + 10, cursorY, legendX + 42, cursorY);
                p.setPen(Qt::black);

                p.drawText(legendX + 50,
                           cursorY - 9,
                           legendW - 58,
                           16,
                           Qt::AlignLeft,
                           "Boltzmann fit");

                p.drawText(legendX + 10,
                           cursorY + 15,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           boltzmannCurve.unavailableReason);
            }

            p.restore();
        };

        if (showGaussian && gaussianStats.valid)
        {
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

        if (showBoltzmann &&
            boltzmannCurve.drawBarrier &&
            boltzmannCurve.barrierMeV >= minX &&
            boltzmannCurve.barrierMeV <= maxX)
        {
            const int barrierX = mapX(boltzmannCurve.barrierMeV);
            p.setPen(barrierPen);
            p.drawLine(barrierX, plotRect.top(), barrierX, plotRect.bottom());

            p.save();
            QFont barrierFont = p.font();
            barrierFont.setPointSize(8);
            p.setFont(barrierFont);
            p.setPen(QColor(70, 70, 70));
            p.translate(barrierX + 8, plotRect.top() + 8);
            p.rotate(-90);
            p.drawText(QRect(-120, 0, 120, 16),
                       Qt::AlignRight | Qt::AlignVCenter,
                       "Coulomb barrier B = " +
                           QString::number(boltzmannCurve.barrierMeV, 'g', 3) + " MeV");
            p.restore();
        }

        if (boltzmannCurve.valid)
        {
            QPolygonF scaledBoltzmannCurve;
            for (const QPointF &point : boltzmannCurve.points)
                scaledBoltzmannCurve << QPointF(mapX(point.x()), mapY(point.y()));

            p.setPen(boltzmannPen);
            p.drawPolyline(scaledBoltzmannCurve);
        }

        drawOneDLegend();

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
    int m_sourceZ = 0;
    int m_sourceA = 0;
    bool m_allowGaussianEnergyOverlay = true;
    bool m_allowBoltzmannEnergyFit = false;
    OverlayParticleKind m_boltzmannKind = OverlayNone;
};

class CMSpectraPlotWidget : public QWidget
{
public:
    explicit CMSpectraPlotWidget(const std::vector<CMSpectraSeries> &series,
                                 QWidget *parent = nullptr)
        : QWidget(parent), m_series(series)
    {
        setMinimumSize(1120, 760);
        resize(1120, 760);
    }

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter p(this);
        p.setRenderHint(QPainter::Antialiasing, true);
        p.fillRect(rect(), Qt::white);

        if (m_series.empty())
        {
            p.drawText(rect(), Qt::AlignCenter, "No C.M. spectra data");
            return;
        }

        QFont titleFont = p.font();
        titleFont.setPointSize(12);
        titleFont.setBold(true);
        p.setFont(titleFont);
        p.setPen(Qt::black);
        p.drawText(0, 10, width(), 26, Qt::AlignCenter,
                   "C.M. Spectra of Emitted Particles");

        const int outerMargin = 18;
        const int titleBottom = 48;
        const int panelGap = 20;
        const int panelW = (width() - 2 * outerMargin - panelGap) / 2;
        const int panelH = (height() - titleBottom - outerMargin - panelGap) / 2;

        for (int i = 0; i < int(m_series.size()); ++i)
        {
            const int row = i / 2;
            const int col = i % 2;
            QRect panelRect(outerMargin + col * (panelW + panelGap),
                            titleBottom + row * (panelH + panelGap),
                            panelW,
                            panelH);
            drawHistogramPanel(p, panelRect, m_series[i]);
        }
    }

private:
    void drawHistogramPanel(QPainter &p,
                            const QRect &panelRect,
                            const CMSpectraSeries &series)
    {
        const int left = 58;
        const int right = 14;
        const int top = 34;
        const int bottom = 58;
        QRect plotRect(panelRect.left() + left,
                       panelRect.top() + top,
                       panelRect.width() - left - right,
                       panelRect.height() - top - bottom);

        int maxCount = 0;
        for (int bin = 1; bin <= kCMSpectraMaxEnergyBins; ++bin)
            maxCount = std::max(maxCount, series.counts[bin]);

        const double maxY = std::max(1.0, std::ceil(double(maxCount) * 1.15));

        auto mapY = [&](double yValue)
        {
            return plotRect.bottom() - int((yValue / maxY) * plotRect.height());
        };

        QFont panelTitleFont = p.font();
        panelTitleFont.setPointSize(10);
        panelTitleFont.setBold(true);
        p.setFont(panelTitleFont);
        p.setPen(Qt::black);
        p.drawText(panelRect.left(),
                   panelRect.top() + 4,
                   panelRect.width(),
                   18,
                   Qt::AlignCenter,
                   series.label);

        QFont axisFont = p.font();
        axisFont.setPointSize(8);
        axisFont.setBold(false);
        p.setFont(axisFont);

        p.setPen(QPen(QColor(150, 150, 150), 1));
        p.drawRect(plotRect);

        const int yTicks = 4;
        for (int i = 0; i <= yTicks; ++i)
        {
            const double frac = double(i) / double(yTicks);
            const double yValue = frac * maxY;
            const int y = plotRect.bottom() - int(frac * plotRect.height());

            p.setPen(QPen(QColor(232, 232, 232), 1));
            p.drawLine(plotRect.left(), y, plotRect.right(), y);

            p.setPen(Qt::black);
            p.drawLine(plotRect.left() - 4, y, plotRect.left(), y);
            p.drawText(panelRect.left() + 2,
                       y - 9,
                       left - 8,
                       18,
                       Qt::AlignRight | Qt::AlignVCenter,
                       QString::number(int(std::round(yValue))));
        }

        const double barW = double(plotRect.width()) / double(kCMSpectraMaxEnergyBins);
        QColor fillColor = series.color;
        fillColor.setAlpha(170);
        QColor edgeColor = series.color.darker(115);
        edgeColor.setAlpha(210);

        for (int bin = 1; bin <= kCMSpectraMaxEnergyBins; ++bin)
        {
            const int count = series.counts[bin];
            const int x = plotRect.left() + int((bin - 1) * barW);
            const int nextX = plotRect.left() + int(bin * barW);
            const int y = mapY(count);
            QRect bar(x + 1,
                      y,
                      std::max(1, nextX - x - 2),
                      std::max(0, plotRect.bottom() - y));

            if (count > 0)
            {
                p.fillRect(bar, fillColor);
                p.setPen(QPen(edgeColor, 1));
                p.drawRect(bar);
            }
        }

        p.setPen(QPen(Qt::black, 1));
        p.drawLine(plotRect.bottomLeft(), plotRect.bottomRight());
        p.drawLine(plotRect.bottomLeft(), plotRect.topLeft());

        for (int binEdge = 0; binEdge <= kCMSpectraMaxEnergyBins; binEdge += 10)
        {
            const int x = plotRect.left() + int((double(binEdge) / kCMSpectraMaxEnergyBins) * plotRect.width());
            p.drawLine(x, plotRect.bottom(), x, plotRect.bottom() + 4);
            p.drawText(x - 18,
                       plotRect.bottom() + 6,
                       36,
                       16,
                       Qt::AlignCenter,
                       QString::number(binEdge));
        }

        const std::array<int, 6> labelBins = {1, 11, 21, 31, 41, 50};
        for (int bin : labelBins)
        {
            const double xFrac = (double(bin) - 0.5) / double(kCMSpectraMaxEnergyBins);
            const int x = plotRect.left() + int(xFrac * plotRect.width());
            const int eLow = bin - 1;
            const int eHigh = (bin == kCMSpectraMaxEnergyBins) ? 99 : bin;

            p.drawText(x - 24,
                       plotRect.bottom() + 25,
                       48,
                       16,
                       Qt::AlignCenter,
                       QString::number(eLow) + "-" + QString::number(eHigh));
        }

        QRect legendRect(plotRect.right() - 148, plotRect.top() + 8, 140, 50);
        p.fillRect(legendRect, QColor(255, 255, 255, 225));
        p.setPen(QPen(QColor(150, 150, 150), 1));
        p.drawRect(legendRect);
        p.fillRect(QRect(legendRect.left() + 8, legendRect.top() + 10, 18, 10), fillColor);
        p.setPen(QPen(edgeColor, 1));
        p.drawRect(QRect(legendRect.left() + 8, legendRect.top() + 10, 18, 10));
        p.setPen(Qt::black);
        p.drawText(legendRect.left() + 32,
                   legendRect.top() + 5,
                   legendRect.width() - 40,
                   16,
                   Qt::AlignLeft,
                   "Counts/bin");
        p.drawText(legendRect.left() + 8,
                   legendRect.top() + 25,
                   legendRect.width() - 16,
                   16,
                   Qt::AlignLeft,
                   "Total " + QString::number(series.counts[kCMSpectraTotalIndex]) +
                       ", Avg " + QString::number(series.averageEnergy, 'f', 2));

        if (series.counts[kCMSpectraTotalIndex] == 0)
        {
            p.setPen(QColor(90, 90, 90));
            p.drawText(plotRect, Qt::AlignCenter, "No counts");
        }

        p.setPen(Qt::black);
        p.drawText(plotRect.left(),
                   panelRect.bottom() - 20,
                   plotRect.width(),
                   16,
                   Qt::AlignCenter,
                   "Energy bin Ex (MeV)");

        p.save();
        p.translate(panelRect.left() + 16, plotRect.top() + plotRect.height() / 2);
        p.rotate(-90);
        p.drawText(QRect(-plotRect.height() / 2, -20, plotRect.height(), 20),
                   Qt::AlignCenter,
                   "Counts");
        p.restore();
    }
    std::vector<CMSpectraSeries> m_series;
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
    std::array<std::array<int, kCMSpectraTotalIndex + 1>, 5> NSPC{};
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
            if (bin > kCMSpectraMaxEnergyBins) bin = kCMSpectraMaxEnergyBins;

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
        for (int bin = 1; bin <= kCMSpectraMaxEnergyBins; ++bin)
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
    html += "<p align=\"center\"><a href=\"gemini://plot_cm_spectra\">Plot C.M. spectra</a></p>";
    html += "<table class=\"cm-spectra\" align=\"center\">";
    html += "<tr><th>Ex(MeV)</th><th>Neut</th><th>Prot</th><th>Alpha</th><th>Gamma</th></tr>";

    for (int bin = 1; bin <= kCMSpectraMaxEnergyBins; ++bin)
    {
        const bool printRow =
            NSPC[1][bin] != 0 ||
            NSPC[2][bin] != 0 ||
            NSPC[3][bin] != 0 ||
            NSPC[4][bin] != 0;

        const double particleEnergy = double(bin) - 0.5;
        for (int mode = 1; mode <= 4; ++mode)
        {
            NSPC[mode][kCMSpectraTotalIndex] += NSPC[mode][bin];
            averageEnergy[mode] += particleEnergy * double(NSPC[mode][bin]);
        }

        if (!printRow) continue;

        const int eLow = bin - 1;
        const int eHigh = (bin == kCMSpectraMaxEnergyBins) ? 99 : bin;

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
        averageEnergy[mode] /= double(NSPC[mode][kCMSpectraTotalIndex]) + 1.0e-9;

    html += "<tr><th>Total</th>";
    for (int mode = 1; mode <= 4; ++mode)
        html += "<th>" + QString::number(NSPC[mode][kCMSpectraTotalIndex]) + "</th>";
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
    const bool isImf = title.contains("IMF", Qt::CaseInsensitive);
    const QString fragmentKind = isImf ? "IMF fragment" : "residual nucleus";
    const QString fragmentKindPlural = isImf ? "IMF fragments" : "residual nuclei";
    const QString velocityLabel = isImf ? "Fragment velocity/c" : "Residual velocity/c";

    html += "<!DOCTYPE html><html><head><meta charset=\"utf-8\">";
    html += "<style>"
            "body{font-family:Sans-Serif;}"
            "table{border-collapse:collapse; margin:10px 0;}"
            "th,td{border:1px solid #909090; padding:4px 6px;}"
            "th{background:#ececec;}"
            ".total-row th,.total-row td{background:#ececec; font-weight:bold;}"
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

    html += isImf
                ? "<p>&nbsp;</p><h2 align=\"center\">Angular distribution results(IMF)</h2>"
                : "<p>&nbsp;</p><h2 align=\"center\">Angular distribution results</h2>";
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

    if (totalAngularCount(allAcc) <= 0
        && !hasAngularSamples(neutronEntry)
        && !hasAngularSamples(protonEntry)
        && !hasAngularSamples(alphaEntry)
        && !hasAngularSamples(gammaEntry))
    {
        html += "<p>No angular-distribution events were found for output.</p></body></html>";
        return html;
    }

    int displayIndex = 1;
    int plotIndex = 0;

    for (const auto &item : selectedAcc)
    {
        const PaceSelectedResidue &r = item.first;
        const PaceAngularAccum &acc = item.second;
        if (totalAngularCount(acc) <= 0) continue;

        html += buildHeaderForResiduePACE(displayIndex, r.z, r.n, fragmentKind);
        appendMainPaceTable(html, acc, r.a, sigmaTotalMb, nCascades, inputMode, false, plotIndex,
                            velocityLabel);
        appendSecondPaceTableIfNeeded(html, acc, plotIndex);
        displayIndex++;
        plotIndex++;
    }

    if (totalAngularCount(allAcc) > 0)
    {
        html += buildHeaderAllPACE(displayIndex, fragmentKindPlural);
        appendMainPaceTable(html, allAcc, compoundA, sigmaTotalMb, nCascades, inputMode, true, plotIndex,
                            velocityLabel);
        appendSecondPaceTableIfNeeded(html, allAcc, plotIndex);
        displayIndex++;
        plotIndex++;
    }

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
        if (totalAngularCount(acc) <= 0) return;

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
                                                     int compoundZ,
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
    compoundZForPlots(compoundZ),
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

    if (url.host() == "plot_cm_spectra")
    {
        openCMSpectraPlotWindow();
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

void AngularDistributionWidget::openCMSpectraPlotWindow()
{
    if (!hasCMSpectraSamples(neutronEntryForPlots,
                             protonEntryForPlots,
                             alphaEntryForPlots,
                             gammaEntryForPlots))
    {
        return;
    }

    const std::vector<CMSpectraSeries> series =
        buildCMSpectraSeries(neutronEntryForPlots,
                             protonEntryForPlots,
                             alphaEntryForPlots,
                             gammaEntryForPlots);

    if (series.empty()) return;

    QDialog *dlg = new QDialog(this);
    dlg->setAttribute(Qt::WA_DeleteOnClose, true);
    dlg->setWindowIcon(QIcon(":/Gemini_logo.png"));
    dlg->setWindowTitle("Gemini: C.M. spectra histograms");
    dlg->setModal(false);

    QVBoxLayout *mainLayout = new QVBoxLayout(dlg);

    QHBoxLayout *topRow = new QHBoxLayout;

    QLabel *topLabel = new QLabel("<h2 style='color:#1d4f91; margin:0;'>C.M. spectra histograms</h2>"
                                  "<p style='margin:2px 0 0 0;'><em>"
                                  + plotTitle.toHtmlEscaped() +
                                  "</em> - one histogram per emitted particle</p>");
    topLabel->setTextFormat(Qt::RichText);
    topRow->addWidget(topLabel);
    topRow->addStretch();

    QPushButton *savePngButton = new QPushButton("Save as PNG");
    savePngButton->setFixedHeight(26);
    savePngButton->setMaximumWidth(110);
    topRow->addWidget(savePngButton, 0, Qt::AlignTop | Qt::AlignRight);

    mainLayout->addLayout(topRow);

    CMSpectraPlotWidget *plot = new CMSpectraPlotWidget(series);
    mainLayout->addWidget(plot);

    QObject::connect(savePngButton, &QPushButton::clicked, dlg, [plot, dlg]()
                     {
                         const QString fileName = QFileDialog::getSaveFileName(
                             dlg,
                             QObject::tr("Save Plot as PNG"),
                             QString(),
                             QObject::tr("PNG Image (*.png)")
                             );

                         if (fileName.isEmpty()) return;

                         QPixmap pixmap(plot->size());
                         pixmap.fill(Qt::white);
                         plot->render(&pixmap);
                         pixmap.save(fileName, "PNG");
                     });

    dlg->setLayout(mainLayout);
    dlg->resize(1200, 880);
    dlg->show();
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
                               mdirForPlots,
                               plotTitle.contains("IMF", Qt::CaseInsensitive));

    if (tables.empty()) return;

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

        const QString displayLabel = tables[i].htmlLabel.isEmpty()
                                         ? tables[i].label.toHtmlEscaped()
                                         : tables[i].htmlLabel;
        QLabel *lbl = new QLabel("<b>" + displayLabel + "</b>");
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
                                      nEventsForPlots,
                                      compoundZForPlots,
                                      compoundAForPlots,
                                      tables[i].allowGaussianEnergyOverlay,
                                      tables[i].allowBoltzmannEnergyFit,
                                      tables[i].boltzmannKind);

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
