// gm_angular_distribution.cpp
#include "gm_angular_distribution.h"

#include <QBoxLayout>
#include <QDesktopServices>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QIcon>
#include <QLabel>
#include <QPainter>
#include <QPen>
#include <QPixmap>
#include <QPolygonF>
#include <QPushButton>
#include <QScrollArea>
#include <QTextBrowser>
#include <QTextEdit>
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
    double count = 0.0;
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

    std::array<std::array<double, 38>, 33> NR{};
    std::array<std::array<double, 38>, 5> NRW{};
};

QString fmt0(double x)
{
    return QString::number(x, 'f', 0);
}

QString fmt1(double x)
{
    return QString::number(x, 'f', 1);
}

QString fmtCount(double x)
{
    if (std::abs(x - std::round(x)) < 1.0e-6)
        return QString::number(x, 'f', 0);
    return QString::number(x, 'g', 4);
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

double isum(const PaceAngularAccum &acc, int energyBin, int l1, int l2)
{
    double sum = 0.0;
    for (int l = l1; l <= l2; ++l)
        sum += acc.NR[energyBin][l];
    return sum;
}

double totalAngularCount(const PaceAngularAccum &acc)
{
    return isum(acc, 32, 1, 37);
}

double sampleWeightAt(const AngularDistEntry &entry, int index)
{
    if (index >= 0 && index < int(entry.weight.size()))
        return std::max(0.0, double(entry.weight[index]));
    return 1.0;
}

double entryWeightTotal(const AngularDistEntry &entry)
{
    const int n = int(entry.kineticEnergy.size());
    double total = 0.0;
    for (int i = 0; i < n; ++i)
        total += sampleWeightAt(entry, i);
    return total;
}

double totalEntryWeight(const std::map<std::pair<int, int>, AngularDistEntry> &entries)
{
    double totalCount = 0.0;
    for (const auto &it : entries)
        totalCount += entryWeightTotal(it.second);
    return totalCount;
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

void addToPaceAccum(PaceAngularAccum &acc,
                    double energyLabMeV,
                    double angleLabDeg,
                    double sampleWeight = 1.0)
{
    if (sampleWeight <= 0.0) return;

    int K = 1;
    if (energyLabMeV >= acc.ELOW)
        K = std::min(31, int((energyLabMeV - acc.ELOW) / acc.DELE) + 2);

    int L = std::min(37, int(angleLabDeg / acc.DELANG) + 1);
    if (L < 1) L = 1;

    acc.NR[K][L] += sampleWeight;
    acc.NR[32][L] += sampleWeight;

    int KW = 4;
    if (energyLabMeV <= acc.EWR[1]) KW = 1;
    else if (energyLabMeV <= acc.EWR[2]) KW = 2;
    else if (energyLabMeV <= acc.EWR[3]) KW = 3;

    acc.NRW[KW][L] += sampleWeight;
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
        const double count = entryWeightTotal(e);
        if (count <= 0.0) continue;

        PaceSelectedResidue r;
        r.z = e.z;
        r.n = e.n;
        r.a = e.z + e.n;
        r.count = count;
        r.percent = 100.0 * count / double(std::max(1, nCascades));
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
                              bool showVelocity,
                              bool isAll,
                              const QString &velocityLabel)
{
    if (!showVelocity || isAll) return QString();

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

void appendAngleHeader(QString &html, const PaceAngularAccum &acc, int l1, int l2, bool includeAboveLabel)
{
    html += "<table border=\"1\" cellspacing=\"0\">";
    const int angleColumns = l2 - l1 + 1 + (includeAboveLabel ? 1 : 0);
    html += "<tr><th>Energy Range</th><th colspan=\"" + QString::number(angleColumns) +
            "\">Angular range (deg)</th></tr>";
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
                         bool showVelocity,
                         bool isAll,
                         int plotIndex,
                         const QString &velocityLabel)
{
    html += buildVelocityLinePACE(acc, particleA, showVelocity, isAll, velocityLabel);

    appendAngleHeader(html, acc, 1, 36, true);

    double E2 = 0.0;
    if (acc.ELOW >= 0.01)
    {
        E2 = acc.ELOW;
        html += "<tr><td>Below " + fmt0(E2) + "</td>";
        for (int k = 1; k <= 37; ++k)
            html += (acc.NR[1][k] > 0) ? "<td align=\"center\">" + fmtCount(acc.NR[1][k]) + "</td>" : "<td></td>";
        html += "</tr>";
    }

    for (int K = 2; K <= 30; ++K)
    {
        const double E1 = (double(K) - 2.0) * acc.DELE + acc.ELOW;
        E2 = E1 + acc.DELE;

        if (isum(acc, K, 1, 37) > 0)
        {
            html += "<tr><td>" + fmt1(E1) + " - " + fmt1(E2) + "</td>";
            for (int k = 1; k <= 37; ++k)
                html += (acc.NR[K][k] > 0) ? "<td align=\"center\">" + fmtCount(acc.NR[K][k]) + "</td>" : "<td></td>";
            html += "</tr>";
        }
    }

    html += "<tr><td>Above " + fmt0(E2) + "</td>";
    for (int k = 1; k <= 37; ++k)
        html += (acc.NR[31][k] > 0) ? "<td align=\"center\">" + fmtCount(acc.NR[31][k]) + "</td>" : "<td></td>";
    html += "</tr>";

    std::array<double, 38> DSIG{};
    const double FAC = sigmaTotalMb / (kTwoPiPACE * double(std::max(1, nCascades)));
    for (int i = 1; i <= 36; ++i)
    {
        const double TET = (double(i) - 0.5) * kDegToRad;
        DSIG[i] = FAC * double(acc.NR[32][i]) / (std::sin(TET) * kDegToRad);
    }

    html += "<tr class=\"total-row\"><th>Total</th>";
    for (int k = 1; k <= 37; ++k)
        html += (acc.NR[32][k] > 0) ? "<th align=\"center\">" + fmtCount(acc.NR[32][k]) + "</th>" : "<th></th>";
    html += "</tr>";

    html += "<tr><td>d&sigma;/d&Omega;</td>";
    for (int k = 1; k <= 36; ++k)
        html += (DSIG[k] > 0.0) ? "<td>" + QString::number(DSIG[k], 'g', 2) + "</td>" : "<td>0.00</td>";
    html += "<td></td>";
    html += "</tr>";

    for (int i = 1; i <= 4; ++i)
    {
        html += (i != 4)
                    ? "<tr><td>" + fmt0(acc.EWR1[i]) + " - " + fmt0(acc.EWR2[i]) + "</td>"
                    : "<tr><td>Above " + fmt0(acc.EWR1[4]) + "</td>";

        for (int k = 1; k <= 37; ++k)
            html += (acc.NRW[i][k] > 0) ? "<td align=\"center\">" + fmtCount(acc.NRW[i][k]) + "</td>" : "<td></td>";
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
        addToPaceAccum(acc, entry.kineticEnergy[i], entry.thetaDeg[i], sampleWeightAt(entry, i));
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
    bool isParticlePlot = false;
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
            const double weight = sampleWeightAt(e, i);
            addToPaceAccum(allAcc, e.kineticEnergy[i], e.thetaDeg[i], weight);
            if (selectedIndex >= 0)
                addToPaceAccum(selectedAcc[selectedIndex].second, e.kineticEnergy[i], e.thetaDeg[i], weight);
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
        for (int i = 0; i < int(e.kineticEnergy.size()); ++i)
            all.weight.push_back(float(sampleWeightAt(e, i)));
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
        table.isParticlePlot = true;
        initPaceAccum(table.acc, compoundA, compoundExcitationMeV, recoilBetaCN, mdir);
        setParticleEnergyInterval(table.acc);
        appendEntryToAccum(table.acc, entry);
        if (totalAngularCount(table.acc) <= 0) return;
        if (overlayParticleKind != OverlayNone)
        {
            table.allowGaussianEnergyOverlay = false;
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

        std::vector<std::vector<double>> counts(
            energyBins,
            std::vector<double>(thetaBins, 0.0)
            );

        double maxCount = 0.0;
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

            counts[energyBin][thetaBin] += sampleWeightAt(m_entry, i);
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
                const double count = counts[e][t];

                int x = plotRect.left() + int(t * cellW);
                int y = plotRect.bottom() - int((e + 1) * cellH);

                QRect cell(x, y, int(cellW) + 1, int(cellH) + 1);

                if (count <= 0.0)
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

double poissonDevianceTerm(double observed, double expected)
{
    const double mu = std::max(expected, 1.0e-12);
    const double n = std::max(0.0, observed);

    if (n <= 0.0)
        return 2.0 * mu;

    return 2.0 * (mu - n + n * std::log(n / mu));
}

double inferUniformBinWidth(const std::vector<double> &xCenters)
{
    if (xCenters.size() < 2)
        return 1.0;

    std::vector<double> spacings;
    spacings.reserve(xCenters.size() - 1);
    for (int i = 1; i < int(xCenters.size()); ++i)
    {
        const double dx = xCenters[i] - xCenters[i - 1];
        if (std::isfinite(dx) && dx > 0.0)
            spacings.push_back(dx);
    }

    if (spacings.empty())
        return 1.0;

    std::sort(spacings.begin(), spacings.end());
    return spacings[spacings.size() / 2];
}

void computeGaussianQualityAndArea(GaussianOverlayStats &stats,
                                   const std::vector<double> &xCenters,
                                   const std::vector<double> &yValues)
{
    stats.reducedChiSquare = 0.0;
    stats.area = 0.0;

    if (!stats.valid || stats.sigma <= 1.0e-12) return;

    int usedPoints = 0;
    double pearsonSum = 0.0;

    for (int i = 0; i < int(yValues.size()); ++i)
    {
        const double observed = std::max(0.0, yValues[i]);
        const double expected = gaussianValue(stats.amplitude,
                                              stats.mean,
                                              stats.sigma,
                                              xCenters[i]);

        const double denom = std::max(1.0, expected);
        const double diff = observed - expected;

        pearsonSum += (diff * diff) / denom;
        usedPoints++;
    }

    const int dof = std::max(1, usedPoints - 3);
    stats.reducedChiSquare = pearsonSum / double(dof);

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

    for (int i = 0; i < int(yValues.size()); ++i)
    {
        if (yValues[i] <= 0.0) continue;

        if (firstPopulated < 0)
            firstPopulated = i;
        lastPopulated = i;
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

    const double binWidth = inferUniformBinWidth(fitX);

    auto evaluateGaussian = [&](double testMu,
                                double testSigma,
                                double &bestAForShape,
                                double &deviance)
    {
        bestAForShape = 0.0;
        deviance = 1.0e300;

        if (testSigma <= 1.0e-12 || !std::isfinite(testMu) || !std::isfinite(testSigma))
            return false;

        std::vector<double> shape;
        shape.reserve(fitY.size());
        double shapeSum = 0.0;
        double observedSum = 0.0;

        for (int i = 0; i < int(fitY.size()); ++i)
        {
            const double s = gaussianValue(1.0, testMu, testSigma, fitX[i]);
            shape.push_back(s);
            shapeSum += s;
            observedSum += std::max(0.0, fitY[i]);
        }

        if (shapeSum <= 1.0e-18 || observedSum <= 0.0)
            return false;

        bestAForShape = observedSum / shapeSum;
        if (!std::isfinite(bestAForShape) || bestAForShape <= 0.0)
            return false;

        deviance = 0.0;
        for (int i = 0; i < int(fitY.size()); ++i)
        {
            const double expected = std::max(1.0e-12, bestAForShape * shape[i]);
            const double observed = std::max(0.0, fitY[i]);
            deviance += poissonDevianceTerm(observed, expected);
        }

        return std::isfinite(deviance);
    };

    auto considerGaussian = [&](double testMu,
                                double testSigma,
                                double &bestA,
                                double &bestMu,
                                double &bestSigma,
                                double &bestDeviance)
    {
        double testA = 0.0;
        double testDeviance = 0.0;
        if (!evaluateGaussian(testMu, testSigma, testA, testDeviance))
            return;

        if (testDeviance < bestDeviance)
        {
            bestA = testA;
            bestMu = testMu;
            bestSigma = testSigma;
            bestDeviance = testDeviance;
        }
    };

    double bestA = A;
    double bestMu = mu;
    double bestSigma = sigma;
    double bestDeviance = 1.0e300;

    considerGaussian(mu, sigma, bestA, bestMu, bestSigma, bestDeviance);

    const double muMin = xMin;
    const double muMax = xMax;
    const double sigmaMin = std::max(0.35 * binWidth, fallbackSigma * 0.35);
    const double sigmaMax = std::max(sigmaMin * 1.05, std::max(xRange * 1.5, fallbackSigma));

    const int coarseMuSteps = std::max(1, int(std::min<double>(80.0, std::max(20.0, fitX.size() * 2.0))));
    const int coarseSigmaSteps = 80;
    for (int mi = 0; mi <= coarseMuSteps; ++mi)
    {
        const double testMu =
            muMin + (muMax - muMin) * double(mi) / double(coarseMuSteps);

        for (int si = 0; si <= coarseSigmaSteps; ++si)
        {
            const double frac = double(si) / double(coarseSigmaSteps);
            const double testSigma =
                sigmaMin * std::pow(sigmaMax / sigmaMin, frac);
            considerGaussian(testMu, testSigma, bestA, bestMu, bestSigma, bestDeviance);
        }
    }

    for (int pass = 0; pass < 3; ++pass)
    {
        const double muHalfWindow =
            std::max(0.25 * binWidth, xRange / std::pow(4.0, double(pass + 1)));
        const double sigmaScale = std::pow(1.8, 1.0 / double(pass + 1));
        const double localMuMin = std::max(xMin - binWidth, bestMu - muHalfWindow);
        const double localMuMax = std::min(xMax + binWidth, bestMu + muHalfWindow);
        const double localSigmaMin = std::max(sigmaMin, bestSigma / sigmaScale);
        const double localSigmaMax = std::min(sigmaMax, bestSigma * sigmaScale);

        for (int mi = 0; mi <= 40; ++mi)
        {
            const double testMu =
                localMuMin + (localMuMax - localMuMin) * double(mi) / 40.0;

            for (int si = 0; si <= 40; ++si)
            {
                const double frac = double(si) / 40.0;
                const double testSigma =
                    localSigmaMin * std::pow(localSigmaMax / localSigmaMin, frac);
                considerGaussian(testMu, testSigma, bestA, bestMu, bestSigma, bestDeviance);
            }
        }
    }

    A = bestA;
    mu = bestMu;
    sigma = bestSigma;

    if (!std::isfinite(A) || A <= 0.0 || !std::isfinite(mu) || !std::isfinite(sigma))
    {
        A = maxY;
        mu = sumX / sumW;
        sigma = std::max(sigma, fallbackSigma);
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
    double normalization = 0.0;
    double coldTemperatureMeV = 0.0;
    double effectiveBarrierMeV = 0.0;
    double barrierDiffusenessMeV = 0.0;
    double pearsonReducedChiSquare = 0.0;
    double maxY = 0.0;
    QPolygonF points;
};

struct HistogramMomentStats
{
    bool valid = false;
    double total = 0.0;
    double mean = 0.0;
    double median = 0.0;
    double area = 0.0;
};

struct BoltzmannFitResult
{
    bool valid = false;
    bool unavailable = false;
    QString unavailableReason;
    double normalization = 0.0;
    double coldTemperatureMeV = 0.0;
    double effectiveBarrierMeV = 0.0;
    double barrierDiffusenessMeV = 0.0;
    double pearsonReducedChiSquare = 0.0;
};

HistogramMomentStats computeHistogramMomentStats(const std::vector<double> &xCenters,
                                                 const std::vector<double> &yValues,
                                                 double binWidth)
{
    HistogramMomentStats stats;
    if (xCenters.size() != yValues.size() || xCenters.empty())
        return stats;

    double weightedSum = 0.0;
    for (int i = 0; i < int(yValues.size()); ++i)
    {
        const double y = std::max(0.0, yValues[i]);
        stats.total += y;
        weightedSum += y * xCenters[i];
    }

    if (stats.total <= 0.0)
        return stats;

    stats.valid = true;
    stats.mean = weightedSum / stats.total;
    stats.area = stats.total * std::max(0.0, binWidth);

    const double halfTotal = 0.5 * stats.total;
    double cumulative = 0.0;
    for (int i = 0; i < int(yValues.size()); ++i)
    {
        cumulative += std::max(0.0, yValues[i]);
        if (cumulative >= halfTotal)
        {
            stats.median = xCenters[i];
            break;
        }
    }

    return stats;
}

struct CoulombBarrierInfo
{
    double valueMeV = 0.0;
};

struct EvaporationModelParams
{
    double coldTemperatureMeV = 1.5;
    double effectiveBarrierMeV = 0.0;
    double barrierDiffusenessMeV = 1.0;
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

    // r0 = 1.25 fm is a conventional nuclear-radius constant for this simple
    // touching-spheres Coulomb estimate.
    constexpr double r0 = 1.25;
    constexpr double e2MeVFm = 1.4399764;
    const double radiusSum = r0 * (std::cbrt(double(aParticle)) + std::cbrt(double(aDaughter)));
    if (radiusSum <= 0.0)
        return barrier;

    const double estimated =
        e2MeVFm * double(zParticle) * double(zDaughter) / radiusSum;
    if (!std::isfinite(estimated) || estimated <= 0.0)
        return barrier;

    barrier.valueMeV = std::clamp(estimated, minBarrier, maxBarrier);
    return barrier;
}

double safeExp(double x)
{
    return std::exp(std::clamp(x, -700.0, 700.0));
}

double smoothCoulombPenetrability(double energyMeV,
                                  double effectiveBarrierMeV,
                                  double barrierDiffusenessMeV)
{
    if (effectiveBarrierMeV <= 0.0)
        return 1.0;
    if (barrierDiffusenessMeV <= 0.0)
        return energyMeV >= effectiveBarrierMeV ? 1.0 : 0.0;

    const double exponent = (effectiveBarrierMeV - energyMeV) / barrierDiffusenessMeV;
    return 1.0 / (1.0 + safeExp(exponent));
}

double evaporationDensity(double energyMeV,
                          OverlayParticleKind particleKind,
                          const EvaporationModelParams &params)
{
    if (!std::isfinite(energyMeV) || energyMeV < 0.0 ||
        params.coldTemperatureMeV <= 0.0)
    {
        return 0.0;
    }

    const bool charged = particleKind == OverlayProton || particleKind == OverlayAlpha;
    double phaseSpace = 0.0;
    if (!charged)
    {
        phaseSpace = std::sqrt(std::max(energyMeV, 0.0));
    }
    else
    {
        phaseSpace = std::sqrt(std::max(energyMeV, 0.0));
    }

    const double penetrability =
        charged
            ? smoothCoulombPenetrability(energyMeV,
                                         params.effectiveBarrierMeV,
                                         params.barrierDiffusenessMeV)
            : 1.0;

    const double density =
        phaseSpace * penetrability * safeExp(-energyMeV / params.coldTemperatureMeV);
    return std::isfinite(density) && density >= 0.0 ? density : 0.0;
}

double integrateEvaporationDensity(double lowMeV,
                                   double highMeV,
                                   OverlayParticleKind particleKind,
                                   const EvaporationModelParams &params)
{
    if (highMeV <= lowMeV)
        return 0.0;

    static constexpr double nodes[8] =
        {
            -0.9602898564975363,
            -0.7966664774136267,
            -0.5255324099163290,
            -0.1834346424956498,
             0.1834346424956498,
             0.5255324099163290,
             0.7966664774136267,
             0.9602898564975363
        };
    static constexpr double weights[8] =
        {
            0.1012285362903763,
            0.2223810344533745,
            0.3137066458778873,
            0.3626837833783620,
            0.3626837833783620,
            0.3137066458778873,
            0.2223810344533745,
            0.1012285362903763
        };

    const double center = 0.5 * (lowMeV + highMeV);
    const double halfWidth = 0.5 * (highMeV - lowMeV);
    double integral = 0.0;

    for (int i = 0; i < 8; ++i)
    {
        const double energy = center + halfWidth * nodes[i];
        integral += weights[i] * evaporationDensity(energy, particleKind, params);
    }

    integral *= halfWidth;
    return std::isfinite(integral) && integral >= 0.0 ? integral : 0.0;
}

void fitRangeForHistogram(const std::vector<double> &yValues,
                          int &firstBin,
                          int &lastBin)
{
    firstBin = -1;
    lastBin = -1;

    double total = 0.0;
    int populated = 0;
    for (int i = 0; i < int(yValues.size()); ++i)
    {
        const double y = std::max(0.0, yValues[i]);
        total += y;
        if (y > 0.0)
        {
            populated++;
            if (firstBin < 0)
                firstBin = i;
            lastBin = i;
        }
    }

    if (firstBin < 0 || lastBin < 0)
        return;

    // With enough statistics, keep the evaporation fit focused on the populated
    // spectral body and nearby empty bins. Isolated extreme tail counts are still
    // shown in the histogram, but they should not determine the peak/barrier fit.
    if (total >= 50.0 && populated >= 6)
    {
        const double lowCut = 0.005 * total;
        const double highCut = 0.995 * total;
        double cumulative = 0.0;
        int coreFirst = firstBin;
        int coreLast = lastBin;

        for (int i = 0; i < int(yValues.size()); ++i)
        {
            cumulative += std::max(0.0, yValues[i]);
            if (cumulative >= lowCut)
            {
                coreFirst = i;
                break;
            }
        }

        cumulative = 0.0;
        for (int i = 0; i < int(yValues.size()); ++i)
        {
            cumulative += std::max(0.0, yValues[i]);
            if (cumulative >= highCut)
            {
                coreLast = i;
                break;
            }
        }

        if (coreLast > coreFirst)
        {
            firstBin = coreFirst;
            lastBin = coreLast;
        }
    }

    if (firstBin > 0)
        firstBin--;
    if (lastBin + 1 < int(yValues.size()))
        lastBin++;
}

BoltzmannFitResult fitBoltzmannToHistogram(const std::vector<double> &xCenters,
                                           const std::vector<double> &yValues,
                                           OverlayParticleKind particleKind,
                                           double energyBinWidthMeV,
                                           const CoulombBarrierInfo &geometricBarrier)
{
    // Boltzmann fitting is used only as a diagnostic for emitted neutron/proton/alpha spectra.
    // It is not applied to residual nuclei or IMF fragments. In GEMINI, IMFs are
    // complex-fragment/asymmetric binary-decay products, so a simple emitted-particle
    // Boltzmann+Coulomb spectrum is not appropriate for IMF plots. The Coulomb barrier
    // used here is an approximate plotting barrier, not the full GEMINI
    // transmission-coefficient treatment.
    //
    // The Coulomb part is an effective smooth penetrability based on a simple
    // touching-spheres geometric estimate. Full GEMINI transmission coefficients,
    // emission-stage source evolution, and detailed barrier distributions are not
    // included in this plotting fit.
    BoltzmannFitResult result;

    if (xCenters.size() != yValues.size() || xCenters.empty() ||
        energyBinWidthMeV <= 0.0 || particleKind == OverlayNone)
    {
        result.unavailable = true;
        result.unavailableReason = "Boltzmann fit unavailable";
        return result;
    }

    const bool charged = particleKind == OverlayProton || particleKind == OverlayAlpha;

    int firstFitBin = -1;
    int lastFitBin = -1;
    fitRangeForHistogram(yValues, firstFitBin, lastFitBin);

    if (firstFitBin < 0 || lastFitBin < firstFitBin)
    {
        result.unavailable = true;
        result.unavailableReason = "Boltzmann fit unavailable: too few points";
        return result;
    }

    double totalObserved = 0.0;
    for (int i = firstFitBin; i <= lastFitBin; ++i)
    {
        const double y = std::max(0.0, yValues[i]);
        totalObserved += y;
    }

    if (totalObserved <= 0.0)
    {
        result.unavailable = true;
        result.unavailableReason = "Boltzmann fit unavailable: no counts";
        return result;
    }

    const int fittedBins = lastFitBin - firstFitBin + 1;

    auto evaluateModel = [&](const EvaporationModelParams &params,
                             double &normalization,
                             double &deviance,
                             double &pearsonReducedChi)
    {
        std::vector<double> integrals;
        integrals.reserve(fittedBins);
        double totalModel = 0.0;

        for (int bin = firstFitBin; bin <= lastFitBin; ++bin)
        {
            const double low = std::max(0.0, xCenters[bin] - 0.5 * energyBinWidthMeV);
            const double high = std::max(low, xCenters[bin] + 0.5 * energyBinWidthMeV);
            const double integral = integrateEvaporationDensity(low, high, particleKind, params);
            integrals.push_back(integral);
            totalModel += integral;
        }

        if (totalModel <= 1.0e-18 || !std::isfinite(totalModel))
            return false;

        normalization = totalObserved / totalModel;
        if (!std::isfinite(normalization) || normalization <= 0.0)
            return false;

        deviance = 0.0;
        double pearsonChi = 0.0;

        for (int j = 0; j < int(integrals.size()); ++j)
        {
            const int bin = firstFitBin + j;
            const double expected = std::max(1.0e-12, normalization * integrals[j]);
            const double observed = std::max(0.0, yValues[bin]);
            const double diff = observed - expected;
            deviance += poissonDevianceTerm(observed, expected);
            pearsonChi += (diff * diff) / std::max(1.0, expected);
        }

        int freeParameters = 2; // normalization plus temperature
        if (charged)
            freeParameters += 1; // effective Coulomb barrier

        const int dof = std::max(1, fittedBins - freeParameters);
        pearsonReducedChi = pearsonChi / double(dof);
        return std::isfinite(deviance) && std::isfinite(pearsonReducedChi);
    };

    EvaporationModelParams bestParams;
    double bestNormalization = 0.0;
    double bestDeviance = 1.0e300;
    double bestPearson = 0.0;

    auto evaluateCandidate = [&](const EvaporationModelParams &params,
                                 double &normalization,
                                 double &deviance,
                                 double &pearson)
    {
        if (params.coldTemperatureMeV <= 0.0 ||
            params.effectiveBarrierMeV < 0.0 ||
            params.barrierDiffusenessMeV <= 0.0)
        {
            return false;
        }

        return evaluateModel(params, normalization, deviance, pearson);
    };

    auto considerModel = [&](const EvaporationModelParams &params,
                             double normalization,
                             double deviance,
                             double pearson)
    {
        if (deviance < bestDeviance)
        {
            bestParams = params;
            bestNormalization = normalization;
            bestDeviance = deviance;
            bestPearson = pearson;
        }
    };

    const double fixedDiffusenessMeV =
        charged ? (particleKind == OverlayAlpha ? 1.5 : 1.0) : 1.0;

    auto paramsForTemperatureAndBarrier = [&](double temperatureMeV,
                                              double effectiveBarrierMeV)
    {
        EvaporationModelParams params;
        params.coldTemperatureMeV = temperatureMeV;
        params.effectiveBarrierMeV = charged ? effectiveBarrierMeV : 0.0;
        params.barrierDiffusenessMeV = fixedDiffusenessMeV;
        return params;
    };

    auto considerTemperatureAndBarrier = [&](double temperatureMeV,
                                             double effectiveBarrierMeV)
    {
        double normalization = 0.0;
        double deviance = 0.0;
        double pearson = 0.0;
        const EvaporationModelParams params =
            paramsForTemperatureAndBarrier(temperatureMeV, effectiveBarrierMeV);
        if (evaluateCandidate(params, normalization, deviance, pearson))
            considerModel(params, normalization, deviance, pearson);
    };

    constexpr double tMin = 0.2;
    constexpr double tMax = 10.0;

    double barrierMin = 0.0;
    double barrierMax = 0.0;
    if (charged)
    {
        const double geometricB = geometricBarrier.valueMeV;
        const double lowerScale = particleKind == OverlayAlpha ? 0.50 : 0.55;
        const double upperScale = particleKind == OverlayAlpha ? 1.60 : 1.55;
        const double hardMax = particleKind == OverlayAlpha ? 60.0 : 30.0;
        barrierMin = std::max(0.0, std::min(geometricB * lowerScale, geometricB - 4.0));
        barrierMax = std::min(hardMax, std::max(geometricB * upperScale, geometricB + 4.0));
        if (barrierMax < barrierMin + 0.5)
            barrierMax = std::min(hardMax, barrierMin + 0.5);
    }

    auto scanGrid = [&](double tLow, double tHigh, int tSteps,
                        double bLow, double bHigh, int bSteps,
                        bool logTemperature)
    {
        tLow = std::clamp(tLow, tMin, tMax);
        tHigh = std::clamp(tHigh, tMin, tMax);
        if (tHigh < tLow)
            std::swap(tLow, tHigh);

        if (charged)
        {
            bLow = std::max(0.0, bLow);
            bHigh = std::max(bLow, bHigh);
        }

        const double logLow = std::log(std::max(tMin, tLow));
        const double logHigh = std::log(std::max(tMin, tHigh));

        for (int ti = 0; ti <= tSteps; ++ti)
        {
            const double tFraction = double(ti) / double(std::max(1, tSteps));
            const double temperature =
                logTemperature
                    ? std::exp(logLow + tFraction * (logHigh - logLow))
                    : tLow + tFraction * (tHigh - tLow);

            const int barrierCount = charged ? bSteps : 0;
            for (int bi = 0; bi <= barrierCount; ++bi)
            {
                const double bFraction = double(bi) / double(std::max(1, barrierCount));
                const double barrier =
                    charged ? bLow + bFraction * (bHigh - bLow) : 0.0;
                considerTemperatureAndBarrier(temperature, barrier);
            }
        }
    };

    scanGrid(tMin, tMax, 360,
             barrierMin, barrierMax, charged ? 90 : 0,
             true);

    if (bestNormalization > 0.0)
    {
        double tHalfWindow = std::max(0.12, bestParams.coldTemperatureMeV * 0.18);
        double bHalfWindow = charged ? std::max(0.15, (barrierMax - barrierMin) * 0.12) : 0.0;

        for (int pass = 0; pass < 5; ++pass)
        {
            const double tLow = std::max(tMin, bestParams.coldTemperatureMeV - tHalfWindow);
            const double tHigh = std::min(tMax, bestParams.coldTemperatureMeV + tHalfWindow);
            const double bLow =
                charged ? std::max(barrierMin, bestParams.effectiveBarrierMeV - bHalfWindow) : 0.0;
            const double bHigh =
                charged ? std::min(barrierMax, bestParams.effectiveBarrierMeV + bHalfWindow) : 0.0;

            scanGrid(tLow, tHigh, 120,
                     bLow, bHigh, charged ? 60 : 0,
                     false);

            tHalfWindow *= 0.35;
            bHalfWindow *= 0.35;
        }
    }

    result.valid = bestNormalization > 0.0 &&
                   std::isfinite(bestNormalization) &&
                   std::isfinite(bestDeviance);
    if (!result.valid)
    {
        result.unavailable = true;
        result.unavailableReason = "Boltzmann fit unavailable";
        return result;
    }

    result.normalization = bestNormalization;
    result.coldTemperatureMeV = bestParams.coldTemperatureMeV;
    result.effectiveBarrierMeV = charged ? bestParams.effectiveBarrierMeV : 0.0;
    result.barrierDiffusenessMeV = charged ? fixedDiffusenessMeV : 0.0;
    result.pearsonReducedChiSquare = bestPearson;
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
                   bool isParticlePlot,
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
        m_boltzmannKind(boltzmannKind),
        m_isParticlePlot(isParticlePlot)
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
        double energyBinWidthMeV = 1.0;
        double histogramBinWidth = 1.0;

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

            const double angularBinWidthDeg = (maxX - minX) / double(bins);
            histogramBinWidth = angularBinWidthDeg;

            std::vector<double> counts(bins, 0.0);

            for (int i = 0; i < n; ++i)
            {
                const double theta = m_entry.thetaDeg[i];

                if (theta < minX || theta > maxX)
                    continue;

                int bin = int((theta - minX) / (maxX - minX) * bins);

                if (bin >= bins) bin = bins - 1;
                if (bin >= 0) counts[bin] += sampleWeightAt(m_entry, i);
            }

            values.resize(bins, 0.0);
            xCenters.resize(bins, 0.0);

            for (int i = 0; i < bins; ++i)
                xCenters[i] = minX + (double(i) + 0.5) * angularBinWidthDeg;

            for (int i = 0; i < bins; ++i)
            {
                if (m_plotKind == PlotCountsVsTheta)
                {
                    values[i] = double(counts[i]);
                }
                else
                {
                    values[i] = m_sigmaTotal * double(counts[i]) /
                                (double(m_nEvents) * angularBinWidthDeg);
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
            energyBinWidthMeV = 1.0;

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

                energyBinWidthMeV = (maxX <= 80.0) ? 1.0 : 2.0;
                bins = int(std::ceil(maxX / energyBinWidthMeV));
                bins = std::clamp(bins, 20, 80);
                energyBinWidthMeV = maxX / double(bins);
                histogramBinWidth = energyBinWidthMeV;

                values.assign(bins, 0.0);
                xCenters.assign(bins, 0.0);

                for (int i = 0; i < bins; ++i)
                    xCenters[i] = minX + (double(i) + 0.5) * energyBinWidthMeV;
            }
            else
            {
                maxX = m_acc.ELOW + 29.0 * m_acc.DELE;
                bins = 31;
                energyBinWidthMeV = m_acc.DELE;
                histogramBinWidth = energyBinWidthMeV;
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
                    bin = int((energy - minX) / energyBinWidthMeV);
                else if (energy >= m_acc.ELOW)
                    bin = std::min(30, int((energy - m_acc.ELOW) / m_acc.DELE) + 1);

                if (bin >= bins)
                    bin = bins - 1;
                if (bin >= 0 && bin < bins)
                    values[bin] += sampleWeightAt(m_entry, i);
            }

            title = "N = f(E)";
            xTitle = useCmEnergy ? "C.M. / source energy (MeV)" : "Lab energy (MeV)";
            yTitle = "Counts N";
        }

        double maxY = 0.0;

        for (double v : values)
            maxY = std::max(maxY, v);

        const bool showGaussian =
            m_allowGaussianEnergyOverlay &&
            (m_plotKind == PlotCountsVsTheta || m_plotKind == PlotCountsVsEnergy);

        GaussianOverlayStats gaussianStats;

        if (showGaussian)
            gaussianStats = computeGaussianOverlayStats(xCenters,
                                                        values,
                                                        m_plotKind == PlotCountsVsEnergy);

        if (gaussianStats.valid)
            maxY = std::max(maxY, gaussianStats.amplitude * 1.15);

        const HistogramMomentStats histogramStats =
            computeHistogramMomentStats(xCenters, values, histogramBinWidth);

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
                fitBoltzmannToHistogram(xCenters, values, m_boltzmannKind,
                                        energyBinWidthMeV, barrier);

            boltzmannCurve.drawBarrier = charged;
            boltzmannCurve.effectiveBarrierMeV = fit.valid ? fit.effectiveBarrierMeV : barrier.valueMeV;
            boltzmannCurve.unavailable = fit.unavailable;
            boltzmannCurve.unavailableReason = fit.unavailableReason;

            if (fit.valid)
            {
                boltzmannCurve.valid = minX < maxX;
                boltzmannCurve.normalization = fit.normalization;
                boltzmannCurve.coldTemperatureMeV = fit.coldTemperatureMeV;
                boltzmannCurve.effectiveBarrierMeV = fit.effectiveBarrierMeV;
                boltzmannCurve.barrierDiffusenessMeV = fit.barrierDiffusenessMeV;
                boltzmannCurve.pearsonReducedChiSquare = fit.pearsonReducedChiSquare;

                if (boltzmannCurve.valid)
                {
                    const int samples = 700;
                    EvaporationModelParams params;
                    params.coldTemperatureMeV = fit.coldTemperatureMeV;
                    params.effectiveBarrierMeV = fit.effectiveBarrierMeV;
                    params.barrierDiffusenessMeV = fit.barrierDiffusenessMeV;

                    for (int i = 0; i <= samples; ++i)
                    {
                        const double frac = double(i) / double(samples);
                        const double xValue = minX + frac * (maxX - minX);
                        const double yValue =
                            fit.normalization *
                            evaporationDensity(xValue, m_boltzmannKind, params) *
                            energyBinWidthMeV;
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
            const bool includeHistogramStats =
                histogramStats.valid &&
                !includeGaussian &&
                (m_plotKind == PlotCountsVsEnergy ||
                 m_plotKind == PlotCrossSectionVsTheta ||
                 (m_isParticlePlot && m_plotKind == PlotCountsVsTheta));

            p.save();

            const int legendW = (includeBoltzmann || includeBoltzmannUnavailable) ? 315 : 210;
            int legendH = 36;
            if (includeHistogramStats)
                legendH += 45;
            if (includeGaussian)
                legendH += 105;
            if (includeBoltzmann)
                legendH += boltzmannCurve.drawBarrier ? 90 : 75;
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

            int cursorY = legendY + 33;

            if (includeHistogramStats)
            {
                p.drawText(legendX + 10,
                           cursorY,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Mean = " + formatGaussianNumber(histogramStats.mean));

                p.drawText(legendX + 10,
                           cursorY + 15,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Median = " + formatGaussianNumber(histogramStats.median));

                p.drawText(legendX + 10,
                           cursorY + 30,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Area = " + formatGaussianNumber(histogramStats.area));

                cursorY += 45;
            }

            if (!includeGaussian && !includeBoltzmann && !includeBoltzmannUnavailable)
            {
                p.restore();
                return;
            }

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
                           "Median = " + formatGaussianNumber(histogramStats.median));

                p.drawText(legendX + 10,
                           cursorY + 60,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Chi^2(red) = " + formatGaussianNumber(gaussianStats.reducedChiSquare));

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
                           "A = " + formatGaussianNumber(boltzmannCurve.normalization));

                p.drawText(legendX + 10,
                           cursorY + 30,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "T = " +
                               QString::number(boltzmannCurve.coldTemperatureMeV, 'g', 3) + " MeV");

                int infoY = cursorY + 45;
                if (boltzmannCurve.drawBarrier)
                {
                    p.drawText(legendX + 10,
                               infoY,
                               legendW - 20,
                               14,
                               Qt::AlignLeft,
                               "B = " +
                                   QString::number(boltzmannCurve.effectiveBarrierMeV, 'g', 3) + " MeV");
                    infoY += 15;
                }

                p.drawText(legendX + 10,
                           infoY,
                           legendW - 20,
                           14,
                           Qt::AlignLeft,
                           "Chi^2(red) = " +
                               formatGaussianNumber(boltzmannCurve.pearsonReducedChiSquare));
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
            boltzmannCurve.effectiveBarrierMeV >= minX &&
            boltzmannCurve.effectiveBarrierMeV <= maxX)
        {
            const int barrierX = mapX(boltzmannCurve.effectiveBarrierMeV);
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
                           QString::number(boltzmannCurve.effectiveBarrierMeV, 'g', 3) + " MeV");
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
    bool m_isParticlePlot = false;
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
    addAngularSample(entry, kineticEnergy, thetaDeg, vz, vxy, cmEnergy, 1.0f);
}

void addAngularSample(AngularDistEntry &entry,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy,
                      float cmEnergy,
                      float weight)
{
    entry.kineticEnergy.push_back(kineticEnergy);
    entry.thetaDeg.push_back(thetaDeg);
    entry.vz.push_back(vz);
    entry.vxy.push_back(vxy);
    entry.cmEnergy.push_back(cmEnergy);
    entry.weight.push_back(weight);
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
    addAngularSample(entries, z, n, kineticEnergy, thetaDeg, vz, vxy, cmEnergy, 1.0f);
}

void addAngularSample(std::map<std::pair<int, int>, AngularDistEntry> &entries,
                      int z,
                      int n,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy,
                      float cmEnergy,
                      float weight)
{
    const std::pair<int, int> key(z, n);
    AngularDistEntry &e = entries[key];
    e.z = z;
    e.n = n;
    addAngularSample(e, kineticEnergy, thetaDeg, vz, vxy, cmEnergy, weight);
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
    double yieldSigmaTotalMb,
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
    const double imfSampleCount = isImf ? totalEntryWeight(entries) : 0.0;
    const double paceSigmaTotalMb =
        (isImf && yieldSigmaTotalMb >= 0.0 && imfSampleCount > 0)
            ? yieldSigmaTotalMb * double(std::max(1, nCascades)) / double(imfSampleCount)
            : sigmaTotalMb;

    html += "<!DOCTYPE html><html><head><meta charset=\"utf-8\">";
    html += "<style>"
            "body{font-family:Sans-Serif;}"
            "table{border-collapse:collapse; margin:10px 0;}"
            "th,td{border:1px solid #909090; padding:4px 6px; white-space:nowrap;}"
            "th{background:#ececec;}"
            ".total-row th,.total-row td{background:#ececec; font-weight:bold;}"
            "h2{color:#1d4f91; text-align:center;}"
            "h3{margin-top:25px;}"
            ".angular-content{display:inline-block; width:fit-content; text-align:left;}"
            ".cm-spectra{margin:24px auto; min-width:520px;}"
            ".cm-spectra th{background:#dfeaf7; color:#1d4f91;}"
            ".cm-spectra td,.cm-spectra th{padding:5px 10px; text-align:center;}"
            "</style></head><body>";
    html += "<div class=\"angular-content\">";

    if ((!hasAngularSamples(entries)
         && !hasAngularSamples(neutronEntry)
         && !hasAngularSamples(protonEntry)
         && !hasAngularSamples(alphaEntry)
         && !hasAngularSamples(gammaEntry))
        || nCascades <= 0)
    {
        html += "<p>No angular-distribution events were found for output.</p></div></body></html>";
        return html;
    }

    const std::vector<PaceSelectedResidue> selected =
        buildSelectedResiduesPACE(entries, nCascades, lowLimitPercent, highLimitPercent);

    const QString tabHeader =
        isImf ? "IMF Angular Distributions" : "Residual Angular Distributions";

    html += "<h2 align=center style=\"color:#1f5b9e;\">"
            + tabHeader +
            "</h2>";

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
            const double weight = sampleWeightAt(e, i);
            addToPaceAccum(allAcc, e.kineticEnergy[i], e.thetaDeg[i], weight);
            if (selectedIndex >= 0)
                addToPaceAccum(selectedAcc[selectedIndex].second, e.kineticEnergy[i], e.thetaDeg[i], weight);
        }
    }

    if (totalAngularCount(allAcc) <= 0
        && !hasAngularSamples(neutronEntry)
        && !hasAngularSamples(protonEntry)
        && !hasAngularSamples(alphaEntry)
        && !hasAngularSamples(gammaEntry))
    {
        html += "<p>No angular-distribution events were found for output.</p></div></body></html>";
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
        appendMainPaceTable(html, acc, r.a, paceSigmaTotalMb, nCascades, true, false, plotIndex,
                            velocityLabel);
        displayIndex++;
        plotIndex++;
    }

    if (totalAngularCount(allAcc) > 0)
    {
        html += buildHeaderAllPACE(displayIndex, fragmentKindPlural);
        appendMainPaceTable(html, allAcc, compoundA, paceSigmaTotalMb, nCascades, true, true, plotIndex,
                            velocityLabel);
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
        appendMainPaceTable(html, acc, particleA, paceSigmaTotalMb, nCascades,
                            showVelocity, false, plotIndex, velocityLabel);

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

    html += "</div></body></html>";
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

    QTextBrowser *text = new QTextBrowser;
    text->setOpenLinks(false);
    text->setLineWrapMode(QTextEdit::NoWrap);
    text->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
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
                                      tables[i].boltzmannKind,
                                      tables[i].isParticlePlot);

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
