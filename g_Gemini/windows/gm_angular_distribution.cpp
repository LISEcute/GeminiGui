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

constexpr int kMaxResiduesToPrint = 15;   // PACE: 15 residues + ALL
constexpr double kDegToRad = 0.01745;     // PACE-style constant
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

    // 1..32 energy bins, 1..37 angle bins
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
    line += "<p><em>" + velocityLabel.toHtmlEscaped() + "</em> Vz = ";
    line += QString::number(vzm, 'e', 2);
    line += " (sig = ";
    line += QString::number(dvz, 'e', 2);
    line += ") rms Vxy = ";
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

    html += "<tr><td>dSig/dOmeg</td>";
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
    html += "<a href=\"gemini://plot_table/" + QString::number(plotIndex) + "\">Plot E vs Theta</a> &nbsp; ";
    html += "<a href=\"gemini://plot_ntheta/" + QString::number(plotIndex) + "\">Plot N vs Theta</a> &nbsp; ";
    html += "<a href=\"gemini://plot_ne/" + QString::number(plotIndex) + "\">Plot N vs E</a> &nbsp; ";
    html += "<a href=\"gemini://plot_dcsdtheta/" + QString::number(plotIndex) + "\">Plot dCS/dTheta vs Theta</a>";
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

    html += "<tr><td>dSig/dOmeg</td>";
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
    html += "<a href=\"gemini://plot_table/" + QString::number(plotIndex) + "\">Plot E vs Theta</a> &nbsp; ";
    html += "<a href=\"gemini://plot_ntheta/" + QString::number(plotIndex) + "\">Plot N vs Theta</a> &nbsp; ";
    html += "<a href=\"gemini://plot_ne/" + QString::number(plotIndex) + "\">Plot N vs E</a> &nbsp; ";
    html += "<a href=\"gemini://plot_dcsdtheta/" + QString::number(plotIndex) + "\">Plot dCS/dTheta vs Theta</a>";
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
    tables.reserve(selected.size() + 4);

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
        appendEntryToAccum(table.acc, entry);
        tables.push_back(table);
        idx++;
    };

    appendParticlePlot("neutrons", neutronEntry);
    appendParticlePlot("protons", protonEntry);
    appendParticlePlot("alpha particles", alphaEntry);

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

        const double minE = 0.0;
        const double maxE = m_acc.ELOW + 29.0 * m_acc.DELE;

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
        }

        if (maxCount <= 0)
        {
            p.drawText(plotRect, Qt::AlignCenter, "No histogram data");
            return;
        }

        const double cellW = double(plotRect.width()) / thetaBins;
        const double cellH = double(plotRect.height()) / energyBins;

        // Draw 2D histogram cells
        for (int e = 0; e < energyBins; ++e)
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

        for (int e = 0; e <= energyBins; ++e)
        {
            const double frac = double(e) / double(energyBins);
            double energyValue = 0.0;
            if (e == 0) energyValue = 0.0;
            else if (e == 1) energyValue = m_acc.ELOW;
            else energyValue = m_acc.ELOW + double(e - 1) * m_acc.DELE;

            const int y = plotRect.bottom() - int(frac * plotRect.height());

            if (e == 0 || e == 1 || e == energyBins || (e - 1) % 5 == 0)
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
                   "Theta range (deg)");

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
    double amplitude = 0.0;
    double mean = 0.0;
    double stddev = 0.0;
    bool valid = false;
};

GaussianOverlayStats computeGaussianOverlayStats(const std::vector<double> &xCenters,
                                                 const std::vector<double> &yValues)
{
    GaussianOverlayStats stats;
    if (xCenters.size() != yValues.size() || xCenters.empty()) return stats;

    double sumW = 0.0;
    double sumX = 0.0;
    stats.amplitude = 0.0;

    for (int i = 0; i < int(yValues.size()); ++i)
    {
        const double w = std::max(0.0, yValues[i]);
        if (w <= 0.0) continue;

        sumW += w;
        sumX += w * xCenters[i];
        stats.amplitude = std::max(stats.amplitude, w);
    }

    if (sumW <= 0.0) return stats;

    stats.mean = sumX / sumW;

    double variance = 0.0;
    for (int i = 0; i < int(yValues.size()); ++i)
    {
        const double w = std::max(0.0, yValues[i]);
        if (w <= 0.0) continue;

        const double dx = xCenters[i] - stats.mean;
        variance += w * dx * dx;
    }

    stats.stddev = std::sqrt(variance / sumW);
    stats.valid = (stats.amplitude > 0.0 && stats.stddev > 1.0e-12);
    return stats;
}

QString formatGaussianNumber(double value)
{
    if (std::fabs(value) >= 1000.0 || (std::fabs(value) > 0.0 && std::fabs(value) < 0.01))
        return QString::number(value, 'e', 2);

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
        const int right = 45;
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
            bins = std::max(1, int(maxX));       // 1-degree histogram bins
            const double binWidth = (maxX - minX) / double(bins);

            std::vector<int> counts(bins, 0);
            for (int i = 0; i < n; ++i)
            {
                const double theta = m_entry.thetaDeg[i];
                if (theta < minX || theta > maxX) continue;

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
                    // Histogram version of d cross section / d theta.
                    // Units are mb/degree because the bin width is in degrees.
                    values[i] = m_sigmaTotal * double(counts[i]) /
                                (double(m_nEvents) * binWidth);
                }
            }

            xTitle = "Theta (deg)";
            if (m_plotKind == PlotCountsVsTheta)
            {
                title = "N = f(theta)";
                yTitle = "Counts N";
            }
            else
            {
                title = "d\u03c3/d\u03b8 = f(theta)";
                yTitle = "d\u03c3/d\u03b8 (mb/deg)";
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

                if (bin >= 0 && bin < bins) values[bin] += 1.0;
            }

            title = "N = f(E)";
            xTitle = "Energy (MeV)";
            yTitle = "Counts N";
        }

        if (values.empty())
        {
            p.drawText(plotRect, Qt::AlignCenter, "No histogram data");
            return;
        }

        double maxY = *std::max_element(values.begin(), values.end());
        if (maxY <= 0.0) maxY = 1.0;
        const bool countPlot = (m_plotKind == PlotCountsVsTheta || m_plotKind == PlotCountsVsEnergy);
        if (countPlot)
            maxY = std::ceil(maxY * 1.10);
        else
            maxY *= 1.10;

        auto mapY = [&](double value)
        {
            return plotRect.bottom() - int(value / maxY * plotRect.height());
        };

        // grid + y labels
        p.setPen(QPen(QColor(220, 220, 220), 1));
        for (int i = 0; i <= 5; ++i)
        {
            const double frac = double(i) / 5.0;
            const int yPix = plotRect.bottom() - int(frac * plotRect.height());
            p.drawLine(plotRect.left(), yPix, plotRect.right(), yPix);
            p.setPen(Qt::black);
            const QString yLabel = countPlot
                                       ? QString::number(int(std::round(frac * maxY)))
                                       : QString::number(frac * maxY, 'g', 3);
            p.drawText(5, yPix - 10, left - 15, 20, Qt::AlignRight | Qt::AlignVCenter, yLabel);
            p.setPen(QPen(QColor(220, 220, 220), 1));
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

        const bool showGaussianOverlay = (m_plotKind == PlotCountsVsTheta ||
                                          m_plotKind == PlotCountsVsEnergy);
        const GaussianOverlayStats gaussianStats = showGaussianOverlay
                                                       ? computeGaussianOverlayStats(xCenters, values)
                                                       : GaussianOverlayStats();
        if (showGaussianOverlay && gaussianStats.valid)
        {
            maxY = std::max(maxY, gaussianStats.amplitude * 1.10);

            auto mapYForGaussian = [&](double value)
            {
                return plotRect.bottom() - int(value / maxY * plotRect.height());
            };

            QPolygonF gaussianCurve;
            gaussianCurve.reserve(std::max(80, bins * 8));

            const int curvePoints = std::max(120, bins * 8);
            const double xMinCurve = xCenters.front();
            const double xMaxCurve = xCenters.back();
            const double xSpanCurve = std::max(1.0e-12, xMaxCurve - xMinCurve);

            for (int i = 0; i < curvePoints; ++i)
            {
                const double frac = double(i) / double(curvePoints - 1);
                const double xValue = xMinCurve + frac * xSpanCurve;
                const double z = (xValue - gaussianStats.mean) / gaussianStats.stddev;
                const double yValue = gaussianStats.amplitude * std::exp(-0.5 * z * z);

                const int xPix = plotRect.left() + int(frac * plotRect.width());
                const int yPix = mapYForGaussian(yValue);
                gaussianCurve << QPointF(xPix, yPix);
            }

            QPen gaussianPen(QColor(220, 40, 40, 145), 2, Qt::DotLine);
            gaussianPen.setCapStyle(Qt::RoundCap);
            p.setPen(gaussianPen);
            p.drawPolyline(gaussianCurve);

            const int legendW = 150;
            const int legendH = 66;
            const int legendX = plotRect.right() - legendW - 10;
            const int legendY = plotRect.top() + 10;
            QRect legendRect(legendX, legendY, legendW, legendH);

            p.fillRect(legendRect, QColor(255, 255, 255, 210));
            p.setPen(QPen(QColor(120, 120, 120), 1));
            p.drawRect(legendRect);

            p.setPen(gaussianPen);
            p.drawLine(legendX + 8, legendY + 14, legendX + 32, legendY + 14);

            QFont legendFont = p.font();
            legendFont.setPointSize(8);
            legendFont.setBold(false);
            p.setFont(legendFont);
            p.setPen(Qt::black);
            p.drawText(legendX + 38, legendY + 5, legendW - 44, 16, Qt::AlignLeft, "Gaussian");
            p.drawText(legendX + 8, legendY + 22, legendW - 16, 14, Qt::AlignLeft,
                       "Amp = " + formatGaussianNumber(gaussianStats.amplitude));
            p.drawText(legendX + 8, legendY + 36, legendW - 16, 14, Qt::AlignLeft,
                       "Mean = " + formatGaussianNumber(gaussianStats.mean));
            p.drawText(legendX + 8, legendY + 50, legendW - 16, 14, Qt::AlignLeft,
                       "Std = " + formatGaussianNumber(gaussianStats.stddev));
        }

        p.setPen(Qt::black);
        if (m_plotKind == PlotCountsVsEnergy)
        {
            for (int i = 0; i <= bins; ++i)
            {
                if (i != 0 && i != 1 && i != bins && (i - 1) % 5 != 0) continue;

                const double frac = double(i) / double(bins);
                double value = 0.0;
                if (i == 0) value = 0.0;
                else if (i == 1) value = m_acc.ELOW;
                else value = m_acc.ELOW + double(i - 1) * m_acc.DELE;

                const int xPix = plotRect.left() + int(frac * plotRect.width());
                p.drawLine(xPix, plotRect.bottom(), xPix, plotRect.bottom() + 6);
                p.drawText(xPix - 30, plotRect.bottom() + 24, 60, 18, Qt::AlignCenter,
                           QString::number(value, 'g', 3));
            }
        }
        else
        {
            for (int i = 0; i <= 6; ++i)
            {
                const double frac = double(i) / 6.0;
                const double value = minX + frac * (maxX - minX);
                const int xPix = plotRect.left() + int(frac * plotRect.width());
                p.drawLine(xPix, plotRect.bottom(), xPix, plotRect.bottom() + 6);
                p.drawText(xPix - 30, plotRect.bottom() + 24, 60, 18, Qt::AlignCenter,
                           QString::number(value, 'g', 3));
            }
        }

        // axes
        p.setPen(QPen(Qt::black, 2));
        p.drawLine(plotRect.bottomLeft(), plotRect.bottomRight());
        p.drawLine(plotRect.bottomLeft(), plotRect.topLeft());

        // titles
        QFont titleFont = p.font();
        titleFont.setPointSize(12);
        titleFont.setBold(true);
        p.setFont(titleFont);
        p.setPen(Qt::black);
        p.drawText(0, 8, width(), 24, Qt::AlignCenter, title);

        QFont axisFont = p.font();
        axisFont.setPointSize(10);
        axisFont.setBold(false);
        p.setFont(axisFont);
        p.drawText(plotRect.left(), height() - 35, plotRect.width(), 24, Qt::AlignCenter, xTitle);

        p.save();
        p.translate(30, plotRect.top() + plotRect.height() / 2);
        p.rotate(-90);
        p.drawText(QRect(-plotRect.height() / 2, -20, plotRect.height(), 20),
                   Qt::AlignCenter, yTitle);
        p.restore();
    }

private:
    AngularDistEntry m_entry;
    PaceAngularAccum m_acc;
    int m_plotKind = PlotCountsVsTheta;
    double m_sigmaTotal = 0.0;
    int m_nEvents = 1;
};
}

void addAngularSample(AngularDistEntry &entry,
                      float kineticEnergy,
                      float thetaDeg,
                      float vz,
                      float vxy)
{
    entry.kineticEnergy.push_back(kineticEnergy);
    entry.thetaDeg.push_back(thetaDeg);
    entry.vz.push_back(vz);
    entry.vxy.push_back(vxy);
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
    const AngularDistEntry &alphaEntry)
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
            "</style></head><body>";

    html += "<p class=\"small\">" + title + "</p>";

    if ((!hasAngularSamples(entries)
         && !hasAngularSamples(neutronEntry)
         && !hasAngularSamples(protonEntry)
         && !hasAngularSamples(alphaEntry))
        || nCascades <= 0)
    {
        html += "<p>No angular-distribution events were found for output.</p></body></html>";
        return html;
    }

    const std::vector<PaceSelectedResidue> selected =
        buildSelectedResiduesPACE(entries, nCascades, lowLimitPercent, highLimitPercent);

    html += "<p>&nbsp;</p><h2 align=\"center\">Angular distribution results</h2>";
    html += "<p>";
    html += "<a href=\"gemini://plot_all\">Plot All: E vs Theta</a> &nbsp; ";
    html += "<a href=\"gemini://plot_all_ntheta\">Plot All: N vs Theta</a> &nbsp; ";
    html += "<a href=\"gemini://plot_all_ne\">Plot All: N vs E</a> &nbsp; ";
    html += "<a href=\"gemini://plot_all_dcsdtheta\">Plot All: dCS/dTheta vs Theta</a>";
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
                                     int particleA)
    {
        if (!hasAngularSamples(entry)) return;

        PaceAngularAccum acc;
        initPaceAccum(acc, compoundA, compoundExcitationMeV, recoilBetaCN, mdir);
        appendEntryToAccum(acc, entry);

        html += buildHeaderForParticlePACE(displayIndex, particleLabel);
        appendMainPaceTable(html, acc, particleA, sigmaTotalMb, nCascades, inputMode, false, plotIndex,
                            velocityLabel);
        appendSecondPaceTableIfNeeded(html, acc, plotIndex);
        displayIndex++;
        plotIndex++;
    };

    appendParticleSection("neutrons", "Neutron velocity/c", neutronEntry, 1);
    appendParticleSection("protons", "Proton velocity/c", protonEntry, 1);
    appendParticleSection("alpha particles", "Alpha velocity/c", alphaEntry, 4);

    html += "</body></html>";
    return html;
}

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
         && !hasAngularSamples(alphaEntryForPlots))
        || nEventsForPlots <= 0) return;

    const std::vector<PlotTableEntry> tables =
        buildPlotTableListPACE(entriesForPlots,
                               neutronEntryForPlots,
                               protonEntryForPlots,
                               alphaEntryForPlots,
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
    QString plotKindLabel = "E vs Theta";
    if (plotKind == PlotCountsVsTheta) plotKindLabel = "N vs Theta";
    else if (plotKind == PlotCountsVsEnergy) plotKindLabel = "N vs E";
    else if (plotKind == PlotCrossSectionVsTheta) plotKindLabel = "dCS/dTheta vs Theta";

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
