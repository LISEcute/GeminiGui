#include "gm_mainwindow.h"
#include "ui_gm_mainwindow.h"

#include <QElapsedTimer>
#include <QEventLoop>
#include <QFuture>
#include <QFutureWatcher>
#include <QProgressDialog>
#include <QThread>
#include <QtConcurrent>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <map>
#include <vector>

#include "g_Gemini/source/CNucleus.h"
#include "g_Gemini/source/CFus.h"
#include "gm_results.h"
#include "gm_about.h"
#include "gm_angular_distribution.h"
// this is an example of using GEMINI CNucleus class to give the
//statistical decay of a compound nucleus
extern bool _useAME;
extern bool _useIMF;
extern bool _useIMFenh;
extern int  _optEvap;
double _LowLimit  = 0.0;
double _HighLimit = 100.0;
extern FILE *mfopen(const QString& filename, const char* operand);
bool Residual_compare(Residual a, Residual b) {return a.index() > b.index();};

namespace {

struct YieldRow
{
    int z = 0;
    int a = 0;
    QString name;
    double events = 0.0;
    double percent = 0.0;
    double sigma = 0.0;
    double err = 0.0;
};

using ImfYieldCounts = std::map<std::pair<int, int>, double>;

QString formatYieldCount(double value)
{
    if (std::abs(value - std::round(value)) < 1.0e-6)
        return QString::number(value, 'f', 0);
    return QString::number(value, 'g', 4);
}

QString cellStyleForYield(bool zSeparator, bool imfSeparator, bool bold = false)
{
    QStringList styles;
    if (bold) styles << "font-weight:bold";
    if (zSeparator) styles << "border-top:1px dotted #b0b0b0";
    if (imfSeparator) styles << "border-left:1px dotted #b0b0b0";
    return styles.isEmpty()
               ? QString()
               : " style=\"" + styles.join("; ") + ";\"";
}

QString yieldCells(const YieldRow *row, bool zSeparator = false, bool imfSeparator = false)
{
    const QString firstCellStyle = cellStyleForYield(zSeparator, imfSeparator, row != nullptr);
    const QString cellStyle = cellStyleForYield(zSeparator, false);

    if (!row)
        return "<td" + firstCellStyle + ">-</td>"
               "<td" + cellStyle + ">-</td>"
               "<td" + cellStyle + ">-</td>"
               "<td" + cellStyle + ">-</td>"
               "<td" + cellStyle + ">-</td>";

    return "<td" + firstCellStyle + ">" + row->name + "</td>"
           "<td" + cellStyle + ">" + formatYieldCount(row->events) + "</td>"
           "<td" + cellStyle + ">" + QString::number(row->percent, 'f', 1) + "%</td>"
           "<td" + cellStyle + ">" + QString::number(row->sigma, 'g', 4) + "</td>"
           "<td" + cellStyle + ">" + QString::number(row->err, 'g', 4) + "</td>";
}

double angularEntryWeightTotal(const AngularDistEntry &entry)
{
    double total = 0.0;
    for (int i = 0; i < int(entry.kineticEnergy.size()); ++i)
    {
        const double weight =
            (i < int(entry.weight.size()))
                ? std::max(0.0, double(entry.weight[i]))
                : 1.0;
        total += weight;
    }
    return total;
}

double angularEntryWeightTotalForMap(const std::map<std::pair<int, int>, AngularDistEntry> &entries)
{
    double total = 0.0;
    for (const auto &it : entries)
        total += angularEntryWeightTotal(it.second);
    return total;
}

double imfYieldCountTotal(const ImfYieldCounts &imfYields)
{
    double total = 0.0;
    for (const auto &it : imfYields)
        total += it.second;
    return total;
}

std::map<int, double> residualYieldCountsByZ(Residual *resid, int length)
{
    std::map<int, double> countsByZ;
    for (int i = 0; i < length; i++)
    {
        if (resid[i].count != 0)
            countsByZ[resid[i].Z] += resid[i].count;
    }
    return countsByZ;
}

YieldPlotData buildYieldPlotDataFromCounts(const std::map<int, double> &residualByZ,
                                           const std::map<int, double> &imfByZ)
{
    std::vector<int> zValues;
    for (const auto &it : residualByZ)
        zValues.push_back(it.first);
    for (const auto &it : imfByZ)
    {
        if (residualByZ.find(it.first) == residualByZ.end())
            zValues.push_back(it.first);
    }
    std::sort(zValues.begin(), zValues.end());

    YieldPlotData plotData;
    for (int z : zValues)
    {
        YieldPlotPoint point;
        point.z = z;

        const auto resIt = residualByZ.find(z);
        if (resIt != residualByZ.end())
            point.residualEvents = resIt->second;

        const auto imfIt = imfByZ.find(z);
        if (imfIt != imfByZ.end())
            point.imfEvents = imfIt->second;

        plotData.points.push_back(point);
    }

    return plotData;
}

YieldPlotData buildYieldPlotDataFromImfCounts(Residual *resid,
                                              int length,
                                              const ImfYieldCounts &imfYields)
{
    std::map<int, double> imfByZ;
    for (const auto &it : imfYields)
        imfByZ[it.first.first] += it.second;

    return buildYieldPlotDataFromCounts(residualYieldCountsByZ(resid, length), imfByZ);
}

std::map<int, std::vector<YieldRow>> buildImfYieldRowsByZ(
    const std::map<std::pair<int, int>, AngularDistEntry> &imfEntries,
    double imfSigmaTotalMb)
{
    const double imfTotal = angularEntryWeightTotalForMap(imfEntries);

    std::map<int, std::vector<YieldRow>> rowsByZ;
    if (imfTotal <= 0.0) return rowsByZ;

    for (const auto &it : imfEntries)
    {
        const AngularDistEntry &entry = it.second;
        const double events = angularEntryWeightTotal(entry);
        if (events <= 0.0) continue;

        CNucleus nuc(entry.z, entry.z + entry.n);

        YieldRow row;
        row.z = entry.z;
        row.a = entry.z + entry.n;
        row.name = nuc.getGName();
        row.events = events;
        row.percent = 100.0 * events / imfTotal;
        row.sigma = imfSigmaTotalMb * events / imfTotal;
        row.err = row.sigma / std::sqrt(events);
        rowsByZ[row.z].push_back(row);
    }

    for (auto &it : rowsByZ)
    {
        std::sort(it.second.begin(), it.second.end(),
                  [](const YieldRow &lhs, const YieldRow &rhs)
                  {
                      if (lhs.a != rhs.a) return lhs.a > rhs.a;
                      return lhs.z > rhs.z;
                  });
    }

    return rowsByZ;
}

std::map<int, std::vector<YieldRow>> buildImfYieldRowsByZ(
    const ImfYieldCounts &imfYields,
    double imfSigmaTotalMb)
{
    const double imfTotal = imfYieldCountTotal(imfYields);

    std::map<int, std::vector<YieldRow>> rowsByZ;
    if (imfTotal <= 0.0) return rowsByZ;

    for (const auto &it : imfYields)
    {
        const int z = it.first.first;
        const int n = it.first.second;
        const double events = it.second;
        if (events <= 0.0) continue;

        CNucleus nuc(z, z + n);

        YieldRow row;
        row.z = z;
        row.a = z + n;
        row.name = nuc.getGName();
        row.events = events;
        row.percent = 100.0 * events / imfTotal;
        row.sigma = imfSigmaTotalMb * events / imfTotal;
        row.err = row.sigma / std::sqrt(events);
        rowsByZ[row.z].push_back(row);
    }

    for (auto &it : rowsByZ)
    {
        std::sort(it.second.begin(), it.second.end(),
                  [](const YieldRow &lhs, const YieldRow &rhs)
                  {
                      if (lhs.a != rhs.a) return lhs.a > rhs.a;
                      return lhs.z > rhs.z;
                  });
    }

    return rowsByZ;
}

QString buildMergedYieldTableHtmlFromImfRows(const QString &title,
                                             Residual *resid,
                                             int length,
                                             int countResidue,
                                             double residualSigmaTotalMb,
                                             const std::map<int, std::vector<YieldRow>> &imfRowsByZ,
                                             double imfTotal,
                                             double imfSigmaTotalMb,
                                             FILE *file_cs)
{
    std::map<int, std::vector<YieldRow>> residueRowsByZ;

    std::sort(resid, resid+length, Residual_compare);
    for(int i=0;i<length;i++)
    {
        if(resid[i].count != 0)
        {
            YieldRow row;
            row.z = resid[i].Z;
            row.a = resid[i].A;
            CNucleus nuc(row.z, row.a);
            row.name = nuc.getGName();
            row.events = resid[i].count;
            row.percent = countResidue > 0 ? 100.0 * double(resid[i].count) / double(countResidue) : 0.0;
            row.sigma = countResidue > 0 ? residualSigmaTotalMb * double(resid[i].count) / double(countResidue) : 0.0;
            row.err = row.sigma / std::sqrt(double(resid[i].count));
            residueRowsByZ[row.z].push_back(row);

            if(file_cs) fprintf(file_cs,"\n%d %d %10.3g %10.3g",resid[i].Z,resid[i].A-resid[i].Z,row.sigma,row.err);
        }
    }

    std::vector<int> zValues;
    for (const auto &it : residueRowsByZ)
        zValues.push_back(it.first);
    for (const auto &it : imfRowsByZ)
    {
        if (residueRowsByZ.find(it.first) == residueRowsByZ.end())
            zValues.push_back(it.first);
    }
    std::sort(zValues.begin(), zValues.end(), std::greater<int>());

    QString results;
    results += "<h3 align=\"center\" style=\"color: blue\"> " + title +
               "(<a href=\"gemini://yield_plot\">plot</a>) </h3><br>";
    results += "<table cellpadding=\"5\" align=\"center\">";
    results += "<tr style=\"color: green\">"
               "<th></th>"
               "<th colspan=\"5\">Residual Nuclei</th>"
               "<th colspan=\"5\" style=\"border-left:1px dotted #b0b0b0;\">IMF Particles</th>"
               "</tr>";
    results += "<tr style=\"color: green\">"
               "<th>Z</th>"
               "<th>Name</th><th>Events</th><th>Percent</th><th>x-section (mb)</th><th>err(mb)</th>"
               "<th style=\"border-left:1px dotted #b0b0b0;\">Name</th><th>Events</th><th>Percent</th><th>x-section (mb)</th><th>err(mb)</th>"
               "</tr>";

    bool firstZGroup = true;
    for (int z : zValues)
    {
        const auto resIt = residueRowsByZ.find(z);
        const auto imfIt = imfRowsByZ.find(z);
        const std::vector<YieldRow> empty;
        const std::vector<YieldRow> &resRows = resIt != residueRowsByZ.end() ? resIt->second : empty;
        const std::vector<YieldRow> &imfRows = imfIt != imfRowsByZ.end() ? imfIt->second : empty;
        const int rowCount = int(std::max(resRows.size(), imfRows.size()));

        for (int rowIndex = 0; rowIndex < rowCount; ++rowIndex)
        {
            const YieldRow *resRow = rowIndex < int(resRows.size()) ? &resRows[rowIndex] : nullptr;
            const YieldRow *imfRow = rowIndex < int(imfRows.size()) ? &imfRows[rowIndex] : nullptr;
            const bool zSeparator = !firstZGroup && rowIndex == 0;

            results += "<tr>";
            results += rowIndex == 0
                           ? "<td" + cellStyleForYield(zSeparator, false) + ">" + QString::number(z) + "</td>"
                           : "<td></td>";
            results += yieldCells(resRow, zSeparator);
            results += yieldCells(imfRow, zSeparator, true);
            results += "</tr>";
        }
        firstZGroup = false;
    }

    results += "<tr>"
               "<td></td>"
               "<td style=\"font-weight:bold; color:green;\">Total</td><td>" + QString::number(countResidue) + "</td>"
               "<td></td><td style=\"font-weight:bold; color:green;\">" + QString::number(residualSigmaTotalMb,'f',2) + "</td><td></td>"
               "<td style=\"font-weight:bold; color:green; border-left:1px dotted #b0b0b0;\">Total</td><td>" + formatYieldCount(imfTotal) + "</td>"
               "<td></td><td style=\"font-weight:bold; color:green;\">" + QString::number(imfSigmaTotalMb,'f',2) + "</td><td></td>"
               "</tr>";
    results += "</table><br>";

    return results;
}

}

QString buildMergedYieldTableHtml(const QString &title,
                                  Residual *resid,
                                  int length,
                                  int countResidue,
                                  double residualSigmaTotalMb,
                                  const std::map<std::pair<int, int>, AngularDistEntry> &imfEntries,
                                  double imfSigmaTotalMb,
                                  FILE *file_cs)
{
    const std::map<int, std::vector<YieldRow>> imfRowsByZ =
        buildImfYieldRowsByZ(imfEntries, imfSigmaTotalMb);
    return buildMergedYieldTableHtmlFromImfRows(title,
                                                resid,
                                                length,
                                                countResidue,
                                                residualSigmaTotalMb,
                                                imfRowsByZ,
                                                angularEntryWeightTotalForMap(imfEntries),
                                                imfSigmaTotalMb,
                                                file_cs);
}

YieldPlotData buildYieldPlotData(Residual *resid,
                                 int length,
                                 const std::map<std::pair<int, int>, AngularDistEntry> &imfEntries)
{
    std::map<int, double> imfByZ;
    for (const auto &it : imfEntries)
        imfByZ[it.first.first] += angularEntryWeightTotal(it.second);

    return buildYieldPlotDataFromCounts(residualYieldCountsByZ(resid, length), imfByZ);
}

namespace {

struct CompoundRunParams
{
    int zCN = 0;
    int aCN = 0;
    double excitationMeV = 0.0;
    double spin = 0.0;
    bool useIMF = false;
    bool useIMFEnhancement = false;
};

struct CompoundChunkJob
{
    int begin = 0;
    int end = 0;
    unsigned int seed = 1;
};

struct CompoundSummary
{
    int attemptedEvents = 0;
    int acceptedEvents = 0;
    double total = 0.0;
    double nFission = 0.0;
    double nImf = 0.0;
    double neutPreSad = 0.0;
    double neutSaddleToScission = 0.0;
    double neutHeavy = 0.0;
    double neutLight = 0.0;
    double imfTotal = 0.0;
    double symTotal = 0.0;
    double aRes = 0.0;
    double zRes = 0.0;
    double resTotal = 0.0;
    double neutMultEv = 0.0;
    double protMultEv = 0.0;
    double alpMultEv = 0.0;
    double gammaEnergy = 0.0;
    int countResidue = 0;
    std::array<Residual, 101> residues;
    ImfYieldCounts imfYieldByZN;

    CompoundSummary()
    {
        for (Residual &residue : residues)
        {
            residue.Z = 0;
            residue.A = 0;
            residue.count = 0;
            residue.name.clear();
        }
    }
};

int hashCompoundResidueName(const std::string &name)
{
    int sum = 0;
    for (char c : name) sum += c;
    return (sum * int(name.length())) % 101;
}

void mergeCompoundSummary(CompoundSummary &target, const CompoundSummary &source)
{
    target.attemptedEvents += source.attemptedEvents;
    target.acceptedEvents += source.acceptedEvents;
    target.total += source.total;
    target.nFission += source.nFission;
    target.nImf += source.nImf;
    target.neutPreSad += source.neutPreSad;
    target.neutSaddleToScission += source.neutSaddleToScission;
    target.neutHeavy += source.neutHeavy;
    target.neutLight += source.neutLight;
    target.imfTotal += source.imfTotal;
    target.symTotal += source.symTotal;
    target.aRes += source.aRes;
    target.zRes += source.zRes;
    target.resTotal += source.resTotal;
    target.neutMultEv += source.neutMultEv;
    target.protMultEv += source.protMultEv;
    target.alpMultEv += source.alpMultEv;
    target.gammaEnergy += source.gammaEnergy;
    target.countResidue += source.countResidue;

    for (int i = 0; i < int(target.residues.size()); ++i)
    {
        if (source.residues[i].count == 0) continue;
        target.residues[i].count += source.residues[i].count;
        target.residues[i].name = source.residues[i].name;
        target.residues[i].Z = source.residues[i].Z;
        target.residues[i].A = source.residues[i].A;
    }

    for (const auto &it : source.imfYieldByZN)
        target.imfYieldByZN[it.first] += it.second;
}

CompoundSummary simulateCompoundChunk(const CompoundRunParams &params,
                                      const CompoundChunkJob &job,
                                      const std::atomic_bool *cancelled)
{
    CRandom::seedCurrentThread(job.seed);

    CompoundSummary summary;
    CNucleus cn(params.zCN, params.aCN);
    cn.setVelocityCartesian();
    CAngle spin(CNucleus::pi / 2, double(0.));
    cn.setSpinAxis(spin);

    for (int event = job.begin; event < job.end; ++event)
    {
        if (cancelled && cancelled->load(std::memory_order_relaxed)) break;
        summary.attemptedEvents++;

        cn.setCompoundNucleus(params.excitationMeV, params.spin);

        if (params.useIMFEnhancement)
            cn.setWeightIMF();

        cn.decay();

        if (cn.abortEvent)
        {
            cn.reset();
            continue;
        }

        const int nFrag = cn.getNumberOfProducts();
        CNucleus *products = cn.getProducts(nFrag - 1);
        const double weight = products->getWeightFactor();

        const bool isSymmetricFission = cn.isSymmetricFission();
        const bool isAsymmetricFission = cn.isAsymmetricFission();
        const bool isResidue = cn.isResidue();

        if (isSymmetricFission)
        {
            summary.symTotal++;
            summary.nFission += weight;
            products = cn.getProducts(0);
            for (int i = 0; i < nFrag - 1; i++)
            {
                if (products->iZ == 0 && products->iA == 1)
                {
                    if (products->isSaddleToScission()) summary.neutSaddleToScission += weight;
                    else if (products->origin == 0) summary.neutPreSad += weight;
                    else if (products->origin == 2) summary.neutLight += weight;
                    else if (products->origin == 3) summary.neutHeavy += weight;
                }
                products = cn.getProducts();
            }
        }

        if (isAsymmetricFission)
        {
            summary.nImf += weight;
            summary.imfTotal++;
            const int zMaxEvap = cn.getZmaxEvap();

            CNucleus *imfProducts = cn.getProducts(0);
            for (int j = 0; j < nFrag && imfProducts; j++)
            {
                if (imfProducts->iZ > zMaxEvap)
                {
                    summary.imfYieldByZN[std::make_pair(imfProducts->iZ,
                                                        imfProducts->iA - imfProducts->iZ)] += weight;
                }
                imfProducts = cn.getProducts();
            }
        }

        summary.total += weight;

        products = cn.getProducts(nFrag - 1);
        if (isResidue)
        {
            const std::string productName = products->getName();
            const int residueIndex = hashCompoundResidueName(productName);
            Residual &residue = summary.residues[residueIndex];

            residue.count++;
            summary.countResidue++;
            residue.name = productName;
            residue.Z = products->iZ;
            residue.A = products->iA;
            summary.aRes += products->iA;
            summary.zRes += products->iZ;
            summary.resTotal += weight;
            summary.gammaEnergy += weight * products->getSumGammaEnergy();
        }

        products = cn.getProducts(0);
        for (int i = 0; i < nFrag - 1; i++)
        {
            if (products->iZ == 0 && products->iA == 1 && isResidue)
                summary.neutMultEv += weight;
            else if (products->iZ == 1 && products->iA == 1 && isResidue)
                summary.protMultEv += weight;
            else if (products->iZ == 2 && products->iA == 4 && isResidue)
                summary.alpMultEv += weight;

            products = cn.getProducts();
        }

        cn.reset();
        summary.acceptedEvents++;
    }

    cn.reset();
    return summary;
}

struct CompoundChunkRunner
{
    typedef CompoundSummary result_type;

    CompoundRunParams params;
    const std::atomic_bool *cancelled = nullptr;

    CompoundSummary operator()(const CompoundChunkJob &job) const
    {
        return simulateCompoundChunk(params, job, cancelled);
    }
};

QVector<CompoundChunkJob> makeCompoundJobs(int eventCount)
{
    QVector<CompoundChunkJob> jobs;
    if (eventCount <= 0) return jobs;

    const int idealThreads = std::max(1, QThread::idealThreadCount());
    const int chunkCount = std::min(eventCount, idealThreads * 12);
    const int chunkSize = (eventCount + chunkCount - 1) / chunkCount;
    const unsigned int baseSeed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());

    for (int begin = 0, index = 0; begin < eventCount; begin += chunkSize, ++index)
    {
        CompoundChunkJob job;
        job.begin = begin;
        job.end = std::min(eventCount, begin + chunkSize);
        job.seed = baseSeed + 1543u * static_cast<unsigned int>(index + 1);
        jobs.push_back(job);
    }

    return jobs;
}

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::execute_compound()
{
    QString filename_cs =  FFileName.split(".",Qt::SkipEmptyParts).at(0) + ".cs4";
    FILE *file_cs;
    file_cs=mfopen(filename_cs,"wt");

    QProgressDialog progress("Calculating..","Cancel",0,num_casc);
    progress.setWindowFlags(Qt::CustomizeWindowHint |
                            Qt::WindowCloseButtonHint);

    progress.setWindowIcon(QIcon(":/Gemini_logo.png"));
    progress.setWindowModality(Qt::WindowModal);
    progress.show();
    //----------------------------------------------
    CNucleus CN(iZCN,iACN); //constructor
    CN.setCompoundNucleus(fEx,fJ); //specify the excitation energy and spin
    CN.setVelocityCartesian(); // set initial CN velocity to zero
    CAngle spin(CNucleus::pi/2,(double)0.);
    CN.setSpinAxis(spin); //set the direction of the CN spin vector
    //--------------------------------------------------------------------

    double _SIGMA = 100;

    length = 101; // for hash;
    Residual resid[101];

    for(int i=0;i<length;i++) {  resid[i].count = 0;}

    int countResidue = 0;
    CN.setEvapMode(_optEvap);

    double total = 0.;
    double Nfission = 0.;
    double Nimf = 0.;
    double neutPreSad = 0.;
    double neutSaddleToScission = 0.;
    double neutHeavy = 0.;
    double neutLight = 0.;

    double imf_total = 0;
    double sym_total = 0;

    double Ares = 0.;
    double Zres = 0.;

    double resTotal = 0.;
    double neutMultEv = 0.;
    double protMultEv = 0.;
    double alpMultEv = 0.;
    double gammaEnergy = 0.;
    int counter=0;
    ImfYieldCounts imfYieldByZN;
    if (_useIMF) CN.setYesIMF();
    else CN.setNoIMF();

    CompoundRunParams params;
    params.zCN = iZCN;
    params.aCN = iACN;
    params.excitationMeV = fEx;
    params.spin = fJ;
    params.useIMF = _useIMF;
    params.useIMFEnhancement = _useIMF && _useIMFenh;

    QVector<CompoundChunkJob> jobs = makeCompoundJobs(num_casc);
    progress.setMaximum(std::max(1, int(jobs.size())));
    progress.setLabelText("Calculating");

    std::atomic_bool cancelled(false);
    CompoundChunkRunner runner;
    runner.params = params;
    runner.cancelled = &cancelled;

    QElapsedTimer parallelTimer;
    parallelTimer.start();

    QFuture<CompoundSummary> future = QtConcurrent::mapped(jobs, runner);
    QFutureWatcher<CompoundSummary> watcher;
    watcher.setFuture(future);

    while (!watcher.isFinished())
    {
        if (progress.wasCanceled())
        {
            cancelled.store(true, std::memory_order_relaxed);
            future.cancel();
        }
        progress.setValue(watcher.progressValue());
        qApp->processEvents(QEventLoop::AllEvents, 50);
        QThread::msleep(20);
    }

    watcher.waitForFinished();
    const qint64 parallelElapsedMs = parallelTimer.elapsed();
    progress.setValue(progress.maximum());
    showRuntime(parallelElapsedMs);

    CompoundSummary parallelSummary;
    for (int i = 0; i < future.resultCount(); ++i)
        mergeCompoundSummary(parallelSummary, future.resultAt(i));

    total = parallelSummary.total;
    Nfission = parallelSummary.nFission;
    Nimf = parallelSummary.nImf;
    neutPreSad = parallelSummary.neutPreSad;
    neutSaddleToScission = parallelSummary.neutSaddleToScission;
    neutHeavy = parallelSummary.neutHeavy;
    neutLight = parallelSummary.neutLight;
    imf_total = parallelSummary.imfTotal;
    sym_total = parallelSummary.symTotal;
    Ares = parallelSummary.aRes;
    Zres = parallelSummary.zRes;
    resTotal = parallelSummary.resTotal;
    neutMultEv = parallelSummary.neutMultEv;
    protMultEv = parallelSummary.protMultEv;
    alpMultEv = parallelSummary.alpMultEv;
    gammaEnergy = parallelSummary.gammaEnergy;
    countResidue = parallelSummary.countResidue;
    counter = parallelSummary.acceptedEvents;
    imfYieldByZN = parallelSummary.imfYieldByZN;
    for (int i = 0; i < length; ++i)
        resid[i] = parallelSummary.residues[i];

    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade END
    int _INPUT = 2;
    double Sconst = _SIGMA/double(qMax(1, counter));

    if(file_cs) {
        fprintf(file_cs,"!-title- Zcomp  Ncomp  Energy  CSfus  Input _AJNUC\n");
        fprintf(file_cs,"!Gemini %d %d %.3f %.3g %d %.3g",iZCN,iACN-iZCN,fEx,_SIGMA, _INPUT, (_INPUT==1 ? 0: fJ));
    }

    results =  "<h1 align=center style=\"color:red;\">Gemini <em><small>GUI</small></em></h1>";
    results += "<h2 align=center style=\"color:green;\">Statistical Decay Code</h2><br><br>";

    results += "<h3 align=center style=\"color:blue;\">Starting Conditions</h3>";
    results += "<table cellpadding=5 align=\"center\"><tr><th></th><th>Z</th><th> N</th><th>A</th></tr>";
    results += "<tr><td><em>Compound nucleus</em></td><td align=center>" + QString::number(CN.iZ) +
               "</td><td align=center>" + QString::number(CN.iA-CN.iZ) + " </td><td align=center> " + QString::number(CN.iA) +
               "</td></tr></table><p>&nbsp;</p>";

    results += "<table cellpadding=5 align=\"center\"><tr><td><em> Compound nucleus Excitation energy</em></td><td>  "
               + QString::number(fEx,'f',2) + " </td><td>MeV</td></tr>";
    results += "<tr><td><em> Compound nucleus spin</em></td><td>  " + QString::number(fJ,'f',2) + " </td><td>&#8463;</td></tr>"
               + "</table><br><br>";
    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

    results += "<h3 align=\"center\" style=\"color: blue\"> Decay Product Summary </h3>";
    results += "<table cellpadding=\"5\" align=\"center\">";

    results += "<tr><th style=\"color:green\" align=\"left\">Result</th><th style=\"color:green;\" align=\"center\">Number</th></tr>";
    results += "<tr><th align=\"left\">Intermediate Mass Fragments</th><td align=\"center\">" + QString::number(imf_total) + "</td></tr>";
    results += "<tr><th align=\"left\">Symmetric Fission</th><td align=\"center\">" + QString::number(sym_total) + "</td></tr>";
    results += "<tr><th align=\"left\">Residual Nuclei</th><td align=\"center\">" + QString::number(countResidue) + "</td></tr>";
    results += "<tr><th align=\"left\">TOTAL</th><td align=\"center\"><b>" + QString::number(imf_total+sym_total+countResidue) + "<b></td></tr>";
    results += "</table>;";
    //--------------------------------------------------------------

    results += buildMergedYieldTableHtmlFromImfRows("Yields of Residual Nuclei and IMF Particles ",
                                                    resid,
                                                    length,
                                                    countResidue,
                                                    _SIGMA,
                                                    buildImfYieldRowsByZ(imfYieldByZN, Nimf*Sconst),
                                                    imfYieldCountTotal(imfYieldByZN),
                                                    Nimf*Sconst,
                                                    file_cs);

    if(file_cs) fclose(file_cs);
    //------------------------------------------------------------------------
    if (resTotal > 0.)
    {
        results += "<BR> <center><b style=\"color:blue\">residue (without IMF or Symmetric fission emissions)</b><br>";
        results += "<table cellpadding=\"5\" align=\"center\">";
        results += " <tr><td>  average residue is Z  </td><td>" + QString::number(Zres/resTotal, 'f', 2) +  "</td></tr>";
        results += " <tr><td>  average residue is A  </td><td>" + QString::number(Ares/resTotal, 'f', 2) +  "</td></tr>";
        results += " <tr><td>  average residue is N  </td><td>" + QString::number((Ares-Zres)/resTotal, 'f', 2) +  "</td></tr>";
        results += " <tr><td> average neutron multiplicity  </td><td> " + QString::number(neutMultEv/resTotal, 'f', 2) + "</td></tr>";
        results += " <tr><td> average proton multiplicity </td><td> " + QString::number(protMultEv/resTotal, 'f', 2) + "</td></tr>";
        results += " <tr><td> average alpha multiplicity  </td><td> " + QString::number(alpMultEv/resTotal, 'f', 2) + "</td></tr>";
        results += " <tr><td> average energy in gamma rays </td><td> " + QString::number(gammaEnergy/resTotal, 'f', 2) +
                   " MeV </td></tr></table></center>";
    }

    if (Nfission > 0.)
    {
        results += "<BR> <center><b style=\"color:blue\">neutron multiplicities in fission</b><br>";
        results += "<table cellpadding=\"5\" align=\"center\">";
        results += " <tr><td> presaddle</td><td> " + QString::number(neutPreSad/Nfission, 'f', 2)  +  "</td></tr>";
        results += " <tr><td> saddle-to-scission</td><td> " + QString::number(neutSaddleToScission/Nfission, 'f', 2) +  "</td></tr>";
        results += " <tr><td> post-scission light frag</td><td> " + QString::number(neutLight/Nfission, 'f', 2)  +  "</td></tr>";
        results += " <tr><td> post-scission heavy frag</td><td> " + QString::number(neutHeavy/Nfission, 'f', 2)  +  "</td></tr>";
        results += " </table></center>";
    }

    if(imf_total > 0)
    {
        results += "<br><center> <b style=\"color:blue\">Intermediate Mass Fragment (Z > " + QString::number(CN.getZmaxEvap()) + ") (IMF) </b><br>" ;
        results += " <table cellpadding=\"5\" align=\"center\">"
                   "<tr><td> IMF prob </td><td> " + QString::number(Nimf/total, 'g', 3) + "</td><td> </td></tr>";
        results += "<tr><td> Cross section </td><td>" + QString::number(Nimf*Sconst, 'g', 3) + "</td><td> mb </td></tr>";
        results += "</table><br>";
    }
    //------------------------------------------------------------------------
    printGeminiProperties(results);
    //------------------------------------------------------------------------
    qApp->processEvents();

    const YieldPlotData yieldPlot =
        buildYieldPlotDataFromImfCounts(resid, length, imfYieldByZN);

    Result_Widget *ress = new Result_Widget(results, yieldPlot);
    ress->show();

    CN.reset();
}
