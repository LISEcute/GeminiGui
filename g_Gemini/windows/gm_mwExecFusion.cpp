#include "gm_mainwindow.h"
#include "ui_gm_mainwindow.h"

#include <QDateTime>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFileInfo>
#include <QFuture>
#include <QFutureWatcher>
#include <QProgressDialog>
#include <QThread>
#include <QtConcurrent>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <map>
#include <vector>
#include "g_Gemini/source/CNucleus.h"
#include "g_Gemini/source/CFus.h"
#include "gm_results.h"
#include "gm_ftype.h"
#include "gm_about.h"
#include "gm_angular_distribution.h"

extern bool _useAME;
extern bool _useIMF;
extern bool _useIMFenh;
extern int  _optEvap;

extern bool Residual_compare(Residual a, Residual b);
extern FILE *mfopen(const QString& filename, const char* operand);
extern QString buildMergedYieldTableHtml(const QString &title,
                                         Residual *resid,
                                         int length,
                                         int countResidue,
                                         double residualSigmaTotalMb,
                                         const std::map<std::pair<int, int>, AngularDistEntry> &imfEntries,
                                         double imfSigmaTotalMb,
                                         FILE *file_cs);
extern YieldPlotData buildYieldPlotData(Residual *resid,
                                        int length,
                                        const std::map<std::pair<int, int>, AngularDistEntry> &imfEntries);

namespace {

QChar superscriptDigit(QChar digit)
{
    switch (digit.unicode())
    {
    case '0': return QChar(0x2070);
    case '1': return QChar(0x00B9);
    case '2': return QChar(0x00B2);
    case '3': return QChar(0x00B3);
    case '4': return QChar(0x2074);
    case '5': return QChar(0x2075);
    case '6': return QChar(0x2076);
    case '7': return QChar(0x2077);
    case '8': return QChar(0x2078);
    case '9': return QChar(0x2079);
    default: return digit;
    }
}

QString superscriptLeadingMassNumber(const QString &nuclideName)
{
    QString displayName;
    int index = 0;
    while (index < nuclideName.size() && nuclideName.at(index).isDigit())
    {
        displayName += superscriptDigit(nuclideName.at(index));
        index++;
    }

    displayName += nuclideName.mid(index);
    return displayName;
}

QString resultWindowTitle(const QString &projectileName,
                          const QString &targetName,
                          double energyMeV)
{
    return "Gemini: Results ("
           + superscriptLeadingMassNumber(projectileName) + "_"
           + superscriptLeadingMassNumber(targetName) + "_"
           + QString::number(energyMeV, 'f', 0) + "MeV)";
}

struct FusionRunParams
{
    int zCN = 0;
    int aCN = 0;
    int projectileA = 0;
    int targetA = 0;
    double excitationMeV = 0.0;
    double labEnergyMeV = 0.0;
    double betaCN = 0.0;
    bool useIMF = false;
    bool useIMFEnhancement = false;
    int evapMode = 2;
    std::vector<double> spinCdf;
};

struct FusionChunkJob
{
    int begin = 0;
    int end = 0;
    unsigned int seed = 1;
};

struct FusionSummary
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
    std::map<std::pair<int, int>, AngularDistEntry> angularDistByZN;
    std::map<std::pair<int, int>, AngularDistEntry> imfAngularByZN;
    AngularDistEntry neutronAngular;
    AngularDistEntry protonAngular;
    AngularDistEntry alphaAngular;
    AngularDistEntry gammaAngular;

    FusionSummary()
    {
        for (Residual &residue : residues)
        {
            residue.Z = 0;
            residue.A = 0;
            residue.count = 0;
            residue.name.clear();
        }
        neutronAngular.z = 0;
        neutronAngular.n = 1;
        protonAngular.z = 1;
        protonAngular.n = 0;
        alphaAngular.z = 2;
        alphaAngular.n = 2;
        gammaAngular.z = 0;
        gammaAngular.n = 0;
    }
};

void appendAngularEntry(AngularDistEntry &target, const AngularDistEntry &source)
{
    target.z = source.z;
    target.n = source.n;
    target.kineticEnergy.insert(target.kineticEnergy.end(),
                                source.kineticEnergy.begin(),
                                source.kineticEnergy.end());
    target.thetaDeg.insert(target.thetaDeg.end(),
                           source.thetaDeg.begin(),
                           source.thetaDeg.end());
    target.vz.insert(target.vz.end(), source.vz.begin(), source.vz.end());
    target.vxy.insert(target.vxy.end(), source.vxy.begin(), source.vxy.end());
    target.cmEnergy.insert(target.cmEnergy.end(),
                           source.cmEnergy.begin(),
                           source.cmEnergy.end());
    target.weight.insert(target.weight.end(), source.weight.begin(), source.weight.end());
}

void mergeAngularMap(std::map<std::pair<int, int>, AngularDistEntry> &target,
                     const std::map<std::pair<int, int>, AngularDistEntry> &source)
{
    for (const auto &it : source)
        appendAngularEntry(target[it.first], it.second);
}

void mergeFusionSummary(FusionSummary &target, const FusionSummary &source)
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

    mergeAngularMap(target.angularDistByZN, source.angularDistByZN);
    mergeAngularMap(target.imfAngularByZN, source.imfAngularByZN);
    appendAngularEntry(target.neutronAngular, source.neutronAngular);
    appendAngularEntry(target.protonAngular, source.protonAngular);
    appendAngularEntry(target.alphaAngular, source.alphaAngular);
    appendAngularEntry(target.gammaAngular, source.gammaAngular);
}

int hashResidueName(const std::string &name)
{
    int sum = 0;
    for (char c : name) sum += c;
    return (sum * int(name.length())) % 101;
}

void computeFusionKinematics(CNucleus *particle,
                             double betaCN,
                             double &keLab,
                             double &thetaLabDeg,
                             double &vzLab,
                             double &vxy,
                             double &cmEnergy)
{
    double *vel = particle->getVelocityVector();

    const double vx = vel[0] / 30.f;
    const double vy = vel[1] / 30.f;
    const double vzLocal = vel[2] / 30.f;
    vzLab = vzLocal + betaCN;
    vxy = std::sqrt(vx * vx + vy * vy);
    const double massMeV = particle->iA * 931.49432f;

    const double betaTot = std::sqrt(vx * vx + vy * vy + vzLab * vzLab);
    thetaLabDeg = 0.f;
    if (betaTot > 0.f)
    {
        double c = vzLab / betaTot;
        if (c > 1.f) c = 1.f;
        if (c < -1.f) c = -1.f;
        thetaLabDeg = static_cast<double>(std::acos(c) * 57.29577951308232);
    }

    keLab = 0.5f * massMeV * (vx * vx + vy * vy + vzLab * vzLab);
    cmEnergy = 0.5f * massMeV * (vx * vx + vy * vy + vzLocal * vzLocal);
}

int chooseSpinFromCdf(const std::vector<double> &spinCdf)
{
    const double ran = CRandom::instance()->Rndm();
    int l = 0;
    while (l + 1 < int(spinCdf.size()) && ran >= spinCdf[l])
        ++l;
    return l;
}

FusionSummary simulateFusionChunk(const FusionRunParams &params,
                                  const FusionChunkJob &job,
                                  const std::atomic_bool *cancelled)
{
    CRandom::seedCurrentThread(job.seed);

    FusionSummary summary;
    CNucleus cn(params.zCN, params.aCN);

    for (int event = job.begin; event < job.end; ++event)
    {
        if (cancelled && cancelled->load(std::memory_order_relaxed)) break;
        summary.attemptedEvents++;

        const int l = chooseSpinFromCdf(params.spinCdf);
        cn.setCompoundNucleus(params.excitationMeV, double(l));

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
                    double keLab = 0.f;
                    double thetaLabDeg = 0.f;
                    double vzLab = 0.f;
                    double vxy = 0.f;
                    double cmEnergy = 0.f;
                    computeFusionKinematics(imfProducts, params.betaCN,
                                            keLab, thetaLabDeg, vzLab, vxy, cmEnergy);
                    addAngularSample(summary.imfAngularByZN,
                                     imfProducts->iZ,
                                     imfProducts->iA - imfProducts->iZ,
                                     keLab,
                                     thetaLabDeg,
                                     vzLab,
                                     vxy,
                                     cmEnergy,
                                     weight);
                }
                imfProducts = cn.getProducts();
            }
        }

        summary.total += weight;

        products = cn.getProducts(nFrag - 1);
        if (isResidue)
        {
            const std::string productName = products->getName();
            const int residueIndex = hashResidueName(productName);
            Residual &residue = summary.residues[residueIndex];

            residue.count++;
            summary.countResidue++;
            residue.name = productName;
            residue.Z = products->iZ;
            residue.A = products->iA;
            summary.aRes += products->iA;
            summary.zRes += products->iZ;
            summary.resTotal += weight;
            const double eventGammaEnergy = products->getSumGammaEnergy();
            summary.gammaEnergy += weight * eventGammaEnergy;

            double keLab = 0.f;
            double thetaLabDeg = 0.f;
            double vzLab = 0.f;
            double vxy = 0.f;
            double cmEnergy = 0.f;

            computeFusionKinematics(products, params.betaCN,
                                    keLab, thetaLabDeg, vzLab, vxy, cmEnergy);

            addAngularSample(summary.angularDistByZN,
                             products->iZ,
                             products->iA - products->iZ,
                             keLab,
                             thetaLabDeg,
                             vzLab,
                             vxy);

            if (eventGammaEnergy > 0.f)
            {
                addAngularSample(summary.gammaAngular,
                                 eventGammaEnergy,
                                 thetaLabDeg,
                                 0.0f,
                                 0.0f,
                                 eventGammaEnergy);
            }
        }

        products = cn.getProducts(0);
        for (int i = 0; i < nFrag - 1; i++)
        {
            if (products->iZ == 0 && products->iA == 1 && isResidue)
            {
                summary.neutMultEv += weight;
                double keLab = 0.f;
                double thetaLabDeg = 0.f;
                double vzLab = 0.f;
                double vxy = 0.f;
                double cmEnergy = 0.f;
                computeFusionKinematics(products, params.betaCN,
                                        keLab, thetaLabDeg, vzLab, vxy, cmEnergy);
                addAngularSample(summary.neutronAngular, keLab, thetaLabDeg, vzLab, vxy,
                                 cmEnergy);
            }
            else if (products->iZ == 1 && isResidue)
            {
                if(products->iA == 1 )
                {
                    summary.protMultEv += weight;
                    double keLab = 0.f;
                    double thetaLabDeg = 0.f;
                    double vzLab = 0.f;
                    double vxy = 0.f;
                    double cmEnergy = 0.f;
                    computeFusionKinematics(products, params.betaCN,
                                            keLab, thetaLabDeg, vzLab, vxy, cmEnergy);
                    addAngularSample(summary.protonAngular, keLab, thetaLabDeg, vzLab, vxy,
                                     cmEnergy);
                }
            }
            else if (products->iZ == 2 && isResidue)
            {
                if(products->iA == 4)
                {
                    summary.alpMultEv += weight;
                    double keLab = 0.f;
                    double thetaLabDeg = 0.f;
                    double vzLab = 0.f;
                    double vxy = 0.f;
                    double cmEnergy = 0.f;
                    computeFusionKinematics(products, params.betaCN,
                                            keLab, thetaLabDeg, vzLab, vxy, cmEnergy);
                    addAngularSample(summary.alphaAngular, keLab, thetaLabDeg, vzLab, vxy,
                                     cmEnergy);
                }
            }
            products = cn.getProducts();
        }

        cn.reset();
        summary.acceptedEvents++;
    }

    cn.reset();
    return summary;
}

struct FusionChunkRunner
{
    typedef FusionSummary result_type;

    FusionRunParams params;
    const std::atomic_bool *cancelled = nullptr;

    FusionSummary operator()(const FusionChunkJob &job) const
    {
        return simulateFusionChunk(params, job, cancelled);
    }
};

QVector<FusionChunkJob> makeFusionJobs(int eventCount)
{
    QVector<FusionChunkJob> jobs;
    if (eventCount <= 0) return jobs;

    const int idealThreads = std::max(1, QThread::idealThreadCount());
    const int chunkCount = std::min(eventCount, idealThreads * 12);
    const int chunkSize = (eventCount + chunkCount - 1) / chunkCount;
    const unsigned int baseSeed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());

    for (int begin = 0, index = 0; begin < eventCount; begin += chunkSize, ++index)
    {
        FusionChunkJob job;
        job.begin = begin;
        job.end = std::min(eventCount, begin + chunkSize);
        job.seed = baseSeed + 7919u * static_cast<unsigned int>(index + 1);
        jobs.push_back(job);
    }

    return jobs;
}

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::execute_fusion()
{
    QString filename_cs =  FFileName.split(".",Qt::SkipEmptyParts).at(0) + ".cs4";
    FILE *file_cs;
    file_cs=mfopen(filename_cs,"wt");

    QProgressDialog progress("Calculating..","Cancel",0,num_events);
    progress.setWindowFlags(Qt::CustomizeWindowHint |
                            Qt::WindowCloseButtonHint);
    progress.setWindowIcon(QIcon(":/Gemini_logo.png"));
    progress.setWindowModality(Qt::WindowModal);
    progress.show();
    //--------------------------------------------

    CFus fus(Zp,Ap,Zt,At,Elab,dif);

    CNucleus CN(fus.iZcn,fus.iAcn); //constructor
    double ffEx = fus.Ex;             //excitation energy of compound nucleus

    if(spinOption==1) l0 = fus.getBassL();  // get maximum spin from Bass model

    //construct fusion spin distribution

    int lmax = (int)l0 +5;  // maxium spin considered
    double prob[301];
    lmax = qMin(300,lmax);
    double sum = 0.;

    for (int l=0;l<=lmax;l++)
    {
        prob[l] =  (double)(2*l+1);

        if (dif > 0.) prob[l] /= (1.+exp(((double)l-l0)/dif));
        else if ( l > l0)  prob[l] = 0.;

        sum += prob[l];
    }

    for (int l=0;l<=lmax;l++)
    {
        prob[l] /= sum;
        if (l > 0) prob[l] += prob[l-1];
    }

    double xfus = sum*fus.plb;

    CN.setEvapMode(_optEvap);  // force a fully Hauser-Feshbach calculation

    // if IMF probablities are small and you are not interested
    // save time and turn them off.
    //CNucleus::setNoIMF();

    // don't like the default evaporation parameters
    //CLevelDensity::setAfAn(1.036);  // change ratio of saddle to ground state
    //level density parameters
    //CLevelDensity::setLittleA(8.); //change level density parameter to A/8.

    //write out statistical model parameters
    //cout << "CN.printParameters():" << std::endl;
    //CN.printParameters();
    //cout << "cone printing parameters" << std::endl;
    //define and zero events parameters

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

    length = 101; // for hash;
    Residual resid[101];
    for(int i=0;i<length;i++) resid[i].count = 0;

    int countResidue = 0;
    int counter=0;
    std::map<std::pair<int, int>, AngularDistEntry> angularDistByZN;
    std::map<std::pair<int, int>, AngularDistEntry> imfAngularByZN;
    AngularDistEntry neutronAngular;
    neutronAngular.z = 0;
    neutronAngular.n = 1;
    AngularDistEntry protonAngular;
    protonAngular.z = 1;
    protonAngular.n = 0;
    AngularDistEntry alphaAngular;
    alphaAngular.z = 2;
    alphaAngular.n = 2;
    AngularDistEntry gammaAngular;
    gammaAngular.z = 0;
    gammaAngular.n = 0;

    const double betaCN = std::sqrt(
        2.0f * (Elab * Ap / (Ap + At)) / (931.49432f * (Ap + At))
        );

    CN.setEvapMode(_optEvap);
    if (_useIMF) CN.setYesIMF();
    else CN.setNoIMF();

    FusionRunParams params;
    params.zCN = fus.iZcn;
    params.aCN = fus.iAcn;
    params.projectileA = Ap;
    params.targetA = At;
    params.excitationMeV = ffEx;
    params.labEnergyMeV = Elab;
    params.betaCN = betaCN;
    params.useIMF = _useIMF;
    params.useIMFEnhancement = _useIMF && _useIMFenh;
    params.evapMode = _optEvap;
    params.spinCdf.assign(prob, prob + lmax + 1);

    QVector<FusionChunkJob> jobs = makeFusionJobs(num_events);
    progress.setMaximum(std::max(1, int(jobs.size())));
    progress.setLabelText("Calculating");

    std::atomic_bool cancelled(false);
    FusionChunkRunner runner;
    runner.params = params;
    runner.cancelled = &cancelled;

    QElapsedTimer parallelTimer;
    parallelTimer.start();

    QFuture<FusionSummary> future = QtConcurrent::mapped(jobs, runner);

    QFutureWatcher<FusionSummary> watcher;
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

    FusionSummary parallelSummary;
    for (int i = 0; i < future.resultCount(); ++i)
        mergeFusionSummary(parallelSummary, future.resultAt(i));

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
    angularDistByZN = parallelSummary.angularDistByZN;
    imfAngularByZN = parallelSummary.imfAngularByZN;
    neutronAngular = parallelSummary.neutronAngular;
    protonAngular = parallelSummary.protonAngular;
    alphaAngular = parallelSummary.alphaAngular;
    gammaAngular = parallelSummary.gammaAngular;
    for (int i = 0; i < length; ++i)
        resid[i] = parallelSummary.residues[i];

    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade BEGIN
    int _INPUT = 1;
    double Sconst = xfus/double(qMax(1, counter));

    if(file_cs) {
        fprintf(file_cs,"!-title- Zcomp  Ncomp  Energy  CSfus  Input _AJNUC\n");
        fprintf(file_cs,"!Gemini %d %d %.3f %.3g %d %.3g",CN.iZ,CN.iN,CN.fEx, xfus, _INPUT, (_INPUT==1 ? l0: l0));
    }


    results = "<h1 align=center style=\"color:red;\">Gemini <em><small>GUI</small></em></h1>";
    results += "<h2 align=center style=\"color:green;\">Statistical Decay Code</h2><br><br>";

    results += "<h3 align=center style=\"color:blue;\">Starting Conditions</h3>";
    //-------------------------------------------------------------------------
    results += "<table align=\"center\" cellpadding=5><tr><th></th><th>Z</th><th> N</th><th>A</th><th style=\"font-weight:bold\"><sup>A </sup>El</th></tr>";

    results += " <tr><td><em>Projectile</em></td><td align=center>" + QString::number(Zp) +
               " </td><td align=center>" + QString::number(Ap-Zp) + " </td><td align=center>" + QString::number(Ap) +
               " </td><td align=center>" + fus.p.getGName()  + "</td></tr>";

    results += "<tr><td><em>Target</em></td><td align=center>" + QString::number(Zt) +
               "</td><td align=center>" + QString::number(At-Zt) + "</td><td align=center>" + QString::number(At) +
               " </td><td align=center>" + fus.t.getGName()  + "</td></tr>";

    results += "<tr><td><em>Compound nucleus</em></td><td align=center>" + QString::number(fus.iZcn) +
               " </td><td align=center>" + QString::number(CN.iN) + " </td><td align=center> " + QString::number(fus.iAcn) +
               " </td><td align=center>" +  fus.c.getGName()  + "</td></tr>" +
               " </table><br><br>";
    //-------------------------------------------------------------------------

    results += "<table align=\"center\" cellpadding=5><tr><td><em> Bombarding energy (MeV)</em></td><td>  " + QString::number(Elab,'f',2) + " </td></tr> " +
               "<tr><td><em> Center of Mass energy (MeV)</em></td><td>  " + QString::number(fus.Ecm,'f',3) + " </td></tr>"
               +  "<tr><td><em> Compound nucleus Excitation energy (MeV)</em></td><td>  " + QString::number(fus.Ex,'f',2) + " </td></tr>"
               +  "<tr><td><em> Q-value of reaction (MeV)</em></td><td>  " + QString::number(fus.Qval,'f',3) + " </td></tr>"
               +  "<tr><td><em> Compound nucleus recoil energy (MeV)</em></td><td>  " + QString::number(Elab*Ap/(Ap+At),'f',3) + " </td></tr>"
               +  "<tr><td><em> Compound nucleus recoil velocity (cm/ns) </em></td><td>  " + QString::number(30.*sqrt(2*(Elab*Ap/(Ap+At))/(931.49432 *(Ap+At))),'e',3) + " </td></tr>"
               +  "<tr><td><em> Compound nucleus recoil (&beta;)</em> </td><td>  " + QString::number(sqrt(2*(Elab*Ap/(Ap+At))/(931.49432 *(Ap+At))),'e',3) + " </td></tr>"
               +  "<tr><td><em> Beam velocity (cm/ns) </em></td><td>  " + QString::number(30.*sqrt(2.0*Elab/(931.49432*Ap)),'e',3) + " </td></tr>"
               +  "<tr><td><em> Beam velocity (&beta;) </em></td><td>  " + QString::number(sqrt(2.0*Elab/(931.49432*Ap)),'e',3) + " </td></tr>"
               +  "</table><p>&nbsp;</p>";

    results += "<table align=\"center\" cellpadding=5>"
               "<tr><td><em> diffuseness </em></td><td>  " + QString::number(fus.dif,'f',2) + "</td><td> &#8463;</td></tr> "
               +  "<tr><td><em> Fusion cross section</em></td><td>  " + QString::number(sum*fus.plb,'f',2) + "</td><td> mb</td></tr>"
               +  "<tr><td><em> Bass L</em></td><td>  " + QString::number(fus.getBassL(),'f',2) + "</td><td> &#8463; </td></tr>"
               +  "<tr><td><em> L0 </em></td><td>  " + QString::number(fus.getL0(sum*fus.plb),'f',2) + "<td> &#8463; </td></tr>"
               +  "<tr><td><em> Bass cross section</em></td><td>  " + QString::number(fus.getBassXsec(),'f',2) + "<td> mb </td></tr>"
               +  "<tr><td><em> Excitation energy</em></td><td>  " + QString::number(ffEx,'f',2) + "<td> MeV</td></tr>"
               +  "<tr><td><em> Critical spin</em></td><td>  " + QString::number(l0,'f',1) + "<td> &#8463;</td></tr>"
               +  "</table><br><br>";

    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

    results += "<h3 align=\"center\" style=\"color: blue\"> Fusion Product Summary </h3>";
    results += "<table cellpadding=\"5\" align=\"center\">";

    results += "<tr><th style=\"color:green\" align=\"left\">Result</th><th style=\"color:green;\" align=\"center\">Number</th></tr>";
    results += "<tr><th align=\"left\">Intermediate Mass Fragments</th><td align=\"center\">" + QString::number(imf_total) + "</td></tr>";
    results += "<tr><th align=\"left\">Symmetric Fission</th><td align=\"center\">" + QString::number(sym_total) + "</td></tr>";
    results += "<tr><th align=\"left\">Residual Nuclei</th><td align=\"center\">" + QString::number(countResidue) + "</td></tr>";
    results += "<tr><th align=\"left\">TOTAL</th><td align=\"center\"><b>" + QString::number(imf_total+sym_total+countResidue) + "<b></td></tr>";
    results += "</table>;";
    //------------------------------------------------------
    results += buildMergedYieldTableHtml("Yields of Residual Nuclei and IMF Particles ",
                                         resid,
                                         length,
                                         countResidue,
                                         xfus,
                                         imfAngularByZN,
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
    extern double _LowLimit, _HighLimit;

    const double recoilBetaCN =
        std::sqrt(2.0 * (Elab * Ap / (Ap + At)) /
                  (931.49432 * (Ap + At)));
    const double imfYieldSigmaMb = double(Nimf*Sconst);

    AngularResultTab residualAngular;
    residualAngular.label = "Residual Angular Distributions";
    residualAngular.html = buildAngularDistributionHtmlPACEStyle(
        angularDistByZN,
        xfus,
        qMax(1, counter),
        _LowLimit,
        _HighLimit,
        fus.Ex,
        fus.iAcn,
        recoilBetaCN,
        "Fusion mode",
        0,
        -1.0,
        neutronAngular,
        protonAngular,
        alphaAngular,
        gammaAngular
        );
    residualAngular.entries = angularDistByZN;
    residualAngular.sigmaTotal = xfus;
    residualAngular.nEvents = qMax(1, counter);
    residualAngular.lowLimitPercent = _LowLimit;
    residualAngular.highLimitPercent = _HighLimit;
    residualAngular.title = "Fusion mode";
    residualAngular.neutronEntry = neutronAngular;
    residualAngular.protonEntry = protonAngular;
    residualAngular.alphaEntry = alphaAngular;
    residualAngular.gammaEntry = gammaAngular;
    residualAngular.compoundExcitationMeV = fus.Ex;
    residualAngular.compoundA = fus.iAcn;
    residualAngular.compoundZ = fus.iZcn;
    residualAngular.recoilBetaCN = recoilBetaCN;
    residualAngular.mdir = 0;

    AngularResultTab imfAngular;
    imfAngular.label = "IMF Angular Distributions";
    imfAngular.html = buildAngularDistributionHtmlPACEStyle(
        imfAngularByZN,
        xfus,
        qMax(1, counter),
        _LowLimit,
        _HighLimit,
        fus.Ex,
        fus.iAcn,
        recoilBetaCN,
        "Fusion mode - IMF",
        0,
        _useIMFenh ? imfYieldSigmaMb : -1.0
        );
    imfAngular.entries = imfAngularByZN;
    imfAngular.sigmaTotal = xfus;
    imfAngular.nEvents = qMax(1, counter);
    imfAngular.lowLimitPercent = _LowLimit;
    imfAngular.highLimitPercent = _HighLimit;
    imfAngular.title = "Fusion mode - IMF";
    imfAngular.compoundExcitationMeV = fus.Ex;
    imfAngular.compoundA = fus.iAcn;
    imfAngular.compoundZ = fus.iZcn;
    imfAngular.recoilBetaCN = recoilBetaCN;
    imfAngular.mdir = 0;

    const YieldPlotData yieldPlot =
        buildYieldPlotData(resid, length, imfAngularByZN);

    Result_Widget *ress = new Result_Widget(results, residualAngular, imfAngular, yieldPlot);
    ress->setWindowTitle(resultWindowTitle(QString::fromStdString(fus.p.getName()),
                                           QString::fromStdString(fus.t.getName()),
                                           Elab));
    ress->show();
    CN.reset();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::printGeminiProperties(QString &results)
{
    QDateTime mydate=QDateTime::currentDateTime();
    QFileInfo fI(FFileName);
    QString windowName = fI.baseName();
    extern bool _useAME;

    results += "<br><center> <b style=\"color:blue\">GEMINI-GUI code properties </b><br>" ;
    results += " <table cellpadding=\"5\" align=\"center\">";
    results += "<tr><td>" + QString(Gemini_version) + "</td><td>" + QString(Gemini_date) + "</td></tr>";
    results += "<tr><td>File name</td><td>" + windowName + "</td></tr>";
    results += "<tr><td>File created</td><td>" + mydate.toString("MM/dd/yyyy HH:mm:ss") + "</td></tr>";
    results += "<tr><td>Mass table</td><td>" + QString(_useAME ? "AME2016" : "GEMINI traditional") + "</td></tr>";
    results += "<tr><td>Evaporation mode</td><td>" + QString::number(_optEvap) + "</td></tr>";
    results += "<tr><td>IMF emission</td><td>" + QString(_useIMF ? "yes" : "no") + "</td></tr>";
    results += "<tr><td>IMF emission enhanced</td><td>" + QString(_useIMFenh ? "yes" : "no") + "</td></tr>";
    results += "<tr><td>spin option</td><td>" + QString::number(spinOption) + "</td></tr>";

    results += "</table></center>";

}
