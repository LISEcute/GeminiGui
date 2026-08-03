#include "gm_mainwindow.h"
#include "ui_gm_mainwindow.h"

#include <QProgressDialog>
#include <algorithm>
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
               "<th rowspan=\"2\">Z</th>"
               "<th colspan=\"5\">Residual Nuclei</th>"
               "<th colspan=\"5\" style=\"border-left:1px dotted #b0b0b0;\">IMF Particles</th>"
               "</tr>";
    results += "<tr style=\"color: green\">"
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
    CAngle spin(CNucleus::pi/2,(float)0.);
    CN.setSpinAxis(spin); //set the direction of the CN spin vector
    //--------------------------------------------------------------------

    float _SIGMA = 100;

    length = 101; // for hash;
    Residual resid[101];

    for(int i=0;i<length;i++) {  resid[i].count = 0;}

    int countResidue = 0;
    CN.setEvapMode(_optEvap);

    float total = 0.;
    float Nfission = 0.;
    float Nimf = 0.;
    float neutPreSad = 0.;
    float neutSaddleToScission = 0.;
    float neutHeavy = 0.;
    float neutLight = 0.;

    float imf_total = 0;
    float sym_total = 0;

    float Ares = 0.;
    float Zres = 0.;

    float resTotal = 0.;
    float neutMultEv = 0.;
    float protMultEv = 0.;
    float alpMultEv = 0.;
    float gammaEnergy = 0.;
    int counter=0;
    ImfYieldCounts imfYieldByZN;
    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade BEGIN
    const int progressStep = std::max(1, num_casc / 1000);
    for (int i=0;i<num_casc;i++)
    {
        if(progress.wasCanceled()){break;}
        if (i % progressStep == 0 || i == num_casc - 1) {
            progress.setValue(i);
            qApp->processEvents();
        }

        // ---------------------  here place of probability random to get L begin
        // ---------------------  here place of probability random to get L end

        CN.setCompoundNucleus(fEx,fJ);

        //---------------------------------------------------------------------- IMF begin
        //if you are interested in IMF emission at low excitation energy
        //then turn IMF weighting on
        if(_useIMF)  {
            CN.setYesIMF();
            if(_useIMFenh) CN.setWeightIMF();// turn on enhanced IMF emission
        }
        else CN.setNoIMF();
        //---------------------------------------------------------------------- IMF end

        CN.decay(); //decay the compound nucleus

        if (CN.abortEvent)
        {
            cout << "abort event\n";
            CN.reset();
            continue;
        }

        int Nfrag = CN.getNumberOfProducts();

        // fragments produced in decay
        CNucleus *products = CN.getProducts(Nfrag-1);

        // the weight will be unity unless setWeightIMF is called
        float weight = products->getWeightFactor();

        const bool isSymmetricFission = CN.isSymmetricFission();
        const bool isAsymmetricFission = CN.isAsymmetricFission();
        const bool isResidue = CN.isResidue();

        if (isSymmetricFission)
        {
            sym_total++;
            Nfission += weight;  //fission event
            products = CN.getProducts(0);  // go to first evaporated particle
            for (int i=0; i<Nfrag-1; i++)
            {
                if (products->iZ == 0 && products->iA == 1) // look for neutrons
                {
                    if (products->isSaddleToScission()) neutSaddleToScission += weight;
                    else if (products->origin == 0) neutPreSad += weight;
                    else if (products->origin == 2) neutLight += weight;
                    else if (products->origin == 3) neutHeavy += weight;
                }
                // go to next particle
                products = CN.getProducts();
            }
        }

        //intermediate mass fragment
        if (isAsymmetricFission)
        {
            Nimf += weight;
            imf_total++;
            const int zMaxEvap = CN.getZmaxEvap();

            CNucleus *imfProducts = CN.getProducts(0);
            for (int j = 0; j < Nfrag && imfProducts; j++)
            {
                if (imfProducts->iZ > zMaxEvap)
                {
                    imfYieldByZN[std::make_pair(imfProducts->iZ,
                                                imfProducts->iA - imfProducts->iZ)] += weight;
                }
                imfProducts = CN.getProducts();
            }
        }

        total += weight;

        if (isResidue)
        {
            const std::string productName = products->getName();
            const int residueIndex = hash(productName);
            Residual &residue = resid[residueIndex];

            residue.count++;
            countResidue++;
            residue.name = productName;
            residue.Z = products->iZ;
            residue.A = products->iA;
            Ares += products->iA;
            Zres += products->iZ;
            resTotal += weight;
            gammaEnergy += weight*products->getSumGammaEnergy();
        }
        //============================================================= analysis BEGIN
        products = CN.getProducts(0);  // go to first evaporated particle

        for (int i=0;i<Nfrag-1;i++)
        {
            if (products->iZ == 0 && products->iA == 1 && isResidue)  //neutrons
            {
                neutMultEv += weight;
            }
            else if (products->iZ == 1 && isResidue) //protons
            {
                if(products->iA == 1 )
                {
                    protMultEv += weight;
                }
            }
            else if (products->iZ == 2  && isResidue)//alpha particles
            {
                if( products->iA == 4)
                {
                    alpMultEv += weight;
                }
            }
            // go to next particle
            products = CN.getProducts();
        }
        //============================================================= analysis END

        //reset the compound nucleus for a new decay
        CN.reset();
        counter++;
    }

    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade END
    int _INPUT = 2;
    float Sconst = _SIGMA/float(qMax(1, counter));

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

    results += buildMergedYieldTableHtmlFromImfRows("Yields of Residual Nuclei and IMF Particles",
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
