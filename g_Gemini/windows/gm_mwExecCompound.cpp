#include "gm_mainwindow.h"
#include "ui_gm_mainwindow.h"

#include <QFile>
#include <QDir>
#include <QDebug>
#include <QMessageBox>
#include <QProgressDialog>
#include <QFileDialog>
#include <QGuiApplication>
#include <QScreen>
#include <array>
#include <map>

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
    float SIG = 0;
    float SIG_ER = 0;


    CN.setEvapMode(_optEvap);

    float total = 0.;
    float Nfission = 0.;
    float Nimf = 0.;
    float neutPreSad = 0.;
    float hell = 0.;
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
    float H2MultEv = 0.;
    float H3MultEv = 0.;
    float He3MultEv = 0.;
    float gammaEnergy = 0.;
    float sum3 = 0.;
    float sum4 = 0.;
    float sum5 = 0.;
    float sum6 = 0.;

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
    const bool showAngDist = _showangdist;
    const bool showAngDistimf = _showangdistimf;

    auto computeLabKinematics = [&](CNucleus *particle,
                                    float &keLab,
                                    float &thetaLabDeg,
                                    float &vzLab,
                                    float &vxy)
    {
        float *vel = particle->getVelocityVector();

        const float vx = vel[0] / 30.f;
        const float vy = vel[1] / 30.f;
        vzLab = vel[2] / 30.f;
        vxy = std::sqrt(vx * vx + vy * vy);

        const float betaTot = std::sqrt(vx * vx + vy * vy + vzLab * vzLab);
        thetaLabDeg = 0.f;
        if (betaTot > 0.f)
        {
            float c = vzLab / betaTot;
            if (c > 1.f) c = 1.f;
            if (c < -1.f) c = -1.f;
            thetaLabDeg = static_cast<float>(std::acos(c) * 57.29577951308232);
        }

        keLab = 0.5f * (particle->iA * 931.49432f) *
                (vx * vx + vy * vy + vzLab * vzLab);
    };

    auto computeCMKineticEnergy = [](CNucleus *particle)
    {
        float *vel = particle->getVelocityVector();

        const float vx = vel[0] / 30.f;
        const float vy = vel[1] / 30.f;
        const float vz = vel[2] / 30.f;

        return 0.5f * (particle->iA * 931.49432f) *
               (vx * vx + vy * vy + vz * vz);
    };
    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade BEGIN
    for (int i=0;i<num_casc;i++)
    {
        if(progress.wasCanceled()){break;}
        if (i % 10 == 0 || i == num_casc - 1) {
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


        if (CN.isSymmetricFission())
        {
            sym_total++;
            Nfission += weight;  //fission event
            products = CN.getProducts(0);  // go to first evaporated particle
            for (int i=0; i<Nfrag-1; i++)
            {
                if (products->iZ == 0 && products->iA == 1) // look for neutrons
                {
                    if (products->isSaddleToScission()) neutSaddleToScission += weight;
                    else if (products->origin == 1) hell += weight; // cout << "hell " << std::endl;
                    else if (products->origin == 0) neutPreSad += weight;
                    else if (products->origin == 2) neutLight += weight;
                    else if (products->origin == 3) neutHeavy += weight;
                }
                // go to next particle
                products = CN.getProducts();
            }
        }

        //intermediate mass fragment
        if (CN.isAsymmetricFission())
        {
            Nimf += weight;
            imf_total++;

            CNucleus *imfProducts = CN.getProducts(0);
            for (int j = 0; j < Nfrag && imfProducts; j++)
            {
                if (imfProducts->iZ > CN.getZmaxEvap())
                {
                    float keLab = 0.f;
                    float thetaLabDeg = 0.f;
                    float vzLab = 0.f;
                    float vxy = 0.f;
                    computeLabKinematics(imfProducts, keLab, thetaLabDeg, vzLab, vxy);
                    addAngularSample(imfAngularByZN,
                                     imfProducts->iZ,
                                     imfProducts->iA - imfProducts->iZ,
                                     keLab,
                                     thetaLabDeg,
                                     vzLab,
                                     vxy,
                                     computeCMKineticEnergy(imfProducts));
                }
                imfProducts = CN.getProducts();
            }
        }

        total += weight;

        if (CN.isResidue())
        {
            resid[hash(products->getName())].count++;
            countResidue++;
            resid[hash(products->getName())].name = products->getName();
            resid[hash(products->getName())].Z = products->iZ;
            resid[hash(products->getName())].A = products->iA;
            Ares += products->iA;
            Zres += products->iZ;
            resTotal += weight;
            gammaEnergy += weight*products->getSumGammaEnergy();

            float keLab = 0.f;
            float thetaLabDeg = 0.f;
            float vzLab = 0.f;
            float vxy = 0.f;
            computeLabKinematics(products, keLab, thetaLabDeg, vzLab, vxy);
            addAngularSample(angularDistByZN,
                             products->iZ,
                             products->iA - products->iZ,
                             keLab,
                             thetaLabDeg,
                             vzLab,
                             vxy);
        }
        //============================================================= analysis BEGIN
        products = CN.getProducts(0);  // go to first evaporated particle

        for (int i=0;i<Nfrag-1;i++)
        {
            if (products->iZ == 0 && products->iA == 1 && CN.isResidue())  //neutrons
            {
                neutMultEv += weight;
                float keLab = 0.f;
                float thetaLabDeg = 0.f;
                float vzLab = 0.f;
                float vxy = 0.f;
                computeLabKinematics(products, keLab, thetaLabDeg, vzLab, vxy);
                addAngularSample(neutronAngular, keLab, thetaLabDeg, vzLab, vxy,
                                 computeCMKineticEnergy(products));
            }
            else if (products->iZ == 1 && CN.isResidue()) //protons
            {
                if(products->iA == 1 )
                {
                    protMultEv += weight;
                    float keLab = 0.f;
                    float thetaLabDeg = 0.f;
                    float vzLab = 0.f;
                    float vxy = 0.f;
                    computeLabKinematics(products, keLab, thetaLabDeg, vzLab, vxy);
                    addAngularSample(protonAngular, keLab, thetaLabDeg, vzLab, vxy,
                                     computeCMKineticEnergy(products));
                }
                else if(products->iA == 2 ) H2MultEv   += weight;
                else if(products->iA == 3 ) H3MultEv   += weight;
            }
            else if (products->iZ == 2  && CN.isResidue())//alpha particles
            {
                if( products->iA == 4)
                {
                    alpMultEv += weight;
                    float keLab = 0.f;
                    float thetaLabDeg = 0.f;
                    float vzLab = 0.f;
                    float vxy = 0.f;
                    computeLabKinematics(products, keLab, thetaLabDeg, vzLab, vxy);
                    addAngularSample(alphaAngular, keLab, thetaLabDeg, vzLab, vxy,
                                     computeCMKineticEnergy(products));
                }
                else if( products->iA == 3)  He3MultEv += weight;
            }

            else if (products->iZ == 3) sum3 += weight;
            else if (products->iZ == 4) sum4 += weight;
            else if (products->iZ == 5) sum5 += weight;
            else if (products->iZ == 6) sum6 += weight;
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

    results += "<h3 align=\"center\" style=\"color: blue\"> Yields of Decay Products </h3>";
    results += "<table cellpadding=\"5\" align=\"center\">";
    results += "<tr style=\"color: green\"> <th>Z</th><th>Name</th><th>Events</th><th>Percent</th><th>x-section (mb)</th><th> err(mb)</th></tr>";

    sort(resid, resid+length, Residual_compare);

    for(int i=0;i<length;i++)
    {
        if(resid[i].count != 0)
        {
            SIG     = _SIGMA *      (float)resid[i].count  / countResidue;
            SIG_ER =  SIG /  sqrt((float)resid[i].count);
            results += "<tr>"
                       "<td>" + QString::number(resid[i].Z) + "</td>"
                                                       "<td style=\"font-weight:bold\">" + QString::fromStdString(resid[i].name) + "</td>"
                                                                 "<td>" + QString::number(resid[i].count) + "</td>"
                                                           "<td>" + QString::number(100*(float)resid[i].count/countResidue,'f',1) +"%</td>"
                                                                                                "<td>" + QString::number(SIG   ,'g',4) + "</td>"
                                                        "<td>" + QString::number(SIG_ER,'g',4) + "</td>"
                                                           "</tr>";

            if(file_cs) fprintf(file_cs,"\n%d %d %10.3g %10.3g",resid[i].Z,resid[i].A-resid[i].Z,SIG,SIG_ER);
        }
    }

    results += "<tr><td></td><td style=\"font-weight:bold; color:green;\">Total</td><td>"+ QString::number(countResidue) + "</td>"
                                                                                                                            "<td> </td>" +           "<td style=\"font-weight:bold; color:green;\">"+ QString::number(_SIGMA,'f',2) +  "</td><td></td></tr>"
                                                                                                     "</table>";

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

        results += " <br><center><b>IMF yields: </b><br><table cellpadding=\"5\" align=\"center\">";
        results += "<tr><td> xsection of Z=3</td><td>" + QString::number(sum3*Sconst, 'g', 3) + " mb </td></tr>";
        results += "<tr><td> xsection of Z=4</td><td>" + QString::number(sum4*Sconst, 'g', 3) + " mb </td></tr>";
        results += "<tr><td> xsection of Z=5</td><td>" + QString::number(sum5*Sconst, 'g', 3) + " mb </td></tr>";
        results += "<tr><td> xsection of Z=6</td><td>" + QString::number(sum6*Sconst, 'g', 3) + " mb </td></tr></table></center>";
    }
    //------------------------------------------------------------------------
    printGeminiProperties(results);
    //------------------------------------------------------------------------
    qApp->processEvents();
    Result_Widget *ress = new Result_Widget(results);
    ress->show();

    if (showAngDist)
    {
        const QString angularHtml = buildAngularDistributionHtmlPACEStyle(
            angularDistByZN,
            _SIGMA,
            qMax(1, counter),
            _LowLimit,
            _HighLimit,
            fEx,
            iACN,
            0.0,
            "Compound mode",
            0,
            2,
            neutronAngular,
            protonAngular,
            alphaAngular
            );

        AngularDistributionWidget *angularWindow =
            new AngularDistributionWidget(
                angularHtml,
                angularDistByZN,
                _SIGMA,
                qMax(1, counter),
                _LowLimit,
                _HighLimit,
                "Compound mode",
                neutronAngular,
                protonAngular,
                alphaAngular,
                AngularDistEntry(),
                fEx,
                iACN,
                iZCN,
                0.0,
                0,
                this
                );
        QScreen *screen = QGuiApplication::screenAt(ress->geometry().center());
        if (!screen) screen = QGuiApplication::primaryScreen();
        if (screen)
        {
            const QRect available = screen->availableGeometry();
            const int gap = 12;
            const int pairWidth = available.width() - gap;
            const int resultWidth = pairWidth / 2;
            const int angularWidth = pairWidth - resultWidth;
            const int pairHeight = qMin(760, available.height());
            const int top = available.top() + (available.height() - pairHeight) / 2;

            ress->setMinimumSize(qMin(640, resultWidth), qMin(560, pairHeight));
            angularWindow->setMinimumSize(qMin(640, angularWidth), qMin(560, pairHeight));
            ress->setGeometry(available.left(), top, resultWidth, pairHeight);
            angularWindow->setGeometry(available.left() + resultWidth + gap, top, angularWidth, pairHeight);
        }

        angularWindow->show();
    }

    if (showAngDistimf)
    {
        const QString imfAngularHtml = buildAngularDistributionHtmlPACEStyle(
            imfAngularByZN,
            _SIGMA,
            qMax(1, counter),
            _LowLimit,
            _HighLimit,
            fEx,
            iACN,
            0.0,
            "Compound mode - IMF",
            0,
            2
            );

        AngularDistributionWidget *imfAngularWindow =
            new AngularDistributionWidget(
                imfAngularHtml,
                imfAngularByZN,
                _SIGMA,
                qMax(1, counter),
                _LowLimit,
                _HighLimit,
                "Compound mode - IMF",
                AngularDistEntry(),
                AngularDistEntry(),
                AngularDistEntry(),
                AngularDistEntry(),
                fEx,
                iACN,
                iZCN,
                0.0,
                0,
                this
                );
        imfAngularWindow->show();
    }

    CN.reset();
}
