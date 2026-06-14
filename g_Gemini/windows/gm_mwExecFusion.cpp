#include "gm_mainwindow.h"
#include "ui_gm_mainwindow.h"

#include <QFile>
#include <QDir>
#include <QDateTime>
#include <QDebug>
#include <QMessageBox>
#include <QProgressDialog>
#include <QFileDialog>
#include <map>
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
    float ffEx = fus.Ex;             //excitation energy of compound nucleus

    if(spinOption==1) l0 = fus.getBassL();  // get maximum spin from Bass model

    //construct fusion spin distribution

    int lmax = (int)l0 +5;  // maxium spin considered
    float prob[301];
    lmax = qMin(300,lmax);
    float sum = 0.;

    for (int l=0;l<=lmax;l++)
    {
        prob[l] =  (float)(2*l+1);

        if (dif > 0.) prob[l] /= (1.+exp(((float)l-l0)/dif));
        else if ( l > l0)  prob[l] = 0.;

        sum += prob[l];
    }

    for (int l=0;l<=lmax;l++)
    {
        prob[l] /= sum;
        if (l > 0) prob[l] += prob[l-1];
    }

    float xfus = sum*fus.plb;

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

    float SIG = 0;
    float SIG_ER = 0;

    float total = 0.;
    float Nfission = 0.;
    float Nimf = 0.;
    float hell = 0.;
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
    float H2MultEv = 0.;
    float H3MultEv = 0.;
    float He3MultEv = 0.;

    float gammaEnergy = 0.;

    float sum3 = 0.;
    float sum4 = 0.;
    float sum5 = 0.;
    float sum6 = 0.;

    length = 101; // for hash;
    Residual resid[101];
    for(int i=0;i<length;i++) resid[i].count = 0;

    int countResidue = 0;
    int counter=0;
    std::map<std::pair<int, int>, AngularDistEntry> angularDistByZN;
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

    const float betaCN = std::sqrt(
        2.0f * (Elab * Ap / (Ap + At)) / (931.49432f * (Ap + At))
        );

    auto computeLabKinematics = [&](CNucleus *particle,
                                    float &keLab,
                                    float &thetaLabDeg,
                                    float &vzLab,
                                    float &vxy)
    {
        float *vel = particle->getVelocityVector();

        const float vx = vel[0] / 30.f;
        const float vy = vel[1] / 30.f;
        const float vzLocal = vel[2] / 30.f;
        vzLab = vzLocal + betaCN;
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
    for (int i=0;i<num_events;i++)
    {
        if(progress.wasCanceled()){break;}
        if (i % 10 == 0 || i == num_events - 1) {
            progress.setValue(i);
            qApp->processEvents();
        }

        //--------------------------------------------------------------
        //choose the spin of the CN from the determed spin distribution
        float ran = CN.ran->Rndm();
        int l = 0;
        for (;;)
        {
            if (ran < prob[l]) break;
            l++;
        }
        //-------------------------------------------------------------------

        //specify the excitation energy and spin
        CN.setCompoundNucleus(ffEx,(float)l);

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

        if (CN.abortEvent)  //gemini had trouble with this event
        {
            CN.reset();
            continue;
        }


        // number of stable fragments produced in the decay
        int Nfrag = CN.getNumberOfProducts();

        //set pointer to last fragment which is the evaporation residue
        // if isResidue is true, otherwise the heavy fission fragment
        // this should be the heaviest fragment produced
        CNucleus *products = CN.getProducts(Nfrag-1);

        // the weight will be unity unless setWeightIMF is called
        float weight = products->getWeightFactor();


        if (CN.isSymmetricFission())
        {
            sym_total++;
            Nfission += weight;             //fission event
            products = CN.getProducts(0);   // go to first evaporated particle
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
        if (CN.isAsymmetricFission()) {Nimf += weight; imf_total++;} //imf event

        total += weight;

        // cout << "fJ = " << products->fJ << std::endl;
        // CAngle ang = products->getAngleDegrees();
        // cout << "theta (deg) = " << products->getThetaDegrees() << std::endl;
        // cout << "angle (deg): phi = " << ang.phi << " pi = " << ang.pi << " theta = " << ang.theta << std::endl;
        // cout << "plb = " << fus.plb << std::endl;

        //evaporation resiudes

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
            const float eventGammaEnergy = products->getSumGammaEnergy();
            gammaEnergy += weight * eventGammaEnergy;

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

            // Gamma angular-distribution sample.
            // Gemini gives total gamma energy from the residue.
            // It does not give individual gamma-ray direction here,
            // so we attach the gamma energy to the residue lab angle.
            if (eventGammaEnergy > 0.f)
            {
                addAngularSample(gammaAngular,
                                 eventGammaEnergy,
                                 thetaLabDeg,
                                 0.0f,
                                 0.0f,
                                 eventGammaEnergy);
            }
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

    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade BEGIN
    int _INPUT = 1;
    float Sconst = xfus/float(qMax(1, counter));

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
    results += "<h3 align=\"center\" style=\"color: blue\"> Yields of Residual Nuclei </h3>";
    results += "<table cellpadding=\"5\" align=\"center\">";
    results += "<tr style=\"color: green\"> <th>Z</th><th>Name</th><th>Events</th><th>Percent</th><th>x-section (mb)</th><th> err(mb)</th></tr>";

    sort(resid,resid+length,Residual_compare);
    for(int i=0;i<length;i++)
    {
        if(resid[i].count != 0)
        {
            SIG     = xfus *      (float)resid[i].count  / countResidue;
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
                                                                                                                            "<td> </td>" +           "<td style=\"font-weight:bold; color:green;\">"+ QString::number(xfus,'f',2) +  "</td><td></td></tr>"
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
    extern double _LowLimit, _HighLimit;

    if (ui->DistAng->isChecked())
    {
        const double recoilBetaCN =
            std::sqrt(2.0 * (Elab * Ap / (Ap + At)) /
                      (931.49432 * (Ap + At)));

        const QString angularHtml = buildAngularDistributionHtmlPACEStyle(
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
            1,
            neutronAngular,
            protonAngular,
            alphaAngular,
            gammaAngular
            );

        AngularDistributionWidget *angularWindow =
            new AngularDistributionWidget(
                angularHtml,
                angularDistByZN,
                xfus,
                qMax(1, counter),
                _LowLimit,
                _HighLimit,
                "Fusion mode",
                neutronAngular,
                protonAngular,
                alphaAngular,
                gammaAngular,
                fus.Ex,
                fus.iAcn,
                recoilBetaCN,
                0,
                this
                );
        angularWindow->show();
    }
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
