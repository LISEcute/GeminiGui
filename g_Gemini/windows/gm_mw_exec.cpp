#include "gm_mainwindow.h"
#include "ui_gm_mainwindow.h"

#include <QFile>
#include <QDir>
#include <QDebug>
#include <QMessageBox>
#include <QProgressDialog>
#include <QFileDialog>

#include "g_Gemini/source/CNucleus.h"
#include "g_Gemini/source/CFus.h"
#include "gm_results.h"
#include "gm_about.h"

// this is an example of using GEMINI CNucleus class to give the
//statistical decay of a compound nucleus

bool Residual_compare(Residual a, Residual b) {return a.Z > b.Z;};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::execute_compound()
{
    QString filename_cs =  FFileName.split(".",Qt::SkipEmptyParts).at(0) + ".cs4";

    float _SIGMA = 100;
    FILE *file_cs;
    file_cs=fopen(filename_cs.toStdString().c_str(),"wt");
    int _INPUT = 2;


    QProgressDialog progress("Calculating..","Cancel",0,num_casc);
    progress.setWindowFlags(Qt::CustomizeWindowHint |
                           Qt::WindowCloseButtonHint);

    progress.setWindowIcon(QIcon(":/pace.png"));
    progress.setWindowModality(Qt::WindowModal);
    progress.show();

    results =  "<h1 align=center style=\"color:red;\">Gemini</h1>";
    results += "<h3 align=center style=\"color:green;\">Statistical Decay Code</h3>";

    iZCN = ui->edit_ZCN->text().toInt();
    iACN = ui->edit_ACN->text().toInt();
    CNucleus CN(iZCN,iACN); //constructor
    fEx = ui->edit_Ex->text().toFloat();
    fJ =  ui->edit_J->text().toFloat();
    CN.setCompoundNucleus(fEx,fJ); //specify the excitation energy and spin
    CN.setVelocityCartesian(); // set initial CN velocity to zero
    CAngle spin(CNucleus::pi/2,(float)0.);
    CN.setSpinAxis(spin); //set the direction of the CN spin vector

    if(file_cs) {
       fprintf(file_cs,"!-title- Zcomp  Ncomp  Energy  CSfus  Input _AJNUC\n");
       fprintf(file_cs,"!Gemini %d %d %.3f %.3g %d %.3g",iZCN,iACN-iZCN,fEx,_SIGMA, _INPUT, (_INPUT==1 ? 0: fJ));
        }

    results += "<h3 style=\"color:blue;\">Starting Conditions</h3>";
    results += "<table cellpadding=5><tr><th></th><th>Z</th><th> N</th><th>A</th></tr>";
    results += "<tr><td><em>Compound nucleus</em></td><td align=center>" + QString::number(CN.iZ) + " </td><td align=center>" + QString::number(CN.iA-CN.iZ) + " </td><td align=center> " + QString::number(CN.iA) +  "</td></tr></table><p>&nbsp;</p>";


    results += "<table cellpadding=5><tr><td><em> Compound nucleus Excitation energy (MeV)</em></td><td>  " + QString::number(fEx,'f',2) + " </td></tr>";
    results += "<tr><td><em> Compound nucleus spin</em></td><td>  " + QString::number(fJ,'f',2) + " </td></tr>"
            +"</table><p>&nbsp;</p>";

    length = 101; // for hash;
    Residual resid[101];

    //int residualCount[length];
   // string residualName[length];
   // int residualZ[length];
   // int residualN[length];
    for(int i=0;i<length;i++)
        {
        resid[i].count = 0;}
        int countResidue = 0;
        float SIG = 0;
        float SIG_ER = 0;
        CN.setEvapMode(1);

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

//        float Ares = 0.;
//        float Zres = 0.;
        float resTotal = 0.;
        float neutMultEv = 0.;
        float protMultEv = 0.;
        float alpMultEv = 0.;
        float gammaEnergy = 0.;
        float sum3 = 0.;
        float sum4 = 0.;
        float sum5 = 0.;
        float sum6 = 0.;

        //float Sconst = xfus/(float)num_casc;
        //qDebug() << "num_casc = " << num_casc;
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade BEGIN
        for (int i=0;i<num_casc;i++)
            {
            //qDebug() << i ;
            if(progress.wasCanceled()){break;}
            progress.setValue(i);
            qApp->processEvents();
            //results += "<h4>event " + QString::number(i) + ":</h4>";
            //ui->textBrowser->setText(results);//.toHtml());
            CN.setCompoundNucleus(fEx,fJ);
            CN.setWeightIMF();// turn on enhanced IMF emission
            CN.decay(); //decay the compound nucleus

            if (CN.abortEvent)
                {
                cout << "abort event\n";
                CN.reset();
                continue;
                }

            int Nfrag = CN.getNumberOfProducts();
            //qDebug() << "NFrag = " << Nfrag;
            // print of number of stable

           // results += "<p>number of products = " + QString::number(CN.getNumberOfProducts()) + "</p>";
           // ui->textBrowser->setText(results);
            // fragments produced in decay
            CNucleus * products = CN.getProducts(Nfrag-1);
            //products->useAME = useAMEmass;
            float weight = products->getWeightFactor();

        //CNucleus * products = CN.getProducts(0); //set pointer to first
        //stable product
        /***********************************************************/
//        for(;;)
 //      {
//            CNucleus * parent;
//            parent = products->getParent();
//            if (products->isResidue())
//            {
//           // if(products->getName() != "n" && products->getName() != "p"){
//                countResidue++;

//               // residualCount[hash(products->getName())]++;
//                resid[hash(products->getName())].count++;
//                //residualName[hash(products->getName())] = products->getName();
//                 resid[hash(products->getName())].name = products->getName();
//               // residualZ[hash(products->getName())] = products->iZ;
//                resid[hash(products->getName())].Z = products->iZ;
//                 //residualN[hash(products->getName())] = products->iN;
//                 resid[hash(products->getName())].A = products->iA;
//          // }
//            }
//            if (parent == NULL) {
//                results += "<p>stable fragment = "+ QString::fromStdString(products->getName()) + "</p>";
//              // cout << "stable fragment = " << products->getName() << "\n";
//               products->getName();
//              //  ui->textBrowser->setText(results);

//           }else {
//                 //parent->print();
//              // results += "<p>stable fragment= " + QString::fromStdString(products->getName()) +        " parent = " + QString::fromStdString(parent->getName()) + "</p>";
//                  //cout <<  "<p>stable fragment= " << products->getName() << " parent = " << parent->getName() << "\n";
//                  products->getName();
//                  parent->getName();
//               // results += parent->printHtml();
//              //  results += products->printHtml();
//                // ui->textBrowser->insertHtml("<p>stable fragment= " + QString::fromStdString(products->getName()) +        " parent = " + QString::fromStdString(parent->getName()) + "</p><br/>");
//            }
/*******************************************************/
        //for(;;) {
            if (CN.isSymmetricFission())
                {
                //qDebug() << "is Symmetric fission!";
                sym_total++;
                Nfission += weight;  //fission event
                products = CN.getProducts(0);  // go to first evaporated particle
                for (int i=0;i<Nfrag-1;i++)
                    {
                    //  results += "<p> symm Z = " + QString::number(products->iZ) + " A = " + QString::number(products->iA) + "</p>";
                    //  ui->textBrowser->setText(results);
                    if (products->iZ == 0 && products->iA == 1) // look for neutrons
                        {
                        if (products->isSaddleToScission())
                            neutSaddleToScission += weight;
                        else if (products->origin == 1) hell += weight; // cout << "hell " << Qt::endl;
                        else if (products->origin == 0) neutPreSad += weight;
                        else if (products->origin == 2) neutLight += weight;
                        else if (products->origin == 3) neutHeavy += weight;
                        }
                    // go to next particle
                    products = CN.getProducts();
                  }
                }
            if (CN.isAsymmetricFission()) {
                Nimf += weight; imf_total++;
                //qDebug() << "is Asymmetric fission!";
                }
            total += weight;
            if (CN.isResidue())
                {
                resid[hash(products->getName())].count++;
                countResidue++;
                resid[hash(products->getName())].name = products->getName();
                resid[hash(products->getName())].Z = products->iZ;
                resid[hash(products->getName())].A = products->iA;
                //Ares += products->iA;
                //Zres += products->iZ;
                //gammaEnergy += weight*products->getSumGammaEnergy();
            }

        products = CN.getProducts();  // go to the next product
           //if (products == NULL) break;
        CN.reset();
        }
   //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade END

    results += "<h3 align=\"center\" style=\"color: blue\"> Decay Product Summary </h3>";
    results += "<table cellpadding=\"5\" align=\"center\">";

    results += "<tr><th style=\"color:green;\">Result</th><th style=\"color:green;\">Number</th></tr>";// + QString::number(imf_total) + "</td></tr>";
    results += "<tr><th>Intermediate Mass Fragments</th><td>" + QString::number(imf_total) + "</td></tr>";
    results += "<tr><th>Symmetric Fission</th><td>" + QString::number(sym_total) + "</td></tr>";
    results += "<tr><th>Residual Nuclei</th><td>" + QString::number(countResidue) + "</td></tr>";

/*    results += "<tr><th style=\"color:green;\">Result</th><th style=\"color:green;\">Number</th><th style=\"color:green;\">Cross section</th></tr>";// + QString::number(imf_total) + "</td></tr>";
    results += "<tr><th>Intermediate Mass Fragments</th><td>" + QString::number(imf_total) + "</td><td>" + QString::number(imf_total/num_casc)+ "</td></tr>";
    results += "<tr><th>Symmetric Fission</th><td>" + QString::number(sym_total) + "</td><td>" + QString::number(sym_total/num_casc)+ "</td></tr>";
    results += "<tr><th>Residual Nuclei</th><td>" + QString::number(countResidue) + "</td><td>" + QString::number((float)countResidue/(float)num_casc)+ "</td></tr>";
*/
    //float _SIGMA = resTotal*Sconst;
    results += "</table>;";

        //CN.reset();
    //}
    results += "<h3 align=\"center\" style=\"color: blue\"> Yields of Decay Products </h3>";
    results += "<table cellpadding=\"5\" align=\"center\">";
    results += "<tr style=\"color: green\"><th>Z</th><th>Name</th><th>Events</th><th>Percent</th><th> </th></tr>";

    sort(resid, resid+length, Residual_compare);

    for(int i=0;i<length;i++)
        {
        //if(residualCount[i] != 0){
            if(resid[i].count != 0){
           // SIG = 100*(float)residualCount[i]/countResidue;
                SIG = 100*(float)resid[i].count/countResidue;
            SIG_ER = sqrt(SIG);
           // results += "<tr><td>" + QString::number(residualZ[i])+ "</td><td style=\"font-weight:bold\">" + QString::fromStdString(residualName[i]) + "</td><td>" + QString::number(residualCount[i]) + "</td><td>" +
             //       QString::number(100*(float)residualCount[i]/countResidue,'f',1) +"%</td><td>" +QString::number(100*(float)residualCount[i]/countResidue,'f',1) + "</td></tr>";
            results += "<tr><td>" + QString::number(resid[i].Z)+ "</td><td style=\"font-weight:bold\">" + QString::fromStdString(resid[i].name) + "</td><td>" + QString::number(resid[i].count) + "</td><td>" +
                              QString::number(100*(float)resid[i].count/countResidue,'f',1) +"%</td><td>" +  "</td></tr>";

            /** add code for .cs4 here OR move file writing to cs_file.cpp or other file (better option filereadwrite.cpp) **/
            //fprintf(file_cs,"\n%d %d %10.3g %10.3g",residualZ[i],residualN[i],SIG,SIG_ER);
            fprintf(file_cs,"\n%d %d %10.3g %10.3g",resid[i].Z,resid[i].A-resid[i].Z,SIG,SIG_ER);
        }
    }
    results += "<tr><td></td><td colspan=\"1\" style=\"font-weight:bold; color:green;\">Total</td><td>"+ QString::number(countResidue) + "</td><td></td></tr></table>";// + QString::number(sum*fus.plb*(float)countResidue)+ "</td></tr></table>";

 //   results += "</table>";

    fclose(file_cs);

    if (Nfission > 0.)
        {
        results += " neutron multiplicities in fission <br>";
        results += "  presaddle = " + QString::number(neutPreSad/Nfission) + "<br>";
        results +=  " saddle-to-scission  = " + QString::number(neutSaddleToScission/Nfission) + "<br>";
        results += "  post-scission light frag = " + QString::number(neutLight/Nfission) + "<br>";
        results += "  post-scission heavy frag = " + QString::number(neutHeavy/Nfission) + "<br><br>";
        // results += Qt::endl;
        }

    results += "<br>";
//Qt-Oleg    results += "<center><b style=\"color:blue\">residue (no IMF or Symmetric fission emissions)</b><br>";

   // results += " residue xsection = " + QString::number(resTotal*Sconst) + " mb <br>";
    // qDebug() << "resTotal = " << resTotal << "Ares = " << Ares;
    if (resTotal > 0.)
        {
     //   results += "  average residue is Z = " + QString::number(Zres/resTotal) + " A = "
            //    + QString::number(Ares/resTotal) + "<br>";

        results += "<table cellpadding=\"5\" align=\"center\"><tr><td>  average neutron multiplicity </td><td>" + QString::number(neutMultEv/resTotal) + "</td></tr>";
        results += " <tr><td> average proton multiplicity </td><td> " + QString::number(protMultEv/resTotal) + "</td></tr>"; //endl;
        results += " <tr><td> average alpha multiplicity </td><td> " + QString::number(alpMultEv/resTotal ) + "</td></tr>";
        results += " <tr><td> average energy in gamma rays </td><td> " + QString::number(gammaEnergy/resTotal) +
                " MeV </td></tr></table></center>";

        results += "<center><br> <b style=\"color:blue\">Intermediate Mass Fragment (Z > " + QString::number(CN.getZmaxEvap()) + ") (IMF) </b><br>" ;//+ Qt::endl;
        results += " <table cellpadding=\"5\" align=\"center\"><tr><td> IMF prob = " + QString::number(Nimf/total)
                + "</td><td> xsection = " + QString::number(0) + " mb </td></tr></table><br>";// <<endl;
        results += " <p>&nbsp;</p><b>IMF yields: </b><br><table cellpadding=\"5\" align=\"center\">";
        results += "<tr><td> xsection of Z=3</td><td>  " + QString::number(sum3/num_casc) + " mb </td></tr>";
        results += "<tr><td> xsection of Z=4</td><td>  " + QString::number(sum4/num_casc) + " mb</td></tr>";
        results += "<tr><td> xsection of Z=5</td><td>  " + QString::number(sum5/num_casc) + " mb </td></tr>";
        results += "<tr><td> xsection of Z=6</td><td>  " + QString::number(sum6/num_casc) + " mb </td></tr></table></center>";
        }

    qApp->processEvents();
    Result_Widget *ress = new Result_Widget(results);
    ress->show();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::execute_fusion()
{

    results = "<h1 align=center style=\"color:red;\">Gemini</h1>";
    results += "<h3 align=center style=\"color:green;\">Statistical Decay Code</h3>";
    Zp = ui->edit_ZP->text().toInt();
    Zt = ui->edit_ZT->text().toInt();
    Ap = ui->edit_AP->text().toInt();
    At = ui->edit_AT->text().toInt();
    Elab = ui->edit_ExP->text().toFloat();
    dif = ui->edit_dif->text().toFloat();
    CFus fus(Zp,Ap,Zt,At,Elab,dif);

    if(ui->radioButton->isChecked()){
        l0 = ui->edit_max_spin->text().toFloat(); //get maximum spin from paper?
        }
    else {
        l0 = fus.getBassL();  // get maximum spin from Bass model
        }
    //construct fusion spin distribution

    int lmax = (int)l0 +5;  // maxium spin considered
    // const int llmax = lmax;
    //float prob[lmax+1];

    float prob[101];
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

    CNucleus CNN(fus.iZcn,fus.iAcn); //constructor
    float ffEx = fus.Ex;             //excitation energy of compound nucleus

    results += "<h3 style=\"color:blue;\">Starting Conditions</h3>";
    results += "<table cellpadding=5><tr><th></th><th>Z</th><th> N</th><th>A</th><th style=\"font-weight:bold\"><sup>A </sup>El</th></tr>";

    results += "<tr><td><em>Projectile</em></td><td alnign=center>" + QString::number(Zp) +
               " </td><td align=center>" + QString::number(Ap-Zp) + " </td><td align=center>" + QString::number(Ap) +
               " </td><td align=center>" + fus.p.getGName()  + "</td></tr>";

    results += "<tr><td><em>Target</em></td><td align=center>" + QString::number(Zt) +
               "</td><td align=center>" + QString::number(At-Zt) + "</td><td align=center>" + QString::number(At) +
            " </td><td align=center>" + fus.t.getGName()  + "</td></tr>";

    results += "<tr><td><em>Compound nucleus</em></td><td align=center>" + QString::number(fus.iZcn) +
               " </td><td align=center>" + QString::number(CNN.iN) + " </td><td align=center> " + QString::number(fus.iAcn) +
            " </td><td align=center>" +  fus.c.getGName()  + "</td></tr>" +
               "</table><p>&nbsp;</p>";

    results += "<table cellpadding=5><tr><td><em> Bombarding energy (MeV)</em></td><td>  " + QString::number(Elab,'f',2) + " </td></tr> " +
            "<tr><td><em> Center of Mass energy (MeV)</em></td><td>  " + QString::number(fus.Ecm,'f',3) + " </td></tr>"
            +  "<tr><td><em> Compound nucleus Excitation energy (MeV)</em></td><td>  " + QString::number(fus.Ex,'f',2) + " </td></tr>"
            +  "<tr><td><em> Q-value of reaction (MeV)</em></td><td>  " + QString::number(fus.Qval,'f',3) + " </td></tr>"
            +  "<tr><td><em> Compound nucleus recoil energy (MeV)</em></td><td>  " + QString::number(Elab*Ap/(Ap+At),'f',3) + " </td></tr>"
            +  "<tr><td><em> Compound nucleus recoil velocity (cm/ns) </em></td><td>  " + QString::number(30.*sqrt(2*(Elab*Ap/(Ap+At))/(931.49432 *(Ap+At))),'e',3) + " </td></tr>"
            +  "<tr><td><em> Compound nucleus recoil velocity/c </em></td><td>  " + QString::number(sqrt(2*(Elab*Ap/(Ap+At))/(931.49432 *(Ap+At))),'e',3) + " </td></tr>"
            +  "<tr><td><em> Beam velocity (cm/ns) </em></td><td>  " + QString::number(30.*sqrt(2.0*Elab/(931.49432*Ap)),'e',3) + " </td></tr>"
            +  "<tr><td><em> Beam velocity/c </em></td><td>  " + QString::number(sqrt(2.0*Elab/(931.49432*Ap)),'e',3) + " </td></tr>"
            +"</table><p>&nbsp;</p>";
    results += "<table cellpadding=5><tr><td><em> diffuseness </em></td><td>  " + QString::number(fus.dif,'f',2) + " </td></tr> "
            +  "<tr><td><em> Fusion cross section (mb) </em></td><td>  " + QString::number(sum*fus.plb,'e',3) + " </td></tr>"
            +  "<tr><td><em> Bass L </em></td><td>  " + QString::number(fus.getBassL(),'f',2) + " </td></tr>"
            +  "<tr><td><em> L0 </em></td><td>  " + QString::number(fus.getL0(sum*fus.plb),'f',2) + " </td></tr>"
            +  "<tr><td><em> Bass cross section </em></td><td>  " + QString::number(fus.getBassXsec(),'f',2) + " </td></tr>"
            +"</table><p>&nbsp;</p>";

    results += " E<sub>X</sub> = " + QString::number(ffEx) + " &nbsp; Critical spin = " + QString::number(l0, 'f' ,1)  + " &#8463;<br>";
    // results += " fusion xsection = " + QString::number(sum*fus.plb) + " mb <br>";

    float xfus = sum*fus.plb;

    CNN.setEvapMode(1);  // force a fully Hauser-Feshbach calculation
   // ui->textBrowser->setText(results);

    // if IMF probablities are small and you are not interested
    // save time and turn them off.
    //CNucleus::setNoIMF();

    // don't like the default evaporation parameters
    //CLevelDensity::setAfAn(1.036);  // change ratio of saddle to ground state
    //level density parameters
    //CLevelDensity::setLittleA(8.); //change level density parameter to A/8.

    //write out statistical model parameters
    //cout << "CNN.printParameters():" << Qt::endl;
    //CNN.printParameters();
    //cout << "cone printing parameters" << Qt::endl;
    //define and zero events parameters

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
    float gammaEnergy = 0.;


    float sum3 = 0.;
    float sum4 = 0.;
    float sum5 = 0.;
    float sum6 = 0.;

    length = 101; // for hash;

    int Nevents = ui->num_Events->text().toInt();
    float Sconst = xfus/(float)Nevents;
    //int residualCount[length];

    Residual resid[101];

    for(int i=0;i<length;i++) resid[i].count = 0;

    int countResidue = 0;
    QProgressDialog progress("Calculating..","Cancel",0,Nevents);
    progress.setWindowFlags(Qt::CustomizeWindowHint |
                           Qt::WindowCloseButtonHint);
    progress.setWindowIcon(QIcon(":/pace.png"));
    progress.setWindowModality(Qt::WindowModal);
    progress.show();

    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade BEGIN
    for (int i=0;i<Nevents;i++)
        {
        if(progress.wasCanceled()){break;}
        progress.setValue(i);
        qApp->processEvents();
        //choose the spin of the CN from the determed spin distribution
        float ran = CNN.ran.Rndm();
        int l = 0;
        for (;;)
            {
            if (ran < prob[l]) break;
            l++;
            }

        //specify the excitation energy and spin
        CNN.setCompoundNucleus(ffEx,(float)l);

        //if you are interested in IMF emission at low excitation energy
        //then turn IMF weighting on
        CNN.setWeightIMF();// turn on enhanced IMF emission

        CNN.decay(); //decay the compound nucleus

        if (CNN.abortEvent)  //gemini had trouble with this event
            {
            CNN.reset();
            continue;
            }


        // number of stable fragments produced in the decay
        int Nfrag = CNN.getNumberOfProducts();

        //set pointer to last fragment which is the evaporation residue
        // if isResidue is true, otherwise the heavy fission fragment
        // this should be the heaviest fragment produced
        CNucleus * products = CNN.getProducts(Nfrag-1);
        //products->useAME = useAMEmass;
        // the weight will be unity unless setWeightIMF is called
        float weight = products->getWeightFactor();

        //cout << "weight factor = " << weight << Qt::endl;

        if (CNN.isSymmetricFission())
            {
            sym_total++;
            Nfission += weight;  //fission event
            products = CNN.getProducts(0);  // go to first evaporated particle
            for (int i=0;i<Nfrag-1;i++)
                {
               // results += "<p> symm Z = " + QString::number(products->iZ) + " A = " + QString::number(products->iA) + "</p>";
                //  ui->textBrowser->setText(results);
                if (products->iZ == 0 && products->iA == 1) // look for neutrons
                    {
                    if (products->isSaddleToScission())
                        neutSaddleToScission += weight;
                    else if (products->origin == 1) hell += weight; // cout << "hell " << Qt::endl;
                    else if (products->origin == 0) neutPreSad += weight;
                    else if (products->origin == 2) neutLight += weight;
                    else if (products->origin == 3) neutHeavy += weight;
                    }
                // go to next particle
                products = CNN.getProducts();
                }
            }

        //intermediate mass fragment
        if (CNN.isAsymmetricFission()) {Nimf += weight; imf_total++;} //imf event
        total += weight;
      //  cout << "fJ = " << products->fJ << Qt::endl;
        // CAngle ang = products->getAngleDegrees();
      //  cout << "theta (deg) = " << products->getThetaDegrees() << Qt::endl;
        // cout << "angle (deg): phi = " << ang.phi << " pi = " << ang.pi << " theta = " << ang.theta << Qt::endl;
        //cout << "plb = " << fus.plb << Qt::endl;
        //evaporation resiudes
        if (CNN.isResidue())
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
            }

        products = CNN.getProducts(0);  // go to first evaporated particle
        for (int i=0;i<Nfrag-1;i++)
            {
            // results += "<p> Z = " + QString::number(products->iZ) + " A = " + QString::number(products->iA) + "</p>";

            if (products->iZ == 0 && products->iA == 1 && CNN.isResidue())  //neutrons
                {
                neutMultEv += weight;
                }
            else if (products->iZ == 1 && products->iA == 1 && CNN.isResidue()) //protons
                {
                protMultEv += weight;
                }
            else if (products->iZ == 2 && products->iA == 4 && CNN.isResidue())//alpha particles
                {
                alpMultEv += weight;
                }

            else if (products->iZ == 3) sum3 += weight;
            else if (products->iZ == 4) sum4 += weight;
            else if (products->iZ == 5) sum5 += weight;
            else if (products->iZ == 6) sum6 += weight;
            // go to next particle
            products = CNN.getProducts();
        }

        //reset the compound nucleus for a new decay
        CNN.reset();
    }
    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  cascade BEGIN

    results += "<h3 align=\"center\" style=\"color: blue\"> Fusion Product Summary </h3>";
    results += "<table cellpadding=\"5\" align=\"center\">";
    results += "<tr><th style=\"color:green;\">Result</th><th style=\"color:green;\">Number</th></tr>";// + QString::number(imf_total) + "</td></tr>";
    results += "<tr><th>Intermediate Mass Fragments</th><td>" + QString::number(imf_total) + "</td></tr>";
    results += "<tr><th>Symmetric Fission</th><td>" + QString::number(sym_total) + "</td></tr>";
    results += "<tr><th>Residual Nuclei</th><td>" + QString::number(countResidue) + "</td></tr>";
/*
    results += "<tr><th style=\"color:green;\">Result</th><th style=\"color:green;\">Number</th><th style=\"color:green;\">Cross section</th></tr>";// + QString::number(imf_total) + "</td></tr>";
    results += "<tr><th>Intermediate Mass Fragments</th><td>" + QString::number(imf_total) + "</td><td>" + QString::number(Nimf*Sconst)+ "</td></tr>";
    results += "<tr><th>Symmetric Fission</th><td>" + QString::number(sym_total) + "</td><td>" + QString::number(Nfission*Sconst)+ "</td></tr>";
    results += "<tr><th>Residual Nuclei</th><td>" + QString::number(countResidue) + "</td><td>" + QString::number(resTotal*Sconst)+ "</td></tr>";

 */
    float _SIGMA = resTotal*Sconst;

    float SIG = 0;
    float SIG_ER = 0;
    results += "<tr><th>TOTAL</th><td><b>" + QString::number(imf_total+sym_total+countResidue) + "<b></td></tr>";
    results += "</table><p>&nbsp;</p>";
    results += "<h3 align=\"center\" style=\"color: blue\"> Yields of Residual Nuclei </h3>";
    results += "<table cellpadding=\"6\" align=\"center\">";
    results += "<tr style=\"color: green\"><th>Z</th><th>Name</th><th>Events</th><th>Percent</th><th>x-section (mb)</th><th> err(mb)</th></tr>";
   // setFileName("test");
    QString filename_cs =  FFileName.split(".",Qt::SkipEmptyParts).at(0) + ".cs4";
    FILE *file_cs;
    file_cs=fopen(filename_cs.toStdString().c_str(),"wt");
    int _INPUT = 2;
    /*
    fJ
*/
    _INPUT = 1;
    if(file_cs) {
       fprintf(file_cs,"!-title- Zcomp  Ncomp  Energy  CSfus  Input _AJNUC\n");
       fprintf(file_cs,"!Gemini %d %d %.3f %.3g %d %.3g",CNN.iZ,CNN.iN,CNN.fEx, _SIGMA, _INPUT, (_INPUT==1 ? l0: l0));
        }
    else { return; }

    sort(resid,resid+length,Residual_compare);
    for(int i=0;i<length;i++)
         {
         if(resid[i].count != 0)
            {
            //SIG = sum*fus.plb*(float)residualCount[i]/countResidue;
            SIG     = xfus *      (float)resid[i].count  / countResidue;
            SIG_ER =  SIG /  sqrt((float)resid[i].count);
            results += "<tr>"
                       "<td>" + QString::number(resid[i].Z) + "</td>"
                       "<td style=\"font-weight:bold\">" + QString::fromStdString(resid[i].name) + "</td>"
                       "<td>" + QString::number(resid[i].count) + "</td>"
                       "<td>" + QString::number(100*(float)resid[i].count/countResidue,'f',1) +"%</td>"
                       "<td>" + QString::number(SIG   ,'f',2) + "</td>"
                       "<td>" + QString::number(SIG_ER,'f',2) + "</td>"
                       "</tr>";

            fprintf(file_cs,"\n%d %d %10.3g %10.3g",resid[i].Z,resid[i].A-resid[i].Z,SIG,SIG_ER);
            }
    }

    results += "<tr><td colspan=\"2\" style=\"font-weight:bold; color:green;\">Total</td><td>"+ QString::number(countResidue) + "</td>"
               "<td> </td>" +           "<td style=\"font-weight:bold; color:green;\">"+ QString::number(xfus,'f',2) +  "</td><td></td></tr>"
            "</table>";


    fclose(file_cs);
   // results += "</table>";
   // cout << "resTotal = " << resTotal << "  countResidue = " << countResidue << Qt::endl;
    // cout << Qt::endl;
   // results += "<p>fission probability = " + QString::number(Nfission/total) + "  xsection = " +
           // QString::number(Nfission*Sconst) + " mb</p>";
   // ui->textBrowser->setText(results);

    if (Nfission > 0.)
        {
        results += " neutron multiplicities in fission <br>";
        results += "  presaddle = " + QString::number(neutPreSad/Nfission) + "<br>";
        results +=  "saddle-to-scission  = " + QString::number(neutSaddleToScission/Nfission) + "<br>";
        results += "  post-scission light frag = " + QString::number(neutLight/Nfission) + "<br>";
        results += "  post-scission heavy frag = " + QString::number(neutHeavy/Nfission) + "<br><br>";
        }

    results += "<br>";


    results += "<center><b style=\"color:blue\">residue (no IMF or Symmetric fission emissions)</b><br>";
   // results += " residue xsection = " + QString::number(resTotal*Sconst) + " mb <br>";
    // qDebug() << "resTotal = " << resTotal << "Ares = " << Ares;
    if (resTotal > 0.)
        {
     //   results += "  average residue is Z = " + QString::number(Zres/resTotal) + " A = "
            //    + QString::number(Ares/resTotal) + "<br>";

        results += "<table cellpadding=\"5\" align=\"center\"><tr><td>  average neutron multiplicity </td><td>" + QString::number(neutMultEv/resTotal) + "</td></tr>";
        results += " <tr><td> average proton multiplicity </td><td> " + QString::number(protMultEv/resTotal) + "</td></tr>"; //endl;
        results += " <tr><td> average alpha multiplicity </td><td> " + QString::number(alpMultEv/resTotal ) + "</td></tr>";
        results += " <tr><td> average energy in gamma rays </td><td> " + QString::number(gammaEnergy/resTotal) +
                " MeV </td></tr></table></center>";

        results += "<center><br> <b style=\"color:blue\">Intermediate Mass Fragment (Z > " + QString::number(CNN.getZmaxEvap()) + ") (IMF) </b><br>" ;//+ Qt::endl;
        results += " <table cellpadding=\"5\" align=\"center\"><tr><td> IMF prob = " + QString::number(Nimf/total)
                + "</td><td> xsection = " + QString::number(Nimf*Sconst) + " mb </td></tr></table><br>";// <<endl;
        results += " <p>&nbsp;</p><b>IMF yields: </b><br><table cellpadding=\"5\" align=\"center\">";
        results += "<tr><td> xsection of Z=3</td><td> " + QString::number(sum3*Sconst) + " mb </td></tr>";
        results += "<tr><td> xsection of Z=4</td><td>  " + QString::number(sum4*Sconst) + " mb</td></tr>";
        results += "<tr><td> xsection of Z=5</td><td>  " + QString::number(sum5*Sconst) + " mb </td></tr>";
        results += "<tr><td> xsection of Z=6</td><td>  " + QString::number(sum6*Sconst) + " mb </td></tr></table></center>";
        }

    qApp->processEvents();
    Result_Widget *ress = new Result_Widget(results);
    ress->show();
    CNN.reset();

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
