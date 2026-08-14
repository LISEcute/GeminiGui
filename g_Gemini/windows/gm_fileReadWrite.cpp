#include "gm_mainwindow.h"
#include "ui_gm_mainwindow.h"
#include <QApplication>
#include <QDir>
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <iostream>
#include <QStandardPaths>
#include <QSettings>

const char *dir_files="/files";
const char *FileNameAbsent = "/gUntitled";
const char *LISEini="/lisepp.ini";

QString FFileNameCS;
QString FFileNameHtml;
QString LISErootPATH;
QString MyDocCompPATH;
QString localPATH;
QString basePATH;
//QString basePATH2;

//OT DBF *s_DBF=nullptr;

extern bool _useAME;
extern bool _useIMF;
extern bool _useIMFenh;
extern int  _optEvap;

extern int fontsizeGlobal;
extern int useHighDpiScaling;
FILE *mfopen(const QString& filename, const char* operand);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void getDPIscaling(void)
{
      MyDocCompPATH = QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).constFirst();
      QString FN1 = MyDocCompPATH + "/LISEcute";  FN1 += LISEini;
      QSettings myLiseIni1(FN1,QSettings::IniFormat);
      useHighDpiScaling  = myLiseIni1.value("font/scaling",  useHighDpiScaling).toInt();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void getInitialDir(void)
{
//---------------------------------------------------------  paths begin

    MyDocCompPATH = QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).constFirst();

    QString FileCheck(LISErootPATH); FileCheck += LISEini;
    FILE *fcheck=mfopen(FileCheck, "at");
    int work_in_LISEroot_main = 0;

    if(fcheck) {                    // work in root directory
          fclose(fcheck);
          QSettings myLiseIni0(FileCheck,QSettings::IniFormat);
          work_in_LISEroot_main = myLiseIni0.value("Version/WorkInROOT",0).toInt();
          if(work_in_LISEroot_main) localPATH = LISErootPATH;
          }

    if(work_in_LISEroot_main==0)
          {
          localPATH = MyDocCompPATH;
          localPATH += "/LISEcute";
          }

 basePATH = localPATH + "/lisecfg/";

 //--------------------------------------------------------- lise.ini  begin
   QString FN1=localPATH;  FN1 += LISEini;
   QSettings myLiseIni1(FN1,QSettings::IniFormat);
   fontsizeGlobal     = myLiseIni1.value("font/size",     fontsizeGlobal).toInt();
   useHighDpiScaling  = myLiseIni1.value("font/scaling",  useHighDpiScaling).toInt();
  //--------------------------------------------------------- lise.ini  end

 localPATH += dir_files;
 QDir pathDir(localPATH);
 if(!pathDir.exists()) pathDir.mkdir(localPATH);

//---------------------------------------------------------  paths end

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::setFileName(const QString& FileName)
{
    FFileName = FileName;

    QFileInfo fI(FFileName);
    QString windowName = fI.baseName();
    if(!windowName.contains(&FileNameAbsent[1]))  this->setWindowTitle("Gemini - " + windowName);
    else                                          this->setWindowTitle("Gemini");

    int posDot   =  FFileName.lastIndexOf('.');
    int posSlash =  FFileName.lastIndexOf('/');
    QString base;

//    if(posDot > 0 && posSlash < posDot) base = FFileName.split(".",Qt::SkipEmptyParts).at(0);
    if(posDot > 0 && posSlash < posDot) base = FFileName.left(posDot);
    else                                base = FFileName;


    FFileNameHtml = base + ".html";
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_action_Open_triggered()
{
QString file = QFileDialog::getOpenFileName(this,tr("Open File"),
                                                FFileName,
                                                tr("GEMINI files (*.gemini)"));
if(file.size()<=0) return;

setFileName(file);
readFile(file);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::readFile(const QString& file)
{
    FILE *f=mfopen(file, "rt");

    if(!f){
        QMessageBox MB;
        MB.setWindowTitle("Warning!");
        MB.setText(tr("Could not open file.\n%1").arg(file));
        MB.setIcon(QMessageBox::Warning);
        MB.exec();
        return;
        }
    else fclose(f);

    QSettings myIni(file,QSettings::IniFormat);


    myIni.beginGroup("Common");
         tab_mode   = myIni.value("mode",1).toInt();
        _useAME     = myIni.value("AME",true).toBool();
        _useIMF     = myIni.value("IMF",true).toBool();
        _useIMFenh  = myIni.value("IMFenh",false).toBool();
        _optEvap    =  myIni.value("optEvap",1).toInt();
    myIni.endGroup();

    myIni.beginGroup("Projectile");
        Ap    = myIni.value("A",16).toInt();
        Zp    = myIni.value("Z",8).toInt();
        Elab  = myIni.value("E",200).toDouble();
        num_events = myIni.value("num_events",50).toInt();
    myIni.endGroup();

    myIni.beginGroup("Target");
        At    = myIni.value("A",9).toInt();
        Zt    = myIni.value("Z",4).toInt();
    myIni.endGroup();

    myIni.beginGroup("Spin");
        spinOption  = myIni.value("S",1).toInt();
        l0          = myIni.value("l0",23).toInt();
        dif         = myIni.value("dif",2).toDouble();
    myIni.endGroup();

    myIni.beginGroup("Compound");
        iACN    = myIni.value("A",198).toInt();
        iZCN    = myIni.value("Z",82).toInt();
        fEx      = myIni.value("Ex",67).toDouble();
        fJ       = myIni.value("J",40).toDouble();
        num_casc = myIni.value("num_casc",300).toInt();
    myIni.endGroup();

writePage();
tabWidget->setCurrentIndex(tab_mode);
//qDebug() << "mode RW" << mode;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionSave_As_triggered()
{
QString file = QFileDialog::getSaveFileName(this,tr("Save File"),
                                                FFileName,tr("GEMINI files (*.gemini)"));

if(file.size()<=0) return;
FFileName = file;
on_action_Save_triggered();
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_action_Save_triggered()
{
QFile f(FFileName);

if(!f.open(QIODevice::ReadWrite | QIODevice::Text)){
        QMessageBox MB;
        MB.setWindowTitle("Warning!");
        QString txt = tr("Could not Save file.\n%1").arg(FFileName);
        MB.setText(txt);
        MB.setIcon(QMessageBox::Warning);
        MB.exec();
        return;
        }
    else f.close();

    setFileName(FFileName);
    readPage();

    QSettings myIni(FFileName,QSettings::IniFormat);

    myIni.beginGroup("Common");
        myIni.setValue("mode",tab_mode);
        myIni.setValue("AME",_useAME);
        myIni.setValue("IMF",_useIMF);
        myIni.setValue("IMFenh",_useIMFenh);
        myIni.setValue("optEvap",_optEvap);
    myIni.endGroup();

    myIni.beginGroup("Projectile");
        myIni.setValue("A",Ap);
        myIni.setValue("Z",Zp);
        myIni.setValue("E",(double)Elab);
        myIni.setValue("num_events",num_events);
    myIni.endGroup();

    myIni.beginGroup("Target");
        myIni.setValue("A",At);
        myIni.setValue("Z",Zt);
    myIni.endGroup();

    myIni.beginGroup("Spin");
        myIni.setValue("S",spinOption);
        myIni.setValue("l0",(double)l0);
        myIni.setValue("dif",(double)dif);
    myIni.endGroup();

    myIni.beginGroup("Compound");
        myIni.setValue("A",iACN);
        myIni.setValue("Z",iZCN);
        myIni.setValue("Ex",(double)fEx);
        myIni.setValue("J" ,(double)fJ);
        myIni.setValue("num_casc",num_casc);
    myIni.endGroup();

writePage();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
FILE *mfopen(const QString& filename, const char* operand)
{
#if !defined(__CYGWIN__) && !defined(_WIN32) && !defined(_WIN64)

  FILE *f =   fopen(filename.toStdString().c_str(),operand);

#else
  QString woperand(operand);
  FILE *f = _wfopen(filename.toStdWString().c_str(), woperand.toStdWString().c_str() );
#endif

return f;
}


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
