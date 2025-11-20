#include "gm_mainwindow.h"
#include "ui_gm_mainwindow.h"

#include <QFile>
#include <QDir>
#include <QDebug>
#include <QMessageBox>
#include <QProgressDialog>
#include <QFileDialog>
#include <QDesktopServices>
#include <QUrl>
#include <QSqlDatabase>


#include "g_Gemini/source/CNucleus.h"
#include "g_Gemini/source/CFus.h"
#include "gm_results.h"
#include "gm_about.h"

// this is an example of using GEMINI CNucleus class to give the
//statistical decay of a compound nucleus


bool _useAME=false;
bool _useIMF=true;
bool _useIMFenh=false;
int  _optEvap=1;

extern QString LISErootPATH;
extern QString localPATH;
extern const char *FileNameAbsent;
extern QString FileArg;
void initchart(std::vector<std::vector<QVariant>>& Array);
std::vector<std::vector<QVariant>> dakChart;
QSqlDatabase AmeDB;
QSqlDatabase Geminidb;


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);
  ui->statusBar->hide();

  //-------------------------------------------------------- dialog init2 begin

  setWindowFlags( Qt::Window | Qt::CustomizeWindowHint |
                  Qt::WindowTitleHint | Qt::WindowSystemMenuHint |
                  Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint );

  //--------------------------------------------------------
  connect(ui->rb_spinMax,  SIGNAL(toggled(bool)), this, SLOT(checked(bool)));
  //--------------------------------------------------------
  tab_mode=0;
  tabWidget = ui->tabWidget;
  tabWidget->setCurrentIndex(tab_mode);
 // qDebug() << tabWidget << tabWidget->count() << tabWidget->currentIndex();

  num_casc = 300;
  iZCN = 82; // proton number of compound nucleus
  iACN = 198; // mass number of compound nucleus
  fEx = 67.; //excitation energy of compound nucleus
  fJ = 40; // spin of compound nucleus

  //--------------------------------------------------------

  num_events = 500;
  Zp = 8; // proton number of projectile
  Ap = 16; // mass number of projectile
  Zt = 6; // proton number of target
  At = 12; //mass number of target
  Elab = 160.; // lab energy in MeV
  dif = 2; //diffuseness of fusion spin distribution in hbar
  l0 = 23;
  spinOption = 1; // Bass
  checked(spinOption==0);

  //-------------------------------------------------------- dialog init2 end


  //------------------------------------------------------- gemini_db begin

  QString Gemini_DbName = "GEMINI.sqlite";
  QString GeminiDbPath = LISErootPATH + "/lisecfg" + "/" + Gemini_DbName;

  QString connectionName = "GeminiConnection";

  Geminidb = QSqlDatabase::addDatabase("QSQLITE","GeminiConnection");
  Geminidb.setDatabaseName(GeminiDbPath);

  if (!Geminidb.open())
      QMessageBox::warning(this, "Database file didn't open", GeminiDbPath + " can't be located or an Error occured");

  //------------------------------------------------------- ame_db begin

  QString AME_DbName = "AME_DB.sqlite"; // Name of the AME Database File
  QString AMEDbPath = LISErootPATH + "/lisecfg" + "/" + AME_DbName; // File Path of the AME Database file

  QString connectionName2 = "AMEConnection";

  AmeDB = QSqlDatabase::addDatabase("QSQLITE", connectionName2);
  AmeDB.setDatabaseName(AMEDbPath);

  if (!AmeDB.open())
      QMessageBox::warning(this, "Database file didn't open", AMEDbPath + " can't be located or an Error occured");
  //------------------------------------------------------- database end
  else
    {
      FFileName  = FileArg.size()>0? FileArg : localPATH + FileNameAbsent;
      setFileName(FFileName);
      if(FileArg.size()>0) readFile(FFileName);

      writePage();
    }

  QSize size = this->minimumSizeHint();
  this->resize(size);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
MainWindow::~MainWindow()
{
  delete ui;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionExecute_triggered()
{
  readPage();

  if(tab_mode == 0) {execute_compound();}
  else              {execute_fusion();  }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::checked(bool b){    ui->edit_max_spin->setEnabled(b);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::makeFusion()
{
  ui->edit_AC->setText(QString::number(Ap+At));
  ui->edit_ZC->setText(QString::number(Zp+Zt));
  ui->edit_NC->setText(QString::number(Ap-Zp+At-Zt));

  CNucleus *comp = new CNucleus(Zp+Zt,Ap+At);
  QString symc = comp->getGName();
  ui->sym_C->setText(symc);
  ui->ME_C->setText(QString::number(comp->getExcessMass()));

  CFus fus(Zp,Ap,Zt,At,Elab,dif);
  ui->edit_ExP->setText(QString::number(fus.Ex));
  ui->edit_QCN->setText(QString::number(fus.Qval));
  ui->edit_E_CM->setText(QString::number(fus.Ecm));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_AP_textEdited(const QString &arg1)
{
  Ap = arg1.toInt(); makeProjectile();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_ZP_textChanged(const QString &arg1)
{
  Zp =  arg1.toInt(); makeProjectile();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::makeProjectile()
{
  int np =  Ap - Zp;
  ui->edit_NP->setText(QString::number(np));
  CNucleus proj(Zp,ui->edit_AP->text().toInt());

  ui->sym_P->setText(proj.getGName());
  ui->ME_P->setText(QString::number(proj.getExcessMass()));

  makeFusion();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_AT_textChanged(const QString &arg1)
{
  At = arg1.toInt();    makeTarget();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_ZT_textChanged(const QString &arg1)
{
  Zt =  arg1.toInt();   makeTarget();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::makeTarget()
{
  int nt = At - Zt;
  ui->edit_NT->setText(QString::number(nt));

  CNucleus targ(Zt,At);
  ui->sym_T->setText(targ.getGName());
  ui->ME_T ->setText(QString::number(targ.getExcessMass()));

  makeFusion();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_ACN_textEdited(const QString &arg1)
{
  iACN = arg1.toInt(); makeCN();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_ZCN_textEdited(const QString &arg1)
{
  iZCN = arg1.toInt(); makeCN();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::makeCN()
{
  int iNCN = iACN - iZCN;
  ui->edit_NCN->setText(QString::number(iNCN));
  CNucleus nuc(iZCN,iACN);
  ui->ME_CN->setText(QString::number(nuc.getExcessMass()));
  ui->sym_CN->setText(nuc.getGName());
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int MainWindow::hash(string key)
{
  int value = 0;

  for(int i=0; i<(int)key.length(); i++)
    {
      value += key[i];
    }

  return (value*key.length()%length);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionAbout_Gemini_triggered() { (new About())->show(); }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_EP_textChanged(const QString &arg1)
{
  Elab = arg1.toFloat();
  CFus fus(Zp,Ap,Zt,At,Elab,dif);
  ui->edit_ExP->setText(QString::number(fus.Ex));
  ui->edit_QCN->setText(QString::number(fus.Qval));
  ui->edit_E_CM->setText(QString::number(fus.Ecm));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_trad_mass_button_clicked()
{
  ::_useAME = false;
  makeProjectile(); makeTarget(); makeCN();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_ame_mass_button_clicked()
{
  ::_useAME = true;
  makeProjectile(); makeTarget(); makeCN();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionGemini_page_triggered()
{
  QDesktopServices::openUrl(QUrl("http://lise.nscl.msu.edu/gemini.html"));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::readPage()
{
  tab_mode = tabWidget->currentIndex();

  iZCN = ui->edit_ZCN->text().toInt();
  iACN = ui->edit_ACN->text().toInt();
  fEx = ui->edit_Ex->text().toFloat();
  fJ =  ui->edit_J ->text().toFloat();

  num_casc   = ui->edit_numEvts->text().toInt();
  //---------------------------------------

  Zp = ui->edit_ZP->text().toInt();
  Zt = ui->edit_ZT->text().toInt();
  Ap = ui->edit_AP->text().toInt();
  At = ui->edit_AT->text().toInt();
  Elab = ui->edit_EP->text().toFloat();
  dif = ui->edit_dif->text().toFloat();

  l0 = ui->edit_max_spin->text().toInt();

  if(ui->rb_spinMax->isChecked()) spinOption=0;
  else                            spinOption=1;

  num_events = ui->num_Events->text().toInt();
  //---------------------------------------

  _useAME    = ui->ame_mass_button->isChecked();
  _useIMF    = ui->cb_IMF->isChecked();
  _useIMFenh = ui->cb_IMFenh->isChecked();

  if(ui->rb_EM0->isChecked()) _optEvap=0;
  else if(ui->rb_EM1->isChecked()) _optEvap=1;
  else                             _optEvap=2;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::writePage()
{
  CNucleus nuc(iZCN,iACN);

  ui->sym_CN->setText(nuc.getGName());
  ui->ME_CN->setText(QString::number(nuc.getExcessMass()));


  ui->edit_ACN->setText(QString::number(iACN));
  ui->edit_ZCN->setText(QString::number(iZCN));
  int iNCN = iACN - iZCN;
  ui->edit_NCN->setText(QString::number(iNCN));
  ui->edit_Ex->setText(QString::number(fEx));
  ui->edit_J->setText(QString::number(fJ));
  ui->edit_numEvts->setText(QString::number(num_casc));

  //---------------------------------------------------------
  ui->edit_ZP->setText(QString::number(Zp));
  ui->edit_ZT->setText(QString::number(Zt));
  ui->edit_AP->setText(QString::number(Ap));
  ui->edit_AT->setText(QString::number(At));
  ui->edit_EP->setText(QString::number(Elab));
  ui->edit_NP->setText(QString::number(Ap-Zp));
  ui->edit_NT->setText(QString::number(At-Zt));


  CNucleus proj(Zp,Ap);
  CNucleus targ(Zt,At);

  ui->ME_P->setText(QString::number(proj.getExcessMass()));
  ui->ME_T->setText(QString::number(targ.getExcessMass()));
  ui->sym_P->setText(proj.getGName());
  ui->sym_T->setText(targ.getGName());

  ui->edit_dif->setText(QString::number(dif));
  ui->edit_max_spin->setText(QString::number(l0));

  ui->rb_spinMax->setChecked(spinOption==0);
  ui->rb_spinBass->setChecked(spinOption==1);

  ui->num_Events->setText(QString::number(num_events));

  makeFusion();

  //---------------------------------------------------------
  tab_mode = tabWidget->currentIndex();
  ui->cb_IMF   ->setChecked(_useIMF);
  ui->cb_IMFenh->setChecked(_useIMFenh);
  ui->rb_EM0->setChecked(_optEvap==0);
  ui->rb_EM1->setChecked(_optEvap==1);
  ui->rb_EM2->setChecked(_optEvap==2);

  if(_useAME) ui->ame_mass_button->setChecked(true);
  else        ui->trad_mass_button->setChecked(true);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void MainWindow::on_action_Exit_triggered()
{
  exit(1);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_pb_testDecay_clicked()
{
  extern int testDecay();

  testDecay();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_pb_testWidth_clicked()
{
  extern int testTheWidth();
  testTheWidth();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_pb_testFusion_clicked()
{
  extern int testFusion();
  testFusion();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
