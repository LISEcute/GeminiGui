#include "gm_about.h"
#include "ui_gm_about.h"

#include "gm_ftype.h"
#include <QDesktopServices>
#include <QUrl>

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

About::About(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::About)
{
    ui->setupUi(this);

    this->setWindowFlags(Qt::WindowCloseButtonHint | Qt::WindowTitleHint);

    ui->label_Version->setText(Gemini_version);
    ui->label_Date->setText(Gemini_date);

    connect(ui->label_LISE, SIGNAL(clicked()), this, SLOT(CmLISE()));
    connect(ui->label_Charge, SIGNAL(clicked()), this, SLOT(CmCharge()));
    connect(ui->label_mail, SIGNAL(clicked()), this, SLOT(CmMail()));

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
About::~About(){    delete ui;}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::CmLISE() {    QDesktopServices::openUrl(QUrl(ui->label_LISE->text()));}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::CmCharge() {    QDesktopServices::openUrl(QUrl(ui->label_Charge->text()));}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::CmMail()
{
QString ss("mailto:tarasov@frib.msu.edu?subject=Gemini++ ");
ss.append(ui->label_Version->text());
QDesktopServices::openUrl(QUrl(ss));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
