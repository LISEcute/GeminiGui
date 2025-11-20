#include "gm_mainwindow.h"
#include <QApplication>
#include <QStyleFactory>

#include <QDir>
#include <QSettings>
#include <QStandardPaths>
//#include <QSqlDatabase>

extern void getInitialDir(void);
extern void getDPIscaling(void);

extern QString LISErootPATH;
extern QString MyDocCompPATH;
extern const char *LISEini;

QString FileArg="";
int fontsizeGlobal=9;
int useHighDpiScaling=1;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

int main(int argc, char *argv[])
{
//--------------------------------------------------------------------------------------
    QApplication app(argc, argv);
    LISErootPATH = QCoreApplication::applicationDirPath();
    if(LISErootPATH.size()<=0)
        {
        QString arg0 = QDir::fromNativeSeparators(argv[0]);
        QDir dir(arg0);
        LISErootPATH = dir.currentPath();
        }

    if(argc>1) FileArg = QDir::fromNativeSeparators(argv[1]);

    //---------------------------------------------------------------------------------- paths + read font size start
    getInitialDir();
    //---------------------------------------------------------------------------------- paths + read font size stop

    //------------------------------------------------------------------ style start

    app.setFont(QFont("Arial", fontsizeGlobal, QFont::Normal));
    app.setOrganizationName("FRIB/MSU");
    app.setStyle(QStyleFactory::create("Fusion"));
    QFile FileStyle(":/w_Main/mainstyle.qss");
    FileStyle.open(QFile::ReadOnly);
    QString StyleSheet = QLatin1String(FileStyle.readAll());
    qApp->setStyleSheet(StyleSheet);
     //------------------------------------------------------------------ style end

    QPalette palette;
    palette.setColor(QPalette::Window, QColor(240, 245, 240)); // set window background color
    palette.setColor(QPalette::WindowText, QColor(0, 0, 0)); // set text color
    palette.setColor(QPalette::Button, QColor(235, 240, 235)); // set button background color
    palette.setColor(QPalette::ButtonText, QColor(0, 10, 0)); // set button text color
    palette.setColor(QPalette::Text, QColor(0, 0, 0)); // set text color for QLineEdit
    palette.setColor(QPalette::Base, QColor(254, 255, 254));
    palette.setColor(QPalette::Light, QColor(255, 255, 255)); //set the light part of the sunken border palette
    palette.setColor(QPalette::Dark, QColor(159, 159, 159)); //set the dark part of the sunken border palette

    QColor bg_disabled(215, 220, 215);
    QColor tc_disabled(145,150,145);

    palette.setColor(QPalette::Disabled, QPalette::Button, bg_disabled); // Set disabled QLineEdit background color
    //  palette.setColor(QPalette::Disabled, QPalette::Base, bg_disabled); // Set disabled QLineEdit background color
    palette.setColor(QPalette::Disabled, QPalette::Text, tc_disabled); // Set disabled QLineEdit text color
    palette.setColor(QPalette::Disabled, QPalette::WindowText, tc_disabled); // Set disabled QLineEdit text color
    palette.setColor(QPalette::Disabled,  QPalette::ButtonText, tc_disabled); // set button text color
    app.setPalette(palette);

//====================================================================================================================


    MainWindow w;
    w.show();

    return app.exec();
}
