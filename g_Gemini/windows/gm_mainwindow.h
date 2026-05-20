#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTabWidget>


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
namespace Ui {
class MainWindow;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    int tab_mode;
    QTabWidget *tabWidget;

    int iZCN;// proton number of compound nucleus
    int iACN;//  mass number of compound nucleus
    float fEx;//excitation energy of compound nucleus
    float fJ;//  spin of compound nucleus
    int num_casc;

    int Zp ; // proton number of projectile
    int Ap; // mass number of projectile
    int Zt; // proton number of target
    int At; //mass number of target
    float Elab; // lab energy in MeV
    float dif; //diffuseness of fusion spin distribution in hbar
    float l0;
    int length;
    int num_events;
    int spinOption;


    QString FFileName;
    QString results;

    int hash(std::string);

    void execute_compound();
    void execute_fusion();

    void makeProjectile();
    void makeTarget();
    void makeFusion();
    void makeCN();

    void setFileName(const QString& FileName);
    void write_cs_file(const QString& filename_cs);

//    bool init_base();
//    void getInitialDir(void);
    void readPage();
    void writePage();

    void readFile(const QString& file);

private:
    void printGeminiProperties(QString &results);

private slots:
    void on_actionExecute_triggered();

    void checked(bool b);

    void on_edit_AP_textEdited(const QString &arg1);
    void on_edit_ZP_textChanged(const QString &arg1);
    void on_edit_AT_textChanged(const QString &arg1);
    void on_edit_ZT_textChanged(const QString &arg1);
    void on_edit_ACN_textEdited(const QString &arg1);
    void on_edit_ZCN_textEdited(const QString &arg1);
    void on_action_Open_triggered();
    void on_action_Save_triggered();
    void on_actionAbout_Gemini_triggered();
    void on_edit_EP_textChanged(const QString &arg1);

    void on_trad_mass_button_clicked();
    void on_ame_mass_button_clicked();

    void on_actionGemini_page_triggered();

    void on_actionSave_As_triggered();

    void on_action_Exit_triggered();

    void on_pb_testDecay_clicked();

    void on_pb_testWidth_clicked();

    void on_pb_testFusion_clicked();

private:
    Ui::MainWindow *ui;
    bool _showangdist = false;

};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

#endif // MAINWINDOW_H
