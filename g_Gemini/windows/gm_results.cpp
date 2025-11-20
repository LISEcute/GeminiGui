#include "gm_results.h"

#include <QBoxLayout>
#include <QIcon>
#include <QToolBar>
#include <QAction>
#include <QFileDialog>
#include <QTextStream>
#include <QPrintDialog>
#include <QPrinter>

#include <QDebug>


extern QString FFileNameHtml;
extern QString localPATH;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
Result_Widget::Result_Widget(QString result_string, QWidget *parent) : QDialog(parent)
{
    QVBoxLayout *layout = new QVBoxLayout;
    result = result_string;

    this->setWindowIcon(QIcon(":/Gemini_logo.png"));
    this->setModal(false);
    QToolBar *toolbar = new QToolBar;

    toolbar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    QAction *saveAction = new QAction(QIcon(":/save29.png"),"&Save",this);
    toolbar->addAction(saveAction);
    connect(saveAction, SIGNAL(triggered()),this, SLOT(save_clicked()));
    saveAction->setIcon(QIcon(":/save29.png"));
    saveAction->setIconText("Save");
    QAction *printAction = new QAction(tr("&Print"), this);
    printAction->setIcon(QIcon(":/printer70.png"));
    printAction->setIconText("Print");
    toolbar->addAction(printAction);
    connect(printAction, SIGNAL(triggered()),this, SLOT(print_clicked()));

    layout->addWidget(toolbar);

    QTextEdit *text = new QTextEdit(result_string);

    text->setReadOnly(true);
    layout->addWidget(text);
    this->setMinimumSize(1000,600);
    layout->setContentsMargins(0,0,0,0);
    layout->setSpacing(0);
    setLayout(layout);

    QFile f(FFileNameHtml);
    f.open(QIODevice::WriteOnly);
    if(!f.isOpen()){
        qDebug() << "Error, File could not be opened.";
        return;
        }
    QTextStream stream(&f);
    stream << result;
    f.close();

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void Result_Widget::save_clicked()
{
    QString filename = QFileDialog::getSaveFileName(this, tr("Save File"),
                                                    FFileNameHtml,tr("Gemini results (*.html)"));
    if(filename.size() <= 0)
        {
        QFile f(filename);
        f.open(QIODevice::WriteOnly);
        if(!f.isOpen()){
            qDebug() << "Error, File could not be opened.";
            return;
            }
        QTextStream stream(&f);
        stream << result;
        f.close();
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void Result_Widget::print_clicked()
{
     QPrinter printer;
     QPrintDialog *printDialog = new QPrintDialog(&printer, this);
     printDialog->setWindowTitle(tr("Print Results File"));

     if (printDialog->exec() != QDialog::Accepted) return;

     QString htmlToPrint(result);
     printer.setFullPage(true);
     QTextDocument textDoc;
     textDoc.setHtml(htmlToPrint);
     textDoc.print(&printer);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

