#ifndef RESULTS_H
#define RESULTS_H

#include <QDialog>
#include <QString>
#include <QTextEdit>
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class Result_Widget : public QDialog
{
    Q_OBJECT
public:
    explicit Result_Widget(QString result_string, QWidget *parent = 0);

signals:

public slots:
    void save_clicked();
    void print_clicked();
private:
    QString result;
};

#endif // RESULTS_H
