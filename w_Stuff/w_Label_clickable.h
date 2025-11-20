#ifndef CLICKABLELABEL_H
#define CLICKABLELABEL_H
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#include <QLabel>
#include <QWidget>
//#include <Qt>
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class ClickableLabel : public QLabel { 
    Q_OBJECT 

public:
    explicit ClickableLabel(QWidget* parent = Q_NULLPTR, Qt::WindowFlags f = Qt::WindowFlags());
    ~ClickableLabel();

signals:
    void clicked();
    void clickedRight();

protected:
    void mousePressEvent(QMouseEvent* event);

};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#endif // CLICKABLELABEL_H
