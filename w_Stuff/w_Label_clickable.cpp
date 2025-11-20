#include "w_Stuff/w_Label_clickable.h"

#include <QMouseEvent>

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
ClickableLabel::ClickableLabel(QWidget* parent, Qt::WindowFlags f)
    : QLabel(parent, f)
{
    QCursor cursorb(QPixmap(":/cursor/finger.cur"),0,0);   //  position left top
    setCursor(cursorb);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
ClickableLabel::~ClickableLabel() {}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ClickableLabel::mousePressEvent(QMouseEvent* event)
{
     if( (event->buttons() & Qt::LeftButton) )  emit clicked();
else if( (event->buttons() & Qt::RightButton) ) emit clickedRight();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
