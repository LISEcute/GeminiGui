#include "gm_results.h"

#include <QBoxLayout>
#include <QIcon>
#include <QToolBar>
#include <QAction>
#include <QFileDialog>
#include <QTextStream>
#include <QPrintDialog>
#include <QPrinter>
#include <QTabWidget>
#include <QTextBrowser>
#include <QPainter>
#include <QDialog>
#include <QUrl>
#include <algorithm>
#include <cmath>

#include <QDebug>


extern QString FFileNameHtml;

namespace {

class YieldTrendPlotWidget : public QWidget
{
public:
    explicit YieldTrendPlotWidget(const YieldPlotData &data, QWidget *parent = nullptr)
        : QWidget(parent), m_data(data)
    {
        setMinimumSize(980, 560);
        resize(980, 560);
    }

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter p(this);
        p.setRenderHint(QPainter::Antialiasing, true);
        p.fillRect(rect(), Qt::white);

        const int left = 120;
        const int right = 120;
        const int top = 90;
        const int bottom = 110;

        QRect plotRect(left, top, width() - left - right, height() - top - bottom);

        p.setPen(QPen(Qt::black, 1));
        p.drawRect(plotRect);

        if (m_data.points.isEmpty())
        {
            p.drawText(plotRect, Qt::AlignCenter, "No yield data");
            return;
        }

        int minZ = m_data.points.first().z;
        int maxZ = m_data.points.first().z;
        double maxY = 0.0;

        for (const YieldPlotPoint &point : m_data.points)
        {
            minZ = std::min(minZ, point.z);
            maxZ = std::max(maxZ, point.z);
            maxY = std::max(maxY, point.residualEvents);
            maxY = std::max(maxY, point.imfEvents);
        }

        if (maxZ == minZ)
        {
            minZ -= 1;
            maxZ += 1;
        }

        if (maxY <= 0.0)
        {
            p.drawText(plotRect, Qt::AlignCenter, "No yield data");
            return;
        }

        const double xMax = double(maxZ) + std::max(1.0, double(maxZ - minZ) * 0.12);
        const int yTicks = 5;
        const int yStep = std::max(1, int(std::ceil((maxY * 1.25) / double(yTicks))));
        const int yAxisMax = yStep * yTicks;

        auto mapX = [&](double z)
        {
            const double frac = (z - double(minZ)) / (xMax - double(minZ));
            return plotRect.left() + int(frac * plotRect.width());
        };

        auto mapY = [&](double events)
        {
            return plotRect.bottom() - int((events / double(yAxisMax)) * plotRect.height());
        };

        QFont titleFont = p.font();
        titleFont.setPointSize(14);
        titleFont.setBold(true);
        p.setFont(titleFont);
        p.setPen(Qt::black);
        p.drawText(0, 30, width(), 30, Qt::AlignCenter,
                   "Yield Events vs Z");

        QFont tickFont = p.font();
        tickFont.setPointSize(11);
        tickFont.setBold(false);
        p.setFont(tickFont);

        p.setPen(QPen(Qt::black, 2));
        p.drawLine(plotRect.bottomLeft(), plotRect.bottomRight());
        p.drawLine(plotRect.bottomLeft(), plotRect.topLeft());

        for (int i = 0; i <= yTicks; ++i)
        {
            const double frac = double(i) / double(yTicks);
            const int yValue = i * yStep;
            const int y = plotRect.bottom() - int(frac * plotRect.height());

            p.setPen(QPen(QColor(225, 225, 225), 1));
            p.drawLine(plotRect.left(), y, plotRect.right(), y);

            p.setPen(Qt::black);
            p.drawLine(plotRect.left() - 6, y, plotRect.left(), y);
            p.drawText(5, y - 10, left - 15, 20,
                       Qt::AlignRight | Qt::AlignVCenter,
                       QString::number(yValue));
        }

        const int zRange = maxZ - minZ;
        int zStep = 1;
        if (zRange > 18) zStep = 2;
        if (zRange > 36) zStep = 5;

        for (int z = minZ; z <= maxZ; z += zStep)
        {
            const int x = mapX(z);
            p.setPen(Qt::black);
            p.drawLine(x, plotRect.bottom(), x, plotRect.bottom() + 6);
            p.drawText(x - 20, plotRect.bottom() + 10, 40, 20,
                       Qt::AlignCenter, QString::number(z));
        }

        auto drawSeries = [&](double YieldPlotPoint::*member,
                              const QColor &color)
        {
            QPolygonF polyline;
            for (const YieldPlotPoint &point : m_data.points)
                polyline << QPointF(mapX(point.z), mapY(point.*member));

            QPen pen(color, 3);
            pen.setCapStyle(Qt::RoundCap);
            pen.setJoinStyle(Qt::RoundJoin);
            p.setPen(pen);
            p.drawPolyline(polyline);

            p.setBrush(color);
            for (const QPointF &point : polyline)
                p.drawEllipse(point, 4, 4);
            p.setBrush(Qt::NoBrush);
        };

        drawSeries(&YieldPlotPoint::residualEvents, QColor(115, 185, 255));
        drawSeries(&YieldPlotPoint::imfEvents, QColor(255, 150, 150));

        const int legendW = 210;
        const int legendH = 78;
        const int legendX = plotRect.right() - legendW - 18;
        const int legendY = plotRect.top() + 18;
        QRect legendRect(legendX, legendY, legendW, legendH);
        p.fillRect(legendRect, QColor(255, 255, 255, 230));
        p.setPen(QPen(QColor(130, 130, 130), 1));
        p.drawRect(legendRect);

        QFont legendFont = p.font();
        legendFont.setPointSize(12);
        p.setFont(legendFont);

        p.setPen(QPen(QColor(115, 185, 255), 3));
        p.drawLine(legendX + 18, legendY + 26, legendX + 58, legendY + 26);
        p.setPen(Qt::black);
        p.drawText(legendX + 70, legendY + 15, legendW - 82, 24,
                   Qt::AlignLeft | Qt::AlignVCenter, "Residuals");

        p.setPen(QPen(QColor(255, 150, 150), 3));
        p.drawLine(legendX + 18, legendY + 54, legendX + 58, legendY + 54);
        p.setPen(Qt::black);
        p.drawText(legendX + 70, legendY + 43, legendW - 82, 24,
                   Qt::AlignLeft | Qt::AlignVCenter, "IMF");

        QFont axisTitleFont = p.font();
        axisTitleFont.setPointSize(12);
        axisTitleFont.setBold(true);
        p.setFont(axisTitleFont);

        p.drawText(plotRect.left(), height() - 45, plotRect.width(), 28,
                   Qt::AlignCenter, "Z");

        p.save();
        p.translate(42, plotRect.top() + plotRect.height() / 2);
        p.rotate(-90);
        p.drawText(QRect(-plotRect.height() / 2, -20,
                         plotRect.height(), 28),
                   Qt::AlignCenter, "Number of events");
        p.restore();
    }

private:
    YieldPlotData m_data;
};

}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
Result_Widget::Result_Widget(QString result_string, QWidget *parent) : QDialog(parent)
{
    init(result_string, nullptr, nullptr, nullptr);
}

Result_Widget::Result_Widget(QString result_string,
                             const YieldPlotData &yieldPlot,
                             QWidget *parent) : QDialog(parent)
{
    init(result_string, nullptr, nullptr, &yieldPlot);
}

Result_Widget::Result_Widget(QString result_string,
                             const AngularResultTab &residualAngular,
                             const AngularResultTab &imfAngular,
                             QWidget *parent) : QDialog(parent)
{
    init(result_string, &residualAngular, &imfAngular, nullptr);
}

Result_Widget::Result_Widget(QString result_string,
                             const AngularResultTab &residualAngular,
                             const AngularResultTab &imfAngular,
                             const YieldPlotData &yieldPlot,
                             QWidget *parent) : QDialog(parent)
{
    init(result_string, &residualAngular, &imfAngular, &yieldPlot);
}

void Result_Widget::init(QString result_string,
                         const AngularResultTab *residualAngular,
                         const AngularResultTab *imfAngular,
                         const YieldPlotData *yieldPlot)
{
    QVBoxLayout *layout = new QVBoxLayout;
    result = result_string;
    if (yieldPlot) yieldPlotData = *yieldPlot;

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

    QTextBrowser *text = new QTextBrowser;
    text->setOpenLinks(false);
    text->setHtml(result_string);
    tabHtml << result_string;
    connect(text, SIGNAL(anchorClicked(QUrl)), this, SLOT(link_clicked(QUrl)));

    const bool hasAngularTabs = residualAngular || imfAngular;
    if (hasAngularTabs)
    {
        tabs = new QTabWidget;
        tabs->addTab(text, "Results");
    }

    auto addAngularTab = [&](const AngularResultTab *tab)
    {
        if (!tab) return;

        AngularDistributionWidget *angular =
            new AngularDistributionWidget(tab->html,
                                          tab->entries,
                                          tab->sigmaTotal,
                                          tab->nEvents,
                                          tab->lowLimitPercent,
                                          tab->highLimitPercent,
                                          tab->title,
                                          tab->neutronEntry,
                                          tab->protonEntry,
                                          tab->alphaEntry,
                                          tab->gammaEntry,
                                          tab->compoundExcitationMeV,
                                          tab->compoundA,
                                          tab->compoundZ,
                                          tab->recoilBetaCN,
                                          tab->mdir,
                                          this);
        tabs->addTab(angular, tab->label);
        tabHtml << tab->html;
    };

    addAngularTab(residualAngular);
    addAngularTab(imfAngular);

    if (hasAngularTabs)
        layout->addWidget(tabs);
    else
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
    if(filename.size() > 0)
        {
        QFile f(filename);
        f.open(QIODevice::WriteOnly);
        if(!f.isOpen()){
            qDebug() << "Error, File could not be opened.";
            return;
            }
        QTextStream stream(&f);
        stream << currentHtml();
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

     QString htmlToPrint(currentHtml());
     printer.setFullPage(true);
     QTextDocument textDoc;
     textDoc.setHtml(htmlToPrint);
     textDoc.print(&printer);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void Result_Widget::link_clicked(QUrl url)
{
    if (url.scheme() != "gemini") return;

    if (url.host() == "yield_plot")
        openYieldPlotWindow();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

QString Result_Widget::currentHtml() const
{
    if (!tabs) return result;

    const int index = tabs->currentIndex();
    if (index >= 0 && index < tabHtml.size())
        return tabHtml[index];

    return result;
}

void Result_Widget::openYieldPlotWindow()
{
    QDialog *dlg = new QDialog(this);
    dlg->setAttribute(Qt::WA_DeleteOnClose);
    dlg->setWindowIcon(QIcon(":/Gemini_logo.png"));
    dlg->setWindowTitle("Gemini: Yield events vs Z");
    dlg->resize(1040, 640);

    QVBoxLayout *layout = new QVBoxLayout(dlg);
    YieldTrendPlotWidget *plot = new YieldTrendPlotWidget(yieldPlotData, dlg);
    layout->addWidget(plot);

    dlg->show();
}
