#include "clustercentroid.h"
#include "clustervertex.h"

#include <QFontMetrics>
#include <QPen>
#include <QPainter>

ClusterCentroid::ClusterCentroid(qreal x, qreal y)
{
    myBackgroundColour = Qt::red;
    myOutlineColour = Qt::black;
    myTextColour = Qt::white;
    myWeight = 1;
    rect_size = 100;
    this->setPos(x,y);
    this->setFlags(QGraphicsItem::ItemIsSelectable);
}

QRectF ClusterCentroid::outlineRect() const
{ // defines shape of the ClusterCentroid
    const int Padding = rect_size;
    QFontMetricsF metrics = QApplication::fontMetrics();
    QRectF rect = metrics.boundingRect(QString(""));
    rect.adjust(-Padding, -Padding, +Padding, +Padding);
    rect.translate(-rect.center());
    return rect;
}

QRectF ClusterCentroid::boundingRect() const
{
    const int Margin = 1;
    return outlineRect().adjusted(-Margin, -Margin, +Margin, +Margin);
}

QPainterPath ClusterCentroid::shape() const
{
    QRectF rect = outlineRect();
    QPainterPath path;
    path.addRoundRect(rect, roundness(rect.width()), roundness(rect.height()));
    return path;
}

void ClusterCentroid::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    QPen pen(myOutlineColour);
    QRectF rect = outlineRect();

    painter->setPen(pen);
    painter->setBrush(myBackgroundColour);
    painter->drawRoundRect(rect, 100,100);


    painter->setPen(QPen(myTextColour, 5));
    QFont font;
    font.setPixelSize(rect_size*3/4);
    painter->setFont(font);
 //   painter->drawText(rect, Qt::AlignCenter, QString::number(myIndex));
}


void ClusterCentroid::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    this->reSize(rect_size*4);
    foreach (ClusterVertex * cv, myCluster)
    {
        cv->reSize(cv->getSize()*4);
        cv->highlight_edge();
    }
}

void ClusterCentroid::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
    this->reSize(rect_size/4);
    foreach (ClusterVertex * cv, myCluster)
    {
        cv->reSize(cv->getSize()/4);
        cv->de_highlight_edge();
    }
}

void ClusterCentroid::addCluster(QList<ClusterVertex*> cluster)
{
    myCluster = cluster;
    foreach (ClusterVertex * cv, myCluster)
        cv->setCentroid(this);
}

QList<ClusterVertex *> ClusterCentroid::getCluster()
{
    return myCluster;
}

