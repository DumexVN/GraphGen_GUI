#include "clustervertex.h"
#include "clustercentroid.h"
#include "arrow.h"
#include "typeinfo"

#include <QFontMetrics>
#include <QPen>
#include <QPainter>
#include <QDebug>

ClusterVertex::ClusterVertex(qreal x, qreal y)
{
    myBackgroundColour = Qt::green;
    myOutlineColour = Qt::black;
    myTextColour = Qt::white;
    myWeight = 1;
    myCentroid = 0;
    rect_size = 50;
    this->setPos(x,y);
    this->setFlags(QGraphicsItem::ItemIsSelectable);
}

QRectF ClusterVertex::outlineRect() const
{ // defines shape of the ClusterVertex
    const int Padding = rect_size;
    QFontMetricsF metrics = QApplication::fontMetrics();
    QRectF rect = metrics.boundingRect(QString(""));
    rect.adjust(-Padding, -Padding, +Padding, +Padding);
    rect.translate(-rect.center());
    return rect;
}

QRectF ClusterVertex::boundingRect() const
{
    const int Margin = 1;
    return outlineRect().adjusted(-Margin, -Margin, +Margin, +Margin);
}

QPainterPath ClusterVertex::shape() const
{
    QRectF rect = outlineRect();
    QPainterPath path;
    path.addRoundRect(rect, roundness(rect.width()), roundness(rect.height()));
    return path;
}

void ClusterVertex::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
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
//    painter->drawText(rect, Qt::AlignCenter, QString::number(myIndex));
}

void ClusterVertex::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    myCentroid->mousePressEvent(event);
}

void ClusterVertex::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
    myCentroid->mouseReleaseEvent(event);
}

void ClusterVertex::setCentroid(ClusterCentroid *C)
{
    myCentroid = C;
}

ClusterCentroid *ClusterVertex::getCentroid() const
{
    if (myCentroid == 0)
        qDebug() << "Has Not Been Clustered";
    else
        return myCentroid;
}

void ClusterVertex::highlight_edge()
{
    foreach (Edge * edge , myEdge)
    {
        Arrow * arr = dynamic_cast<Arrow*>(edge);
        if (arr==0)
            edge->setColour(Qt::red);
        else
            arr->changePen(Qt::red,20);
    }
}

void ClusterVertex::de_highlight_edge()
{
    foreach (Edge * edge , myEdge)
    {
        Arrow * arr = dynamic_cast<Arrow*>(edge);
        if (arr==0)
            edge->setColour(Qt::black);
        else
            arr->changePen(Qt::black,10);
    }
}
