#ifndef CLUSTERVERTEX_H
#define CLUSTERVERTEX_H

#include "vertex.h"

class ClusterCentroid;
class Arrow;

class ClusterVertex : public Vertex
{
public:
    ClusterVertex(qreal x, qreal y);

    QRectF outlineRect() const;
    QRectF boundingRect() const;
    QPainterPath shape() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

    void setCentroid(ClusterCentroid * C);
    ClusterCentroid* getCentroid() const;

    void highlight_edge();
    void de_highlight_edge();

private:
    ClusterCentroid * myCentroid;
};

#endif // CLUSTERVERTEX_H
