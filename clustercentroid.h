#ifndef CLUSTERCENTROID_H
#define CLUSTERCENTROID_H

#include "vertex.h"


class ClusterVertex;
class Arrow;

class ClusterCentroid : public Vertex
{
public:
    ClusterCentroid(qreal x, qreal y);

    QRectF outlineRect() const;
    QRectF boundingRect() const;
    QPainterPath shape() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

    void addCluster(QList<ClusterVertex*> cluster);
    QList<ClusterVertex*> getCluster();

private:
    QList<ClusterVertex*> myCluster;

};

#endif // CLUSTERCENTROID_H
