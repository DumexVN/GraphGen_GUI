#ifndef LINEANIMATOR_H
#define LINEANIMATOR_H

#include <QGraphicsLineItem>
#include <QObject>
#include <QPen>

class LineAnimator : public QObject, public QGraphicsLineItem
{
    Q_OBJECT
    Q_PROPERTY(QPointF endPoint READ endPoint WRITE setEndPoint)
public:
    LineAnimator(QPointF p1, QPointF p2, double size)
    {
        this->setLine(QLineF(p1,p2));
        setPen(QPen(Qt::black, size, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        this->setZValue(-1);
    }

    void setSize(double p)
    {
        setPen(QPen(Qt::black, p, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        update();
    }

    void setEndPoint(const QPointF &p)
    {
        QPointF p1 = this->line().p1();
        this->setLine(QLineF(p1, p));
    }

    QPointF endPoint() const
    {
        return this->line().p2();
    }

    void drawDeletedEdge()
    {
        setPen(QPen(Qt::black, 10, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin));
        update();
    }

public slots:
    void highlight_edge()
    {
        setPen(QPen(Qt::red, 20, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        update();
    }

    void dehighlight_edge()
    {
        setPen(QPen(Qt::black, 10, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        update();
    }
};

#endif // LINEANIMATOR_H
