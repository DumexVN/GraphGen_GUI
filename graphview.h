#ifndef GRAPHVIEW_H
#define GRAPHVIEW_H

#include <QGraphicsView>

class GraphView : public QGraphicsView
{
    Q_OBJECT

public:
    GraphView(QWidget *parent = 0);

protected:
    void wheelEvent(QWheelEvent *event);
};

#endif // GRAPHVIEW_H
