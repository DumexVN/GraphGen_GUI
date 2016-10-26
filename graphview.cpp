#include <QtGui>
#include <cmath>

#include "GraphView.h"

GraphView::GraphView(QWidget *parent) : QGraphicsView(parent)
{
    setDragMode(ScrollHandDrag);
}

void GraphView::wheelEvent(QWheelEvent *event)
{
    double numDegrees = -event->delta() / 8.0;
    double numSteps = numDegrees / 15.0;
    double factor = std::pow(1.125, numSteps);
    scale(factor, factor);
}

