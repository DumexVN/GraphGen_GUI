#include <QtGui>

#include "edge.h"
#include "vertex.h"

Edge::Edge(Vertex *fromVertex, Vertex *toVertex, int index)
{
    myFromVertex = fromVertex;
    myToVertex = toVertex;

    myFromVertex->addEdge(this);
    myToVertex->addEdge(this);

    myFromVertex->addAdj(toVertex->getIndex());
    myToVertex->addAdj(fromVertex->getIndex());

    setZValue(-1);
    setPen(QPen(Qt::black, 5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

    trackVertex();
    this->index = index;

    this->setVisible(false);
}

Edge::~Edge()
{
    myFromVertex->removeEdge(this);
    myToVertex->removeEdge(this);
    myFromVertex->removeAdj(myToVertex->getIndex());
    myToVertex->removeAdj(myFromVertex->getIndex());
}

Vertex *Edge::fromVertex() const
{
    return myFromVertex;
}

Vertex *Edge::toVertex() const
{
    return myToVertex;
}

void Edge::setColour(const QColor &colour)
{
    setPen(QPen(colour, 10.0));
    update();
}

QColor Edge::getcolour() const
{
    return pen().color();
}

void Edge::trackVertex()
{
    setLine(QLineF(myFromVertex->pos(), myToVertex->pos()));
    update();
}

void Edge::removeAll()
{
    myFromVertex->removeEdge(this);
    myToVertex->removeEdge(this);
    myFromVertex->removeAdj(myToVertex->getIndex());
    myToVertex->removeAdj(myFromVertex->getIndex());
}

int Edge::getIndex() const
{
    return index;
}

int Edge::getWeight() const
{
    return Weight;
}

void Edge::setWeight(const int &w)
{
    Weight = w;
}



