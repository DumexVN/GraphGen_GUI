#ifndef EDGE_H
#define EDGE_H

#include <QGraphicsLineItem>

class Vertex;

class Edge :  public QGraphicsLineItem
{

public:
    Edge(Vertex *fromVertex, Vertex *toVertex, int index);
    ~Edge();

    Vertex *fromVertex() const;
    Vertex *toVertex() const;

    void setColour(const QColor &colour);
    QColor getcolour() const;

    void trackVertex();
    void removeAll();

    int getIndex() const;

    int getWeight() const;
    void setWeight(const int &w);

protected:
    Vertex *myFromVertex;
    Vertex *myToVertex;
    int index;
    int Weight;
};

#endif
