#ifndef VERTEX_H
#define VERTEX_H

#include <QColor>
#include <QApplication>
#include <QGraphicsItem>
#include <QSet>
#include <QDebug>

#include "edge.h"

class Vertex : public QObject, public QGraphicsItem
{
    Q_OBJECT
    Q_PROPERTY(QPointF position READ pos WRITE setPos)
public:
    Vertex();
    ~Vertex();
    void setIndex(const int &number);
    int getIndex() const;
    void setName(const QString &name);
    void addAdj(const int &index);
    int getNumAdj() const;
    int getOneNeighbourIndex(const int &index);
    void removeAdj(const int &index);
    void removeAll();

    void setWeight(const int &w);
    void setWeightAsNumberOfAbsorbed();
    int getWeight() const;

    void setBackgroundColour(const QColor &colour);
    QColor getbackgroundColour() const;
    void setOutlineColour(const QColor &colour);
    QColor getOutLineColour() const;
    void setTextColour(QColor textColor);
    QString getName() const;

    QRectF outlineRect() const;
    QRectF boundingRect() const;
    QPainterPath shape() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    QVariant itemChange(GraphicsItemChange change, const QVariant &value);

    void addEdge(Edge *edge);
    void removeEdge(Edge *edge);
    int getNumberEdge() const;
    void remove_all_edges();

    QList<int> getEdgesIndexForRemovalAnimation(Edge * mainEdge);
    void absorb_removeEdge(int edge_index);
    void absorb_removeEdge(Edge * e);
    void absorb_removeVertex_retainEdge(Edge * e);
    void absorb_retainEdge(Edge * e);
    void absorb_retainEdge_setParentPointer(Edge * e);
    void absorb_singleton(Vertex * v);

    Edge * getEdge(int edgeIndex) const;
    Edge * getHighestWeightEdge();
    Edge * getWeightedProbabilisticEdge();
    Edge * getDegreeProbabilisticEdge();
    Edge * getEdgeFromVertex(Vertex * v2);
    Edge * getSmallestCurrentDegreeNeighbour();
    Edge * getSmallestCurrentWeightNeighbour();
    Edge * getHighestDegreeNeighbour();
    QList<Edge*> getAllEdge() const;

    Vertex * aggregate_get_degree_biased_neighbour();
    Vertex *get_neighbour_fromEdge(int edge_index);
    Vertex *get_neighbour_fromEdge(Edge * e);


    void addAdjVertex(Vertex *adjVertex, QString direction);

    void reSize(int size);
    int getSize() const;

    void setParent(Vertex * v);
    void setParentPointerOnly(Vertex * v);
    Vertex * getParent() const;

    void incrementNoChild();
    void setExtraWeight(const int &w);
    int getExtraWeight() const;
    int getNoChild() const;

    QList<Vertex*> getAbsorbedList();
    QList<int> getNeighbourIndexes();

    void disable_item_selection();

    void setOriginPos(const QPointF &origin);
    QPointF getOriginPos() const;

    void loser_drag_vertex_with_degree_one(Edge * exception);
    void set_vertex_as_dragged_along(bool val);
    void set_vertex_as_absorbed(bool val);

    bool is_vertex_absorbed() const;
    bool is_vertex_dragged_along() const;

    Edge * getMostMutualVertex();
    Edge * correct_getMostMutualVertex();
    Edge * getHighestTriangulateCluster();
    Edge * getProbabilisticTriangulationCoeffVertex();
    Edge * getProbabilisticTriangulationAndWeightVertex();

    void setAffection(const double &p);
    double getAffection() const;

    int getNoOfTriangles(Vertex * v);
    QList<Vertex*> getMyCluster();
    void addMemberToCluster(Vertex * v);
    void addMemberToCluster(QList<Vertex*> v);
    void clearCluster();
    void clearAbsorbed();

    void setDeselectedColour(const QColor &col);
    QColor getDeselectedColour() const;

    void setTruthCommunity(const int &p);
    int getTruthCommunity() const;

    void resetClusterRelevant();
signals:
    void vertex_selected(QList<int> info);
    void vertex_deselected();

public slots:
    void disappear();
    void highlight_v();
    void de_highlight_v();
    void highligh_absorber();
    void de_highlight_absorber();

private slots:
    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

private:
    Vertex * parent;

    QList<int> myNeighbours;
    QList<Vertex*> absorbed;

protected:
    QList<Edge *> myEdge;
    QList<Vertex*> myCluster;
    int roundness(double size) const;

    QColor myBackgroundColour;
    QColor myOutlineColour;
    QColor myTextColour;
    QColor SelectedColour;
    QColor DeSelectedColour;
    QPointF originPos;
    int myIndex;
    QString myName;
    int myWeight;
    int rect_size;
    double affection;
    bool isDraggedAlong;
    bool isAbsorbed;
    int noOfChild;
    int ExtraWeight;
    int myRealCommunity;
};

#endif // VERTEX_H
