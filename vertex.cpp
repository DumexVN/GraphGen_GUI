#include "vertex.h"
#include "edge.h"

#include <QPainter>
#include <QFontMetrics>
#include <QTime>
#include <QDebug>
#include <QTime>

std::default_random_engine gen;


Vertex::Vertex()
{
    myBackgroundColour = Qt::green;
    myOutlineColour = Qt::black;
    myTextColour = Qt::red;
    SelectedColour = Qt::red;
    DeSelectedColour = Qt::green;
    myWeight = 1;
    rect_size = 10;
    parent = 0;
    this->setPos(0,0);
    this->setFlags(QGraphicsItem::ItemIsSelectable);
    isDraggedAlong = false;
    isAbsorbed = false;
    myName = "";
    affection = 0.0;
    noOfChild = 0;
    ExtraWeight = 0;
    myRealCommunity = -1;
    gen.seed(QTime::currentTime().msec());

}

Vertex::~Vertex()
{
    foreach (Edge *edge, myEdge)
        delete edge;
}



void Vertex::setIndex(const int &number)
{
    myIndex = number;
}

int Vertex::getIndex() const
{
    return myIndex;
}

void Vertex::setName(const QString &name)
{
    myName = name;
    update();
}

void Vertex::addAdj(const int &index)
{
    if (!myNeighbours.contains(index))
        myNeighbours.append(index);
    else
    {
        //qDebug() << "NEIGHBOUR ALREADY EXISTS, SKIPPING";
    }
}

int Vertex::getNumAdj() const
{
    return myNeighbours.size();
}

int Vertex::getOneNeighbourIndex(const int &index)
{
    if (myNeighbours.size()+1 > index)
        return myNeighbours.at(index);
    else
        qDebug() << "Out of Bound While Getting A Neighbour";
}

void Vertex::removeAdj(const int &index)
{
    if (myNeighbours.contains(index))
        myNeighbours.removeAll(index);
    else{}
      //  qDebug() << "Out Of Bound While Trying to Remove A Neighbour";
}

void Vertex::removeAll()
{
    foreach (Edge *edge, myEdge)
       edge->removeAll();
}

void Vertex::setWeight(const int &w)
{
    myWeight = w;
}

void Vertex::setWeightAsNumberOfAbsorbed()
{
    if (absorbed.size() == 0)
        return;
    else
        myWeight = absorbed.size();
}

int Vertex::getWeight() const
{
    return myWeight;
}

void Vertex::setBackgroundColour(const QColor &colour)
{
    myBackgroundColour = colour;
    update();
}

QColor Vertex::getbackgroundColour() const
{
    return myBackgroundColour;
}

void Vertex::setOutlineColour(const QColor &colour)
{
    myOutlineColour = colour;
    update();
}

QColor Vertex::getOutLineColour() const
{
    return myOutlineColour;
}

void Vertex::setTextColour(QColor textColor)
{
    myTextColour = textColor;
    update();
}

QString Vertex::getName() const
{
    return myName;
}

QRectF Vertex::outlineRect() const
{ // defines shape of the vertex
    const int Padding = rect_size;
    QFontMetricsF metrics = QApplication::fontMetrics();
    QRectF rect;
  //  if (myName == "")
        rect = metrics.boundingRect(QString(""));
  //  else
 //       rect = metrics.boundingRect(myName);
    rect.adjust(-Padding, -Padding, +Padding, +Padding);
    rect.translate(-rect.center());
    return rect;
}

QRectF Vertex::boundingRect() const
{
    const int Margin = 1;
    return outlineRect().adjusted(-Margin, -Margin, +Margin, +Margin);
}

QPainterPath Vertex::shape() const
{
    QRectF rect = outlineRect();

    QPainterPath path;
    path.addRoundRect(rect, roundness(rect.width()), roundness(rect.height()));
    return path;
}

void Vertex::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    QPen pen(myOutlineColour);
    QRectF rect = outlineRect();

    //draw a dummy boundingRect
    painter->setBrush(Qt::white);
    painter->drawRoundedRect(rect, 100, 100);

    painter->setPen(pen);
    painter->setBrush(myBackgroundColour);


 //   painter->drawRoundRect(rect, roundness(rect.width()),roundness(rect.height()));
    painter->drawRoundRect(rect, 100,100);


    painter->setPen(QPen(myTextColour, 5));
    QFont font;
    font.setPixelSize(rect_size/2);
    painter->setFont(font);
 //   if (myName == "")
 //       painter->drawText(rect, Qt::AlignCenter, QString::number(myIndex));
 //   else
        painter->drawText(rect, Qt::AlignCenter, myName);
}


void Vertex::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    QColor text_changed(Qt::white);
    int origin_size = rect_size;
    this->setBackgroundColour(SelectedColour);
    this->reSize(origin_size*2);
    this->setTextColour(text_changed);
    foreach (Edge *edge, myEdge)
    {
        Vertex * other;
        if (edge->fromVertex() == this)
            other = edge->toVertex();
        else
            other = edge->fromVertex();
        other->setBackgroundColour(SelectedColour);
        other->setTextColour(text_changed);
        other->reSize(origin_size*2);
        edge->setColour(SelectedColour);
        edge->setVisible(true);
    }
    //info: {index, degree, weight}
    QList<int> info;
    info << this->myIndex << myEdge.size() << myWeight;
    emit vertex_selected(info);
}

void Vertex::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
    QColor text_changed(Qt::red);
    int origin_size = rect_size;
    this->setBackgroundColour(DeSelectedColour);
    this->reSize(origin_size/2);
    this->setTextColour(text_changed);
    foreach (Edge *edge, myEdge)
    {
        Vertex * other;
        if (edge->fromVertex() == this)
            other = edge->toVertex();
        else
            other = edge->fromVertex();
        other->setBackgroundColour(other->getDeselectedColour());
        other->setTextColour(text_changed);
        other->reSize(origin_size/2);
        edge->setVisible(false);
    }
    emit vertex_deselected();

}

QVariant Vertex::itemChange(GraphicsItemChange change,
                          const QVariant &value)
{
    if (change == ItemPositionHasChanged) {
        foreach (Edge *link, myEdge)
            link->trackVertex();
            qDebug() << "HERE";
    }
    return QGraphicsItem::itemChange(change, value);
}

int Vertex::roundness(double size) const
{
    const int Diameter = rect_size;
    return 100 * Diameter / int(size);
}

void Vertex::addAdjVertex(Vertex *adjVertex, QString direction)
{

}

void Vertex::reSize(int size)
{
    rect_size = size;
    prepareGeometryChange();
    update();
}

int Vertex::getSize() const
{
    return rect_size;
}

void Vertex::setParent(Vertex *v)
{
    if (v == 0)
        qDebug() << "NULL POINTER PARENT";
    else
    {
        parent = v;
        isAbsorbed = true;
        if (v == this)
            return;
        v->incrementNoChild();
        v->setExtraWeight(this->getWeight());

        //readjustcluter
        v->addMemberToCluster(this);
        v->addMemberToCluster(myCluster);
        myCluster.clear();
    }
}

void Vertex::setParentPointerOnly(Vertex *v)
{
    if (v == 0)
        qDebug() << "NULL POINTER PARENT";
    else
    {
        parent = v;
        isAbsorbed = true;
        if (v == this)
            return;
        v->incrementNoChild();
        v->setExtraWeight(this->getWeight());

        //readjustcluter
    }
}

void Vertex::incrementNoChild()
{
    noOfChild++;
}

void Vertex::setExtraWeight(const int &w)
{
    ExtraWeight += w;
}

int Vertex::getExtraWeight() const
{
    return ExtraWeight;
}

int Vertex::getNoChild() const
{
    return noOfChild;
}

Vertex *Vertex::getParent() const
{
    return parent;
}

QList<Vertex *> Vertex::getAbsorbedList()
{
    return absorbed;
}

QList<int> Vertex::getNeighbourIndexes()
{
    QList<int> indexes;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Vertex * neighbour;
        Edge * e = myEdge.at(i);
        if (e->fromVertex() == this)
            neighbour = e->toVertex();
        else
            neighbour = e->fromVertex();
        indexes.append(neighbour->getIndex());
    }
    return indexes;
}

void Vertex::disable_item_selection()
{
    this->setFlag(QGraphicsItem::ItemIsSelectable, false);
}

void Vertex::setOriginPos(const QPointF &origin)
{
    originPos = origin;
}

QPointF Vertex::getOriginPos() const
{
    return originPos;
}



void Vertex::disappear()
{
    this->setPos(-99999,-99999);
}

void Vertex::highlight_v()
{
    this->setBackgroundColour(Qt::red);
    this->reSize(rect_size*2);
    qDebug() << "HERE";
    update();
}

void Vertex::de_highlight_v()
{
    this->setBackgroundColour(Qt::green);
    this->reSize(rect_size/2);
    update();
}

void Vertex::highligh_absorber()
{
    this->setBackgroundColour(Qt::blue);
    this->reSize(rect_size*2);
    update();
}

void Vertex::de_highlight_absorber()
{
    this->setBackgroundColour(Qt::green);
    this->reSize(rect_size/4);
    update();
}

void Vertex::addEdge(Edge *edge)
{
    if (myEdge.contains(edge))
    {
        //DUP
    }
    else
        myEdge.append(edge);
}

void Vertex::removeEdge(Edge *edge)
{
    myEdge.removeOne(edge);
    Vertex * neighbour = this->get_neighbour_fromEdge(edge);
    myNeighbours.removeOne(neighbour->getIndex());
}

int Vertex::getNumberEdge() const
{
    return myEdge.size();
}

void Vertex::remove_all_edges()
{
    foreach (Edge *edge, myEdge)
        delete edge;
}

Edge *Vertex::getEdgeFromVertex(Vertex * v2)
{
    for (int i = 0; i < myEdge.size(); i++)
    {
        Edge * e = myEdge.at(i);
        if (e->fromVertex() == v2 || e->toVertex() == v2)
            return e;
    }

    qDebug() << "WARNING ! CANNOT FIND NEIGHBOUR FROM AN EDGE! TERMINATING";
    return 0;
}

Edge *Vertex::getSmallestCurrentDegreeNeighbour()
{
    QList<int> indexes;
    int smallest = 9999;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Edge * e = myEdge.at(i);
        Vertex * neighbour = this->get_neighbour_fromEdge(e);
        int w = neighbour->getNumberEdge();
        if (w < smallest)
        {
            indexes.clear();
            indexes.append(i);
            smallest = w;
        }
        else if (w == smallest)
        {
            indexes.append(i);
        }
    }

    if (indexes.size() == 1)
        return myEdge[indexes[0]];
    else
    {
        std::uniform_int_distribution<int> distribution(0,indexes.size()-1);
     //   std::default_random_engine gen;
     //   gen.seed(QTime::currentTime().msec());
        int ran = distribution(gen);
        return myEdge.at(indexes[ran]);
    }
}

Edge *Vertex::getSmallestCurrentWeightNeighbour()
{
    QList<int> indexes;
    int smallest = 9999;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Edge * e = myEdge.at(i);
        Vertex * neighbour = this->get_neighbour_fromEdge(e);
        int w = neighbour->getWeight();
        if (w < smallest)
        {
            indexes.clear();
            indexes.append(i);
            smallest = w;
        }
        else if (w == smallest)
        {
            indexes.append(i);
        }
    }

    if (indexes.size() == 1)
        return myEdge[indexes[0]];
    else
    {
        std::uniform_int_distribution<int> distribution(0,indexes.size()-1);
     //   std::default_random_engine gen;
     //   gen.seed(QTime::currentTime().msec());
        int ran = distribution(gen);
        return myEdge.at(indexes[ran]);
    }
}

QList<int> Vertex::getEdgesIndexForRemovalAnimation(Edge * mainEdge)
{
    QList<int> anim_mater;
    anim_mater.append(mainEdge->getIndex());

    for (int i = 0 ; i < myEdge.size(); i++)
    {
        if (mainEdge != myEdge.at(i))
            anim_mater.append(myEdge.at(i)->getIndex());
    }
    return anim_mater;
}

void Vertex::absorb_removeEdge(int edge_index)
{
    Vertex * neighbour = 0;
    Edge * edge = myEdge.at(edge_index);
    if (edge->fromVertex() == this)
        neighbour = edge->toVertex();
    else
        neighbour = edge->fromVertex();
    neighbour->remove_all_edges();
    absorbed.append(neighbour);
    absorbed.append(neighbour->getAbsorbedList());
    neighbour->setParent(this);
}

void Vertex::absorb_removeEdge(Edge *e)
{
    if (!myEdge.contains(e))
    {
        qDebug() << "Absorb Edge NOT FOUND!";
        return;
    }
    Vertex * neighbour = 0;
    if (e->fromVertex() == this)
        neighbour = e->toVertex();
    else
        neighbour = e->fromVertex();
   // neighbour->loser_drag_vertex_with_degree_one(e);
    neighbour->remove_all_edges();
    absorbed.append(neighbour);
    absorbed.append(neighbour->getAbsorbedList());
    neighbour->setParent(this);
}

void Vertex::absorb_removeVertex_retainEdge(Edge *e)
{
    if (!myEdge.contains(e))
    {
        qDebug() << "Absorb Edge NOT FOUND!";
        return;
    }
    Vertex * neighbour = 0;
    if (e->fromVertex() == this)
        neighbour = e->toVertex();
    else
        neighbour = e->fromVertex();
    absorbed.append(neighbour);
    absorbed.append(neighbour->getAbsorbedList());
    neighbour->setParent(this);
}

void Vertex::absorb_retainEdge(Edge *e)
{
    if (!myEdge.contains(e))
    {
        qDebug() << "Absorb Edge NOT FOUND!";
        return;
    }
    Vertex * neighbour = 0;
    if (e->fromVertex() == this)
        neighbour = e->toVertex();
    else
        neighbour = e->fromVertex();

    absorbed.append(neighbour);
    absorbed.append(neighbour->getAbsorbedList());
    neighbour->setParent(this);
}

void Vertex::absorb_retainEdge_setParentPointer(Edge *e)
{
    if (!myEdge.contains(e))
    {
        qDebug() << "Absorb Edge NOT FOUND!";
        return;
    }
    Vertex * neighbour = 0;
    if (e->fromVertex() == this)
        neighbour = e->toVertex();
    else
        neighbour = e->fromVertex();

    absorbed.append(neighbour);
    absorbed.append(neighbour->getAbsorbedList());
    neighbour->setParentPointerOnly(this);
}

void Vertex::absorb_singleton(Vertex *v)
{
    absorbed.append(v->getAbsorbedList());
    absorbed.append(v);
    v->setParent(this);
    v->remove_all_edges();
}

Vertex *Vertex::get_neighbour_fromEdge(int edge_index)
{
    Vertex * neighbour = 0;
    Edge * edge = myEdge.at(edge_index);
    if (edge->fromVertex() == this)
        neighbour = edge->toVertex();
    else
        neighbour = edge->fromVertex();
    return neighbour;
}

Vertex *Vertex::get_neighbour_fromEdge(Edge *edge)
{
    Vertex * neighbour = 0;
    if (edge->fromVertex() == this)
        neighbour = edge->toVertex();
    else
        neighbour = edge->fromVertex();
    return neighbour;
}

Edge *Vertex::getHighestDegreeNeighbour()
{
    QList<Edge*> edge;
    Edge * final;
    int highest = 0;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Edge * e = myEdge.at(i);
        Vertex * v = get_neighbour_fromEdge(e);
        int d = v->getNumberEdge();
        if (d > highest)
        {
            highest = d;
            edge.clear();
            edge.append(e);
        }
        else if (d == highest)
        {
            edge.append(e);
        }
    }
    if (edge.size() > 1)
    {
        std::uniform_int_distribution<int> distribution(0,edge.size()-1);
        int ran = distribution(gen);
        final = edge.at(ran);
    }
    else
        final = edge.at(0);
    return final;
}

QList<Edge *> Vertex::getAllEdge() const
{
    return myEdge;
}

Edge *Vertex::getEdge(int edgeIndex) const
{
    return myEdge.at(edgeIndex);
}

Edge *Vertex::getHighestWeightEdge()
{
    int highest_w = 0, index = -1;
    for (int i = 0; i < myEdge.size(); i++)
    {
        int cur_w = myEdge.at(i)->getWeight();
        if (cur_w > highest_w)
            index = i;
    }
    return myEdge.at(index);
}

/**  Return a neigbour vertex which was selected with uniform selectiong with a given bias to one's weight
 * @brief Vertex::getWeightedProbabilisticEdge
 * @return
 */
Edge *Vertex::getWeightedProbabilisticEdge()
{
    QList<Edge*> edge;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Edge * e = myEdge.at(i);
        Vertex * v = get_neighbour_fromEdge(e);
        int w = v->getWeight();
        for (int j = 0; j < w; j++)
            edge.append(e);
    }
    std::uniform_int_distribution<int> distribution(0,edge.size()-1);
    int ran = distribution(gen);
    return edge.at(ran);
}

Edge *Vertex::getDegreeProbabilisticEdge()
{
    QList<Edge*> edge;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Edge * e = myEdge.at(i);
        Vertex * v = get_neighbour_fromEdge(e);
        int w = v->getNumberEdge();
        for (int j = 0; j < w; j++)
            edge.append(e);
    }

    std::uniform_int_distribution<int> distribution(0,edge.size()-1);

    int ran = distribution(gen);
    return edge.at(ran);
}


/** For aggregate with degree bias:
 * Return a neigbour vertex which was selected with uniform selectiong with a given bias to one's weight
 * @brief Vertex::aggregate_get_degree_biased_neighbour
 * @return
 */
Vertex *Vertex::aggregate_get_degree_biased_neighbour()
{
    QList<Vertex*> neighbours;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Vertex * neighbour = 0;
        Edge * edge = myEdge.at(i);
        if (edge->fromVertex() == this)
            neighbour = edge->toVertex();
        else
            neighbour = edge->fromVertex();
        int weight = neighbour->getWeight();
        for (int j = 0; j < weight; j++)
            neighbours.append(neighbour);
    }
    std::uniform_int_distribution<int> distribution(0,neighbours.size()-1);
  //  std::default_random_engine gen;
  //  gen.seed(QTime::currentTime().msec());
    int ran = distribution(gen);
    return neighbours.at(ran);
}

/** DRAG ALONG ANY ADJACENT VERTICES WITH DEGREE ONE
 * @brief Vertex::drag_vertex_with_degree_one
 */
void Vertex::loser_drag_vertex_with_degree_one(Edge *exception)
{
    for (int i = 0 ;i < myEdge.size(); i++)
    {
        Edge * e = myEdge.at(i);
        if (e == exception)
            continue;

        Vertex * neighbour;
        if (e->toVertex()==this)
            neighbour = e->fromVertex();
        else
            neighbour = e->toVertex();
        if (neighbour->getNumberEdge() == 1)
        {
            this->absorb_singleton(neighbour);
            neighbour->set_vertex_as_dragged_along(true);
        }
    }
}

void Vertex::set_vertex_as_dragged_along(bool val)
{
    isDraggedAlong = val;
}

void Vertex::set_vertex_as_absorbed(bool val)
{
    isAbsorbed = val;
}

bool Vertex::is_vertex_absorbed() const
{
    return isAbsorbed;
}

bool Vertex::is_vertex_dragged_along() const
{
    return isDraggedAlong;
}



/** GET THE NEIGHBOUR VERTEX THAT HAS THE HIGHEST NUMBER OF MUTUAL TRIANGULATION
 * @brief Vertex::getMostMutualVertex
 */
Edge * Vertex::getMostMutualVertex()
{
    int highest = -1, id = -1, abs_id = -1;
    double highest_coeff = -1;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Vertex * neighbour;
        Edge * e = myEdge.at(i);
        if (e->fromVertex() == this)
            neighbour = e->toVertex();
        else
            neighbour = e->fromVertex();

        if (neighbour->getParent() == this)
            continue;

        QList<int> thisAdj = this->getNeighbourIndexes();
        QList<int> neighbourAdj = neighbour->getNeighbourIndexes();
        QList<int> combine = thisAdj;
        int similar = 0;
        for (int j = 0; j < thisAdj.size(); j++)
        {
            int adj = thisAdj.at(j);
            if (neighbourAdj.contains(adj))
                similar++;
            else
                combine.append(adj);
        }
        double size = combine.size();
        double coeff = (similar +1.0) / size;
        if (coeff > highest_coeff)
        {
            highest_coeff = coeff;
            id = i;
            abs_id = neighbour->getIndex();
        }
        /*
        if (similar > highest)
        {
            highest = similar;
            id = i;
            abs_id = neighbour->getIndex();
        }
        */
    }
    if (id < 0)
    {
        std::uniform_int_distribution<int> distribution(0,myEdge.size()-1);
        int ran = distribution(gen);
        return myEdge.at(ran);
    }
    return myEdge.at(id);
}

/** GET THE NEIGHBOUR VERTEX THAT HAS THE HIGHEST NUMBER OF MUTUAL TRIANGULATION
 * @brief Vertex::correct_getMostMutualVertex
 * @return
 */
Edge *Vertex::correct_getMostMutualVertex()
{
    int highest = 0;
    QList<Edge*> candidate;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Vertex * neighbour = this->get_neighbour_fromEdge(myEdge[i]);

        if (neighbour->getParent() == this)
            continue;

        QSet<int> thisAdj = this->getNeighbourIndexes().toSet();
        QSet<int> neighbourAdj = neighbour->getNeighbourIndexes().toSet();
        int similar = thisAdj.intersect(neighbourAdj).size();
        if (similar > highest)
        {
            highest = similar;
            candidate.clear();
            candidate.append(myEdge[i]);
        }
        else if (similar == highest)
        {
            candidate.append(myEdge[i]);
        }
    }
    if (candidate.empty())
    {
        std::uniform_int_distribution<int> distribution(0,myEdge.size()-1);
        int ran = distribution(gen);
        return myEdge.at(ran);
    }
    else if (candidate.size() > 1)
    {
        std::uniform_int_distribution<int> distribution(0,candidate.size()-1);
        int ran = distribution(gen);
        return candidate.at(ran);
    }
    else
        return candidate.at(0);
}


/** Get the Highest Triangulate Cluster
 * @brief Vertex::getHighestTriangulateCluster
 * @return
 */
Edge *Vertex::getHighestTriangulateCluster()
{
    //first get all neighbour cluster

    QList<Vertex*> centroids;
    QList<int> queried_edge;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Edge * e = myEdge.at(i);
        Vertex * neighbour = this->get_neighbour_fromEdge(e);
        //get the neighbour cluster
        Vertex * par = neighbour->getParent();
        QList<Vertex*> tree;
        while (par != 0)
        {
            neighbour = par;
            if (tree.contains(neighbour))
                break;
            else
                tree.append(neighbour);
            par = neighbour->getParent();
        }
        if (centroids.contains(neighbour))
            continue;
        else
        {
                centroids.append(neighbour);
                queried_edge.append(i);
        }

    }

    int highest_score = 0;
    QList<int> index;

    for (int i = 0 ; i < queried_edge.size(); i++)
    {
        int edge_index = queried_edge.at(i);
        Edge * e = myEdge.at(edge_index);
        Vertex * adjacent = this->get_neighbour_fromEdge(e);

        Vertex * queried_centroid = centroids[i];
        int score = 0;
        // count number of real triangles between V and This
        QList<int> current_iteration_neighbour = adjacent->getNeighbourIndexes();
        //get current centroid indexes
        QList<Vertex*> current_iteration_cluster = queried_centroid->getMyCluster();
        QList<int> current_iteration_cluster_indexes;
        for (int j = 0; j < current_iteration_cluster.size(); j++ )
        {
            Vertex * v = current_iteration_cluster.at(j);
            current_iteration_cluster_indexes.append(v->getIndex());

        }

        //counting score
        for (int j = 0; j < myNeighbours.size(); j++)
        {
            int adj = myNeighbours[j];
            if (adj == adjacent->getIndex())
                continue;
            if (current_iteration_neighbour.contains(adj))
            {
                score++;
            }
            if (current_iteration_cluster_indexes.contains(adj))
            {
                score++;
            }
        }

        if (queried_centroid->getIndex() == this->getIndex())
            score -= myCluster.size();
        if (score < highest_score)
        {
            if (index.size() == 0)
            {
                highest_score = score;
                index.append(edge_index);
            }
        }
        else if (score == highest_score)
        {
            index.append(edge_index);
        }
        else if (score > highest_score)
        {
            highest_score = score;
            index.clear();
            index.append(edge_index);
        }
    }


    int selected_index = -1;
    if (index.size() > 1)
    {
        std::uniform_int_distribution<int> distribution(0,index.size()-1);
     //   std::default_random_engine gen;
      //  gen.seed(QTime::currentTime().msec());
        int ran = distribution(gen);
        selected_index = index.at(ran);
    }
    else
    {
        selected_index = index.at(0);
    }

    Edge * final = myEdge.at(selected_index);
    Vertex * other = this->get_neighbour_fromEdge(final);
    return final;
}

/** A vertex is selected with a probability toward its number of triangulation
 * Sample using uniform-simulate
 * @brief Vertex::getProbabilisticTriangulationCoeffVertex
 * @return selected Edge
 */
Edge *Vertex::getProbabilisticTriangulationCoeffVertex()
{
    QList<Edge*> sample;
    for (int i = 0; i < myEdge.size(); i++)
    {
        Vertex * neighbour;
        Edge * e = myEdge.at(i);
        if (e->fromVertex() == this)
            neighbour = e->toVertex();
        else
            neighbour = e->fromVertex();

        if (neighbour->getParent() == this)
            continue;

        QList<int> thisAdj = this->getNeighbourIndexes();
        QList<int> neighbourAdj = neighbour->getNeighbourIndexes();

        int similar = 0;
        for (int j = 0; j < thisAdj.size(); j++)
        {
            int adj = thisAdj.at(j);
            if (neighbourAdj.contains(adj))
                similar++;
        }

        sample.append(e);
        for (int j =0; j < similar; j++)
            sample.append(e);
    }
    if(!sample.empty())
    {
        std::uniform_int_distribution<int> distribution(0,sample.size()-1);
      //  std::default_random_engine gen;
      //  gen.seed(QTime::currentTime().msec());
        int ran = distribution(gen);
        return sample.at(ran);
    }
    else
    {
        std::uniform_int_distribution<int> distribution(0,myEdge.size()-1);
      //  std::default_random_engine gen;
      //  gen.seed(QTime::currentTime().msec());
        int ran = distribution(gen);
        return myEdge.at(ran);
    }
}

/** GET A NEIGHBOUR VERTEX FOR ABSORPTION. THE RULE IS AS FOLLOWS:
 * For a vertex v, the probability of selecting a neighbour vertex u is:
 * Pr(u) = ((number of triangulation) * (sum/weight) ) / (sum over all neighbour)
 * @brief Vertex::getProbabilisticTriangulationAndWeightVertex
 * @return the selected vertex
 */
Edge *Vertex::getProbabilisticTriangulationAndWeightVertex()
{
    QList<Edge*> sample;
    if (this->getNumberEdge() == 1)
        return myEdge.at(0);

    for (int i = 0; i < myEdge.size(); i++)
    {
        Vertex * neighbour;
        Edge * e = myEdge.at(i);
        if (e->fromVertex() == this)
            neighbour = e->toVertex();
        else
            neighbour = e->fromVertex();


        QList<int> thisAdj = this->getNeighbourIndexes();
        QList<int> neighbourAdj = neighbour->getNeighbourIndexes();

        int similar = 0;
        for (int j = 0; j < thisAdj.size(); j++)
        {
            int adj = thisAdj.at(j);
            if (neighbourAdj.contains(adj))
                similar++;
        }
        int normalise_w = 0;
        if (neighbour->getNoChild() > 0)
            normalise_w = neighbour->getExtraWeight() / neighbour->getNoChild();
        int new_w = (similar*2) * (neighbour->getWeight() + normalise_w);
        for (int j =0; j < new_w; j++)
            sample.append(e);
    }
    if (sample.size() == 0)
        return 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    gen.seed(QTime::currentTime().msec());
    std::uniform_int_distribution<int> distribution(0, sample.size() - 1);
    int ran = distribution(gen);
    return sample.at(ran);
}

void Vertex::setAffection(const double &p)
{
    affection = p;
}

double Vertex::getAffection() const
{
    return affection;
}

int Vertex::getNoOfTriangles(Vertex *v)
{
    QList<int> thisAdj = this->getNeighbourIndexes();
    QList<int> neighbourAdj = v->getNeighbourIndexes();
    qSort(thisAdj);
    qSort(neighbourAdj);
    int similar = 0;
    for (int j = 0; j < thisAdj.size(); j++)
    {
        int adj = thisAdj.at(j);
        if (neighbourAdj.contains(adj))
        {
            similar++;
        }
    }
    return similar;
}

QList<Vertex *> Vertex::getMyCluster()
{
    return myCluster;
}

void Vertex::addMemberToCluster(Vertex *v)
{
    if (myCluster.contains(v))
        qDebug() << "ERR: MEMBER ALREADY IN CLUSTER";
    else
        myCluster.append(v);
}

void Vertex::addMemberToCluster(QList<Vertex *> v)
{
    if (v.size() == 0)
        return;
    for (int i = 0; i < v.size(); i++)
    {
        if (!myCluster.contains(v[i]))
            myCluster.append(v[i]);
    }
}

void Vertex::clearCluster()
{
    myCluster.clear();
}

void Vertex::clearAbsorbed()
{
    absorbed.clear();
}

void Vertex::setDeselectedColour(const QColor &col)
{
    DeSelectedColour = col;
}

QColor Vertex::getDeselectedColour() const
{
    return DeSelectedColour;
}

void Vertex::setTruthCommunity(const int &p)
{
    myRealCommunity = p;
}

int Vertex::getTruthCommunity() const
{
    if (myRealCommunity == -1)
        qDebug() << "Have No Truth Community";
    return myRealCommunity;
}

void Vertex::resetClusterRelevant()
{
    myCluster.clear();
    absorbed.clear();
    myWeight = 1;
    parent = 0;
    isDraggedAlong = false;
    isAbsorbed = false;
    noOfChild = 0;
    ExtraWeight = 0;
    myNeighbours.clear();
    myEdge.clear();
}
