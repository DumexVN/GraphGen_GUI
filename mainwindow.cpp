#include "mainwindow.h"
#include "ui_mainwindow.h"


#include "GraphView.h"
#include "arrow.h"
#include "clustercentroid.h"
#include "clustervertex.h"
#include "infowidget.h"
#include <limits>
#include <random>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/topology.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/circle_layout.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/bc_clustering.hpp>

#include <QInputDialog>
#include <QDialogButtonBox>
#include <QMessageBox>
#include <QColor>
#include <QApplication>
#include <QFileDialog>
#include <QTime>
#include <QPropertyAnimation>
#include <QParallelAnimationGroup>
#include <QGraphicsScene>
#include <QSequentialAnimationGroup>
#include <QDockWidget>
#include <QMainWindow>
#include <QPushButton>
#include <QXmlStreamReader>
#include <QFormLayout>

int VERTEX_GLOBAL_SIZE = 50;
double EDGE_GLOBAL_SIZE = 1.0;
bool dense_graph = true;
bool GRAPHICS = true;
std::random_device rd;
std::mt19937 generator(rd());
MainWindow::MainWindow()
{   //set up graphic scenes to display all kinds of stuff
    myScene = new QGraphicsScene(0, 0, 30000, 30000);
    myView = new GraphView;

    myView->setScene(myScene);
    myView->setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    myView->setContextMenuPolicy(Qt::ActionsContextMenu);
    setCentralWidget(myView);
    setWindowTitle(QString("Graph Generator"));
    //set up Actions and Menus bars/options
    createActions();
    createMenus();
    myInfoWidget = 0;
    originGraph = 0;
    graphIsReady = false;
    generator.seed(sqrt(QTime::currentTime().msec()*QTime::currentTime().msec()));
}




/*
 *  BASIC I/O
 *
 */

/** READ GML FILE
 * @brief MainWindow::parseGMLfile
 */
void MainWindow::parseGMLfile()
{
    QString filePath = QFileDialog::getOpenFileName(this,
                                                     ("Load Vertex")
                                                     , "C:/Users/Dumex/Desktop/SocialNetworksCollection"
                                                     , ("GML File (*.gml)"));
    QFile file(filePath);
    QFileInfo fileInfo(file);
    QDir dir = fileInfo.absoluteDir();
    if (!file.exists())
    {
        qDebug() << "FILE NOT FOUND";
        return;
    }

    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);
    QString all = in.readAll();
    QStringList split = all.split("\n");
    QString directed = split[3];
    if (!directed.contains("0"))
    {
        qDebug() << "DIRECTED GRAPH";
        return; //directed graph
    }

    NormalGraph g;
    QStringListIterator iter(split);
    while(iter.hasNext())
    {
        QString str = iter.next();
        if (str.contains("node") || str.contains("edge"))
        {
            bool node = str.contains(" node");
            bool edge = str.contains(" edge");
            if (!node && !edge)
            {
                qDebug() << "FILE FORMATE ERROR";
                return;
            }
            else
                iter.next();
            if (node)
            {
                while (iter.peekNext() != "  ]")
                {
                    QStringList cur = iter.next().split(" ");
                    qDebug() << cur.last();
                    if (cur.contains("id"))
                    {
                        boost::add_vertex(g);
                    }
                    else if (cur.contains("label"))
                    {
                        int index = cur.indexOf("label");
                        index++;
                        QString label;
                        for (index; index < cur.size(); index++)
                        {
                            label.append(cur.at(index));
                        }
                        vertex_label.append(label);
                    }
                    else if (cur.contains("value"))
                    {
                        bool ok;
                        vertex_weight.append(cur.last().toInt(&ok));
                        if (!ok)
                        {
                            QString c = cur.last();
                            QString name = vertex_label.last();
                            QString parsed_name = "";
                            for (int i = 1; i < name.size()-1;i++) parsed_name += name[i];
                            qDebug() << parsed_name;
                            if (vertex_value.contains(c))
                            {
                                vertex_value.insertMulti(c, parsed_name);
                            }
                            else
                                vertex_value.insert(cur.last(), parsed_name);
                        }
                    }
                }
            }
            else if (edge)
            {
                int source = -1, des = -1;
                while (iter.peekNext() != "  ]")
                {
                    QStringList cur = iter.next().split(" ");
                    if (cur.contains("source"))
                    {
                        source = cur.last().toInt();
                    }
                    else if (cur.contains("target"))
                    {
                        des = cur.last().toInt();
                    }
                    else if (cur.contains("value"))
                    {
                        edge_weight.append(cur.last().toInt());
                    }
                }
                if (source != -1 && des != -1)
                {
                    boost::add_edge(source,des,g);
                }
                else
                {
                    qDebug() << "EDGE ERROR" << source << des;
                    return;
                }
            }
            else
            {
                qDebug() << "ERROR: UNKNOWN ERROR WHILE PARSING! RETURNING";
                return;
            }

        }
        else
            continue;
    }
  //  qDebug() << boost::num_vertices(g) << boost::num_edges(g);
    myGraph = g;
    debugOriginalGraph(myGraph);
    if (boost::num_vertices(g) > 1000)
    {
        add_components_without_graphic();
    }
    else
        loadGephiLayout(dir.absolutePath());
}


void MainWindow::add_components_without_graphic()
{
    NormalGraph g = myGraph;
    for (int i = 0; i < boost::num_vertices(g); i++)
    {
        Vertex * v = new Vertex;
        v->setIndex(i);
        myVertexList.append(v);
    }

    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);

        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);

        Edge *edge = new Edge(v1, v2, num_edge);
        LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE);
        myLine.append(line);
        myEdgeList.append(edge);
        myScene->addItem(line);
        num_edge++;
    }
    graphIsReady = true;
}


void MainWindow::loadGephiLayout(QString dirPath)
{
    QList<QPair<double,double> > position;
    for (int i = 0; i < boost::num_vertices(myGraph);i++)
        position.append(qMakePair(0.0,0.0));

    QList<QString> label;
    for (int i = 0; i < boost::num_vertices(myGraph);i++)
        label.append(QString(""));

    QDir dir(dirPath);
    QStringList list = dir.entryList();
    QString filePath = "";
    for (int i = 0; i < list.size(); i++)
    {
        QString filename = list.at(i);
        QString filetype = filename.split(".").at(1);
        if (filetype == "gexf")
        {
            filePath = dirPath + "/" + filename;
        }
    }
    if (filePath == "")
        qDebug() << "ERR: LAYOUT FILE NOT FOUND";
    else
    {
        QFile file(filePath);
        file.open(QFile::ReadOnly | QFile::Text);
        QXmlStreamReader xml(&file);
        while(!xml.atEnd())
        {
            xml.readNextStartElement();
            /* SAMPLE XML
             * <node id="0.0" label="BrighamYoung">
                <attvalues>
                  <attvalue for="value" value="7.0"></attvalue>
                  <attvalue for="componentnumber" value="0"></attvalue>
                  <attvalue for="degree" value="12"></attvalue>
                </attvalues>
                <viz:size value="10.0"></viz:size>
                <viz:position x="-138.91212" y="-88.019615"></viz:position>
                <viz:color r="153" g="153" b="153"></viz:color>
              </node>
                */
            if (xml.isStartElement() && xml.name() == "node")
            {   //Process a node
                QString name;
                int id = -1;
                double x = 0.0 ,y = 0.0;

                QXmlStreamAttributes atts = xml.attributes();
                name = atts.at(1).value().toString();
                QString str = atts.at(0).value().toString().split('.')[0];
                bool ok;
                id = str.toInt(&ok, 10);
                label.replace(id, name);
                while(xml.readNext() != 0)
                {
                    if (xml.name() == "node" && xml.tokenType() == QXmlStreamReader::EndElement)
                        break;
                    else
                    {
                        if (xml.name() == "attvalue")
                        {
                            // do stuff wit attvalue
                        }
                        else if (xml.name() == "position" && xml.tokenType() == QXmlStreamReader::StartElement)
                        {
                             QXmlStreamAttributes pos = xml.attributes();
                             x = pos.at(0).value().toString().toDouble();
                             y = pos.at(1).value().toString().toDouble();
                        }
                    }
                }

                qDebug() << "A NODE" << id << name << "Pos: " << x << "," << y;
                position.replace(id, qMakePair(x,y));
            }

        }
    }
    VERTEX_GLOBAL_SIZE = 45;
    EDGE_GLOBAL_SIZE = 1;

    for (int i = 0; i < position.size(); i++)
    {
        if (position[i] == qMakePair(0.0,0.0))
        {
            boost::remove_vertex(i, myGraph);
        }
        else
        {
            Vertex *v = new Vertex;
            v->reSize(VERTEX_GLOBAL_SIZE*3/4);
            v->setIndex(i);
            v->setPos(position[i].first + 15000, position[i].second + 15000);
            v->setOriginPos(v->pos());
            v->setName(label.at(i));
            originVertexPos.append(v->pos());
            myVertexList.push_back(v);
            myScene->addItem(v);
        }
        //sort out index
    }
    NormalGraph g = myGraph;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);

        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);

        Edge *edge = new Edge(v1, v2, num_edge);
        LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE);
        myLine.append(line);
        myEdgeList.append(edge);
        myScene->addItem(line);
        num_edge++;
    }
    //load community if available
    if (!vertex_value.empty())
    {
       QList<QString> keys = vertex_value.uniqueKeys();
       for (int i = 0; i < keys.size(); i++)
       {
           QList<int> comm_id;
           QList<QString> comm = vertex_value.values(keys[i]);
           for (int j = 0; j < comm.size(); j++)
           {
                QString name = comm[j];
                for (int k = 0; k < myVertexList.size(); k++)
                {
                    QString v_name = myVertexList[k]->getName();
                    v_name.replace(" ","");
                    if (v_name == name)
                        comm_id.append(myVertexList[k]->getIndex());
                }
           }
           ground_truth_communities.append(comm_id);
       }
       colouring_ground_truth_communities();
    }
    export_string_ground_truth_to_id_ground_truth(dirPath);
 //   setUpGraphInfoWidget();
  //  createSeperateScene();
    graphIsReady = true;
}

void MainWindow::read_seperated_graph_input()
/*
 * READ GRAPH WITH *EDGE* FILE ONLY
 */
{
    QString filePath = QFileDialog::getOpenFileName(this,
                                                     ("Load Vertex")
                                                     , "C:/Users/Dumex/Desktop"
                                                     , ("Text Files (*.txt *csv)"));
    QFile file(filePath);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);

    QList<QPair<int,int> > edge;
    while (!in.atEnd())
    {
        QStringList str = in.readLine().split('\t');
        if (str[0].startsWith("#")) continue;
        else
        {
            int v1 = str[0].toInt(), v2 = str[1].toInt();
            edge.append(qMakePair(v1,v2));
        }
    }
    file.close();
    NormalGraph g;

    for (int i = 0; i < edge.size(); i++)
    {
        QPair<int,int> e = edge.at(i);
        boost::add_edge(e.first, e.second, g);
    }


    myGraph = g;
    layout_graph(g);
}


/** Investigating Bridges for small graph (access by index, extreme slow for large graph)
 * Given the Ground Truth is Loaded, Inspect the relation between:
 * a. vertex's degree and bridge i.e. how many bridges a high-degree vertex has vs small
 * @brief MainWindow::examining_bridges
 */
void MainWindow::examining_bridges()
{
    //first get the community size: easy
    qDebug() << "Community and its size:";
    for (int i = 0; i < ground_truth_communities.size(); i++)
    {
        qDebug() << i << ground_truth_communities[i].size();
    }
    //second get the community sum degree: easy
    qDebug() << "Community and its sum degree:";
    for (int i = 0; i < ground_truth_communities.size(); i++)
    {
        QList<int> c = ground_truth_communities[i];
        int sum_degree = 0;
        for (int j = 0; j < c.size(); j++)
        {
            int vi = c[j];
            for (int x = 0; x < myVertexList.size(); x++)
            {
                Vertex * v = myVertexList[x];
                if (v->getIndex() == vi)
                {
                    sum_degree += v->getNumberEdge();
                }
            }
        }
        qDebug() << i << sum_degree;
    }
    //third get the number of bridges for each community: little troublesome
    qDebug() << "Community and bridges (note: bridge are counted twice)";
    QSet<int> bridge; //list of bridges with two end point degree
    for (int i = 0; i < ground_truth_communities.size(); i++)
    {
        QList<int> c = ground_truth_communities[i]; //community i
        int intra = 0 , inter = 0;
        for (int j = 0; j < c.size(); j++)
        {
            int vi = c[j];
            for (int x = 0; x < myVertexList.size(); x++)
            {
                Vertex * v = myVertexList[x];
                if (v->getIndex() == vi) //confirming identity
                {
                    QList<Edge*> e = v->getAllEdge();
                    for (int y = 0; y < e.size(); y++)
                    {
                        Vertex * u = v->get_neighbour_fromEdge(e[y]);
                        if (c.contains(u->getIndex()))
                            intra++;
                        else
                        {
                            inter++;
                            int e_id = e[y]->getIndex();
                            if (!bridge.contains(e_id))
                                bridge.insert(e_id);
                        }
                    }
                }
            }
        }
        qDebug() << i << "Intra:" << intra << "; Inter:" << inter;
    }
    qDebug() << "Total Number of Bridge: " << bridge.size();
    qDebug() << "Edge Index ; v_degree ; u_degree ";
    QSetIterator<int> i(bridge);
    while(i.hasNext())
    {
        int e_id = i.next();
        Edge * e = myEdgeList.at(e_id);
        qDebug() << e_id << ";" << e->fromVertex()->getNumberEdge() << ";" << e->toVertex()->getNumberEdge();
    }
}



/** Load Ground Truth From Here
 * @brief MainWindow::load_ground_truth_communities
 */
void MainWindow::load_ground_truth_communities()
{
    QStringList items;
    items << QString("Community Defined Using String") << QString("Community Defined By Indices");

    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select:"),
                                         tr("Method:"), items, 0, false, &ok);
    if (ok && !item.isEmpty())
    {
        if (item == "Community Defined Using String")
            read_string_ground_truth();
        else if (item == "Community Defined By Indices")
            read_ground_truth_communities();
    }
}

/** LOAD GROUND TRUTH COMMUNITIES FILES
 * COMMUNITY ARE SEPERATED BY CONNECTED COMPONENTS ?
 * @brief MainWindow::read_ground_truth_communities
 */
void MainWindow::read_ground_truth_communities()
{
    qDebug() << "PARSING GROUND TRUTH COMMUNITIES";
    ground_truth_communities.clear();
    QString filePath = QFileDialog::getOpenFileName(this,
                                                     ("Load Vertex")
                                                     , "C:/Users/Dumex/Desktop"
                                                     , ("Text Files (*.txt *csv)"));
    QFile file(filePath);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);

    QList<bool> check_list;
    for (int i = 0; i < myVertexList.size(); i++)
        check_list.append(false);
    int smallest_index = 99999;
    while (!in.atEnd())
    {
        QStringList str = in.readLine().split('\n'); //a community
        QStringList sub_str = str[0].split('\t'); //community split by '\t'
        QList<int> community;
        for (int i = 0; i < sub_str.size(); i++)
        {
            int indx = sub_str[i].toInt(); //get vertex index
            if (indx < smallest_index)
                smallest_index = indx;
            community.append(indx);
            if (check_list.at(indx))
            {
                qDebug() << "ERROR: DUPLICATION WHILE PARSING GROUND TRUTH!";
                return;
            }
            else
                check_list.replace(indx, true);
        }
        ground_truth_communities.append(community);
    }
    //readjust if index does not start from 0
    if (smallest_index > 0)
    {
        for(int i = 0; i < ground_truth_communities.size(); i++)
        {
            QList<int> c = ground_truth_communities[i];
            for (int j = 0; j < c.size(); j++)
            {
                int id = c[j];
                c.replace(j,id-smallest_index);
            }
        }
    }

    qDebug() << "FINISHED! Number of Comm: " << ground_truth_communities.size();
    colouring_ground_truth_communities();
    file.close();
}


/** LOAD GROUND TRUTH COMMUNITIES STRING/NAME
 * COMMUNITY ARE SEPERATED BY CONNECTED COMPONENTS ?
 * @brief MainWindow::read_ground_truth_communities
 */
void MainWindow::read_string_ground_truth()
{
    qDebug() << "PARSING GROUND TRUTH COMMUNITIES";
    ground_truth_communities.clear();
    QString filePath = QFileDialog::getOpenFileName(this,
                                                     ("Load Vertex")
                                                     , "C:/Users/Dumex/Desktop"
                                                     , ("Text Files (*.txt *csv)"));
    QFile file(filePath);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);
    QList<bool> check_list;
    for (int i = 0; i < myVertexList.size(); i++)
        check_list.append(false);
    while (!in.atEnd())
    {
        QStringList str = in.readLine().split('\n'); //a community
        QStringList sub_str = str[0].split('\t'); //community split by '\t'
        QList<int> community;
        for (int i = 0; i < sub_str.size(); i++)
        {
            QString v_name = sub_str[i];
            bool found = false;
            for (int j = 0; j < myVertexList.size(); j++)
            {
                QString label = myVertexList.at(j)->getName();
                if (label == v_name)
                {
                    community.append(myVertexList.at(j)->getIndex());
                    if (check_list.at(j))
                    {
                        qDebug() << "ERROR: DUPLICATION WHILE PARSING GROUND TRUTH!";
                        return;
                    }
                    else
                    {
                        check_list.replace(j, true);
                        found = true;
                        break;
                    }
                }
            }
            if (!found)
                qDebug() << v_name << "NOT FOUND";
        }
        ground_truth_communities.append(community);
    }
    //check sum
    int n = 0;
    for (int i = 0; i < ground_truth_communities.size(); i++)
        n += ground_truth_communities[i].size();

    if (n != myVertexList.size())
        qDebug() << "CHECK SUM: NUMBER OF ELEMENTS DOES NOT MATCH, RESTART!";

    qDebug() << "FINISHED! Number of GROUND TRUTH Comm: " << ground_truth_communities.size();
    colouring_ground_truth_communities();

    for(int i = 0; i < ground_truth_communities.size(); i++)
        qDebug() << ground_truth_communities[i];
    file.close();
    //write the reindexed file
    QFileInfo info(file);
    export_string_ground_truth_to_id_ground_truth(info.absolutePath());
}

/** Change Ground Truth from STRING to INDEX
 * File is written in the same dir as the string truth
 * @brief MainWindow::export_string_ground_truth_to_id_ground_truth
 */
void MainWindow::export_string_ground_truth_to_id_ground_truth(const QString &fileDir)
{
    qDebug() << "Rewritting truth file";
    QString filePath(fileDir+"/truth_file.txt");
    QFile file(filePath);
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream in(&file);
    for (int i = 0; i < ground_truth_communities.size(); i++)
    {
        QList<int> com = ground_truth_communities.at(i);
        for(int j = 0; j < com.size(); j++)
        {
            in << com.at(j);
            if (j != (com.size()-1)){ in << '\t';}
        }
        in << '\n';
    }
    file.close();
}


/** Colour ground truth communities
 * @brief MainWindow::colouring_ground_truth_communities
 */
void MainWindow::colouring_ground_truth_communities()
{
    //default colour set
    QList<QColor> color;
    color << Qt::GlobalColor::red << Qt::GlobalColor::blue << Qt::GlobalColor::green << Qt::GlobalColor::black <<
             Qt::GlobalColor::cyan << Qt::GlobalColor::white << Qt::GlobalColor::yellow << Qt::GlobalColor::gray <<
             Qt::GlobalColor::darkBlue << Qt::GlobalColor::darkGreen << Qt::GlobalColor::darkRed << Qt::GlobalColor::lightGray
          << Qt::GlobalColor::darkMagenta;

    if (ground_truth_communities.empty() || ground_truth_communities.size() > color.size())
    {
        qDebug() << "Colouring Ground Truth Not APPLICABLE!";
        return;
    }

    for (int i = 0; i < ground_truth_communities.size(); i++)
    {
        QList<int> c = ground_truth_communities[i];
        for (int j = 0; j < c.size(); j++)
        {
            Vertex * v = myVertexList[c[j]];
            v->setBackgroundColour(color[i]);
        }
    }
}

// ----------------------- GRAPH GENERATOR -------------------------------------------
/** Generate Hierarchical Layered Gnp
 * from level 0 to k:
 * p0 > p1 > p2 > ... > pk
 * let m be the number of cluster at level 0, m = k, then k = log_2(m);
 * @brief MainWindow::generateLayerGnp
 */

void MainWindow::generateLayerGnp()
{
    int n = 512, m = 16, v_per_c = n/m;
    QList<QPair<int,int> > e;
    std::uniform_real_distribution<double> dis(0,1);
    QList<double> p;
    p << 0.5 << 0.1 << 0.025 << 0.005;
    int count = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i+1; j < n; j++)
        {
            int clus = -1;
            for (int k = 1; k <= 3; k++)
            {
                int exp = qPow(2,k-1);
                if (i/(v_per_c*exp) == j/(v_per_c*exp))
                {
                    clus+=k;
                    break;
                }
            }
            if (clus == -1)
                clus = p.size()-1;
            double ran = dis(generator), edge_p = p[clus];
            bool edge = (ran <= edge_p);
            if (edge)
            {
                count++;
                QPair<int,int> edge = qMakePair(i,j);
                e.append(edge);
            }
        }
    }
    NormalGraph g;
    for (int i = 0; i < n; i++)
        boost::add_vertex(g);
    for (int i = 0; i < e.size(); i++)
        boost::add_edge(e.at(i).first, e.at(i).second, g);
    myGraph = g;
    layout_artificial_network(g,m);

    QList<QList<int> > C;
    for (int i = 0; i < m; i++)
    {
        QList<int> c;
        C.append(c);
    }

    for (int i = 0; i < n;i++)
    {
        int c_id = i/(v_per_c);
        C[c_id].append(i);
    }
    ground_truth_communities = C;
    colouring_ground_truth_communities();
}


/** GENERATE SIMPLE CLUSTER GNP
 * @brief MainWindow::generateSimpleClusterGnp
 */
void MainWindow::generateSimpleClusterGnp()
{
    bool ok;
    double z_in  = QInputDialog::getDouble(this, tr("Enter Value for intra edge:"),
                                          tr("Z_IN:"), 0.333, 0, 1, 5, &ok);
    double z_out = QInputDialog::getDouble(this, tr("Enter Value for inter edge:"),
                                          tr("Z_OUT:"), 0.033, 0, 1, 5, &ok);
    int n = 128, m = 4, v_per_c = n/m;
    QList<QPair<int,int> > e;
    std::uniform_real_distribution<double> dis(0,1);

    for (int i = 0; i < n; i++)
    {
        for (int j = i+1; j < n; j++)
        {
            double ran = dis(generator);
            bool edge = false;
            if (i == j)
                continue;
            else if (i/v_per_c == j/v_per_c)
            {
                if (ran <= z_in)
                {
                    edge = true;
                }
            }
            else
            {
                if (ran <= z_out)
                {
                    edge = true;
                }
            }

            if (edge)
            {
                QPair<int,int> edge = qMakePair(i,j);
                e.append(edge);
            }
        }
    }

    NormalGraph g;
    for (int i = 0; i < n; i++)
        boost::add_vertex(g);
    for (int i = 0; i < e.size(); i++)
        boost::add_edge(e.at(i).first, e.at(i).second, g);
    myGraph = g;
    layout_artificial_network(g,m);

    QList<QList<int> > C;
    for (int i = 0; i < m; i++)
    {
        QList<int> c;
        C.append(c);
    }

    for (int i = 0; i < n;i++)
    {
        int c_id = i/(v_per_c);
        C[c_id].append(i);
    }
    ground_truth_communities = C;
    colouring_ground_truth_communities();

}


/** GENERATING A SIMPLE CYCLE
 * @brief MainWindow::generateCycle
 */
void MainWindow::generateCycle()
{
    bool ok;
    int n = QInputDialog::getInt(this, tr("Number of INTRA edge for each v!"),
                                 tr("z_in: "), 100, 0, 20000, 1, &ok);
    if (!ok)
    {
        qDebug() << "- Error While Generating Cycle";
        return;
    }
    NormalGraph g;
    for (int i = 0; i < n; i++) //v
        boost::add_vertex(g);

    for (int i = 0; i <= n-2; i++) //e
        boost::add_edge(i, i+1, g);

    boost::add_edge(0, n-1, g);

    QList<QPointF> pos;
    for (int i = 0; i < n; i++) //pos
    {
        QPointF p;
        double step = 360.0/n;
        double radius = (n+VERTEX_GLOBAL_SIZE*5)*10;
        double angle = (double) qDegreesToRadians(step*i);
        double x = qCos(angle) * radius + 5000;
        double y = qSin(angle) * radius + 5000;
        p.setX(x);
        p.setY(y);
        pos.append(p);
    }

    myGraph = g;
    typedef boost::graph_traits<NormalGraph>::vertex_iterator ver_iter;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    //get index map
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    for (std::pair<ver_iter, ver_iter> vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(index[*vp.first]);
        v->setPos(pos.at(v->getIndex()));
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
        myScene->addItem(v);
    }

    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE);
        myLine.append(line);
        myEdgeList.append(edge);
        myScene->addItem(line);
        num_edge++;
    }
    //add ground truthar
    dense_graph = true;
    graphIsReady = true;

}


/** GENERATING ARTIFICIAL CLUSTER
 * Let m be number of clusters and n be number of vertices, then:
 * Each cluster consists of n/m vertices,
 * Then for each pair of vertices, an edge is generated independently and uniformly at random:
 * - for vertices in the same cluster: with Pi
 * - for vertices in different clusters: with Pj
 * where: Pj < Pi
 * @brief MainWindow::generateArtificialSocialNetwork
 */
void MainWindow::generateArtificialSocialNetwork()
{
    //set global for now
    //get user input for levels and pi and pj
    QList<double> p;
    bool ok;
    double p0  = QInputDialog::getDouble(this, tr("Enter Value for p0:"),
                                          tr("Amount:"), 0.8, 0, 1, 5, &ok);
    double p1 = QInputDialog::getDouble(this, tr("Enter Value for p1:"),
                                          tr("Amount:"), 0.15, 0, 1, 5, &ok);
    double p2 = QInputDialog::getDouble(this, tr("Enter Value for p2:"),
                                          tr("Amount:"), 0.01, 0, 1, 5, &ok);
    if (!ok)
    {
        qDebug() << "- Please Enter p0 > p1 > p2";
        return;
    }
    p << p0 << p1 << p2;
    //now generate
    int n = 128, m = 4, v_per_c = n/m, levels = 2;
    QList<QPair<int,int> > edges;
    std::uniform_real_distribution<> dis(0, 1);
    int sum_intra = 0, sum_inter = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i+1; j < n; j++)
        {
            if (i == j)
                continue;
            double ran = dis(generator);
            bool edge = false;

            for (int k = 1; k <= levels; k++) //if in the same nested cluster
            {
                if ((j/(v_per_c*k)) == (i/(v_per_c*k)))
                {
                    if (ran <= p[k-1])
                    {
                        edge = true;
                        sum_intra++;
                        break;
                    }
                }
            }
            if (!edge) //else connect a global
            {
                if (ran <= p.last())
                {
                    edge = true;
                    sum_inter++;
                }
            }

            if (edge)
            {
                QPair<int,int> edge = qMakePair(i,j);
                edges.append(edge);
            }
        }
    }
    qDebug() << "- Artificial Network Generated:";
    qDebug() << "- Nested Clusters: " << levels;
    qDebug() << "- Intra Edges: " << sum_intra << "; Inter Edges: " << sum_inter;
    qDebug() << "- Intra/Inter Ratio:" << sum_intra/sum_inter;
    qDebug() << "- Average Degree: " << (sum_intra + sum_inter)/n;
    NormalGraph g;
    for (int i = 0; i < n; i++)
        boost::add_vertex(g);
    for (int i = 0; i < edges.size(); i+=2)
        boost::add_edge(edges.at(i).first, edges.at(i).second, g);

    qDebug() << boost::num_vertices(g) << boost::num_edges(g);
    myGraph = g;
    layout_artificial_network(g, m); //layout graphically
    // add the ground truth here
    int no_cluster = n/(v_per_c*levels);
    QList<QList<int> > C;
    for (int i = 0; i < no_cluster; i++)
    {
        QList<int> c;
        C.append(c);
    }

    for (int i = 0; i < n;i++)
    {
        int c_id = i/(v_per_c*levels);
        C[c_id].append(i);
    }
    ground_truth_communities = C;
    colouring_ground_truth_communities();
}

/** GN Network which consists of fixed z_out and z_in
 * @brief MainWindow::generateGNArtificialNetwork
 */
void MainWindow::generateGNArtificialNetwork()
{
    int const_d = 16;
    bool ok;
    int z_in = QInputDialog::getInt(this, tr("Number of INTRA edge for each v!"),
                                 tr("z_in: "), 1, 0, 16, 1, &ok);
    if (!ok)
    {
        qDebug() << "- Error While Generating GN Network";
        return;
    }
    int z_out = const_d - z_in;
    if (!ok)
    {
        qDebug() << "- Error While Generating GN Network";
        return;
    }
    int n = 128, m = 4, v_per_c = n/m;
    QList<QPair<int,int> > edges;
    for (int i = 0; i < n; i++)
    {
        int clus = i/v_per_c;
        QList<int> in,out;
        for (int j = i+1; j < n; j++)
        {
            if (j == i)
                continue;
            else if (j/v_per_c == clus)
                in << j;
            else
                out << j;
        }
        if (in.size() != 31 || out.size() != 96)
        {
            qDebug() << "- Error While Generating Clusters";
        }
        QList<int> intra = generate_random_edge(z_in, in);
        QList<int> inter = generate_random_edge(z_out, out);
        //create intra edges
        for (int j = 0; j < intra.size(); j++)
        {
            int neigh = intra[j];
            QPair<int,int> edge = qMakePair(i,neigh);
            edges.append(edge);
        }
        //create inter edges
        for (int j = 0; j < inter.size(); j++)
        {
            int neigh = inter[j];
            QPair<int,int> edge = qMakePair(i,neigh);
            edges.append(edge);
        }
    }
    NormalGraph g;
    for (int i = 0; i < n; i++)
        boost::add_vertex(g);
    for (int i = 0; i < edges.size(); i++)
        boost::add_edge(edges.at(i).first, edges.at(i).second, g);

    qDebug() << boost::num_vertices(g) << boost::num_edges(g);
    myGraph = g;
    layout_artificial_network(g, m); //layout graphically
    // add the ground truth here
    int no_cluster = n/(v_per_c);
    QList<QList<int> > C;
    for (int i = 0; i < no_cluster; i++)
    {
        QList<int> c;
        C.append(c);
    }

    for (int i = 0; i < n;i++)
    {
        int c_id = i/(v_per_c);
        C[c_id].append(i);
    }
    ground_truth_communities = C;
    colouring_ground_truth_communities();
}


/*
 * SOME BASIC GRAPH GENERATORS
 */
void MainWindow::generateErdosReyni(double p)
{
    GraphGenerator gen;
    NormalGraph g;
    int no_v = 200;
    g = gen.erdos_reyni_generator(no_v, p);
    myGraph = g;
    typedef boost::graph_traits<NormalGraph>::vertex_iterator ver_iter;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    //get index map
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    for (std::pair<ver_iter, ver_iter> vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(index[*vp.first]);
        v->setPos(0, 0);
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
    }

    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        myEdgeList.append(edge);
        num_edge++;
    }

    setVertexDetails();
    setEdgeWeight();
}

void MainWindow::debugOriginalGraph(NormalGraph g)
{
    typedef boost::graph_traits<NormalGraph>::vertex_iterator ver_iter;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    //get index map
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    for (std::pair<ver_iter, ver_iter> vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        qDebug() << *vp.first <<";" << boost::degree(*vp.first, g);
    }
}

void MainWindow::generateErdosReyni()
/*
 * SIMPLE G(n,p) ERDOS - REYNI RANDOM GRAPH
 */
{
    GraphGenerator gen;
    NormalGraph G;

    int no_v = 20;
    double p = 0.4;
    G = gen.erdos_reyni_generator(no_v, p);
    std::vector<int> component(boost::num_vertices(G));
    int num = boost::connected_components(G, &component[0]);
    qDebug() << "Gnp Generated, No of Connected Components:" << num;
    myGraph = G;
    layout_graph(G);
/*
// GENERATE 2 CLIQUES
    QList<QPair<int,int> > edges;
    for (int i = 0; i < 20; i++)
        boost::add_vertex(G);

    for (int i = 0; i < 20; i++)
    {
        if (i < 10)
        {
            for (int j = 0; j < 10; j++)
            {
                if ( i != j)
                {
                    QPair<int,int> p = qMakePair(i,j);
                    QPair<int,int> r_p = qMakePair(j,i);
                    if (edges.contains(p) || edges.contains(r_p))
                    {}
                    else
                    {
                        edges.append(p);
                        edges.append(r_p);
                    }
                }
            }
        }
        else
        {
            for (int j = 10; j < 20; j++)
            {
                if (i != j)
                {
                    QPair<int,int> p = qMakePair(i,j);
                    QPair<int,int> r_p = qMakePair(j,i);
                    if (edges.contains(p) || edges.contains(r_p))
                    {}
                    else
                    {
                        edges.append(p);
                        edges.append(r_p);
                    }
                }
            }
        }
    }

    int bridge_per_v = 1;
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < bridge_per_v; j++)
        {
            if ((10+j+i) > 19)
            {
                QPair<int,int> p = qMakePair(i,j+i);
                QPair<int,int> r_p = qMakePair(j+i,i);
                if (edges.contains(p) || edges.contains(r_p))
                {}
                else
                {
                    edges.append(p);
                    edges.append(r_p);
                }
            }
            else
            {
                QPair<int,int> p = qMakePair(i,i+j+10);
                QPair<int,int> r_p = qMakePair(i+j+10,i);
                if (edges.contains(p) || edges.contains(r_p))
                {}
                else
                {
                    edges.append(p);
                    edges.append(r_p);
                }
                break;
            }
        }
        break;
    }

    for (int i = 0; i < edges.size(); i+=2)
        boost::add_edge(edges[i].first, edges[i].second, G);

    myGraph = G;

    typedef boost::graph_traits<NormalGraph>::vertex_iterator ver_iter;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    //get index map
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, G);
    int n = 10;
    qreal angle = 0, rad = 500, step = 360/n;

    for (int i = 0; i < 10; i++)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(i);
        qreal this_angle = (i%10)*step;
        this_angle = qDegreesToRadians(this_angle);
        qreal x = rad*qCos(this_angle), y = rad*qSin(this_angle);
        v->setPos(x+5000,y+5000);
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
        myScene->addItem(v);
    }

    for (int i = 10; i < 20; i++)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(i);
        qreal this_angle = (i%10)*step;
        this_angle = qDegreesToRadians(this_angle);
        qreal x = rad*qCos(this_angle), y = rad*qSin(this_angle);
        v->setPos(6500-x,5000+y);
        v->setOriginPos(v->pos());
        originVertexPos.push_back(v->pos());
        myVertexList.push_back(v);
        myScene->addItem(v);
    }

    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = boost::edges(G); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, G);
        v = target(e_desc, G);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE);
        myLine.append(line);
        myEdgeList.append(edge);
        myScene->addItem(line);
        myScene->addItem(edge);
        num_edge++;
    }
    setVertexDetails();
    graphIsReady = true;

   // setUpGraphInfoWidget();
   // createSeperateScene();
    /*
    int no_partition = 15, parition_size = 20;
    QList<double> inter;
    QList<double> intra;

    for (int i = 0; i < no_partition; i++)
    {
        std::uniform_real_distribution<double> out(0, 0.01);
        std::uniform_real_distribution<double> within(0.5, 0.7);
        double pout = out(generator);
        qDebug() << pout;
        double pwithin = within(generator);
        inter.append(pout);
        intra.append(pwithin);
    }

    GraphGenerator gen;
    NormalGraph G;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    for (int i = 0; i < no_partition; i++)
    {
        double p = intra.at(i);
        NormalGraph g = gen.erdos_reyni_generator(parition_size,p);
        for (int j = 0; j < parition_size; j++)
            boost::add_vertex(G);

        IndexMap index = get(boost::vertex_index, g);
        for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
        {
            boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
            e_desc = *ep.first;
            boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
            u = source(e_desc, g);
            v = target(e_desc, g);
            int ui = index[u] + parition_size*i, vi = index[v] + parition_size*i;
            boost::add_edge(ui,vi,G);
        }
    }
    //add random inter cluster edges
    QList<QPair<int,int> > extra;
    for (int i = 0 ; i < boost::num_vertices(G); i++)
    {
        int Ci = i / parition_size;
        double Pi = inter.at(Ci);
        for (int j = 0; j < boost::num_vertices(G);j++)
        {
            int Cj = j / parition_size;
            if (Ci == Cj)
                continue;
            else
            {
                std::uniform_real_distribution<double> roll(0, 1);
                double p = roll(generator);
                if ( p <= Pi)
                {
                    QPair<int,int> e = qMakePair(i,j);
                    QPair<int,int> reverse = qMakePair(j,i);
                    if (!extra.contains(e) && !extra.contains(reverse))
                    {
                        extra.append(e);
                        extra.append(reverse);
                    }
                }
            }
        }
    }

    for (int i = 0; i < extra.size(); i+=2)
    {
        QPair<int,int> pair = extra.at(i);
        boost::add_edge(pair.first, pair.second, G);
    }

    myGraph = G;
    // add random edge between 2 cluster
    //write to file
    QString filepath = "C:/Users/Dumex/Desktop/Gnp/Gnp.csv";
    QFile file(filepath);
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    IndexMap index = get(boost::vertex_index, G);
    for (std::pair<edge_iter, edge_iter> ep = edges(G); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, G);
        v = target(e_desc, G);
        int origin_u = index[u], origin_v = index[v];
        out << origin_u << ";" << origin_v << endl;
    }
  //  layout_graph(G);
  */
}

/** GENERATE ZANARCHY KARATE CLUB
 * @brief MainWindow::generateKarateClub
 */
void MainWindow::generateKarateClub()
{
    int edge[155][2] = {1,	2,
                       1,	3,
                       1,	4,
                       1,	5,
                       1,	6,
                       1,	7,
                       1,	8,
                       1,	9,
                       1,	11,
                       1,	12,
                       1,	13,
                       1,	14,
                       1,	18,
                       1,	20,
                       1,	22,
                       1,	32,
                       2,	1,
                       2,	3,
                       2,	4,
                       2,	8,
                       2,	14,
                       2,	18,
                       2,	20,
                       2,	22,
                       2,	31,
                       3,	1,
                       3,	2,
                       3,	4,
                       3,	8,
                       3,	9,
                       3,	10,
                       3,	14,
                       3,	28,
                       3,	29,
                       3,	33,
                       4,	1,
                       4,	2,
                       4,	3,
                       4,	8,
                       4,	13,
                       4,	14,
                       5,	1,
                       5,	7,
                       5,	11,
                       6,	1,
                       6,	7,
                       6,	11,
                       6,	17,
                       7,	1,
                       7,	5,
                       7,	6,
                       7,	17,
                       8,	1,
                       8,	2,
                       8,	3,
                       8,	4,
                       9,	1,
                       9,	3,
                       9,	31,
                       9,	33,
                       9,	34,
                       10,	3,
                       10,	34,
                       11,	1,
                       11,	5,
                       11,	6,
                       12,	1,
                       13,	1,
                       13,	4,
                       14,	1,
                       14,	2,
                       14,	3,
                       14,	4,
                       14,	34,
                       15,	33,
                       15,	34,
                       16,	33,
                       16,	34,
                       17,	6,
                       17,	7,
                       18,	1,
                       18,	2,
                       19,	33,
                       19,	34,
                       20,	1,
                       20,	2,
                       20,	34,
                       21,	33,
                       21,	34,
                       22,	1,
                       22,	2,
                       23,	33,
                       24,	26,
                       24,	28,
                       24,	30,
                       24,	33,
                       24,	34,
                       25,	26,
                       25,	28,
                       25,	32,
                       26,	24,
                       26,	25,
                       26,	32,
                       27,	30,
                       27,	34,
                       28,	3,
                       28,	24,
                       28,	25,
                       28,	34,
                       29,	3,
                       29,	32,
                       29,	34,
                       30,	24,
                       30,	27,
                       30,	33,
                       30,	34,
                       31,	2,
                       31,	9,
                       31,	33,
                       31,	34,
                       32,	1,
                       32,	25,
                       32,	26,
                       32,	29,
                       32,	33,
                       32,	34,
                       33,	3,
                       33,	9,
                       33,	15,
                       33,	16,
                       33,	19,
                       33,	21,
                       33,	23,
                       33,	24,
                       33,	30,
                       33,	31,
                       33,	32,
                       33,	34,
                       34,	9,
                       34,	10,
                       34,	14,
                       34,	15,
                       34,	16,
                       34,	19,
                       34,	20,
                       34,	21,
                       34,	23,
                       34,	24,
                       34,	27,
                       34,	28,
                       34,	29,
                       34,	30,
                       34,	31,
                       34,	32,
                       34,	33,
    };

    NormalGraph g;
    for (int i = 0; i < 34; i++)
        boost::add_vertex(g);

    QList<QPair<int,int> > edges;
    for (int i = 0 ; i < 155; i++)
    {
        QPair<int,int> p = qMakePair(edge[i][0]-1, edge[i][1]-1);
        QPair<int,int> reverse_p = qMakePair(edge[i][1]-1, edge[i][0]-1);
        if (edges.contains(p) || edges.contains(reverse_p)){}
        else
        {
            edges.append(p);
            edges.append(reverse_p);
        }
    }

    for (int i = 0; i < edges.size(); i+=2)
    {
        QPair<int,int> p = edges.at(i);
        boost::add_edge(p.first, p.second, g);
    }


    qDebug() << boost::num_vertices(g) << boost::num_edges(g);
    myGraph = g;
    //load ground truth
    QList<int> c1,c2;
    //  9	14	20	13	3	12	1	2	4	11	5	6	7	17	22	18	8
    c1 << 8 << 13 << 19 << 12 << 2 << 11 << 0 << 1 << 3 << 10 << 4 << 5 << 6 << 16 << 21 << 17 << 7;
    // 31	33	34	32	29	28	25	10	19	21	23	16	15	27	30	24	26
    c2 << 30 << 32 << 33 << 31 << 28 << 27 << 24 << 9 << 18 << 20 << 22 << 15 << 14 << 26 << 29 << 23 << 25;
    ground_truth_communities << c1 << c2;
    layout_KarateClub();

}

void MainWindow::generateWeightedKarateClub()
{
    int edge[155][3] = {1,2,4,
                        1,3,5,
                        1,4,3,
                        1,5,3,
                        1,6,3,
                        1,7,3,
                        1,8,2,
                        1,9,2,
                        1,11,2,
                        1,12,3,
                        1,13,2,
                        1,14,3,
                        1,18,2,
                        1,20,2,
                        1,22,2,
                        1,32,2,
                        2,1,4,
                        2,3,6,
                        2,4,3,
                        2,8,4,
                        2,14,5,
                        2,18,1,
                        2,20,2,
                        2,22,2,
                        2,31,2,
                        3,1,5,
                        3,2,6,
                        3,4,3,
                        3,8,4,
                        3,9,5,
                        3,10,1,
                        3,14,3,
                        3,28,2,
                        3,29,2,
                        3,33,3,
                        4,1,3,
                        4,2,3,
                        4,3,3,
                        4,8,3,
                        4,13,3,
                        4,14,3,
                        5,1,3,
                        5,7,2,
                        5,11,3,
                        6,1,3,
                        6,7,5,
                        6,11,3,
                        6,17,3,
                        7,1,3,
                        7,5,2,
                        7,6,5,
                        7,17,3,
                        8,1,2,
                        8,2,4,
                        8,3,4,
                        8,4,3,
                        9,1,2,
                        9,3,5,
                        9,31,3,
                        9,33,4,
                        9,34,3,
                        10,3,1,
                        10,34,2,
                        11,1,2,
                        11,5,3,
                        11,6,3,
                        12,1,3,
                        13,1,1,
                        13,4,3,
                        14,1,3,
                        14,2,5,
                        14,3,3,
                        14,4,3,
                        14,34,3,
                        15,33,3,
                        15,34,2,
                        16,33,3,
                        16,34,4,
                        17,6,2,
                        17,7,1,
                        18,1,1,
                        18,2,2,
                        19,33,1,
                        19,34,2,
                        20,1,2,
                        20,2,2,
                        20,34,1,
                        21,33,3,
                        21,34,1,
                        22,1,2,
                        22,2,2,
                        23,33,2,
                        24,26,5,
                        24,28,4,
                        24,30,2,
                        24,33,5,
                        24,34,4,
                        25,26,2,
                        25,28,3,
                        25,32,2,
                        26,24,5,
                        26,25,2,
                        26,32,7,
                        27,30,4,
                        27,34,2,
                        28,3,2,
                        28,24,4,
                        28,25,3,
                        28,34,4,
                        29,3,2,
                        29,32,2,
                        29,34,2,
                        30,24,3,
                        30,27,4,
                        30,33,3,
                        30,34,2,
                        31,2,2,
                        31,9,3,
                        31,33,3,
                        31,34,3,
                        32,1,2,
                        32,25,2,
                        32,26,7,
                        32,29,2,
                        32,33,4,
                        32,34,4,
                        33,3,2,
                        33,9,3,
                        33,15,3,
                        33,16,3,
                        33,19,1,
                        33,21,3,
                        33,23,2,
                        33,24,5,
                        33,30,4,
                        33,31,3,
                        33,32,4,
                        33,34,5,
                        34,9,4,
                        34,10,2,
                        34,14,3,
                        34,15,2,
                        34,16,4,
                        34,19,2,
                        34,20,1,
                        34,21,1,
                        34,23,3,
                        34,24,4,
                        34,27,2,
                        34,28,4,
                        34,29,2,
                        34,30,2,
                        34,31,3,
                        34,32,4,
                        34,33,5,
                       };


    NormalGraph g;
    for (int i = 0; i < 34; i++)
        boost::add_vertex(g);


    QList<QPair<int,int> > edges;
    QList<int> w;
    for (int i = 0 ; i < 155; i++)
    {
        QPair<int,int> p = qMakePair(edge[i][0]-1, edge[i][1]-1);
        QPair<int,int> reverse_p = qMakePair(edge[i][1]-1, edge[i][0]-1);
        if (edges.contains(p) || edges.contains(reverse_p)){}
        else
        {
            edges.append(p);
            edges.append(reverse_p);
            w.append(edge[i][2]);
            w.append(edge[i][2]);
        }
    }

    for (int i = 0; i < edges.size(); i+=2)
    {
        QPair<int,int> p = edges.at(i);
        boost::add_edge(p.first, p.second, g);
        edge_weight.append(w.at(i));
    }
    qDebug() << boost::num_vertices(g) << boost::num_edges(g);
    myGraph = g;
    layout_KarateClub();
    setEdgeWeight();
}

void MainWindow::generateSimpleSample()
{
    int no_vertex = 8;
    int edge [13][2] = {1,2,
                      1,3,
                      1,4,
                      2,3,
                      2,4,
                      3,4,
                      3,5,
                      5,6,
                      5,7,
                      5,8,
                      6,7,
                      6,8,
                      7,8

    };

    QList<QPointF > pos;
    pos << QPointF(500,500) << QPointF(500,1000) <<  QPointF(1000,1000) << QPointF(1000,500)
        << QPointF(1500,1500) << QPointF(1500,2000) << QPointF(2000,1500) << QPointF(2000,2000);

    NormalGraph g;
    for (int i = 0; i < no_vertex; i++)
        boost::add_vertex(g);

    for (int i = 0; i < 13; i++)
        boost::add_edge(edge[i][0]-1, edge[i][1]-1, g);

    myGraph = g;
    for (int i = 0; i < no_vertex; i++)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(i);
        v->setPos(pos.at(i));
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
        myScene->addItem(v);
    }

    int num_edge = 0;
    for (int i =0; i < 13; i++)
    {
        Vertex *v1 = myVertexList.at(edge[i][0]-1);
        Vertex *v2 = myVertexList.at(edge[i][1]-1);
        Edge *edge = new Edge(v1, v2, num_edge);
        LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE);
        myLine.append(line);
        myEdgeList.append(edge);
        myScene->addItem(line);
        myScene->addItem(edge);
        num_edge++;
    }

    graphIsReady = true;
}

/** LAYOUT KARATE CLUB
 * @brief MainWindow::layout_KarateClub
 */
void MainWindow::layout_KarateClub()
{
    QList<QPointF> pos;
    pos << QPointF(3526.49, 1403.28)
    << QPointF(3731.3, 1934.58)
    << QPointF(2602.3, 1469.33)
    << QPointF(4068.63, 2272.44)
    << QPointF(3630.72, 1048.68)
    << QPointF(3682.95, 705.598)
    << QPointF(3956.71, 1060.28)
    << QPointF(4507.42, 1996.25)
    << QPointF(2818.4, 1063.78)
    << QPointF(1947.13, 726.904)
    << QPointF(3438.85, 739.917)
    << QPointF(3255.65, 834.49)
    << QPointF(3319.78, 2394.32)
    << QPointF(2927.11, 1750.42)
    << QPointF(1247.82, 1779.51)
    << QPointF(1201.93, 1612.81)
    << QPointF(3946.25, 648.159)
    << QPointF(4507.83, 1677.61)
    << QPointF(1696.61, 885.266)
    << QPointF(2794.17, 2167.45)
    << QPointF(1437.77, 1080.42)
    << QPointF(4196.78, 1267.48)
    << QPointF(1255.5, 1290.76)
    << QPointF(1246.83, 2387.25)
    << QPointF(1917.64, 2607.94)
    << QPointF(1551.92, 2551.56)
    << QPointF(1110.12, 1966.81)
    << QPointF(1885.88, 2030.16)
    << QPointF(2271.27, 2410.43)
    << QPointF(1237.82, 2109.73)
    << QPointF(2658.84, 2380.23)
    << QPointF(2423.11, 1781.9)
    << QPointF(2416.93, 981.651)
    << QPointF(2119.25, 1646.93);

    VERTEX_GLOBAL_SIZE = 75;

    for (int i = 0 ; i < pos.size(); i++)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(i);
        v->setPos(pos.at(i));
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
        myScene->addItem(v);
    }
    NormalGraph g = myGraph;
    int num_edge = 0;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE);
        myLine.append(line);
        myEdgeList.append(edge);
        myScene->addItem(line);
        myScene->addItem(edge);
        num_edge++;
    }
    setVertexDetails();

   // setUpGraphInfoWidget();
    graphIsReady = true;
}


/** Draw the graph using Fructherman_Reingold Force-Directed Algorithm
 * @brief MainWindow::layout_graph
 * @param g
 */
void MainWindow::layout_graph(NormalGraph &g)
{
 //   dense_graph  = (boost::num_vertices(g) > 1000 && boost::num_edges(g) > 1000);
    typedef boost::graph_traits<NormalGraph>::vertex_iterator ver_iter;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    //get index map
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);

    //get position map
    typedef boost::rectangle_topology<> topology_type;
    typedef topology_type::point_type point_type;
    //position vertices
    typedef std::vector<point_type> PositionVec;
    PositionVec position_vec(num_vertices(g));
    typedef boost::iterator_property_map<PositionVec::iterator,
                                          boost::property_map<NormalGraph, boost::vertex_index_t>::type>
                                                        PositionMap;
    PositionMap position(position_vec.begin(), get(boost::vertex_index, g));
    boost::minstd_rand rand(QTime::currentTime().msec());

    /* LAYOUT
     */
    QMessageBox msgBox;
    msgBox.setInformativeText(QString("CHOOSE A LAYOUT!"));
    QPushButton *circleLayoutButton = msgBox.addButton(QString("Circle Layout"), QMessageBox::ActionRole);
    QPushButton *fruchtermanButton = msgBox.addButton(QString("Fruchterman Reingold"), QMessageBox::ActionRole);
    msgBox.exec();

    dense_graph = false;
    int layout = -1;
    if (dense_graph)
    {
        topology_type topo(rand, 500, 500, 10000, 8000);
        VERTEX_GLOBAL_SIZE = 30;
        EDGE_GLOBAL_SIZE = 0.5;
        random_graph_layout(g, position, topo);
        layout = 2;
    }
    else if (msgBox.clickedButton() == fruchtermanButton) {
        // circle layout

        topology_type topo(rand, 500, 500, 5000, 5000);
        random_graph_layout(g, position, topo);
        fruchterman_reingold_force_directed_layout(g, position, topo);
        layout = 0;
    } else if (msgBox.clickedButton() == circleLayoutButton) {
        //fructherman
        VERTEX_GLOBAL_SIZE = 100;
        boost::circle_graph_layout(g, position, myScene->width()/5);
        layout = 1;
    }
    else
    {
        qDebug() << "LAYOUT NOT DEFINED, RETURNING!";
        return;
    }


    for (std::pair<ver_iter, ver_iter> vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(index[*vp.first]);
        if (layout == 1)
            v->setPos(position[*vp.first][0] + 10000, position[*vp.first][1] + 10000 );
        else
            v->setPos(position[*vp.first][0], position[*vp.first][1] );
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
        myScene->addItem(v);
    }

    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE);
        myLine.append(line);
        myEdgeList.append(edge);
        myScene->addItem(line);
        myScene->addItem(edge);
        num_edge++;
    }

    setVertexDetails();
    setEdgeWeight();

   // setUpGraphInfoWidget();
   // createSeperateScene();
    graphIsReady = true;
}

void MainWindow::layout_artificial_network(NormalGraph &g, int no_cluster)
{
    /* LAYOUT
     */
    qDebug() << "G: " << boost::num_vertices(g) << boost::num_edges(g);
    int n = boost::num_vertices(g);
    int v_per_cluster = n/no_cluster;
    QList<QPointF> pos;

    QList<QPointF> cen;
    for (int i = 0; i < no_cluster; i++) //pos
    {
        QPointF p;
        double step = 360.0/no_cluster;
        double radius = 10000;
        double angle = (double) qDegreesToRadians(step*i);
        double x = qCos(angle) * radius + 15000;
        double y = qSin(angle) * radius + 15000;
        p.setX(x);
        p.setY(y);
        cen.append(p);
    }
    for (int i = 0 ; i < n; i++)
    {
        QPointF p;
        std::uniform_real_distribution<double> rad(0, 1);
        double r = rad(generator);
        std::uniform_real_distribution<double> angle(0, 360);
        double phi = angle(generator);
        double x = sqrt(r)*qCos(phi)*1000; //rescale and relocate to middle of screen
        double y = sqrt(r)*qSin(phi)*1000; //rescale and relocate to middle of screen
        if (no_cluster > 4) //spread
        {
            int k = i/v_per_cluster;
            p.setX(x);
            p.setY(y);
            p+=cen[k];
        }
        else // just normal 4 cluster
        {
            if ( i/v_per_cluster == 0)
            {
                x+=3000;
                y+=3000;
                p.setX(x);
                p.setY(y);
            }
            else if (i/v_per_cluster == 1)
            {
                x+=7000;
                y+=3000;
                p.setX(x);
                p.setY(y);
            }
            else if (i/v_per_cluster == 2)
            {
                x+=3000;
                y+=8000;
                p.setX(x);
                p.setY(y);
            }
            else
            {
                x+=7000;
                y+=8000;
                p.setX(x);
                p.setY(y);
            }
        }
        pos.append(p);
    }

    typedef boost::graph_traits<NormalGraph>::vertex_iterator ver_iter;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    //get index map
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);

    for (std::pair<ver_iter, ver_iter> vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(index[*vp.first]);
        v->setPos(pos.at(v->getIndex()));
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
        myScene->addItem(v);
    }

    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE);
        myLine.append(line);
        myEdgeList.append(edge);
        myScene->addItem(line);
        num_edge++;
    }
    //add ground truthar
    dense_graph = true;
    graphIsReady = true;
}

void MainWindow::setVertexDetails()
{
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        if(!vertex_label.empty())
        {
            v->setName(vertex_label.at(i));
        }
        if (!vertex_weight.empty())
            v->setWeight(vertex_weight.at(i));
    }
}

void MainWindow::setEdgeWeight()
{
    if (edge_weight.empty())
        return;
    for (int i = 0; i < myEdgeList.size(); i++)
    {
        Edge * e = myEdgeList.at(i);
        e->setWeight(edge_weight.at(i));
    }
}


/*
 * CHECK IF THE GRAPH IS IN GOOD CONDITION
 */
bool MainWindow::checkGraphCondition()
{
    if (myVertexList.empty())
    {
        qDebug() << "V is empty, GENERATE A GRAPH FIRST!";
        return false;
    }
    if (!graphIsReady)
    {
       // qDebug() << "Graph is DISCONNECTED!";
        return false;
    }
    else return true;
}

/*
 * RECONNECT THE GRAPH AFTER AN AGGREGATION HAS BEEN DONE
 */
void MainWindow::reConnectGraph()
{
    resetGraphics();
    hierarchy.clear();
    centroids.clear();
    for (int i = 0; i < originVertexPos.size(); i++)
    {
        Vertex * v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(i);
        v->setPos(originVertexPos.at(i));
        v->setOriginPos(v->pos());
        myVertexList.push_back(v);
        myScene->addItem(v);
      }

    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, myGraph);

    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = edges(myGraph); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, myGraph);
        v = target(e_desc, myGraph);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        myEdgeList.append(edge);
        if (GRAPHICS)
        {
            LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE);
            myLine.append(line);
            myScene->addItem(line);
            myScene->addItem(edge);
        }
        num_edge++;
    }
 //   setVertexDetails();
    setEdgeWeight();
 //   setUpGraphInfoWidget();
    graphIsReady = true;
}

// --------------------------- RANDOM REMOVAL OF EDGES -----------------------------
// ---------------------------------------------------------------------------------

/** Randomly Remove Edges of a Vertex
 * @brief MainWindow::random_edge_removal
 */
void MainWindow::random_edge_removal()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays

    int step = 3, final_step = 1, min_d = 2;
    while (step != 0)
    {
        QList<Vertex*> players = myVertexList;
        while(!players.empty()) //start
        {
            //select a vertex uniformly at random
            int size = players.size();
            std::uniform_int_distribution<int> distribution(0,size-1);
            int selected_index = distribution(generator);
            Vertex * selected = players.at(selected_index);
            //get a neighbour
            int no_neighbour = selected->getNumberEdge();
            if (no_neighbour == min_d) // if there is no neighbour, declare a winner
            {
                players.removeOne(selected);
            }
            else // else absorb
            {
                int before = no_neighbour;
                //number of edges to be removed
                int to_be_removed = -1;
                if (step == final_step)
                    to_be_removed = no_neighbour - min_d;
                else
                    to_be_removed = no_neighbour / step;
                for (int i = 0; i < to_be_removed; i++)
                {
                    int no_adj = selected->getNumberEdge();
                    std::uniform_int_distribution<int> distribution2(0,no_adj-1);
                    int selected_edge_index = distribution2(generator);
                    Edge * e = selected->getEdge(selected_edge_index);
                    e->removeAll();
                    myScene->removeItem(myLine.at(e->getIndex()));
                }
                qDebug() << "BEFORE REMOVAL: " << before << "; AFTER REMOVAL: " << selected->getNumberEdge()<< "; Step: " << step;
                players.removeOne(selected);
            }
        }
        step--;
    }
    graphIsReady = false;
}



// -----------------------------RANDOM AGGREGATE CLUSTERING -------------------------
// ----------------------------------------------------------------------------------
/** Type I.a - Uniform (Everything is Uniformly at Random)
 * @brief MainWindow::random_aggregate
 */

void MainWindow::random_aggregate()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
  //  QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0, count_sucess = -1, self_loop = 0, fail = 0;
    QTime t0;
    t0.start();
    qDebug() << "STARTING...";
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random

        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        if (selected->is_vertex_absorbed() || selected->getParent() != 0)
        {
            qDebug() << "Candidate Has Been Clustered!!!!! ";
            return;
        }
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
            qDebug() << t << "*\t"  << selected->getIndex() << selected->getIndex();
            self_loop++;
        }
        else // else absorb
        {
            std::uniform_int_distribution<int> distribution2(0,no_neighbour-1);
            int selected_edge_index = distribution2(generator);
            Edge * e = selected->getEdge(selected_edge_index);
            Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
            int dv = selected->getNumAdj(), du = neighbour->getNumAdj();
            hierarchy.append(qMakePair(neighbour->getIndex(), selected->getIndex()));
            selected->absorb_removeEdge(e);
            players.removeOne(neighbour);
            t++;

            if (dv == 2 && du == 2){ qDebug() << t << "S\t"; count_sucess++;}
            else if (dv == 1 && du == 2){ qDebug() << t  << "S\t"; count_sucess++;}
            else if (dv == 2 && du == 1){ qDebug() << t << "F\t"; fail++;}
            else{ qDebug() << t << "F\t";fail++;}
        }
    }
    centroids = winners;
    qDebug() << "SUCESS:" << count_sucess;
    qDebug() << "Self Loop:" << self_loop;
    qDebug("I.a - Time elapsed: %d ms", t0.elapsed());

    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));

  //  draw_hierarchy_tree();
    draw_aggregation_result();
 //   draw_aggregation_result_by_shape();
}

void MainWindow::random_aggregate_n_steps_walk()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    //input n number of steps
    bool ok;
    int n = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"),
                                 tr("Number of Steps:"), 1, 0, 100, 1, &ok);
    if (!ok)
        return;

    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            // do the walk
            for (int i = 0; i < n; i++)
            {
                std::uniform_int_distribution<int> distribution2(0,no_neighbour-1);
                int selected_edge_index = distribution2(generator);
                Edge * e = selected->getEdge(selected_edge_index);
                Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
                hierarchy.append(qMakePair(selected->getIndex(), neighbour->getIndex()));
            }
            players.removeOne(selected);
            t++;
        }
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
    // group_anim->start();
    // connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    //draw_aggregation_result();
    draw_aggregation_retain_vertex_result();
}

void MainWindow::random_walking_constant_restart_remove_edges()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0, steps = 2;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else
        {
            for (int i = 0; i < steps; i++)
            {
                int ran_size = selected->getNumberEdge();
                if (ran_size == 0)
                    break;
                std::uniform_int_distribution<int> distribution2(0,ran_size-1);
                int selected_edge_index = distribution2(generator);
                Edge * e = selected->getEdge(selected_edge_index);
                Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
                e->removeAll();
                myScene->removeItem(e);
                //single disaapearing animation
                LineAnimator * line = myLine.at(e->getIndex());
                QPointF from = line->line().p1(), to = line->line().p2();
                QPointF animFrom, animTo;
                if (selected->x() == from.x() || selected->y() == from.y())
                {
                    line->setLine(QLineF(to, from));
                    animFrom = from;
                    animTo = to;
                }
                else
                {
                    animFrom = to;
                    animTo = from;
                }
                QPropertyAnimation * anim = new QPropertyAnimation(line, "endPoint");
                anim->setStartValue(animFrom);
                anim->setEndValue(animTo);
                anim->setDuration(1000);
                group_anim->addAnimation(anim);
                connect(anim, SIGNAL(stateChanged(QAbstractAnimation::State,QAbstractAnimation::State)),
                        line, SLOT(highlight_edge()));
                connect(anim, SIGNAL(finished()),
                        line, SLOT(dehighlight_edge()));
                hierarchy.append(qMakePair(selected->getIndex(), neighbour->getIndex()));
                //
            }
            players.removeOne(selected);
            t++;
        }
    }

    //group_anim->start()
    connect(group_anim, SIGNAL(finished()), this, SLOT(draw_walked_directed_graph()));
    draw_walked_directed_graph();
}

void MainWindow::random_walking_normal_retain_edges()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0, steps = 2;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else
        {
            Vertex * current = selected;
            for (int i = 0; i < steps; i++)
            {
                int ran_size = current->getNumberEdge();
                int from = -1, to = -1;
                if (ran_size == 0)
                    break;
                std::uniform_int_distribution<int> distribution2(0,ran_size-1);
                int selected_edge_index = distribution2(generator);
                Edge * e = current->getEdge(selected_edge_index);
                from = current->getIndex();
                current = current->get_neighbour_fromEdge(e);
                to = current->getIndex();

                //single disaapearing animation
                /*
                LineAnimator * line = myLine.at(e->getIndex());
                QPointF from = line->line().p1(), to = line->line().p2();
                QPointF animFrom, animTo;
                if (selected->x() == from.x() || selected->y() == from.y())
                {
                    line->setLine(QLineF(to, from));
                    animFrom = from;
                    animTo = to;
                }
                else
                {
                    animFrom = to;
                    animTo = from;
                }
                QPropertyAnimation * anim = new QPropertyAnimation(line, "endPoint");
                anim->setStartValue(animFrom);
                anim->setEndValue(animTo);
                anim->setDuration(1000);
                group_anim->addAnimation(anim);
                connect(anim, SIGNAL(stateChanged(QAbstractAnimation::State,QAbstractAnimation::State)),
                        line, SLOT(highlight_edge()));
                connect(anim, SIGNAL(finished()),
                        line, SLOT(dehighlight_edge()));
                        */
                hierarchy.append(qMakePair(from,to));
                //
            }
            players.removeOne(selected);
            t++;
        }
    }
    connect(group_anim, SIGNAL(finished()), this, SLOT(draw_walked_directed_graph()));
    draw_walked_directed_graph();
}


void MainWindow::random_walking_testing()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    int a,b;
    bool ok;
    a = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"),
                                 tr("First Index:"), 1, 0, 116, 1, &ok);

    b = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"),
                                 tr("Second Index:"), 1, 0, 116, 1, &ok);


    int max_steps = 10, repeat = 10;
    QPair<int,int> pair = qMakePair(a,b);
    Vertex * v1 = myVertexList.at(pair.first);
    Vertex * v2 = myVertexList.at(pair.second);
    int selection = 1;
    if (selection == 1)
    {
       Edge * e = v1->getEdgeFromVertex(v2);
       e->removeAll();
    }
    else
    {}

    QList<int> firstw;
    for (int i = 0; i < 2; ++i)
    {
        Vertex * startV, * nextV;
        if (i == 0)
        {
            startV = myVertexList.at(pair.first);
            startV->setBackgroundColour(Qt::red);
        }
        else
        {
            startV = myVertexList.at(pair.second);
            startV->setBackgroundColour(Qt::blue);
        }
        int t = 0;
        QString str = "";
        qDebug() << "Walk For V" << startV->getIndex() << ":";
        for (int m = 0; m < max_steps; m++)
        {
            int adj = startV->getNumberEdge();
            std::uniform_int_distribution<int> distribution2(0,adj-1);
            int selected_edge_index = distribution2(generator);
            Edge * e = startV->getEdge(selected_edge_index);
            nextV = startV->get_neighbour_fromEdge(e);
            if (i == 0)
                nextV->setBackgroundColour(Qt::red);
            else
                nextV->setBackgroundColour(Qt::blue);

            if (i == 0)
                firstw << nextV->getIndex();
            else
            {
                if (firstw.contains(nextV->getIndex()))
                    nextV->setBackgroundColour(Qt::black);
            }

            startV = nextV;
            str += QString::number(nextV->getIndex()) + ",";
            t++;
        }
        qDebug() << str;
    }
    graphIsReady = false;
}


/** Type I.b - Uniformly and Comparing the CURRENT DEGREE
 * Pr(v) = u.a.r
 * Pr(u) = u.a.r
 * Graph Type: Destructive
 * @brief MainWindow::random_aggregate_with_degree_comparison
 */
void MainWindow::random_aggregate_with_degree_comparison()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }

    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            std::uniform_int_distribution<int> distribution2(0,no_neighbour-1);
            int selected_edge_index = distribution2(generator);
            Edge * e = selected->getEdge(selected_edge_index);
            Vertex * neighbour = selected->get_neighbour_fromEdge(selected_edge_index); //get the neighbour (not clean)
            Vertex * winner, * loser;
            int selected_d = selected->getNumberEdge(), neighbour_d = neighbour->getNumberEdge();
            if (selected_d >= neighbour_d)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }

            //abosbr
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            winner->absorb_removeEdge(e);
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
    qDebug("I.b - Time elapsed: %d ms", t0.elapsed());
    // draw_dense_graph_aggregation_result();
    // group_anim->start();
    // connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}


/** Type I.c - Uniformly and Comparing the ORIGINAL DEGREE
 * Pr(v) = u.a.r
 * Pr(u) = u.a.r
 * Graph Type: Destructive
 * @brief MainWindow::random_aggregate_with_weight_comparison
 */
void MainWindow::random_aggregate_with_weight_comparison()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }

    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
 //   QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {

            std::uniform_int_distribution<int> distribution2(0,no_neighbour-1);
            int selected_edge_index = distribution2(generator);
            Edge * e = selected->getEdge(selected_edge_index);
            Vertex * neighbour = selected->get_neighbour_fromEdge(selected_edge_index); //get the neighbour (not clean)
            Vertex * winner, * loser;
            int selected_w = selected->getWeight(), neighbour_w = neighbour->getWeight();
            if (selected_w >= neighbour_w)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            //create the animation
           // QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
           // QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
           // group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
    qDebug("I.c - Time elapsed: %d ms", t0.elapsed());
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}

/** Type II.a - Select Neighbour With the ORIGINAL DEGREE BIAS
 * Pr(v) = u.a.r
 * Pr(u) = w(u) / w(i) forall i in adj(v)
 * @brief MainWindow::random_aggregate_with_neighbour_degree_bias
 */
void MainWindow::random_aggregate_with_neighbour_initial_degree_bias()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
  //  QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Vertex * neighbour = selected->aggregate_get_degree_biased_neighbour();
            Edge * e = selected->getEdgeFromVertex(neighbour);
            //create the animation
         //   QList<int> anim_mater = neighbour->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
         //   QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, neighbour, selected);
          //  group_anim->addAnimation(single_anim);
            selected->absorb_removeEdge(e);
            hierarchy.append(qMakePair(neighbour->getIndex(), selected->getIndex()));
            players.removeOne(neighbour);
            t++;
        }
    }
    centroids = winners;
    qDebug("II.a - Time elapsed: %d ms", t0.elapsed());
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}

/** Type II.a(i) - Select Neighbour With the ORIGINAL DEGREE BIAS
 * Pr(v) = u.a.r
 * Pr(u) = w(u) / w(i) forall i in adj(v)
 * if w(v) > w(u) then u->v;
 * else vice versa
 * @brief MainWindow::random_aggregate_with_neighbour_degree_bias
 */
void MainWindow::random_aggregate_with_neighbour_initial_degree_bias_with_comparison()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;

    int t = 0;
    QTime t0;
    t0.start();
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Vertex * neighbour = selected->aggregate_get_degree_biased_neighbour();
            Edge * e = selected->getEdgeFromVertex(neighbour);
            Vertex * winner, * loser;
            if (selected->getWeight() >= neighbour->getWeight())
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }
            winner->absorb_removeEdge(e);
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
    qDebug("II.a(i) - Time elapsed: %d ms", t0.elapsed());
    draw_aggregation_result();
}

/** Type II.b - Select Neighbour With the CURRENT Degree Bias
 * Pr(v) = u.a.r
 * Pr(u) = d(u) / sum d(i) forall i in adj(v)
 * u -> v
 * @brief MainWindow::random_aggregate_with_neighbour_CURRENT_degree_bias
 */
void MainWindow::random_aggregate_with_neighbour_CURRENT_degree_bias()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;

    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->getDegreeProbabilisticEdge();
            Vertex * neighbour = selected->get_neighbour_fromEdge(e);
            selected->absorb_removeEdge(e);
            hierarchy.append(qMakePair(neighbour->getIndex(), selected->getIndex()));
            players.removeOne(neighbour);
            t++;
        }
    }
    centroids = winners;
    qDebug("II.b - Time elapsed: %d ms", t0.elapsed());
    draw_aggregation_result();
}



/** Type II.b(i) - Select Neighbour With the CURRENT Degree Bias
 * Pr(v) = u.a.r
 * Pr(u) = d(u) / sum d(i) forall i in adj(v)
 * u -> v
 * @brief MainWindow::random_aggregate_with_neighbour_CURRENT_degree_bias
 */
void MainWindow::random_aggregate_with_neighbour_CURRENT_degree_bias_with_comparison()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;

    int t = 0;
    QTime t0;
    t0.start();
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->getDegreeProbabilisticEdge();
            Vertex * neighbour = selected->get_neighbour_fromEdge(e);
            Vertex * winner, * loser;
            if (selected->getNumberEdge() >= neighbour->getNumberEdge())
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }
            winner->absorb_removeEdge(e);
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
    qDebug("II.b(i) - Time elapsed: %d ms", t0.elapsed());
    draw_aggregation_result();
}


/** Type II.c - Aggregate HIGHEST DEGREE neighbour
 * Pr(v) = u.a.r
 * Select u: arg max d(u)
 * if d(v) < d(u) ...
 * @brief MainWindow::random_aggregate_highest_CURRENT_degree_neighbour
 */
void MainWindow::random_aggregate_highest_CURRENT_degree_neighbour()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
   // QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->getHighestDegreeNeighbour();
            Vertex * neighbour = selected->get_neighbour_fromEdge(e);
            int dv = selected->getNumberEdge(), du = neighbour->getNumberEdge();
            Vertex * winner, * loser;
            if (dv >= du)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }
            //create the animation
          //  QList<int> anim_mater = neighbour->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
          //  QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, neighbour, selected);
         //   group_anim->addAnimation(single_anim);
            //take snapshot

            winner->absorb_removeEdge(e);
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            players.removeOne(loser);
            t++;
        }
    }
//    qDebug() << "C: " << winners.size();
    centroids = winners;
    qDebug("II.c - Time elapsed: %d ms", t0.elapsed());
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}


/** Type II.d - Select Min DEGREE Neighbour
 * Pr(v) = u.a.r
 * Select u: arg min d(u): u in adj(v)
 * @brief MainWindow::random_aggregate_with_minimum_weight_neighbour
 */
void MainWindow::random_aggregate_with_minimum_weight_neighbour()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
   // QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->getSmallestCurrentDegreeNeighbour();
            Vertex * neighbour = selected->get_neighbour_fromEdge(e);
            //create the animation
        //    QList<int> anim_mater = neighbour->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
        //    QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, neighbour, selected);
        //    group_anim->addAnimation(single_anim);
            selected->absorb_removeEdge(e);
            hierarchy.append(qMakePair(neighbour->getIndex(), selected->getIndex()));
            players.removeOne(neighbour);
            t++;
        }
    }
    centroids = winners;
    qDebug("II.d - Time elapsed: %d ms", t0.elapsed());
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}


/** Type II.e Probabilistic Minimum Degree Neighbour (Destructive, Greedy)
 * Pr(v) = d(v)/ sum d(i) forall
 * Select u: arg min u forall u in adj(v)
 * @brief MainWindow::random_aggregate_probabilistic_lowest_degree_neighbour_destructive
 */
void MainWindow::random_aggregate_probabilistic_lowest_degree_neighbour_destructive()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
   // QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        QList<Vertex*> ran_list;
        for (int i = 0; i < players.size(); i++)
        {
            Vertex * v = players.at(i);
            for (int j = 0; j < v->getNumberEdge(); j++)
                ran_list.append(v);
        }
        if (ran_list.size() == 0)
        {
            int size = players.size();
            std::uniform_int_distribution<int> distribution(0,size-1);
            int selected_index = distribution(generator);
            Vertex * selected = players.at(selected_index);
            winners.append(selected);
            players.removeOne(selected);
        }
        else
        {
            int size = ran_list.size();
            std::uniform_int_distribution<int> distribution(0,size-1);
            int selected_index = distribution(generator);
            Vertex * selected = ran_list.at(selected_index);
            //get a neighbour
            int no_neighbour = selected->getNumberEdge();
            if (no_neighbour == 0) // if there is no neighbour, declare a winner
            {
                winners.append(selected);
                players.removeOne(selected);
            }
            else // else absorb
            {
                Edge * e = selected->getSmallestCurrentDegreeNeighbour();
                Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
                Vertex * winner, * loser;
                winner = selected;
                loser = neighbour;
                //create the animation
               // QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
              //  QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
                hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
               // group_anim->addAnimation(single_anim);
                winner->absorb_removeEdge(e);
                players.removeOne(loser);
            }
        }
        t++;
    }
    centroids = winners;
    qDebug("II.e - Time elapsed: %d ms", t0.elapsed());
    // draw_dense_graph_aggregation_result();
    // group_anim->start();
    // connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}


/** Type II.f - Probabilistic Min DEGREE Neighbour (RETENTIVE)
 * Pr(v) = w(v)/ sum w(i) forall i in V
 * Select u: arg min d(u): u in adj(v)
 * u -> v: w(v) += w(u)
 * @brief MainWindow::random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour
 */
void MainWindow::random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
  //  QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        QList<Vertex*> ran_list;
        for (int i = 0 ;i < players.size(); i++)
        {
            Vertex * v = players.at(i);
            int w = v->getWeight();
            for (int j = 0; j < w; j++)
                ran_list.append(v);
        }
        if (ran_list.size() == 0)
        {
            int size = players.size();
            std::uniform_int_distribution<int> distribution(0,size-1);
            int selected_index = distribution(generator);
            Vertex * selected = players.at(selected_index);
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else
        {
            int size = ran_list.size();
            std::uniform_int_distribution<int> distribution(0,size-1);
            int selected_index = distribution(generator);
            Vertex * selected = ran_list.at(selected_index);
            //get a neighbour
            int no_neighbour = selected->getNumberEdge();
            if (no_neighbour == 0) // if there is no neighbour, declare a winner
            {
                winners.append(selected);
                players.removeOne(selected);
                t++;
            }
            else // else absorb
            {
                Edge * e = selected->getSmallestCurrentDegreeNeighbour();
                Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
                Vertex * winner, * loser;
                winner = selected;
                loser = neighbour;
                hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
                winner->absorb_removeEdge(e);
                winner->setWeight(loser->getWeight() + winner->getWeight());
                players.removeOne(loser);
                t++;
            }
        }
    }
    centroids = winners;
    qDebug("II.f - Time elapsed: %d ms", t0.elapsed());
    draw_aggregation_result();
}


/** Type II.g - Deterministic Max Degree Candidate (Destructive)
 * Select candidate: v: arg max d(v)
 * Select u: arg min d(u)
 * @brief MainWindow::random_aggregate_greedy_max_degree
 */
void MainWindow::random_aggregate_greedy_max_degree()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        QList<Vertex*> ran_list;
        int max_d = -1, id = -1;
        for (int i = 0; i < players.size(); i++)
        {
            Vertex * v = players.at(i);
            int dv = v->getNumberEdge();
            if (dv > max_d)
            {
                max_d = dv;
                id = i;
                ran_list.clear();
                ran_list.append(v);
            }
            else if (dv == max_d)
            {
                ran_list.append(v);
            }
        }

        int size = ran_list.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = ran_list.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
        }
        else // else absorb
        {
            Edge * e = selected->getSmallestCurrentDegreeNeighbour();
            Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
            Vertex * winner, * loser;
            winner = selected;
            loser = neighbour;

            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            winner->absorb_removeEdge(e);
            players.removeOne(loser);
        }

        t++;
    }
    centroids = winners;
    qDebug("II.g - Time elapsed: %d ms", t0.elapsed());
    draw_aggregation_result();
}


/** Type II.h (Retentive) Greedy Max Weight (candidate selection)
 * Select Candidate: v = arg max w(v)
 * Select Neighbour: u = arg min w(u)
 * @brief MainWindow::random_aggregate_greedy_max_weight
 */
void MainWindow::random_aggregate_greedy_max_weight()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    //QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        QList<Vertex*> ran_list;
        int max_w = -1;
        for (int i = 0; i < players.size(); i++)
        {
            Vertex * v = players[i];
            int wv = v->getWeight();
            if (wv > max_w)
            {
                max_w = wv;
                ran_list.clear();
                ran_list.append(v);
            }
            else if (wv == max_w)
            {
                ran_list.append(v);
            }
        }

        int size = ran_list.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = ran_list.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->getSmallestCurrentWeightNeighbour();
            Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
            Vertex * winner, * loser;
            winner = selected;
            loser = neighbour;

            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));

            winner->absorb_removeEdge(e);
            winner->setWeight(loser->getWeight() + winner->getWeight());
            players.removeOne(loser);
            t++;
        }

    }
    centroids = winners;
    qDebug("II.h - Time elapsed: %d ms", t0.elapsed());

    draw_aggregation_result();
}


/** Type:: Special Case, Check Again Later
 * SELECT A VERTEX <U> UNIFORMLY AT RANDOM
 * SELECT AN ASSOCIATE EDGE WITH HIGHEST WEIGHT, HENCE VERTEX <V>
 * VERTEX WITH HIGHER WEIGHT ABSORBS THE OTHER
 * @brief MainWindow::random_aggregate_with_highest_edge_weight_and_weight_comparison
 */
void MainWindow::random_aggregate_with_highest_edge_weight_and_weight_comparison()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->getHighestWeightEdge();
            Vertex * neighbour, * winner, * loser;
            if (e->toVertex() == selected)
                neighbour = e->fromVertex();
            else
                neighbour = e->toVertex();

            if (selected == neighbour)
            {   qDebug() << "BUG CHECK" << "SELECTED POINTER == NEIGHBOUR POINTER";
                return;
            }
            int u_w = selected->getWeight(), v_w = neighbour->getWeight();
            if (u_w >= v_w)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            //create the animation
            QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
            QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
            group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            players.removeOne(loser);
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            t++;
        }
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}


/** Type:: Weighted Graph Check Again LAter
 * SELECT A VERTEX <U> UNIFORMLY AT RANDOM
 * SELECT AN EDGE WITH PROBABILITY GIVEN TO THE WEIGHT:
 * Prob(edge_i) = w(i) / sum(w)
 * VERTEX WITH HIGHER WEIGHT ABSORBS
 * @brief MainWindow::random_aggregate_with_edge_weight_bias
 */
void MainWindow::random_aggregate_with_edge_weight_bias_and_weight_comparison()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->getWeightedProbabilisticEdge();
            Vertex * neighbour, * winner, * loser;
            if (e->toVertex() == selected)
                neighbour = e->fromVertex();
            else
                neighbour = e->toVertex();

            if (selected == neighbour)
            {   qDebug() << "BUG CHECK" << "SELECTED POINTER == NEIGHBOUR POINTER";
                return;
            }
            int u_w = selected->getWeight(), v_w = neighbour->getWeight();
            if (u_w >= v_w)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }
            //create the animation
            QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
            QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
            group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
   // draw_dense_graph_aggregation_result();
   // group_anim->start();
   // connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();

}


/** Type III.c - Select Highest Triangles Neighbour Destructive
 * Pr(v) = u.a.r
 * Select u: arg max tri(u) forall u in adj(v)
 * if w(v) > w(u) then u -> v and vice versa
 * SELECTED VERTEX V ABSORBS THE HIGHEST-TRIANGULATED VERTEX U
 * @brief MainWindow::random_aggregate_with_highest_triangulated_vertex
 */
void MainWindow::random_aggregate_with_highest_triangulated_vertex()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
   // QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->correct_getMostMutualVertex();
            Vertex * neighbour, * winner, * loser;
            if (e->toVertex() == selected)
                neighbour = e->fromVertex();
            else
                neighbour = e->toVertex();

            if (selected == neighbour)
            {   qDebug() << "BUG CHECK" << "SELECTED POINTER == NEIGHBOUR POINTER";
                return;
            }
            int u_w = selected->getWeight(), v_w = neighbour->getWeight();
            if (u_w >= v_w)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }
            //create the animation
          //  QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
          //  QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
           // group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
    qDebug("III.c - Time elapsed: %d ms", t0.elapsed());
   // draw_dense_graph_aggregation_result();
   // group_anim->start();
  //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}

// ---------------------------------------- AGGREGATION WHITE RETAINING VERTEX -----------------------------
/** The absorbed vertices are now retained in the graph.
 * Type III.a - Highest Triangles Neighbour
 * Pr(v) = u.a.r
 * Selet u: arg max tri(u)
 * Stationary
 * @brief MainWindow::random_aggregate_retain_vertex_using_triangulation_and_weight_comparison
 */
void MainWindow::random_aggregate_retain_vertex_using_triangulation()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
       // Edge * e = selected->getProbabilisticTriangulationCoeffVertex();
        Edge * e = selected->correct_getMostMutualVertex();
        Vertex * neighbour, * winner, * loser;
        if (e->toVertex() == selected)
            neighbour = e->fromVertex();
        else
            neighbour = e->toVertex();

        winner = neighbour;
        loser = selected;
        //create the animation
        winner->absorb_retainEdge(e);
        hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
        players.removeOne(loser);
        t++;
    }
    qDebug("III.a - Time elapsed: %d ms", t0.elapsed());
    draw_aggregation_retain_vertex_result();

}


/** Type III.b - Probabilistic Triangle Neighbour
 * Select candidate u.a.r
 * Select neighbour with Pr = tri(u)/sum_tri(u)
 * @brief MainWindow::random_aggregate_retain_vertex_using_probabilistic_triangulation
 */
void MainWindow::random_aggregate_retain_vertex_using_probabilistic_triangulation()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }

    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);

        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
       // Edge * e = selected->getProbabilisticTriangulationCoeffVertex();
        Edge * e = selected->getProbabilisticTriangulationCoeffVertex();
        Vertex * neighbour, * winner, * loser;
        if (e->toVertex() == selected)
            neighbour = e->fromVertex();
        else
            neighbour = e->toVertex();

        winner = neighbour;
        loser = selected;
        //create the animation
        winner->absorb_retainEdge(e);
        hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
        players.removeOne(loser);
        t++;
    }
    qDebug("III.b - Time elapsed: %d ms", t0.elapsed());
    draw_aggregation_retain_vertex_result();
}


/** The absorbed vertices are now retained in the graph.
 * Type III.d - Probabilistic Triangles Emphaised Cluster
 * Pr(v) = u.a.r
 * f(u) = (tri(u)*2) * (extra_w(u)/no_absorbed(u))
 * Pr(u) = f(u) / sum f(i) forall i in adj(v)
 * @brief MainWindow::random_aggregate_retain_vertex_using_triangulation_and_weight_comparison
 */
void MainWindow::random_aggregate_retain_vertex_using_triangulation_times_weight()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
       // v->setWeight(v->getNumberEdge());
        v->setWeight(1);
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);

        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        Edge * e = selected->getProbabilisticTriangulationAndWeightVertex();
        if (e == 0)
        {
            selected->setParent(selected);
            players.removeOne(selected);
            continue;
        }
        Vertex * neighbour, * winner, * loser;
        if (e->toVertex() == selected)
            neighbour = e->fromVertex();
        else
            neighbour = e->toVertex();

        winner = neighbour;
        loser = selected;
        //create the animation
        winner->absorb_retainEdge(e);
        hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
        players.removeOne(loser);
        t++;
    }
    qDebug("III.d - Time elapsed: %d ms", t0.elapsed());
    draw_aggregation_retain_vertex_result();

}

/** Type III.e - Highest Tri(Cluster)
 * Pr(v) = u.a.r
 * LEt C(u) be the set of vertex in u cluster:
 * Select u: arg max tri(C(u)) forall u in adj(v)
 * @brief MainWindow::random_aggregate_retain_vertex_using_triangulation_of_cluster
 */
void MainWindow::random_aggregate_retain_vertex_using_triangulation_of_cluster()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    QTime t0;
    t0.start();

    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        Edge * e = selected->getHighestTriangulateCluster();
        Vertex * neighbour, * winner, * loser;
        if (e->toVertex() == selected)
            neighbour = e->fromVertex();
        else
            neighbour = e->toVertex();

        winner = neighbour;
        loser = selected;
        //create the animation
        winner->absorb_retainEdge(e);
        hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
        players.removeOne(loser);
        t++;
    }

    qDebug("III.e - Time elapsed: %d ms", t0.elapsed());
    draw_aggregation_retain_vertex_result();

}






// ------------------------------------- ANIMATION ---------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------

/*
 * PREPARE ANIMATION
 * @PARAM:
 * edges: first edge is the 'appearing' edge, rest are disappearing edges
 * abs: the vertex being removed.
 */
QSequentialAnimationGroup* MainWindow::prepare_animation_for_one_vertex(QList<int> edges, Vertex * abs, Vertex * absorber)
// FOR A SINGLE VERTEX
{
    QSequentialAnimationGroup * sequence = new QSequentialAnimationGroup;
    QParallelAnimationGroup * dis_anim = new QParallelAnimationGroup;
    QParallelAnimationGroup * last_anim = new QParallelAnimationGroup;
    QPointF abs_v = abs->pos();
    for (int i = 0; i < edges.size(); i++)
    {
        int index = edges.at(i);
        LineAnimator * e = myLine.at(index);
        QPointF from = e->line().p1(), to = e->line().p2();
        QPointF animFrom, animTo;
        if (i == 0) //appearing edge
        {
            if (abs_v.x() == from.x() || abs_v.y() == from.y())
            {
                e->setLine(QLineF(to, from));
                animFrom = from;
                animTo = to;
            }
            else
            {
                animFrom = to;
                animTo = from;
            }
            // moving vertex anim
            QPropertyAnimation * x_anim = new QPropertyAnimation(abs, "position");
            x_anim->setStartValue(animFrom);
            x_anim->setEndValue(animTo);
            x_anim->setDuration(1000);
            connect(x_anim, SIGNAL(finished()), abs, SLOT(disappear()));
            last_anim->addAnimation(x_anim);
            //disaeapring edge anim
            QPropertyAnimation * anim = new QPropertyAnimation(e, "endPoint");
            anim->setStartValue(animFrom);
            anim->setEndValue(animTo);
            anim->setDuration(1000);
            last_anim->addAnimation(anim);
        }
        else
        {
            QPropertyAnimation * anim = new QPropertyAnimation(e, "endPoint");
            //all disaapearing edges
            if (abs_v.x() == from.x() || abs_v.y() == from.y())
            {
                animFrom = to;
                animTo = from;
            }
            else
            {
                e->setLine(QLineF(to, from));
                animFrom = from;
                animTo = to;
            }
            anim->setStartValue(animFrom);
            anim->setEndValue(animTo);
            anim->setDuration(1000);
            connect(anim, SIGNAL(stateChanged(QAbstractAnimation::State,QAbstractAnimation::State)),
                    e, SLOT(highlight_edge()));
            dis_anim->addAnimation(anim);
        }
    }
    sequence->addAnimation(dis_anim);
    sequence->addAnimation(last_anim);
    connect(sequence, SIGNAL(stateChanged(QAbstractAnimation::State,QAbstractAnimation::State)),
            abs, SLOT(highlight_v()));
    connect(sequence, SIGNAL(stateChanged(QAbstractAnimation::State,QAbstractAnimation::State)),
            absorber, SLOT(highligh_absorber()));
    connect(sequence, SIGNAL(finished()), abs, SLOT(de_highlight_v()));
    connect(sequence, SIGNAL(finished()), absorber, SLOT(de_highlight_absorber()));

    sequence->addPause(300);
    return sequence;
}

void MainWindow::draw_hierarchy_tree()
{
    graphIsReady = false;
    if (centroids.empty())
    {
        qDebug() << "WINNER SET IS EMPTY";
        return;
    }
  //  preserveVertexPointer();
 //   myScene->clear();
    QList<ClusterCentroid*> C;
    QList<Arrow*> arrows;
    for (int i = 0; i < centroids.size(); i++)
    {
        Vertex * v = centroids.at(i);
        if (v->is_vertex_dragged_along() || v->getParent() != 0)
            continue;
        QPointF origin = v->getOriginPos();
        ClusterCentroid * centroid = new ClusterCentroid(origin.x(), origin.y());
        myScene->addItem(centroid);
        centroid->reSize(VERTEX_GLOBAL_SIZE);
        centroid->setIndex(v->getIndex());
        C.append(centroid);
        QList<Vertex*> absorbed = v->getAbsorbedList();
        QList<ClusterVertex*> clus;
        QList<int> ClusterIndexes;
        for (int j = 0; j < absorbed.size(); j++)
        {
            Vertex * v2 = absorbed.at(j);
            QPointF originPos = v2->getOriginPos();
            ClusterVertex * c_v = new ClusterVertex(originPos.x(), originPos.y());
            myScene->addItem(c_v);
            c_v->reSize(VERTEX_GLOBAL_SIZE);
            cluster_vertex.append(c_v);
            c_v->setIndex(v2->getIndex());
            clus.append(c_v);
            ClusterIndexes.append(c_v->getIndex());
        }
        centroid->addCluster(clus);
        //connecting parent with child
        for (int j = 0; j < clus.size(); j++)
        {
            Vertex * v2 = absorbed.at(j);
            Vertex * parent = v2->getParent();
            ClusterVertex * cv = clus.at(j);
            int pIndex = parent->getIndex();
            if (ClusterIndexes.contains(pIndex))
            {
                int id = ClusterIndexes.indexOf(pIndex);
                ClusterVertex * pv = clus.at(id);
                Arrow * arr = new Arrow (cv, pv);
                arrows.append(arr);
                myScene->addItem(arr);
            }
            else // points to the centroid
            {
                Arrow * arr = new Arrow(cv, centroid);
                arrows.append(arr);
                myScene->addItem(arr);
            }
        }
    }
}

void MainWindow::draw_walked_directed_graph()
{
    graphIsReady = false;
    preserveVertexPointer();
    myScene->clear();
    for(int i = 0 ; i < myVertexList.size(); i++)
        myScene->addItem(myVertexList[i]);

    for (int i = 0; i < hierarchy.size(); i++)
    {
        QPair<int,int> p = hierarchy.at(i);
        Vertex * from = myVertexList.at(p.first),
                * to = myVertexList.at(p.second);
        Arrow * arr = new Arrow(from, to);
        arr->changePen(Qt::black, 2);
        myScene->addItem(arr);
    }


    QList<QPair<int,int> > d_edges;
    for (int i = 0; i < hierarchy.size(); i++)
    {
        QPair<int,int> e = hierarchy.at(i);
        if (!d_edges.contains(e))
            d_edges.append(e);
    }

    using namespace boost;
    adjacency_list<vecS, vecS, directedS> G;
    adjacency_list<vecS, vecS, undirectedS> g;

    for (int i = 0; i < myVertexList.size(); i++)
    {
        add_vertex(G);
        add_vertex(g);
    }

    for (int i = 0; i < d_edges.size(); i++)
    {
        boost::add_edge(d_edges[i].first, d_edges[i].second, G);
    }

    for (int i = 0; i < hierarchy.size(); i++)
    {
        add_edge(hierarchy[i].first, hierarchy[i].second, g);
    }
    //
    typedef boost::graph_traits<adjacency_list<vecS, vecS, undirectedS> >::vertex_iterator ver_iter;
    for (std::pair<ver_iter, ver_iter> vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        qDebug() << *vp.first << ";" << boost::degree(*vp.first, g);
    }

    typedef graph_traits<adjacency_list<vecS, vecS, directedS> >::vertex_descriptor Vertex;

    std::vector<int> component(num_vertices(G)), discover_time(num_vertices(G));
    std::vector<default_color_type> color(num_vertices(G));
    std::vector<Vertex> root(num_vertices(G));
    int num = strong_components(G, make_iterator_property_map(component.begin(), get(vertex_index, G)),
                                root_map(make_iterator_property_map(root.begin(), get(vertex_index, G))).
                                color_map(make_iterator_property_map(color.begin(), get(vertex_index, G))).
                                discover_time_map(make_iterator_property_map(discover_time.begin(), get(vertex_index, G))));

    std::vector<int> ccomponent(boost::num_vertices(g));
    int numc = connected_components(g, &ccomponent[0]);

    std::cout << "Total number of Strongly Connected components: " << num << std::endl;
    std::cout << "Total number of Weakly Connected Components: " << numc << std::endl;

    /*
    QList<QColor> colorV;
    colorV << Qt::black	 <<  Qt::red
    << Qt::darkRed
    << Qt::green
    << Qt::darkGreen
    << Qt::blue
    << Qt::darkBlue
    << Qt::cyan
    << Qt::darkCyan
    << Qt::magenta
    << Qt::darkMagenta
    << Qt::yellow
    << Qt::darkYellow
    << Qt::gray
    << Qt::darkGray
    << Qt::lightGray
    << Qt::transparent	;

    std::vector<int>::size_type i;
      for (i = 0; i != component.size(); ++i)
      {
        std::cout << "Vertex " << myVertexList[i]->getIndex()
             <<" is in component " << component[i] << std::endl;

        if (num < colorV.size())
        {
            myVertexList[i]->setBackgroundColour(colorV.at(component[i]));
            myVertexList[i]->setDeselectedColour(colorV.at(component[i]));
        }
      }
    */
    graphIsReady = false;
}


/*
 * DRAW RESULT FROM ALL AGGREGATION USING THE WINNER SET
 */
void MainWindow::draw_aggregation_result()
{
    graphIsReady = false;
    if (centroids.empty())
    {
        qDebug() << "WINNER SET IS EMPTY";
        return;
    }
    if (!GRAPHICS)
    {
        parse_aggregation_result();
        return;
    }
    qDebug() << "- Aggregation Finished! Number of C: " << centroids.size();
    preserveVertexPointer();
    myScene->clear();
    QList<ClusterCentroid*> C;
    QList<Arrow*> arrows;
    //qDebug() << "C: " << centroids.size();
    for (int i = 0; i < centroids.size(); i++)
    {
        Vertex * v = centroids.at(i);
        if (v->is_vertex_dragged_along() || v->getParent() != 0)
            continue;
        QPointF origin = v->getOriginPos();
        ClusterCentroid * centroid = new ClusterCentroid(origin.x(), origin.y());
        centroid->reSize(VERTEX_GLOBAL_SIZE);
        centroid->setIndex(v->getIndex());
        C.append(centroid);
        QList<Vertex*> absorbed = v->getAbsorbedList();
        QList<ClusterVertex*> clus;
        QList<int> ClusterIndexes;
        for (int j = 0; j < absorbed.size(); j++)
        {
            Vertex * v2 = absorbed.at(j);
            QPointF originPos = v2->getOriginPos();
            ClusterVertex * c_v = new ClusterVertex(originPos.x(), originPos.y());
            c_v->reSize(VERTEX_GLOBAL_SIZE);
            cluster_vertex.append(c_v);
            c_v->setIndex(v2->getIndex());
            clus.append(c_v);
            ClusterIndexes.append(c_v->getIndex());
        }
        centroid->addCluster(clus);
        //connecting parent with child
        for (int j = 0; j < clus.size(); j++)
        {
            Vertex * v2 = absorbed.at(j);
            Vertex * parent = v2->getParent();
            ClusterVertex * cv = clus.at(j);
            int pIndex = parent->getIndex();
            if (ClusterIndexes.contains(pIndex))
            {
                int id = ClusterIndexes.indexOf(pIndex);
                ClusterVertex * pv = clus.at(id);
                Arrow * arr = new Arrow (cv, pv);
                arrows.append(arr);
                myScene->addItem(arr);
            }
            else // points to the centroid
            {
                Arrow * arr = new Arrow(cv, centroid);
                arrows.append(arr);
                myScene->addItem(arr);
            }
        }
    }

    //alternative to the arrow

    NormalGraph g = myGraph;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        LineAnimator * line = new LineAnimator(v1->pos(),v2->pos(), EDGE_GLOBAL_SIZE/2);
        myScene->addItem(line);
    }

    // begin drawing

    QList<QColor> color;
    color << Qt::GlobalColor::red << Qt::GlobalColor::blue << Qt::GlobalColor::green << Qt::GlobalColor::black <<
             Qt::GlobalColor::cyan << Qt::GlobalColor::white << Qt::GlobalColor::yellow << Qt::GlobalColor::gray <<
             Qt::GlobalColor::darkBlue << Qt::GlobalColor::darkGreen << Qt::GlobalColor::darkRed << Qt::GlobalColor::lightGray
          << Qt::GlobalColor::darkMagenta;
    for (int i = 0; i < C.size(); i++)
    {
        myScene->addItem(C[i]);
        if (C.size() < color.size())
        {
            C[i]->setBackgroundColour(color[i]);
            C[i]->setDeselectedColour(color[i]);
        }
        QList<ClusterVertex*> cluster = C[i]->getCluster();
        for (int j = 0; j < cluster.size(); j++)
        {
            myScene->addItem(cluster[j]);
            if (C.size() < color.size())
            {
                cluster[j]->setBackgroundColour(color[i]);
                cluster[j]->setDeselectedColour(color[i]);
            }
        }
    }
    // prepare for cluster matching;
    QList<QList<int> > clusters;
    for (int i = 0; i < C.size(); i++)
    {
        QList<int> c;
        c.append(C[i]->getIndex());
        QList<ClusterVertex*> cluster = C[i]->getCluster();
        for (int j = 0; j < cluster.size(); j++)
        {
            c.append(cluster[j]->getIndex());
        }
        clusters.append(c);
    }
    result = clusters;
}

void MainWindow::draw_aggregation_result_by_shape()
{
    graphIsReady = false;
    if (centroids.empty())
    {
        qDebug() << "WINNER SET IS EMPTY";
        return;
    }

    QList<ClusterCentroid*> C;
 //   qDebug() << "C: " << centroids.size();
    for (int i = 0; i < centroids.size(); i++)
    {
        Vertex * v = centroids.at(i);
        if (v->is_vertex_dragged_along() || v->getParent() != 0)
            continue;
        QPointF origin = v->getOriginPos();
        ClusterCentroid * centroid = new ClusterCentroid(origin.x(), origin.y());
        centroid->reSize(VERTEX_GLOBAL_SIZE);
        centroid->setIndex(v->getIndex());
        C.append(centroid);
        QList<Vertex*> absorbed = v->getAbsorbedList();
        QList<ClusterVertex*> clus;
        QList<int> ClusterIndexes;
        for (int j = 0; j < absorbed.size(); j++)
        {
            Vertex * v2 = absorbed.at(j);
            QPointF originPos = v2->getOriginPos();
            ClusterVertex * c_v = new ClusterVertex(originPos.x(), originPos.y());
            c_v->reSize(VERTEX_GLOBAL_SIZE);
            cluster_vertex.append(c_v);
            c_v->setIndex(v2->getIndex());
            clus.append(c_v);
            ClusterIndexes.append(c_v->getIndex());
        }
        centroid->addCluster(clus);
    }
    preserveVertexPointer();
    myScene->clear();
    //readding edges
    for (int i = 0; i < myLine.size(); i++)
    {
        myLine[i]->dehighlight_edge();
        myLine[i]->setSize(EDGE_GLOBAL_SIZE);
        myScene->addItem(myLine[i]);
    }
    //drawing shapes
    // 1 for rectangle, 2 triangle, 3 circular (normal)
    QList<QColor> color;
    color << Qt::green << Qt::red << Qt::blue << Qt::magenta;
    if (C.size() < 4)
    {
        for (int i = 0; i < C.size(); i++)
        {
            QGraphicsPolygonItem * rect = new QGraphicsPolygonItem;
            rect->setPen(QPen(Qt::black, 10, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
            rect->setBrush(QBrush(color[i]));
            ClusterCentroid * Cen = C.at(i);
            if ( i==0 )
            {
                C[i]->setBackgroundColour(Qt::green);
                myScene->addItem(C[i]);
                QList<ClusterVertex*> cluster = C[i]->getCluster();
                for (int j = 0; j < cluster.size(); j++)
                {
                    cluster[j]->setBackgroundColour(Qt::green);
                    myScene->addItem(cluster[j]);
                }
            }
            else if (i == 1)
            {

                QPolygonF poly;
                QRectF Crect = Cen->boundingRect();
                Crect.translate(Cen->pos());
                poly << Crect.topLeft() << Crect.topRight() << Crect.bottomRight() << Crect.bottomLeft();
                rect->setPolygon(poly);
                myScene->addItem(rect);
                QList<ClusterVertex*> cluster = C[i]->getCluster();
                for (int j = 0; j < cluster.size(); j++)
                {
                    QGraphicsPolygonItem * sub_rect = new QGraphicsPolygonItem;
                    sub_rect->setPen(QPen(Qt::black, 10, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
                    sub_rect->setBrush(QBrush(color[i]));
                    QPolygonF sub_poly;
                    QRectF Srect = cluster[j]->boundingRect();
                    Srect.translate(cluster[j]->pos());
                    sub_poly << Srect.topLeft() << Srect.topRight() << Srect.bottomRight() << Srect.bottomLeft();
                    sub_rect->setPolygon(sub_poly);
                    myScene->addItem(sub_rect);
                }
            }
            else if ( i == 2)
            {
                QPolygonF poly;
                QRectF Crect = Cen->boundingRect();
                Crect.translate(Cen->pos());
                QPointF trian_top = (Crect.topLeft() + Crect.topRight())/2;
                poly << trian_top << Crect.bottomRight() << Crect.bottomLeft();
                rect->setPolygon(poly);
                myScene->addItem(rect);
                QList<ClusterVertex*> cluster = C[i]->getCluster();
                for (int j = 0; j < cluster.size(); j++)
                {
                    QGraphicsPolygonItem * sub_rect = new QGraphicsPolygonItem;
                    sub_rect->setPen(QPen(Qt::black, 10, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
                    sub_rect->setBrush(QBrush(color[i]));
                    QRectF Srect = cluster[j]->boundingRect();
                    Srect.translate(cluster[j]->pos());
                    QPolygonF sub_poly;
                    QPointF sub_trian_top = (Srect.topLeft() + Srect.topRight())/2;
                    sub_poly << sub_trian_top << Srect.bottomRight() << Srect.bottomLeft();
                    sub_rect->setPolygon(sub_poly);
                    myScene->addItem(sub_rect);
                }
            }
        }
    }

}



void MainWindow::draw_dense_graph_aggregation_result()
{
    graphIsReady = false;
    qDebug() << centroids.size();
    QList<QColor> color;
    color << Qt::GlobalColor::red << Qt::GlobalColor::blue << Qt::GlobalColor::green << Qt::GlobalColor::black <<
             Qt::GlobalColor::cyan << Qt::GlobalColor::white << Qt::GlobalColor::yellow << Qt::GlobalColor::gray <<
             Qt::GlobalColor::darkBlue << Qt::GlobalColor::darkGreen << Qt::GlobalColor::darkRed << Qt::GlobalColor::lightGray
          << Qt::GlobalColor::darkMagenta;

    if (centroids.size() < color.size())
    {
        QList<ClusterCentroid*> C;
        for (int i = 0; i < centroids.size(); i++)
        {
            QColor cor = color.at(i);
            Vertex * v = centroids.at(i);
            QPointF origin = v->pos();
            ClusterCentroid * centroid = new ClusterCentroid(origin.x(), origin.y());
            centroid->setBackgroundColour(cor);
            centroid->setIndex(v->getIndex());
            C.append(centroid);
            QList<Vertex*> absorbed = v->getAbsorbedList();
            QList<ClusterVertex*> clus;
            QList<int> ClusterIndexes;
            for (int j = 0; j < absorbed.size(); j++)
            {
                Vertex * v2 = absorbed.at(j);
                QPointF originPos = v2->pos();
                ClusterVertex * c_v = new ClusterVertex(originPos.x(), originPos.y());
                cluster_vertex.append(c_v);
                c_v->setIndex(v2->getIndex());
                c_v->setBackgroundColour(cor);
                clus.append(c_v);
                ClusterIndexes.append(c_v->getIndex());
            }
            centroid->addCluster(clus);
        }
        // begin drawing
        myScene->clear();
        for (int i = 0; i < C.size(); i++)
        {
            myScene->addItem(C[i]);
            QList<ClusterVertex*> cluster = C[i]->getCluster();
            for (int j = 0; j < cluster.size(); j++)
            {
                myScene->addItem(cluster[j]);
            }
        }
        qDebug() << "HERE";
    }
}

/** JUST DRAW THE HIERARCHY TREE
 * @brief MainWindow::draw_aggregation_retain_vertex_result
 */
void MainWindow::draw_aggregation_retain_vertex_result()
{
    graphIsReady = false;
    if (!GRAPHICS)
    {
        parse_retain_result();
    }
    preserveVertexPointer();
    myScene->clear();
    NormalGraph g;
    for(int i = 0; i < myVertexList.size(); i++)
        boost::add_vertex(g);
    for (int i = 0; i < hierarchy.size(); i++)
    {
        QPair<int,int> p = hierarchy.at(i);
        boost::add_edge(p.first, p.second, g);
    }

    std::vector<int> component(boost::num_vertices(g));
    int num = boost::connected_components(g, &component[0]);
    qDebug() << "- Aggregation Finished! Number of C: " << num;
    std::vector<int>::size_type i;


    QList<QColor> color;
    color << Qt::black	 <<  Qt::red
    << Qt::darkRed
    << Qt::green
    << Qt::darkGreen
    << Qt::blue
    << Qt::darkBlue
    << Qt::cyan
    << Qt::darkCyan
    << Qt::magenta
    << Qt::darkMagenta
    << Qt::yellow
    << Qt::darkYellow
    << Qt::gray
    << Qt::darkGray
    << Qt::lightGray
    << Qt::transparent	;

    if (num < color.size())
    {
        for (i = 0; i != component.size(); ++i)
        {
            Vertex * v = myVertexList.at(i);
          //  myScene->addItem(v);
            v->setBackgroundColour(color.at(component[i]));
            v->setDeselectedColour(color.at(component[i]));
            v->reSize(10);
            myScene->addItem(v);
        }
    }


    for (int i = 0; i < hierarchy.size(); i++)
    {
        QPair<int,int> p = hierarchy.at(i);
        Vertex * v1 = myVertexList.at(p.first);
        Vertex * v2 = myVertexList.at(p.second);
        Arrow * arr = new Arrow(v1,v2);
        myScene->addItem(arr);
    }

    calculate_modularity_for_clusters();

    //compute indices
    QList<QList<int> > clusters;
    for (int i = 0 ; i < num; i++)
    {
        QList<int> c;
        clusters.append(c);
    }

    for (int i = 0; i < component.size(); i++)
    {
        QList<int> c = clusters[component[i]];
        c.append(i);
        clusters.replace(component[i],c);
    }
    result = clusters;
}


void MainWindow::preserveVertexPointer()
{
    QList<QGraphicsItem*> items = myScene->items();
    QMutableListIterator<QGraphicsItem*> i(items);
    while(i.hasNext())
    {
        Vertex * v = dynamic_cast<Vertex*>(i.next());
        if (v)
            myScene->removeItem(v);
    }
}

// ------------------------- CREATE ACTIONS AND MENUS -------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

void MainWindow::createActions()
{
    generateLayerGnpAction = new QAction(QString("Generate Layered Gnp"), this);
    connect(generateLayerGnpAction ,SIGNAL(triggered()), this, SLOT(generateLayerGnp()));

    generateCycleAction = new QAction(QString("Generate Simple Cycle"), this);
    connect(generateCycleAction, SIGNAL(triggered()), this, SLOT(generateCycle()));

    generateSimpleClusterGnpAction = new QAction(QString("Generate Simple Cluster-Gnp"), this);
    connect(generateSimpleClusterGnpAction, SIGNAL(triggered()), this, SLOT(generateSimpleClusterGnp()));

    generateArtificialSocialNetworkAction = new QAction(QString("Generate Artificial Clusters"), this);
    connect(generateArtificialSocialNetworkAction, SIGNAL(triggered()), this, SLOT(generateArtificialSocialNetwork()));

    generateGNArtificialNetworkAction = new QAction(QString("Generate Girvan&Newman Artificial Network"), this);
    connect(generateGNArtificialNetworkAction, SIGNAL(triggered()), this, SLOT(generateGNArtificialNetwork()));

    generateErdosReyniAction = new QAction(QString("Generate Erdos Reyni"), this);
    connect(generateErdosReyniAction,SIGNAL(triggered()), this, SLOT(generateErdosReyni())) ;

    generateKarateClubAction = new QAction(QString("Generate KarateClub"), this);
    connect(generateKarateClubAction, SIGNAL(triggered()), this, SLOT(generateKarateClub()));

    generateWeightedKarateClubAction = new QAction(QString("Generate Weighted Karate club"), this);
    connect(generateWeightedKarateClubAction, SIGNAL(triggered()), this, SLOT(generateWeightedKarateClub()));

    generateSmallSampleAction = new QAction(QString("Generate Small Sample"), this);
    connect(generateSmallSampleAction, SIGNAL(triggered()), this, SLOT(generateSimpleSample()));

    readInputAction = new QAction(QString("Read Input(EDGE File Only)"), this);
    connect(readInputAction, SIGNAL(triggered()), this, SLOT(read_seperated_graph_input()));

    readSNAPAction = new QAction(QString("Read Input for SNAP"), this);
    connect(readSNAPAction, SIGNAL(triggered()), this, SLOT(read_large_graph_with_ground_truth_communities()));

    loadGroundTruthCommAction = new QAction(QString("Load Ground Truth"), this);
    connect(loadGroundTruthCommAction, SIGNAL(triggered()), this, SLOT(load_ground_truth_communities()));

    readGMLFileAction = new QAction(QString("Read GML File"), this);
    connect(readGMLFileAction, SIGNAL(triggered()), this, SLOT(parseGMLfile()));

    comparesKarateClubResultAction = new QAction(QString("Compare Karate Club Result"), this);
    connect(comparesKarateClubResultAction, SIGNAL(triggered()), this, SLOT(compareResultWithKarateClubData()));

    BCAction = new QAction(QString("Betweenness Centrality Clustering"), this);
    connect(BCAction, SIGNAL(triggered()), this, SLOT(betweenness_centrality_clutering()));

    randomAggregateAction = new QAction(QString("Type I.a - Uniform"), this);
    connect(randomAggregateAction, SIGNAL(triggered()), this, SLOT(random_aggregate()));

    randomAggregationNStepsWalkAction = new QAction(QString("Uniformly N Walks Aggregate"), this);
    connect(randomAggregationNStepsWalkAction, SIGNAL(triggered()), this, SLOT(random_aggregate_n_steps_walk()));

    randomAggregationWithDegreeComparisonAction = new QAction(QString("Type I.b - Uniform and Comparing CURRENT Degree"), this);
    connect(randomAggregationWithDegreeComparisonAction, SIGNAL(triggered()), this, SLOT(random_aggregate_with_degree_comparison()));

    randomAggregationWithWeightComparisonAction = new QAction(QString("Type I.c - Uniform and Comparing ORIGINAL Degree"), this);
    connect(randomAggregationWithWeightComparisonAction, SIGNAL(triggered()), this, SLOT(random_aggregate_with_weight_comparison()));

    randomAggregateWithInitialDegreeBiasAction = new QAction(QString("Type II.a - Uniform v, Neighbour W Biased"), this);
    connect(randomAggregateWithInitialDegreeBiasAction, SIGNAL(triggered()), this, SLOT(random_aggregate_with_neighbour_initial_degree_bias()));

    randomAggregationII_a_i_Action = new QAction(QString("Type II.a_(i) - Uniform v, Neighbour's W Biased, Comparing W"), this);
    connect(randomAggregationII_a_i_Action, SIGNAL(triggered()), this, SLOT(random_aggregate_with_neighbour_initial_degree_bias_with_comparison()));

    randomAggregationSelectHighestEdgeWeightWithVertexWeightComparisonAction = new QAction(QString("Edge-Weighted Graph I"), this);
    connect(randomAggregationSelectHighestEdgeWeightWithVertexWeightComparisonAction, SIGNAL(triggered()),
            this, SLOT(random_aggregate_with_highest_edge_weight_and_weight_comparison()));

    randomAggregationWithEdgeWeightBiasAndVertexWeightComparisonAction = new QAction(QString("Edge-Weighted Graph II"), this);
    connect(randomAggregationWithEdgeWeightBiasAndVertexWeightComparisonAction, SIGNAL(triggered()),
            this, SLOT(random_aggregate_with_edge_weight_bias_and_weight_comparison()));

    randomAggregationWithHighestTriangulatedVertexAndWeightComparisonAction = new QAction(QString("Type III.c - (Destructive) Highest Tri Neighbour, Constraint by Weight Comparison"), this);
    connect(randomAggregationWithHighestTriangulatedVertexAndWeightComparisonAction, SIGNAL(triggered()),
            this, SLOT(random_aggregate_with_highest_triangulated_vertex()));

    randomAggregationRetainingVertexWithTriangulationCoeffAndWeightComparisonAction = new QAction(QString("Type III.a - (Stationary) Highest Tri Neighbour"), this);
    connect(randomAggregationRetainingVertexWithTriangulationCoeffAndWeightComparisonAction, SIGNAL(triggered()),
            this, SLOT(random_aggregate_retain_vertex_using_triangulation()));

    randomAggregationProbabilisticTriangulationAction = new QAction(QString("Type III.b - (Stationary) Probabilisitc Triangles Neighbour"), this);
    connect(randomAggregationProbabilisticTriangulationAction, SIGNAL(triggered()),
            this, SLOT(random_aggregate_retain_vertex_using_probabilistic_triangulation()));

    randomAggregationRetainingVertexWithTriangulationxWeightIndexAction = new QAction(QString("Type III.d - (Stationary) PRobabilistic: Tri*2*(average(extra_w)"), this);
    connect(randomAggregationRetainingVertexWithTriangulationxWeightIndexAction, SIGNAL(triggered()),
            this, SLOT(random_aggregate_retain_vertex_using_triangulation_times_weight()));

    randomAggregationRetainingVertexWithHighestTrianglesClusterAction = new QAction(QString("Type III.e - (Stationary) Highest Tri(Cluster)"), this);
    connect(randomAggregationRetainingVertexWithHighestTrianglesClusterAction, SIGNAL(triggered()),
            this, SLOT(random_aggregate_retain_vertex_using_triangulation_of_cluster()));

    randomAggregationProbabilisticWeightSelectionWithMinimumDegreeNeighbourAction = new QAction(QString("Type II.f - (Retentive) Probabilistic Lowest Degree Neighbour"), this);
    connect(randomAggregationProbabilisticWeightSelectionWithMinimumDegreeNeighbourAction, SIGNAL(triggered()),
            this, SLOT(random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour()));

    randomAggregationWithMinimumDegreeNeighbourAction = new QAction(QString("Type II.d - (Destructive) Lowest Degree Neighbour"), this);
    connect(randomAggregationWithMinimumDegreeNeighbourAction, SIGNAL(triggered()), this, SLOT(random_aggregate_with_minimum_weight_neighbour()));

    randomAggregateMaxDegreeAction = new QAction(QString("Type II.g - (Destructive) Max Degree"), this);
    connect(randomAggregateMaxDegreeAction, SIGNAL(triggered()), this, SLOT(random_aggregate_greedy_max_degree()));

    randomAggregateMaxWeightAction = new QAction(QString("Type II.h - (Retentive) Max Weight"), this);
    connect(randomAggregateMaxWeightAction, SIGNAL(triggered()), this, SLOT(random_aggregate_greedy_max_weight()));

    randomWalkingRemovingEdgesAction = new QAction(QString("Random Walking Removing Edges"), this);
    connect(randomWalkingRemovingEdgesAction, SIGNAL(triggered()), this, SLOT(random_walking_constant_restart_remove_edges()));

    randomWalkingNormalAction = new QAction(QString("Random Walking Normal"), this);
    connect(randomWalkingNormalAction, SIGNAL(triggered()), this, SLOT(random_walking_normal_retain_edges()));

    randomWalkingTestAction = new QAction(QString("Random Walking Testing"), this);
    connect(randomWalkingTestAction, SIGNAL(triggered()), this, SLOT(random_walking_testing()));

    randomEdgeRemovalAction = new QAction(QString("Random Removal of Edges"), this);
    connect(randomEdgeRemovalAction, SIGNAL(triggered()), this, SLOT(random_edge_removal()));

    typeIretainAction = new QAction(QString("Type I Aggregation Retain Vertices"), this);
    connect(typeIretainAction, SIGNAL(triggered()), this, SLOT(typeI_retaining_vertex()));

    typeIIretainAction = new QAction(QString("Type II Aggregation Retain Vertices"), this);
    connect(typeIIretainAction, SIGNAL(triggered()), this, SLOT(typeII_retaining_vertex()));

    typeIIIretainAction = new QAction(QString("Type III Aggregation Retain Vertices"), this);
    connect(typeIIIretainAction, SIGNAL(triggered()), this, SLOT(typeIIIc_retaining_vertex()));

    typeIIIdremoveAction= new QAction(QString("Type III(d) Aggregation Remove Vertices, P Disconnect Edge"), this);
    connect(typeIIIdremoveAction, SIGNAL(triggered()), this, SLOT(typeIIId_removing_vertex()));

    typeIVretainAction = new QAction(QString("Type IV Aggregation Retain Vertices"), this);
    connect(typeIVretainAction, SIGNAL(triggered()), this, SLOT(typeIV_retaining_vertex()));
    // ------------------------------- TESTING --------------
    multipleAggregateAction = new QAction(QString("Multirun of Aggregation"), this);
    connect(multipleAggregateAction, SIGNAL(triggered()), this, SLOT(multirun_aggregate()));

    multipleAggregateEdgeWeightedGraphAction = new QAction(QString("Multirun of Aggregation on Edge-Weighted Graph"), this);
    connect(multipleAggregateEdgeWeightedGraphAction, SIGNAL(triggered()), this, SLOT(multirun_aggregate_edge_weight_graph()));

    randomAggregationWithNeighbourCurrentDegreeBiasedAction = new QAction(QString("Type II.b - Aggregate Neighbour with d Biased"), this);
    connect(randomAggregationWithNeighbourCurrentDegreeBiasedAction, SIGNAL(triggered()), this, SLOT(random_aggregate_with_neighbour_CURRENT_degree_bias()));

    randomAggregationII_b_i_Action = new QAction(QString("Type II.b_(i) - Aggregate Neighbour with d Biased, with comparison"), this);
    connect(randomAggregationII_b_i_Action, SIGNAL(triggered()), this, SLOT(random_aggregate_with_neighbour_CURRENT_degree_bias_with_comparison()));

    randomAggregateHighestDegreeNeighbourAction = new QAction(QString("Type II.c - Aggregate Highest Neighbour (compare Degree)"), this);
    connect(randomAggregateHighestDegreeNeighbourAction, SIGNAL(triggered()), this, SLOT(random_aggregate_highest_CURRENT_degree_neighbour()));

    randomAggregateLowestDegreeNeighbourDestructiveAction = new QAction(QString("Type II.e - (Destructive) Probabilisic Lowest Degree Neighbour"), this);
    connect(randomAggregateLowestDegreeNeighbourDestructiveAction, SIGNAL(triggered()), this, SLOT(random_aggregate_probabilistic_lowest_degree_neighbour_destructive()));

    resetAction = new QAction(QString ("Reset"), this);
    connect(resetAction, SIGNAL(triggered()), this, SLOT(reset()));

    saveGMLAction = new QAction(QString("Save File as GML"), this);
    connect(saveGMLAction, SIGNAL(triggered()), this, SLOT(exportGML()));

    examineBridgeAction = new QAction(QString("Examining Bridges"), this);
    connect(examineBridgeAction, SIGNAL(triggered()), this, SLOT(get_bridge_stats()));

    //large graph
    reindexingSNAPAction = new QAction(QString("Import and Reindexing SNAP Graphs"), this);
    connect(reindexingSNAPAction, SIGNAL(triggered()), this, SLOT(read_large_graph_with_ground_truth_communities()));

    //experiment
    GNGraphExperimentAction = new QAction(QString("Girvan Newman Graph Experiment"), this);
    connect(GNGraphExperimentAction, SIGNAL(triggered()), this, SLOT(GN_experiment()));

    calculateIndicesAction = new QAction(QString("Calculate Ground Truth Pairwise Matching"), this);
    connect(calculateIndicesAction, SIGNAL(triggered()), this, SLOT(calculate_indices()));

    calculateQAction = new QAction(QString("Calculate Q"), this);
    connect(calculateQAction, SIGNAL(triggered()), this, SLOT(calculate_Q()));

}

void MainWindow::createMenus()
{
    generateGraphMenu = menuBar()->addMenu(QString("Generating Graphs (DO THIS 1st)!"));
    generateGraphMenu->addAction(generateLayerGnpAction);
    generateGraphMenu->addAction(generateCycleAction);
    generateGraphMenu->addAction(generateSimpleClusterGnpAction);
    generateGraphMenu->addAction(generateArtificialSocialNetworkAction);
    generateGraphMenu->addAction(generateGNArtificialNetworkAction);
    generateGraphMenu->addAction(generateErdosReyniAction);
    generateGraphMenu->addAction(generateKarateClubAction);
    generateGraphMenu->addAction(generateWeightedKarateClubAction);
    generateGraphMenu->addAction(generateSmallSampleAction);
    generateGraphMenu->addSeparator();
    generateGraphMenu->addAction(readInputAction);
    generateGraphMenu->addAction(readSNAPAction);
    generateGraphMenu->addAction(readGMLFileAction);
    generateGraphMenu->addAction(loadGroundTruthCommAction);
    generateGraphMenu->addSeparator();
    generateGraphMenu->addAction(saveGMLAction);

    aggregateMenu = menuBar()->addMenu(QString("Aggregate (DO THIS 2nd)"));
    aggregateMenu->addAction(BCAction);
    aggregateMenu->addSeparator();
    aggregateMenu->addAction(randomAggregateAction);
    aggregateMenu->addAction(randomAggregationWithDegreeComparisonAction);
    aggregateMenu->addAction(randomAggregationWithWeightComparisonAction);
    aggregateMenu->addAction(randomAggregateWithInitialDegreeBiasAction);
    aggregateMenu->addAction(randomAggregationII_a_i_Action);
    aggregateMenu->addAction(randomAggregationWithNeighbourCurrentDegreeBiasedAction);
    aggregateMenu->addAction(randomAggregationII_b_i_Action);
    aggregateMenu->addAction(randomAggregateHighestDegreeNeighbourAction);
    aggregateMenu->addAction(randomAggregationWithMinimumDegreeNeighbourAction);
    aggregateMenu->addAction(randomAggregateLowestDegreeNeighbourDestructiveAction);
    aggregateMenu->addAction(randomAggregationProbabilisticWeightSelectionWithMinimumDegreeNeighbourAction);
    aggregateMenu->addAction(randomAggregateMaxDegreeAction);
    aggregateMenu->addAction(randomAggregateMaxWeightAction);
 //   aggregateMenu->addAction(randomAggregationWithHighestTriangulatedVertexAndWeightComparisonAction);
 //   aggregateMenu->addSeparator(); //edge weighted graph, removed for now
 //   aggregateMenu->addAction(randomAggregationSelectHighestEdgeWeightWithVertexWeightComparisonAction);
 //   aggregateMenu->addAction(randomAggregationWithEdgeWeightBiasAndVertexWeightComparisonAction);
    aggregateMenu->addSeparator();
    aggregateMenu->addAction(randomAggregationRetainingVertexWithTriangulationCoeffAndWeightComparisonAction);
    aggregateMenu->addAction(randomAggregationProbabilisticTriangulationAction);
    aggregateMenu->addAction(randomAggregationWithHighestTriangulatedVertexAndWeightComparisonAction);
    aggregateMenu->addAction(randomAggregationRetainingVertexWithTriangulationxWeightIndexAction);
    aggregateMenu->addAction(randomAggregationRetainingVertexWithHighestTrianglesClusterAction);

    aggregateMenu->addSeparator();
    aggregateMenu->addAction(comparesKarateClubResultAction);

    /*
    QMenu * rwMenu = menuBar()->addMenu(QString("Random Walking Based"));
    rwMenu->addAction(randomWalkingRemovingEdgesAction);
    rwMenu->addAction(randomWalkingNormalAction);
    rwMenu->addSeparator();
    rwMenu->addAction(randomEdgeRemovalAction);
    rwMenu->addSeparator();
    rwMenu->addAction(randomWalkingTestAction);

    QMenu * retainAggregationMenu = menuBar()->addMenu(QString("Aggregation Retain Vertices"));
    retainAggregationMenu->addAction(typeIretainAction);
    retainAggregationMenu->addAction(typeIIretainAction);
    retainAggregationMenu->addAction(typeIIIretainAction);
    retainAggregationMenu->addAction(typeIIIdremoveAction);
    retainAggregationMenu->addAction(typeIVretainAction);*/

    QMenu * testMenu = menuBar()->addMenu(QString("Multi Runs"));
   // testMenu->addAction(multipleAggregateAction);
   // testMenu->addAction(multipleAggregateEdgeWeightedGraphAction);
    testMenu->addAction(resetAction);
    testMenu->addSeparator();
    testMenu->addAction(GNGraphExperimentAction);

    QMenu * experimentMenu = menuBar()->addMenu(QString("Calculate Result"));
    experimentMenu->addAction(calculateIndicesAction);
    experimentMenu->addAction(calculateQAction);

    QMenu * SNAPMenu = menuBar()->addMenu(QString("SNAP"));
    SNAPMenu->addAction(reindexingSNAPAction);
}

void MainWindow::setUpGraphInfoWidget()
{
    if (myInfoWidget != 0) // if one is up already
        return;
    int g_num_v = boost::num_vertices(myGraph);
    int g_num_e = boost::num_edges(myGraph);
    QDockWidget * dock = new QDockWidget(QString("Info:"),this);
    InfoWidget * widget = new InfoWidget(dock, g_num_v, g_num_e);
    widget->setAttribute(Qt::WA_DeleteOnClose, true);
    int w = widget->geometry().width(), h = widget->geometry().height();
    dock->setWidget(widget);
    dock->setAllowedAreas(Qt::RightDockWidgetArea | Qt::LeftDockWidgetArea);
    dock->setMinimumWidth(w);
    dock->setMinimumHeight(h);
    dock->setFloating(false);
    dock->setFeatures(QDockWidget::NoDockWidgetFeatures);
    this->addDockWidget(Qt::LeftDockWidgetArea, dock);
    myInfoWidget = widget;
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        connect(v, SIGNAL(vertex_selected(QList<int>)), widget, SLOT(display_vertex_selected_info(QList<int>)));
    }
}

void MainWindow::createSeperateScene()
{
    QMainWindow * seperate = new QMainWindow;
    seperate->setAttribute(Qt::WA_DeleteOnClose);
    QGraphicsScene * newScene = new QGraphicsScene(0,0,10000,10000);
    GraphView * newView = new GraphView;

    newView->setScene(newScene);
    newView->setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    newView->setContextMenuPolicy(Qt::ActionsContextMenu);
    seperate->setCentralWidget(newView);
    seperate->setWindowTitle(QString("Origin Graph"));
    //set up Actions and Menus bars/options
    int size = 150, text_size = size*3/5;
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        QGraphicsEllipseItem * rect = new QGraphicsEllipseItem(v->x()-size/2,v->y()-size/2,
                                                           size,size);
        rect->setPen(QPen(Qt::black, 10, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        rect->setBrush(QBrush(Qt::green));
        rect->setZValue(1);
        newScene->addItem(rect);

        QGraphicsTextItem * index = new QGraphicsTextItem;
        index->setZValue(2);
        QFont font;
        font.setPixelSize(text_size);
        index->setFont(font);
        index->setPos(v->getOriginPos().x()-text_size*3/5,v->getOriginPos().y()-text_size*3/5);
        index->setDefaultTextColor(Qt::red);
        index->setPlainText(QString::number(v->getIndex()));
        newScene->addItem(index);
    }

    for (int i = 0; i < myLine.size(); i++)
    {
        LineAnimator * line = myLine.at(i);
        QLineF l = line->line();
        QGraphicsLineItem * newline = new QGraphicsLineItem(l);
        newline->setPen(QPen(Qt::black, 20, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        newline->setZValue(0);
        newScene->addItem(newline);
    }

    seperate->show();
    originGraph = seperate;
}

void MainWindow::resetGraphics()
{
    myScene->clear();
    myLine.clear();
    myEdgeList.clear();
    myVertexList.clear();
    centroids.clear();
    cluster_vertex.clear();
}

void MainWindow::saveAggregationResultInSeperateWindow()
{
    QMainWindow * seperate = new QMainWindow;
    seperate->setAttribute(Qt::WA_DeleteOnClose);
    seperate->setWindowTitle("Result");
    QGraphicsScene * newScene = new QGraphicsScene(0,0,10000,10000);
    GraphView * newView = new GraphView;
    newView->setScene(newScene);
    newView->setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    newView->setContextMenuPolicy(Qt::ActionsContextMenu);
    seperate->setCentralWidget(newView);

    QList<Vertex*> newV;
    for(int i = 0; i < myVertexList.size(); i++)
        newV.append(0);

    QList<QGraphicsItem*> items = myScene->items();
    QMutableListIterator<QGraphicsItem *> i(items);
    while (i.hasNext())
    {
        Vertex * v = dynamic_cast<Vertex*>(i.next());
        if (v)
        {
            int index = v->getIndex();
            QPointF p = v->pos();
            Vertex * copy = new Vertex;
            if (dynamic_cast<ClusterCentroid*>(v))
                copy->setBackgroundColour(Qt::red);
            newV.replace(index,copy);
            copy->setPos(p);
            copy->reSize(VERTEX_GLOBAL_SIZE);
            newScene->addItem(copy);
        }
    }
    i.toFront();
    while(i.hasNext())
    {
        Arrow * arr = dynamic_cast<Arrow*>(i.next());
        if(arr)//if object is a vertex
        {
            Vertex * from = arr->fromVertex();
            Vertex * to = arr->toVertex();
            int fromindex = from->getIndex();
            int toindex = to->getIndex();
            Vertex * newFrom = newV.at(fromindex);
            Vertex * newTo = newV.at(toindex);
            Arrow * newArr = new Arrow(newFrom, newTo);
            newScene->addItem(newArr);
        }
    }
    seperate->show();
}

void MainWindow::compareResultWithKarateClubData()
{
    graphIsReady = false;
    QList<int> actual_result;
    actual_result << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 <<33
                  << 0 << 0 << 0 << 0 << 33 << 33 << 0 << 0  << 33 << 0
                  << 33 << 0 << 33 << 33 << 33 << 33 << 33 << 33 << 33 << 33 << 33
                  << 33 << 33 << 33;
    if (centroids.size() == 0)
    {
        qDebug() << "HAVE NOT BEEN CLUSTERED";
        return;
    }

    QList<int> strict;
    for (int i = 0; i < 34; i++)
        strict.append(-1);

    for (int i = 0; i < centroids.size(); i++)
    {
        Vertex * v = centroids.at(i);
        if (v == 0)
            qDebug() << "NULL POINTER CENTROID";
        if (v->is_vertex_dragged_along())
            continue;
        int centroid_index = v->getIndex();

        strict.replace(centroid_index, centroid_index);
        QList<Vertex*> absorbed = v->getAbsorbedList();
        for (int j = 0; j < absorbed.size(); j++)
        {
            Vertex * v2 = absorbed.at(j);
            int cluster_index = v2->getIndex();
            strict.replace(cluster_index, centroid_index);
        }
    }
    //compare
/*
    qDebug() << "INDEX /t ACTUAL DATA CENTROID /t AGGREGATED RESULT CENTROID /t HIT/MISS";
    int hitc = 0, miss = 0;
    for (int i = 0; i < actual_result.size(); i++)
    {
        int actual = actual_result.at(i);
        int aggregated = strict.at(i);
        bool hit = (actual == aggregated);
        if (hit)
        {
            qDebug() << i << "/t" << actual << "/t" << aggregated << "(HIT)";
            hitc++;
        }
        else
        {
            qDebug() << i << "/t" << actual << "/t" << aggregated << "(MISS)";
            miss++;
        }
    }
    qDebug() << "TOTAL: HIT-" << hitc << "; MISS-" << miss;
    */
}


QList<int> hit_arr;
QList<int> miss_arr;

void MainWindow::multirun_aggregate()
{
/*
    int hit = 0, miss = 0, n = 1156;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type I.a: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_degree_comparison();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type I.b: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_weight_comparison();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type I.c: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_neighbour_initial_degree_bias();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type II.a: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_neighbour_CURRENT_degree_bias();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type II.b: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_highest_CURRENT_degree_neighbour();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type II.c: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_minimum_weight_neighbour();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type II.d: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_probabilistic_lowest_degree_neighbour_destructive();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type II.e: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;


    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_triangulation();
        QPair<int,int> pair = persistentGraphCompareWithKarate();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type III.a: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_highest_triangulated_vertex();
        QPair<int,int> pair = persistentGraphCompareWithKarate();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type III.c: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_triangulation_times_weight();
        QPair<int,int> pair = persistentGraphCompareWithKarate();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type III.d: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;

    hit = 0, miss = 0;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_triangulation_of_cluster();
        QPair<int,int> pair = persistentGraphCompareWithKarate();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type III.e: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;
*/
    /*
    double Q = 0.0;
    QList<QPair<double,int> > Qlist;
    int n = 100;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    qDebug() << "Sum Q: " << Q <<"; Average Q: " << Q/n;
    QString fileName = "TypeI.a";
    writeModularityReport(fileName, Qlist);
        // -----------------------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_degree_comparison();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeI.b";
    writeModularityReport(fileName, Qlist);
         // -----------------------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_weight_comparison();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeI.c";
    writeModularityReport(fileName, Qlist);
         // -----------------------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_neighbour_initial_degree_bias();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeII.a";
    writeModularityReport(fileName, Qlist);
         // -----------------------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_neighbour_CURRENT_degree_bias();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeII.b";
    writeModularityReport(fileName, Qlist);

   // -----------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_highest_CURRENT_degree_neighbour();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeII.c";
    writeModularityReport(fileName, Qlist);
    // -----------------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_minimum_weight_neighbour();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeII.d";
    writeModularityReport(fileName, Qlist);
    // ----------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_probabilistic_lowest_degree_neighbour_destructive();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeII.e";
    writeModularityReport(fileName, Qlist);
    */
        //------------------
    double Q = 0.0;
    QList<QPair<double,int> > Qlist;
    int n = 100;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    QString fileName = "TypeII.f";
    writeModularityReport(fileName, Qlist);
    //------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_greedy_max_degree();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeII.g";
    writeModularityReport(fileName, Qlist);
    //------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_greedy_max_weight();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeII.h";
    writeModularityReport(fileName, Qlist);
    // -----------------------
    /*
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_triangulation();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeIII.a";
    writeModularityReport(fileName, Qlist);
    // -----------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_probabilistic_triangulation();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeIII.b";
    writeModularityReport(fileName, Qlist);
            // -------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_highest_triangulated_vertex();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeIII.c";
    writeModularityReport(fileName, Qlist);
    // ---------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_triangulation_times_weight();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeIII.d";
    writeModularityReport(fileName, Qlist);
    // ------------------
    Q = 0.0;
    Qlist.clear();
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_triangulation_of_cluster();
        QPair<double,int> pair = calculate_modularity_for_clusters();
        Qlist.append(pair);
        Q+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        qDebug() << pair.first << ";" << pair.second << ";" << normal_q;
    }
    fileName = "TypeIII.e";
    writeModularityReport(fileName, Qlist);

/*
    // Matching Indicies
    QString filePath = "C:/Users/Dumex/Desktop/ExperimentStatFile/Pairwise_measures.txt";
    QFile file(filePath);
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    out << "RAND Index; Jaccard Index; Adjusted Rand Index; Girvan&Newman Index" << endl;

    int n = 100;
    out << "************* Type I.a *****************" << endl;
    out << "************* START ********************" << endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
        // -----------------------------------
    out<< "************* Type I.b *****************" << endl;
    out<< "************* START ********************" << endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_degree_comparison();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
         // -----------------------------------
    out<< "************* Type I.c *****************" << endl;
    out<< "************* START ********************" << endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_weight_comparison();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
         // -----------------------------------
    out<< "************* Type II.a *****************" << endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_neighbour_initial_degree_bias();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
         // -----------------------------------
    out<< "************* Type II.b *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_neighbour_CURRENT_degree_bias();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
   // -----------------------
    out<< "************* Type II.c *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_highest_CURRENT_degree_neighbour();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
    // -----------------------------
    out<< "************* Type II.d *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_minimum_weight_neighbour();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
    // ----------------------
    out<< "************* Type II.e *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_probabilistic_lowest_degree_neighbour_destructive();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
        //------------------
    out<< "************* Type II.f *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }

    // -----------------------
    out<< "************* Type III.a *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_triangulation();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
    // -----------------------
    out<< "************* Type III.b *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_probabilistic_triangulation();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
            // -------------
    out<< "************* Type III.c *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_with_highest_triangulated_vertex();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
    // ---------------
    out<< "************* Type III.d *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_triangulation_times_weight();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
    // ------------------

    out<< "************* Type III.e *****************"<< endl;
    out<< "************* START ********************"<< endl;
    for (int i = 0 ; i < n; i++)
    {
        random_aggregate_retain_vertex_using_triangulation_of_cluster();
        QList<double> id = compute_pairwise_matching_efficient();
        double newman = compute_Newman_fraction_of_classified(result);
        out << id[0] << ";" << id[1] << ";" << id[2] << ";" << newman << endl;
    }
*/
    qDebug()<< "END";
}

void MainWindow::multirun_aggregate_edge_weight_graph()
{

    int hit = 0, miss = 0;

    for (int i = 0 ; i < 1000; i++)
    {
        typeV_with_weight_update();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type V: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;
    writeDetailRun(5);

    hit = 0, miss = 0;
    for (int i = 0 ; i < 1000; i++)
    {
        typeVI_with_weight_update();
        QPair<int,int> pair = compareResultWithKarateClubDataStats();
        hit+= pair.first;
        miss+=pair.second;
    }
    qDebug() << "Type VI: *** TOTAL HIT:" << hit <<"; TOTAL MISS: " << miss;
    writeDetailRun(6);
}

void MainWindow::reset()
{
    if (GRAPHICS)
        resetGraphics();
    myVertexList.clear();
    myEdgeList.clear();
    hierarchy.clear();
    ground_truth_communities.clear();
    originVertexPos.clear();
    myLine.clear();
    centroids.clear();
    cluster_vertex.clear();
    vertex_label.clear();
    vertex_value.clear();
    vertex_weight.clear();
    edge_weight.clear();
    result.clear();
    graphIsReady = false;
}



QPair<int,int> MainWindow::compareResultWithKarateClubDataStats()
{
    graphIsReady = false;
    QList<int> actual_result;
    actual_result << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 <<33
                  << 0 << 0 << 0 << 0 << 33 << 33 << 0 << 0  << 33 << 0
                  << 33 << 0 << 33 << 33 << 33 << 33 << 33 << 33 << 33 << 33 << 33
                  << 33 << 33 << 33;
    if (centroids.size() == 0)
    {
        qDebug() << "HAVE NOT BEEN CLUSTERED";
        return qMakePair(0,0);
    }
    //interpret result
    QList<int> strict;
    for (int i = 0; i < 34; i++)
        strict.append(-1);
    for (int i = 0; i < centroids.size(); i++)
    {
        Vertex * v = centroids.at(i);
        if (v == 0)
            qDebug() << "NULL POINTER CENTROID";
        if (v->is_vertex_dragged_along())
            continue;
        int centroid_index = v->getIndex();
        strict.replace(centroid_index, centroid_index);
        QList<Vertex*> absorbed = v->getAbsorbedList();
        for (int j = 0; j < absorbed.size(); j++)
        {
            Vertex * v2 = absorbed.at(j);
            int cluster_index = v2->getIndex();
            strict.replace(cluster_index, centroid_index);
        }
    }
    //compare
    int hitc = 0, miss = 0;
    for (int i = 0; i < actual_result.size(); i++)
    {
        int actual = actual_result.at(i);
        int aggregated = strict.at(i);
        bool hit = (actual == aggregated);
        if (hit)
        {
            hitc++;
        }
        else
        {
            miss++;
        }
    }
    hit_arr.append(hitc);
    miss_arr.append(miss);
    return qMakePair(hitc, miss);
}

QPair<int, int> MainWindow::persistentGraphCompareWithKarate()
{
    QList<int> actual_result;
    graphIsReady = false;
    actual_result << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1
                  << 0 << 0 << 0 << 0 << 1 << 1 << 0 << 0  << 1 << 0
                  << 1 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
                  << 1 << 1 << 1;
    NormalGraph g;
    for(int i = 0; i < myVertexList.size(); i++)
        boost::add_vertex(g);

    for (int i = 0; i < hierarchy.size(); i++)
    {
        QPair<int,int> p = hierarchy.at(i);
        boost::add_edge(p.first, p.second, g);
    }

    std::vector<int> component(boost::num_vertices(g));
    int num = boost::connected_components(g, &component[0]);

    std::vector<int>::size_type i;

    int hit = 0, miss = 0;
    for (i = 0; i != component.size(); ++i)
    {
        if (component[i] == actual_result[i])
        {
            hit++;
        }
        else
        {
            miss++;
        }
    }
    hit_arr.append(hit);
    miss_arr.append(miss);
    return qMakePair(hit,miss);

}



/** Compute various indices for cluster matching
 * RAND, Jaccard, etc.
 * @brief MainWindow::compute_cluster_matching
 */
void MainWindow::compute_cluster_matching(QList<QList<int> > clusters)
{
    //check sum
    int n = 0;
    for (int i = 0; i < clusters.size(); i++)
        n += clusters[i].size();
    //checking ground truth
    if (ground_truth_communities.empty())
    {
        qDebug() << "GROUND TRUTH COMMUNITIES HAS NOT BEEN LOADED";
        return;
    }

    if (n == myVertexList.size())
        qDebug() << "Cluster Matching Pre-Check-Sum: OK!";
    else
    {
        qDebug() << "Cluster Matching Pre-Check-Sum: n does not match";
        return;
    }

    double RAND = compute_RAND_index(clusters),
            JAcc = compute_Jaccard_index(clusters),
            Newman = compute_Newman_fraction_of_classified(clusters);
    if (RAND == 1.0 && JAcc == 1.0 && Newman == 0)
    {
        compute_Newman_fraction_of_classified(clusters);
    }
    qDebug() << RAND << ";" << JAcc << ";" << Newman;
    return;
}

/** Compare result using RAND index
 * ground truth is X = {x1, x2 ..., xr }
 * result is Y = {y1, y2, .., ys}
 * with r != s
 * then, count:
    a, the number of pairs of elements in S that are in the same set in X and in the same set in Y
    b, the number of pairs of elements in S that are in different sets in X and in different sets in Y
    c, the number of pairs of elements in S that are in the same set in X and in different sets in Y
    d, the number of pairs of elements in S that are in different sets in X and in the same set in Y
 *  R = {a+b}/{a+b+c+d} = {a+b}/{n \choose 2 }
 * NAIVE IMPLEMENTATION: PREFERABLY FOR SMALL COMMUNITIES SET
 * @brief MainWindow::compare_using_RAND_index
 * @return
 */

double MainWindow::compute_RAND_index(QList<QList<int> > result)
{
    if (ground_truth_communities.empty())
    {
        qDebug() << "GROUND TRUTH COMMUNITIES HAS NOT BEEN LOADED";
        read_ground_truth_communities();
    }

    int n = 0;
    for (int i = 0; i < ground_truth_communities.size(); i++)
        n+= ground_truth_communities.at(i).size();

    int denom = double(n*(n-1))/2.0;
    QList<QPair<int,int> > truth_same, truth_diff;
    truth_same = create_same_set(ground_truth_communities);
    truth_diff = create_diff_set(ground_truth_communities);
    //Result Diff & Same
    QList<QPair<int,int> > result_same, result_diff;
    result_same = create_same_set(result);
    result_diff = create_diff_set(result);
    //compare
    int A = 0, B = 0;
    for (int i = 0; i < truth_same.size(); i++)
    {
        QPair<int,int> p = truth_same[i];
        QPair<int,int> reverse_p = qMakePair(p.second, p.first);
        if (result_same.contains(p) || result_same.contains(reverse_p))
            A++;
    }
    for (int i = 0; i < truth_diff.size(); i++)
    {
        QPair<int,int> p = truth_diff[i];
        QPair<int,int> reverse_p = qMakePair(p.second, p.first);
        if (result_diff.contains(p) || result_diff.contains(reverse_p))
            B++;
    }
    qDebug() << "DONE RAND;";
    double rand = double(A+B)/double(denom);
    return rand;
}


/** Compare result using JACCARD index
 * ground truth is X = {x1, x2 ..., xr }
 * result is Y = {y1, y2, .., ys}
 * with r != s
 * then, count:
    a, the number of pairs of elements in S that are in the same set in X and in the same set in Y
    c, the number of pairs of elements in S that are in the same set in X and in different sets in Y
    d, the number of pairs of elements in S that are in different sets in X and in the same set in Y
 *  R = {a}/{a+c+d}
 * NAIVE IMPLEMENTATION: PREFERABLY FOR SMALL COMMUNITIES SET
 * @brief MainWindow::compute_Jaccard_index
 * @return
 */

double MainWindow::compute_Jaccard_index(QList<QList<int> > result)
{
    if (ground_truth_communities.empty())
    {
        qDebug() << "GROUND TRUTH COMMUNITIES HAS NOT BEEN LOADED";
        read_ground_truth_communities();
    }

    QList<QPair<int,int> > truth_same, truth_diff;
    truth_same = create_same_set(ground_truth_communities);
    truth_diff = create_diff_set(ground_truth_communities);
    //Result Diff & Same
    QList<QPair<int,int> > result_same, result_diff;
    result_same = create_same_set(result);
    result_diff = create_diff_set(result);
    //compare
    int A = 0, C = 0, D = 0;
    for (int i = 0; i < truth_same.size(); i++)
    {
        QPair<int,int> p = truth_same[i];
        QPair<int,int> reverse_p = qMakePair(p.second, p.first);
        if (result_same.contains(p) || result_same.contains(reverse_p))
            A++;
    }
    for (int i = 0; i < truth_same.size(); i++)
    {
        QPair<int,int> p = truth_same[i];
        QPair<int,int> reverse_p = qMakePair(p.second, p.first);
        if (result_diff.contains(p) || result_diff.contains(reverse_p))
            C++;
    }
    for (int i = 0; i < truth_diff.size(); i++)
    {
        QPair<int,int> p = truth_diff[i];
        QPair<int,int> reverse_p = qMakePair(p.second, p.first);
        if (result_same.contains(p) || result_same.contains(reverse_p))
            D++;
    }
    qDebug() << "DONE JACC;";
    double jacc = double(A)/double(A+C+D);
    return jacc;
}


QList<QPair<int, int> > MainWindow::create_same_set(QList<QList<int> > C)
{
    QList<QPair<int, int> > same;
    for (int i = 0; i < C.size(); i++)
    {
        QList<int> ci = C.at(i);
        for (int j = 0; j < ci.size()-1; j++)
        {
            for(int k = j+1; k < ci.size(); k++)
            {
                if (ci[j] == ci[k])
                    qDebug() << "DUPLICATION DETECTED; Result is INCORRECT! " << ci[j] << ci[k];
                same.append(qMakePair(ci[j],ci[k]));
            }
        }
    }
    return same;
}


QList<QPair<int, int> > MainWindow::create_diff_set(QList<QList<int> > C)
{
    QList<QPair<int, int> > diff;
    for (int i = 0; i < C.size()-1; i++)
    {
        QList<int> ci = C.at(i);
        for (int j = i + 1; j < C.size(); j++)
        {
            QList<int> cj = C.at(j);
            for (int x = 0; x < ci.size(); x++)
            {
                for (int y = 0; y < cj.size(); y++)
                {
                    if (ci[x] == cj[y])
                        qDebug() << "DUPLICATION DETECTED; Result is INCORRECT! " << ci[x] << cj[y];
                    diff.append(qMakePair(ci[x],cj[y]));
                }
            }
        }
    }
    return diff;
}



/** Compute Girvan and Newman (2004) Fraction of Correctly Classified Vertices
 * a vertex is correctly classified if it is clustered into a community which has more than half of its natural (or truth)
 * the index is then simply: i = (number of correct)/(number of vertices)
 * @brief MainWindow::compute_Newman_fraction_of_classified
 * @param result_clusters
 * @return
 */

double MainWindow::compute_Newman_fraction_of_classified(QList<QList<int> > result_clusters)
{
    int n = 0;
    for (int i = 0; i < result.size(); i++)
        n+= result_clusters[i].size();
    //indexing ground truth community
    QMap<int,int> truth_c; //mapping vertex_index - truth community
    for (int i = 0; i < ground_truth_communities.size(); i++)
    {
        QList<int> c = ground_truth_communities[i];
        for (int j = 0; j < c.size(); j++)
            truth_c.insert(c[j],i);
    }
    //indexing result community
    QMap<int, int> result_c;
    for (int i = 0; i < result_clusters.size(); i++)
    {
        QList<int> c = result_clusters[i];
        for (int j = 0; j < c.size(); j++)
            result_c.insert(c[j],i);
    }

    QList<int> v;
    for (int i = 0; i < n; i++)
        v.append(-1);
    // -1 = not compared, 0 = incorrect, 1 = correct;
    QList<QList<int> > permutation = generateP(ground_truth_communities.size());
    QList<int> weight;
    for (int i = 0; i < permutation.size(); i++)
    {
        QList<int> p = permutation[i];
        int w = 0;
        for (int j = 0; j < p.size(); j++)
            w += ground_truth_communities[p[j]].size();
        weight.append(w);
    }

    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] == -1)
        {
            int truth_comm_index = truth_c.value(i);
            QList<int> truth = ground_truth_communities[truth_comm_index];
            int result_comm_index = result_c.value(i);
            QList<int> result = result_clusters[result_comm_index];
            //check for merging first
            int k = result.size();
            QList<QList<int> > p;
            for (int i = 0; i < weight.size(); i++)
            {
                if (weight[i] == k)
                    p.append(permutation[i]);
            }
            bool merged_sets = false;
            if (!p.empty())
                merged_sets = check_merge(result, p);

            if(merged_sets)
            {
                for (int i = 0; i < result.size(); i++)
                    v[result[i]] = 0;
            }
            else
            {
            //compare matching
                QList<int> matching = community_matching(truth, result);
                int correct = matching.size() - 1;
                if (correct < truth.size() / 2)
                {
                    for (int j = 0; j < matching.size(); j++)
                    {
                        int index = matching[j];
                        v[index] = 0;
                    }
                }
                else if (correct >= truth.size() / 2)
                {
                    for (int j = 0; j < matching.size(); j++)
                    {
                        int index = matching[j];
                        v[index] = 1;
                    }
                }
            }
        }
    }

    int match = 0;
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] < 0)
            qDebug() << "Newman Index - ERROR Potential Duplication/Un-parsed result";
        else
            match += v[i];
    }
    double newman_index = double(match)/double(n);
    return newman_index;
}

/** Check if a partition is a result of two or more original community merge
 * @brief MainWindow::check_merge
 * @param result_c
 * @return
 */
bool MainWindow::check_merge(QList<int> result_c, QList<QList<int> > permutation)
{
    qSort(result_c.begin(), result_c.end());
    for (int i = 0; i < permutation.size(); i++)
    {
        QList<int> combo = permutation[i];
        QList<int> merge;
        for (int j = 0; j < combo.size(); j++)
            merge.append(ground_truth_communities[combo[j]]);
        qSort(merge.begin(), merge.end());
        if (result_c == merge)
            return true;
    }
    return false;
}

/** Permutate, the HARD way
 * @brief MainWindow::generateP
 * @param size
 * @return
 */
QList<QList<int> > MainWindow::generateP(int n)
{
    QList<QList<int> > p;
    for(int r = 2; r <= n; r++)
    {
        std::vector<bool> v(n);
        std::fill(v.begin(), v.end() - n + r, true);

        do {
            QList<int> c;
            for (int i = 0; i < n; ++i) {
                if (v[i]) {
                    c.append(i);
                }
            }
            p.append(c);
        } while (std::prev_permutation(v.begin(), v.end()));
    }
    return p;
}


/** Sub Problem: given a set {x, y...} and a number k, find the subset that sum up equals k
 * @brief MainWindow::subset_sum
 * @return
 */
void MainWindow::subset_sum(QList<int> numbers, int target)
{
    QList<int> partial, partial_id;
    sum_up_recursive(numbers, target, partial_id);
}


void MainWindow::sum_up_recursive(QList<int> numbers, int target, QList<int> partial)
{
    int s = 0;
    foreach (int x, partial) s += x;

    if (s == target)
        qDebug() << partial;

    if (s >= target)
        return;

    for (int i = 0; i < numbers.size(); i++)
    {
        QList<int> remaining;
        int n = numbers[i];
        for (int j = i + 1; j < numbers.size(); j++) remaining.append(numbers[j]);

        QList<int> partial_rec = partial;
        partial_rec.append(n);
        sum_up_recursive(remaining, target, partial_rec);
    }
}


QList<int> MainWindow::community_matching(QList<int> truth_c, QList<int> result_c)
{
    QList<int> matching;
    for (int i = 0; i < result_c.size(); i++)
    {
        int id = result_c.at(i);
        if (truth_c.contains(id))
            matching.append(id);
    }
    return matching;
}


void MainWindow::writeDetailRun(int type)
{
    QString filePath = "C:/Users/Dumex/Desktop/ExperimentStatFile/Jan_2016/DetailRunofEachType/Type" + QString::number(type) + "_with_w_update.txt";
    QFile file(filePath);
    if (!file.exists())
        qDebug() << "HAHA";
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    out << "HIT/tMISS" << endl;
    for (int i = 0; i < hit_arr.size(); i++)
    {
        out << hit_arr[i] << '/t' << miss_arr[i] << endl;
    }
    file.close();
    hit_arr.clear();
    miss_arr.clear();
}


void MainWindow::writeModularityReport(QString fileName, QList<QPair<double, int> > Q)
{
    QString filePath = "C:/Users/Dumex/Desktop/ExperimentStatFile/May_2016/Modularity/"+ fileName + ".txt";
    QFile file(filePath);
    if (!file.exists()) {}
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    out << "Q;Number of Cluster;Normalized_Q" << endl;
    double sumQ = 0.0;
    for (int i = 0 ; i < Q.size(); i++)
    {
        QPair<double,int> pair = Q.at(i);
        sumQ+= pair.first;
        double normal_q = pair.first / (double)pair.second;
        out << pair.first << ";" << pair.second << ";" << normal_q <<endl;
    }
    out << "Sum Q: " << sumQ <<"; Average Q: " << sumQ/Q.size();
    file.close();
}


void MainWindow::writeCommunityMatchingReport(QString fileName, QList<QList<double> > Indices)
{
    QString filePath = "C:/Users/Dumex/Desktop/ExperimentStatFile/April_2016/Modularity/"+ fileName + ".txt";
    QFile file(filePath);
    if (file.exists()) {
        file.remove();
    }
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    out << "RAND Index; Jaccard Index; Girvan&Newman Index" << endl;
    double sumQ = 0.0;
    for (int i = 0 ; i < Indices.size(); i++)
    {
        QList<double> index;
        out << index[0] << ";" << index[1] << ";" << index[2] <<endl;
    }
    file.close();
}

QList<double> MainWindow::compute_pairwise_matching_efficient()
{
    QList<double> measures;
    quint64 nij = 0, nij_minus = 0, nij_square = 0, nij_choose_2 = 0;
    int row = ground_truth_communities.size(), column = result.size();
  //  qDebug() << "R:" << row << "; C:" << column;
    QList<quint64> ni, nj; //ni: sum row, nj: sum column
    for(int j = 0; j < column; j++)
        nj.append(0);

    for (int i = 0; i < row; i++)
    {
        quint64 sum_row = 0;
        QSet<int> X = ground_truth_communities[i].toSet();
        for (int j = 0; j < column; j++)
        {
            QSet<int> Y = result[j].toSet();
            QSet<int> copy = X;
            quint32 entry = copy.intersect(Y).size();
            nij_minus += (quint64)entry*(entry-1); // nij(nij-1)
            nij_square += (quint64)entry*entry; // nij^2
            nij_choose_2 += (quint64) entry*(entry-1)/2; //(nij choose 2) for adjust rand
            nij += (quint64)entry;
            sum_row += (quint64)entry;
            quint64 sum_col = nj[j];
            sum_col += entry;
            nj.replace(j, sum_col);
        }
        ni.append(sum_row);
    }
    int n = myVertexList.size();
    quint64 n_square = 0,
            n_choose_2 = n*(n-1)/2,
            ni_sum = 0, //sum row
            nj_sum =0, // sum column
            ni_choose_2 = 0, //bionomial row
            nj_choose_2 = 0, //bionomial column
            ni_square = 0,  // sum each row square
            nj_square = 0; // sum each column square

    for (int i = 0; i < ni.size(); i++)
    {
        int entry = ni[i];
        quint64 entry_square = entry*entry;
        quint64 entry_choose_2 = entry*(entry-1)/2;
        ni_square += entry_square;
        ni_choose_2 += entry_choose_2;
        ni_sum+=  entry;
    }
    ni.clear();
    for (int i = 0; i < nj.size(); i++)
    {
        int entry = nj[i];
        quint64 entry_square = entry*entry;
        quint64 entry_choose_2 = entry*(entry-1)/2;
        nj_square += entry_square;
        nj_choose_2 += entry_choose_2;
        nj_sum+=  entry;
    }
    nj.clear();

    n_square = qPow(n,2);

    QList<quint64> param_a,param_b,param_c,param_d;
    param_a << nij_minus;
    param_b << ni_square << nij_square;
    param_c << nj_square << nij_square;
    param_d << n_square << nij_square << ni_square << nj_square;

    quint64 a = calA(param_a);
    quint64 d = calD(param_d);
    quint64 c = calC(param_c); // type iii: same diff
    quint64 b = calB(param_b); //type iv: diff same

    double RAND = (double) (a+d)/(a+b+c+d);
    double Jaccard = (double) a/(a+b+c);
    //
    QList<quint64> param_ARI;
    param_ARI << ni_choose_2 << nj_choose_2 << n_choose_2 << nij_choose_2;
    double ARI = calAdRand(param_ARI);

    measures << RAND << Jaccard << ARI;
    return measures;
}

quint64 MainWindow::calA(QList<quint64> param)
{
    quint64 a = param[0];
    a = a/2;
    return a;
}

quint64 MainWindow::calB(QList<quint64> param)
{
    quint64 ai = param[0], nij_square = param[1];
    quint64 b = (ai-nij_square)/2;
    return b;
}

quint64 MainWindow::calC(QList<quint64> param)
{
    quint64 ni_square = param[0], nij_square = param[1];
    quint64 c = (ni_square - nij_square)/2;
    return c;
}

quint64 MainWindow::calD(QList<quint64> param)
{
    quint64 n_square = param[0],
            nij_square = param[1],
            ni_square = param[2],
            nj_square = param[3];
    quint64 d = n_square + nij_square - ni_square - nj_square; //type ii: diff diff
    d /= 2;
    return d;
}

double MainWindow::calAdRand(QList<quint64> param)
{
    quint64 ni_choose_2 = param[0],
            nj_choose_2 = param[1],
            n_choose_2 = param[2],
            nij_choose_2 = param[3];
    double nc = (double)ni_choose_2*nj_choose_2/n_choose_2;
    double nom = (double)(nij_choose_2 - nc);
    double sum = (double) (ni_choose_2 + nj_choose_2)/2;
    double denom = sum - nc;
    double ARI = nom/denom;
    return ARI;
}

// --------------------------- AGGREGATE WITH WEIGHT UPDATE --------------------------

void MainWindow::typeII_with_weight_update()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            std::uniform_int_distribution<int> distribution2(0,no_neighbour-1);
            int selected_edge_index = distribution2(generator);
            Edge * e = selected->getEdge(selected_edge_index);
            Vertex * neighbour = selected->get_neighbour_fromEdge(selected_edge_index); //get the neighbour (not clean)
            Vertex * winner, * loser;
            int selected_d = selected->getNumberEdge(), neighbour_d = neighbour->getNumberEdge();
            if (selected_d >= neighbour_d)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }

            //create the animation
            QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
            QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
            group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            winner->setWeight(winner->getWeight() + loser->getWeight());
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
    // group_anim->start();
    // connect(group_anim, SIGNAL(finished()), this, SLOT(_()));
    draw_aggregation_result();
}

void MainWindow::typeII_probabilistic_current_degree_uniform_neighbour()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        QList<Vertex*> ran_list;
        for (int i = 0 ;i < players.size(); i++)
        {
            Vertex * v = players.at(i);
            int w = v->getNumberEdge();
            for (int j = 0; j < w; j++)
                ran_list.append(v);
        }
        if (ran_list.size() == 0)
        {
            int size = players.size();
            std::uniform_int_distribution<int> distribution(0,size-1);
            int selected_index = distribution(generator);
            Vertex * selected = players.at(selected_index);
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else
        {
            int size = ran_list.size();
            std::uniform_int_distribution<int> distribution(0,size-1);
            int selected_index = distribution(generator);
            Vertex * selected = ran_list.at(selected_index);
            //get a neighbour
            int no_neighbour = selected->getNumberEdge();
            if (no_neighbour == 0) // if there is no neighbour, declare a winner
            {
                winners.append(selected);
                players.removeOne(selected);
                t++;
            }
            else // else absorb
            {
                std::uniform_int_distribution<int> distribution2(0,no_neighbour-1);
                int selected_edge_index = distribution2(generator);
                Edge * e = selected->getEdge(selected_edge_index);
                Vertex * neighbour = selected->get_neighbour_fromEdge(selected_edge_index); //get the neighbour (not clean)
                Vertex * winner, * loser;
                winner = selected;
                loser = neighbour;

                //create the animation
                QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
                QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
                group_anim->addAnimation(single_anim);
                winner->absorb_removeEdge(e);
                players.removeOne(loser);
                t++;
            }
        }
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
    // group_anim->start();
    // connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}

void MainWindow::typeII_probabilistic_degree_smallest_neighbour()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        QList<Vertex*> ran_list;
        for (int i = 0 ;i < players.size(); i++)
        {
            Vertex * v = players.at(i);
            int w = v->getNumberEdge();
            for (int j = 0; j < w; j++)
                ran_list.append(v);
        }
        if (ran_list.size() == 0)
        {
            int size = players.size();
            std::uniform_int_distribution<int> distribution(0,size-1);
            int selected_index = distribution(generator);
            Vertex * selected = players.at(selected_index);
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else
        {
            int size = ran_list.size();
            std::uniform_int_distribution<int> distribution(0,size-1);
            int selected_index = distribution(generator);
            Vertex * selected = ran_list.at(selected_index);
            //get a neighbour
            Edge * e = selected->getSmallestCurrentDegreeNeighbour();
            Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
            Vertex * winner, * loser;
            winner = selected;
            loser = neighbour;

            //create the animation
            QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
            QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
            group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
    // group_anim->start();
    // connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}


/*
 * RANDOM AGGREGATION WITH WEIGHT COMPARISON
 */
void MainWindow::typeIII_with_weight_update()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            std::uniform_int_distribution<int> distribution2(0,no_neighbour-1);
            int selected_edge_index = distribution2(generator);
            Edge * e = selected->getEdge(selected_edge_index);
            Vertex * neighbour = selected->get_neighbour_fromEdge(selected_edge_index); //get the neighbour (not clean)
            Vertex * winner, * loser;
            int selected_w = selected->getWeight(), neighbour_w = neighbour->getWeight();
            if (selected_w >= neighbour_w)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }

            //create the animation
            QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
            QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
            group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            winner->setWeight(winner->getWeight() + loser->getWeight());
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}

void MainWindow::typeIII_probabilistic_with_weight_update()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        QList<Vertex*> ran_list;
        for (int i = 0; i < players.size(); i++)
        {
            Vertex * player = players[i];
            int w = player->getWeight();
            for (int j = 0; j < w; j++)
                ran_list.append(player);
        }
        std::uniform_int_distribution<int> distribution(0,ran_list.size() - 1);
        int selected_index = distribution(generator);
        Vertex * selected = ran_list.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            std::uniform_int_distribution<int> distribution2(0,no_neighbour-1);
            int selected_edge_index = distribution2(generator);
            Edge * e = selected->getEdge(selected_edge_index);
            Vertex * neighbour = selected->get_neighbour_fromEdge(selected_edge_index); //get the neighbour (not clean)
            Vertex * winner, * loser;
            int selected_w = selected->getWeight(), neighbour_w = neighbour->getWeight();
            if (selected_w >= neighbour_w)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }

            //create the animation
            QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
            QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
            group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            winner->setWeight(winner->getWeight() + loser->getWeight());
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}

/** A Vertex (u) is selected u.a.r
 * Then a neighbour (v) is selected with an INITIAL degree bias (meaning degree does not change), means:
 * Pr(v) = /sum_d(v) = 1, /forall v adj u
 * @brief MainWindow::random_aggregate_with_neighbour_degree_bias
 */
void MainWindow::typeIV_with_weight_update()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Vertex * neighbour = selected->aggregate_get_degree_biased_neighbour();
            Edge * e = selected->getEdgeFromVertex(neighbour);
            //create the animation
            QList<int> anim_mater = neighbour->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
            QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, neighbour, selected);
            group_anim->addAnimation(single_anim);
            selected->absorb_removeEdge(e);
            selected->setWeight(selected->getWeight() + neighbour->getWeight());
            players.removeOne(neighbour);
            t++;
        }
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}

/** SELECT A VERTEX <U> UNIFORMLY AT RANDOM
 * SELECT AN ASSOCIATE EDGE WITH HIGHEST WEIGHT, HENCE VERTEX <V>
 * VERTEX WITH HIGHER WEIGHT ABSORBS THE OTHER
 * @brief MainWindow::random_aggregate_with_highest_edge_weight_and_weight_comparison
 */
void MainWindow::typeV_with_weight_update()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->getHighestWeightEdge();
            Vertex * neighbour, * winner, * loser;
            if (e->toVertex() == selected)
                neighbour = e->fromVertex();
            else
                neighbour = e->toVertex();

            if (selected == neighbour)
            {   qDebug() << "BUG CHECK" << "SELECTED POINTER == NEIGHBOUR POINTER";
                return;
            }
            int u_w = selected->getWeight(), v_w = neighbour->getWeight();
            if (u_w >= v_w)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }
            //create the animation
            QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
            QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
            group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            winner->setWeight(winner->getWeight() + loser->getWeight());
            players.removeOne(loser);
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            t++;
        }
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();
}

/** SELECT A VERTEX <U> UNIFORMLY AT RANDOM
 * SELECT AN EDGE WITH PROBABILITY GIVEN TO THE WEIGHT:
 * Prob(edge_i) = w(i) / sum(w)
 * VERTEX WITH HIGHER WEIGHT ABSORBS
 * @brief MainWindow::random_aggregate_with_edge_weight_bias
 */
void MainWindow::typeVI_with_weight_update()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    QSequentialAnimationGroup * group_anim = new QSequentialAnimationGroup;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Edge * e = selected->getWeightedProbabilisticEdge();
            Vertex * neighbour, * winner, * loser;
            if (e->toVertex() == selected)
                neighbour = e->fromVertex();
            else
                neighbour = e->toVertex();

            if (selected == neighbour)
            {   qDebug() << "BUG CHECK" << "SELECTED POINTER == NEIGHBOUR POINTER";
                return;
            }
            int u_w = selected->getWeight(), v_w = neighbour->getWeight();
            if (u_w >= v_w)
            {
                winner = selected;
                loser = neighbour;
            }
            else
            {
                winner = neighbour;
                loser = selected;
            }
            //create the animation
            QList<int> anim_mater = loser->getEdgesIndexForRemovalAnimation(e); //get all edges going to be removed
            QSequentialAnimationGroup * single_anim = prepare_animation_for_one_vertex(anim_mater, loser, winner);
            group_anim->addAnimation(single_anim);
            winner->absorb_removeEdge(e);
            winner->setWeight(winner->getWeight() + loser->getWeight());
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            players.removeOne(loser);
            t++;
        }
    }
    centroids = winners;
   // draw_dense_graph_aggregation_result();
   // group_anim->start();
   // connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_result();

}


// ---------------------------------------- AGGREGATION WHITE RETAINING VERTEX -----------------------------
/** Vertices are only removed from the playing list instead of completely removed
 * @brief MainWindow::typeI_retaining_vertex
 */
void MainWindow::typeI_retaining_vertex()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        std::uniform_int_distribution<int> distribution(0,players.size() - 1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        std::uniform_int_distribution<int> distribution2(0,no_neighbour-1);
        int selected_edge_index = distribution2(generator);
        Edge * e = selected->getEdge(selected_edge_index);
        Vertex * neighbour = selected->get_neighbour_fromEdge(selected_edge_index); //get the neighbour (not clean)
        Vertex * winner, * loser;
        winner = neighbour;
        loser = selected;

        //create the animation
        winner->absorb_retainEdge_setParentPointer(e);
        hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
        players.removeOne(loser);
        t++;
    }
    centroids = winners;
    // draw_dense_graph_aggregation_result();
   //  group_anim->start();
   //  connect(group_anim, SIGNAL(finished()), this, SLOT(draw_aggregation_result()));
    draw_aggregation_retain_vertex_result();
}

void MainWindow::typeII_retaining_vertex()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }

    //initialise arrays
    QList<Vertex*> players = myVertexList;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        std::uniform_int_distribution<int> distribution(0,players.size()-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        Edge * e = selected->getDegreeProbabilisticEdge();
        Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
        Vertex * winner, * loser;
        if (!selected->is_vertex_absorbed())
        {
            winner = neighbour;
            loser = selected;
        }
        else if (!neighbour->is_vertex_absorbed())
        {
            winner = selected;
            loser = neighbour;
        }
        else if (neighbour->is_vertex_absorbed() && selected->is_vertex_absorbed())
        {
            qDebug() << "ERR: BOTH VERTEX ABSORBED ?";
        }

        //create the animation
        hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
        winner->absorb_retainEdge(e);
        players.removeOne(loser);
        t++;
    }
    draw_aggregation_retain_vertex_result();
}

void MainWindow::typeIIIc_retaining_vertex()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }

    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        std::uniform_int_distribution<int> distribution(0,players.size()-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);

        //get a neighbour
        Edge * e = selected->getWeightedProbabilisticEdge();
        Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
        Vertex * winner, * loser;
        if (!selected->is_vertex_absorbed())
        {
            winner = neighbour;
            loser = selected;
        }
        else if (!neighbour->is_vertex_absorbed())
        {
            winner = selected;
            loser = neighbour;
        }
        else if (neighbour->is_vertex_absorbed() && selected->is_vertex_absorbed())
        {
            qDebug() << "ERR: BOTH VERTEX ABSORBED ?";
        }
        //create the animation
        hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
        winner->absorb_retainEdge(e);
        int new_w = winner->getWeight() + loser->getWeight();
        winner->setWeight(new_w);
        loser->setWeight(new_w);
        players.removeOne(loser);
        t++;
    }
    centroids = winners;
    draw_aggregation_retain_vertex_result();
}

/** Probability of disconnecting an edge
 * @brief MainWindow::typeIIId_removing_vertex
 */
void MainWindow::typeIIId_removing_vertex()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    int sum_w = 0;
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
        sum_w += v->getWeight();
    }


    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        std::uniform_int_distribution<int> distribution(0,players.size()-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);

        //get a neighbour
        Edge * e = selected->getWeightedProbabilisticEdge();
        Vertex * neighbour = selected->get_neighbour_fromEdge(e); //get the neighbour (not clean)
        int total_w = selected->getWeight() + neighbour->getWeight();
        double p_dis = double(total_w) / double(sum_w);
        std::uniform_real_distribution<double> ran(0,1);
        double ran_n = ran(generator);

        if (ran_n < p_dis)
        { // disconnect the edge
            qDebug() << "DIS";
            selected->removeEdge(e);

        }
        else
        {
            Vertex * winner, * loser;
            if (!selected->is_vertex_absorbed())
            {
                winner = neighbour;
                loser = selected;
            }
            else if (!neighbour->is_vertex_absorbed())
            {
                winner = selected;
                loser = neighbour;
            }
            else if (neighbour->is_vertex_absorbed() && selected->is_vertex_absorbed())
            {
                qDebug() << "ERR: BOTH VERTEX ABSORBED ?";
            }
            //create the animation
            hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
            winner->absorb_retainEdge(e);
            int new_w = winner->getWeight() + loser->getWeight();
            winner->setWeight(new_w);
            loser->setWeight(new_w);
            players.removeOne(loser);
        }
        t++;
    }
    centroids = winners;
    draw_aggregation_retain_vertex_result();
}

void MainWindow::typeIV_retaining_vertex()
{
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        //get a neighbour
        int no_neighbour = selected->getNumberEdge();
        if (no_neighbour == 0) // if there is no neighbour, declare a winner
        {
            winners.append(selected);
            players.removeOne(selected);
            t++;
        }
        else // else absorb
        {
            Vertex * neighbour = selected->aggregate_get_degree_biased_neighbour();
            Edge * e = selected->getEdgeFromVertex(neighbour);
            //create the animation
            selected->absorb_removeEdge(e);
            selected->setWeight(selected->getWeight() + neighbour->getWeight());
            players.removeOne(neighbour);
            t++;
        }
    }
    draw_aggregation_retain_vertex_result();
}


/** The absorbed vertices are now retained in the graph.
 *
 * @brief MainWindow::random_aggregate_retain_vertex_using_triangulation_and_weight_comparison
 */
void MainWindow::typeVII_with_weight_update()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        v->setWeight(v->getNumberEdge());
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
       // Edge * e = selected->getProbabilisticTriangulationCoeffVertex();
        Edge * e = selected->correct_getMostMutualVertex();
        Vertex * neighbour, * winner, * loser;
        if (e->toVertex() == selected)
            neighbour = e->fromVertex();
        else
            neighbour = e->toVertex();

        winner = neighbour;
        loser = selected;
        //create the animation
        winner->absorb_retainEdge(e);
        winner->setWeight(winner->getWeight() + loser->getWeight());
        hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
        players.removeOne(loser);
        t++;
    }
 //   draw_aggregation_retain_vertex_result();

}

/** The absorbed vertices are now retained in the graph.
 *
 * @brief MainWindow::random_aggregate_retain_vertex_using_triangulation_and_weight_comparison
 */
void MainWindow::typeVIII_with_weight_update()
{
    hierarchy.clear();
    if (!checkGraphCondition())
    {
        reConnectGraph();
    }
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
       // v->setWeight(v->getNumberEdge());
        v->setWeight(1);
    }
    //initialise arrays
    QList<Vertex*> players = myVertexList;
    QList<Vertex*> winners;
    int t = 0;
    while(!players.empty()) //start
    {
        //select a vertex uniformly at random
        int size = players.size();
        std::uniform_int_distribution<int> distribution(0,size-1);

        int selected_index = distribution(generator);
        Vertex * selected = players.at(selected_index);
        Edge * e = selected->getProbabilisticTriangulationAndWeightVertex();
        if (e == 0)
        {
            selected->setParent(selected);
            players.removeOne(selected);
            continue;
        }
        Vertex * neighbour, * winner, * loser;
        if (e->toVertex() == selected)
            neighbour = e->fromVertex();
        else
            neighbour = e->toVertex();

        winner = neighbour;
        loser = selected;
        //create the animation
        winner->absorb_retainEdge(e);
        winner->setWeight(winner->getWeight() + loser->getWeight());
        hierarchy.append(qMakePair(loser->getIndex(), winner->getIndex()));
        players.removeOne(loser);
        t++;
    }

    // draw_aggregation_retain_vertex_result();
}

/**
 * @brief MainWindow::calculate_modularity_for_clusters
 * @return modulairty for each cluster
 */
QPair<double,int> MainWindow::calculate_modularity_for_clusters()
{
    // ADD A STEP OF SYCHNROSIE A NORMAL GRAPH TYPE
    graphIsReady = false;
    NormalGraph g;
    for(int i = 0; i < myVertexList.size(); i++)
        boost::add_vertex(g);

    for (int i = 0; i < hierarchy.size(); i++)
    {
        QPair<int,int> p = hierarchy.at(i);
        boost::add_edge(p.first, p.second, g);
    }

    std::vector<int> component(boost::num_vertices(g));
    int num = boost::connected_components(g, &component[0]);
    g = myGraph;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    //count intra and inter edge
    int intra [num];
    for (int i = 0; i < num; i++)
        intra[i] = 0;
    int inter [num];
    for (int i = 0; i < num; i++)
        inter[i] = 0;

    int count = 0;
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        int ui = index[u];
        int vi = index[v];
        if (ui > component.size() || vi > component.size() )
            qDebug() << "ERR: OUT OF BOUNDS";
        int c_u = component[ui];
        int c_v = component[vi];
        if (c_u == c_v)
        {
            intra[c_u]++;
        }
        else
        {
            inter[c_u]++;
            inter[c_v]++;
        }
        count++;
    }

    double Q = 0.0;
    int num_edge = boost::num_edges(g);
    int m = num_edge;

    for (int i = 0; i < num; i++)
    {
        double eii = intra[i];
        double ai = inter[i];
        double e = eii/m;
        double a = (2*eii + ai)/(2*m);
        double Qi = e - a*a;

        Q += Qi;
    }
  //  qDebug() << "Modularity: Q" << Q <<"; Average: " << Q/num;
    return qMakePair(Q, num);
}

/** Calculate Modularity for LARGE graph (Ground Truth)
 * @brief MainWindow::calculate_modularity_for_ground_truth
 * @return
 */
double MainWindow::calculate_modularity_for_ground_truth()
{
    quint32 global_e = myEdgeList.size();
    double Q = 0.0;
    for (int i = 0; i < ground_truth_communities.size(); i++)
    {
        QList<int> c = ground_truth_communities[i];
        QSet<int> vi = c.toSet();
        quint32 intra = 0, inter = 0;
        for (int j = 0; j < c.size(); j++)
        {
            int id = c[j];
            Vertex * v = myVertexList.at(id);
            QList<Edge*> adj = v->getAllEdge();
            for (int k = 0; k < adj.size(); k++)
            {
                Vertex * other = v->get_neighbour_fromEdge(adj[k]);
                int other_id = other->getIndex();
                if (vi.contains(other_id))
                    intra++;
                else
                    inter++;
            }
        }

        quint64 m = 2*global_e;
        double e = (double)intra/m;
        double a = (double)(intra + inter)/m;
        double Qi = e - qPow(a,2);
        Q += Qi;
    }
    return Q;
}

/** Investiaget Bridges
 * Get:
 * a: number of bridges
 * b: degree of two end points of bridges
 * c: how many bridges a high degree vertex has
 * d: how many bridges a low degree vertex has
 * e: degree of vertices that have no bridge
 * @brief MainWindow::get_bridge_stats
 */
void MainWindow::get_bridge_stats()
{
    if (!graphIsReady && ground_truth_communities.empty())
    {
        qDebug() << "GRAPH IS NOT READY!";
        return;
    }
    examining_bridges();
}

/** Brandes Betweenness Centrality Clustering
 * @brief MainWindow::betweenness_centrality_clutering
 */
bool pairCompareFirst(const QPair<double,int> &s1, const QPair<double,int> &s2)
{
    return s1.first < s2.first;
}

void MainWindow::betweenness_centrality_clutering()
{
    // Declear Weighted Graph as in Boost
    if (!graphIsReady)
    {
        reConnectGraph();
    }
    typedef boost::adjacency_list< boost::setS, boost::vecS, boost::undirectedS, boost::no_property, boost::property_kind<int> > BGraph;
    BGraph bg;
    typedef BGraph::vertex_descriptor gVertex;
    typedef BGraph::edge_descriptor gEdge;
    //copy graph
    for (int i = 0; i < myVertexList.size(); i++)
        boost::add_vertex(bg);

    typedef std::map<gEdge, int> StdEdgeIndexMap;
    StdEdgeIndexMap my_e_index;
    // associative property map needed for iterator property map-wrapper
    typedef boost::associative_property_map< StdEdgeIndexMap > EdgeIndexMap;
    EdgeIndexMap e_index(my_e_index);

    for(int i = 0; i < myEdgeList.size(); i++)
    {
        Edge * e = myEdgeList.at(i);
        int u = e->fromVertex()->getIndex(), v = e->toVertex()->getIndex();
        gEdge edge = boost::add_edge(u,v,bg).first;
        my_e_index.insert(std::pair< gEdge, int >( edge, i));
    }
    /*
    int i = 0;
    BGL_FORALL_EDGES(edge, bg, BGraph)
    {
      my_e_index.insert(std::pair< gEdge, int >( edge, i));
      ++i;
    }
    */
    std::vector< double > e_centrality_vec(boost::num_edges(bg), 0.0);
    // Create the external property map
    boost::iterator_property_map< std::vector< double >::iterator, EdgeIndexMap >
            e_centrality_map(e_centrality_vec.begin(), e_index);

    //threshold
    double max_centrality = -1.0;
    boost::bc_clustering_threshold< double > terminate(max_centrality, bg, false);
    //invoke
    boost::betweenness_centrality_clustering( bg, terminate, e_centrality_map );
    //parse result, sort first
    QList<QPair<double,int> > edge_c;
    for (int i = 0; i < e_centrality_vec.size(); i++)
    {
        double e_c = e_centrality_vec[i];
        edge_c.append(qMakePair(e_c, i));
    }
    qSort(edge_c.begin(), edge_c.end(), pairCompareFirst);
    //remove from original
    NormalGraph g;
    g = myGraph;
    std::vector<int> component(boost::num_vertices(g));
    int num = boost::connected_components(g, &component[0]);
    qDebug() << num;
    int last = edge_c.size()-1, truth_c = ground_truth_communities.size();
    if (truth_c == 0)
    {
        qDebug() << "- Ground Truth Has Not Been Initialised!";
        return;
    }
    while (num != truth_c && last >= 0)
    {
        int e_i = edge_c.at(last).second;
        Edge * e = myEdgeList.at(e_i);
        int ui = e->fromVertex()->getIndex(), vi = e->toVertex()->getIndex();
     //   qDebug() << "REMOVED E:" << e->getIndex() << "; Centrality: " << edge_c.at(last).first
     //   << e_centrality_vec.at(last)<<"; From: " << ui << "; To: " << vi;
        boost::remove_edge(ui,vi,g);
        last--;
        num = boost::connected_components(g, &component[0]);
    }

    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        int ui = index[u];
        int vi = index[v];
        hierarchy.append(qMakePair(ui,vi));
    }

    qDebug() << "Done";
    draw_aggregation_retain_vertex_result();
}


void MainWindow::createSnapShot(int t)
{
    QRect viewport = myView->viewport()->rect();
    QImage image(viewport.size(), QImage::Format_ARGB32);
    image.scaled(viewport.width()/2, viewport.height()/2, Qt::IgnoreAspectRatio);
    image.fill(Qt::transparent);                                              // Start all pixels transparent

    QPainter painter(&image);
    myView->render(&painter);
    QString defaulpath = "C:/Users/Dumex/Desktop/ExperimentSnapShot/";
    QDir dir(defaulpath);
    QFileInfoList tree = dir.entryInfoList();
    if (tree.size() > 0)
    {

    }
    image.save(defaulpath + "file_name" + QString::number(t) + ".png");
}

void MainWindow::readjust_graph(Edge *mainE, QList<int> subE, Vertex *delete_vertex, int step)
{
    Vertex * v1 = mainE->fromVertex(), * v2 = mainE->toVertex();
    v1->setBackgroundColour(Qt::red);
    v2->setBackgroundColour(Qt::red);
    //change painter
    for (int i = 0; i < subE.size(); i++)
    {
        LineAnimator * line = myLine.at(subE[i]);
        line->drawDeletedEdge();
    }
    LineAnimator * mainl = myLine.at(mainE->getIndex());
    mainl->highlight_edge();
    createSnapShot(step); //take screenshot
    //revert the scene
    v1->setBackgroundColour(Qt::green);
    v2->setBackgroundColour(Qt::green);
    delete_vertex->setBackgroundColour(Qt::black);
    for (int i = 0; i < subE.size(); i++)
    {
        LineAnimator * line = myLine.at(subE[i]);
        myScene->removeItem(line);
    }
    myScene->removeItem(mainl);
}

void MainWindow::parse_aggregation_result()
{
    // prepare for cluster matching;
    QList<QList<int> > clusters;
    for (int i = 0; i < centroids.size(); i++)
    {
        QList<int> c;
        c.append(centroids[i]->getIndex());
        QList<Vertex*> cluster = centroids[i]->getMyCluster();
        for (int j = 0; j < cluster.size(); j++)
        {
            c.append(cluster[j]->getIndex());
        }
        clusters.append(c);
    }
    result = clusters;
}

/** Parse Result from Retain Agg
 * @brief MainWindow::parse_retain_result
 */
void MainWindow::parse_retain_result()
{
    NormalGraph g;
    for(int i = 0; i < myVertexList.size(); i++)
        boost::add_vertex(g);
    for (int i = 0; i < hierarchy.size(); i++)
    {
        QPair<int,int> p = hierarchy.at(i);
        boost::add_edge(p.first, p.second, g);
    }

    std::vector<int> component(boost::num_vertices(g));
    int num = boost::connected_components(g, &component[0]);


    QList<QList<int> > clusters;
    for (int i = 0 ; i < num; i++)
    {
        QList<int> c;
        clusters.append(c);
    }

    for (int i = 0; i < component.size(); i++)
    {
        QList<int> c = clusters[component[i]];
        c.append(i);
        clusters.replace(component[i],c);
    }
    result = clusters;
}



/** SAMPLE
 * @brief MainWindow::exportGML
 */
void MainWindow::exportGML()
{
    QString filepath = "C:/Users/Dumex/Desktop/Test.gml";
    QFile file(filepath);
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    out << "Creator 'Ngoc Vu'" << endl;
    out << "graph" << endl;
    out << "[" << endl;
    out << "  directed 0" << endl;
    //begin writing node
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        QString colourRGB = v->getDeselectedColour().name();
        out << "  node \n"
            << "  [" << endl
            << "    id " << QString::number(v->getIndex()) << endl
            << "    graphics \n"
            << "    [ \n "
            << "      fill \"" << colourRGB << "\"" << endl
            << "    ]" << endl
            << "  ]" << endl;
    }
    //write the edge
    NormalGraph g;
    g = myGraph;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        int ui = index[u];
        int vi = index[v];
        out << "  edge \n"
            << "  [" << endl
            << "    source " + QString::number(ui) + "\n"
            << "    target " + QString::number(vi) + "\n"
            << "  ]" << endl;
    }

    out << "]" << endl;
    file.close();
}

/** Calculate Maching Indices Compare With Ground Truth
 * @brief MainWindow::calculate_indices
 */
void MainWindow::calculate_indices()
{
    if (ground_truth_communities.size() == 0 || result.empty())
    {
        qDebug() << "Result or Ground Truth Empty";
        return;
    }
    QList<double> r =compute_pairwise_matching_efficient();
    qDebug() << "Rand: " << r[0]
             << "\nJaccard: " << r[1]
             << "\nARI: " << r[2];

}

/** Calculate Q
 * @brief MainWindow::calculate_Q
 */
void MainWindow::calculate_Q()
{
    if (hierarchy.empty())
    {
        qDebug() << "Results are empty";
        return;
    }
    qDebug() << "Q: " <<calculate_modularity_for_clusters().first;
}

void MainWindow::exportGML(int take, double p)
{
        QString filepath = "C:/Users/Dumex/Desktop/Gnp/Take" + QString::number(take)+ +"_p" + QString::number(p) + ".gml";
        QFile file(filepath);
        file.open(QFile::WriteOnly | QFile::Text);
        QTextStream out(&file);
        out << "Creator 'Ngoc Vu'" << endl;
        out << "graph" << endl;
        out << "[" << endl;
        out << "  directed 0" << endl;
        //begin writing node
        for (int i = 0; i < myVertexList.size(); i++)
        {
            Vertex * v = myVertexList.at(i);
            QString colourRGB = v->getDeselectedColour().name();
            out << "  node \n"
                << "  [" << endl
                << "    id " << QString::number(v->getIndex()) << endl
                << "    graphics \n"
                << "    [ \n "
                << "      fill \"" << colourRGB << "\"" << endl
                << "    ]" << endl
                << "  ]" << endl;
        }
        //write the edge
        NormalGraph g;
        g = myGraph;
        typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
        typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
        IndexMap index = get(boost::vertex_index, g);
        for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
        {
            boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
            e_desc = *ep.first;
            boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
            u = source(e_desc, g);
            v = target(e_desc, g);
            int ui = index[u];
            int vi = index[v];
            out << "  edge \n"
                << "  [" << endl
                << "    source " + QString::number(ui) + "\n"
                << "    target " + QString::number(vi) + "\n"
                << "  ]" << endl;
        }

        out << "]" << endl;
        file.close();
}

void MainWindow::exportGML(int take, double p, double Q)
{
    QString filepath = "C:/Users/Dumex/Desktop/Gnp/Take" + QString::number(take)+ +"_p" + QString::number(p) +"_Q" +QString::number(Q) + ".gml";
    QFile file(filepath);
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    out << "Creator 'Ngoc Vu'" << endl;
    out << "graph" << endl;
    out << "[" << endl;
    out << "  directed 0" << endl;
    //begin writing node
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        QString colourRGB = v->getDeselectedColour().name();
        out << "  node \n"
            << "  [" << endl
            << "    id " << QString::number(v->getIndex()) << endl
            << "    graphics \n"
            << "    [ \n "
            << "      fill \"" << colourRGB << "\"" << endl
            << "    ]" << endl
            << "  ]" << endl;
    }
    //write the edge
    NormalGraph g;
    g = myGraph;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    for (std::pair<edge_iter, edge_iter> ep = edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        int ui = index[u];
        int vi = index[v];
        out << "  edge \n"
            << "  [" << endl
            << "    source " + QString::number(ui) + "\n"
            << "    target " + QString::number(vi) + "\n"
            << "  ]" << endl;
    }

    out << "]" << endl;
    file.close();
}



// ------------------------- END OF ACTIONS AND MENUS -------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------


// ------------------------- FOR LARGE GRAPH ----------------------------------------
// ------------ TOO LAZY TO SEPERATE TO A DIFFRENT PROJECT --------------------------
/** READ LARGE GRAPH WITH GROUND TRUTH COMMUNITIES
 * PARSE AND REINDEXING SO IT FOLLOWS A CONSISTENT INDICES
 * FILES ARE FROM SNAP - STANFORD
 * @brief MainWindow::read_large_graph_with_ground_truth_communities
 */
void MainWindow::read_large_graph_with_ground_truth_communities()
{
    //testing
    //READING ...
    //parsing, file are edges file with \t seperator
    //NOTE: Vertices are not uniquely label in increasing order
    QString filePath = QFileDialog::getOpenFileName(this,
                                                     ("Load Vertex")
                                                     , "C:/Users/Dumex/Desktop"
                                                     , ("Text Files (*.txt *csv)"));
    QFile file(filePath);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);
    QList<QPair<int,int> > edge;
    while (!in.atEnd())
    {
        QStringList str = in.readLine().split('\t');
        if (str.size() == 2 && !str[0].startsWith("#"))
        {
            bool ok;
            int v1 = str[0].toInt(&ok), v2 = str[1].toInt(&ok);
            if (!ok)
            {
                qDebug() << "ERROR WHEN PARSING: Data not int! Terminating";
                return;
            }
            edge.append(qMakePair(v1,v2));
        }
    }
    file.close();

    qDebug() << "HERE";
    qDebug() << "FILE CLOSED \nPARSING NOW";
    QMap<int, Vertex*> v_list;
    qDebug() << edge.size();
    for (int i = 0; i < edge.size(); i++)
    {
        QPair<int,int> p = edge[i];
        int v1 = p.first, v2 = p.second;
        if (!v_list.contains(v1))
        {
            Vertex * v = new Vertex;
            v->setIndex(v1);
            v_list.insert(v1,v);
            myVertexList.append(v);
        }
        if (!v_list.contains(v2))
        {
            Vertex * v = new Vertex;
            v->setIndex(v2);
            v_list.insert(v2,v);
            myVertexList.append(v);
        }
    }

    for (int i = 0; i < edge.size(); i++)
    {
        QPair<int,int> p = edge[i];
        int v1 = p.first, v2 = p.second;
        Vertex * from = v_list.value(v1);
        Vertex * to = v_list.value(v2);
        Edge *edge = new Edge(from, to, i);
        myEdgeList.append(edge);
    }

    graphIsReady = true;
    qDebug() << "DONE HASHING" << myVertexList.size() << myEdgeList.size();
   // reindexing();
   // reindexing_ground_truth();
    read_large_ground_truth_communities();
//    qDebug() << "Ground Truth Q:" << calculate_modularity_for_ground_truth();
}

/** REINDEXING THE SNAP GRAPH (THE PROBLEM IS THAT THE INDEX IS NOT CONTINOUS)
 *  SO A VERTEX HAS 2 INDICES: SNAP INDICIES and DUMEX INDICES
 *  ALGORITHM RUNS USING DUMEX INDICES
 *  CONTENT MATCHING USING SNAP INDICIES (FOR GROUND TRUTH)
 * @brief MainWindow::reindexing
 */
void MainWindow::reindexing()
{
    qDebug() << "Reindexing Started ...";
    qDebug() << "Writing Edge!";
    QString edgePath = "C:/Users/Dumex/Desktop/SocialNetworksCollection/SNAP_DumexTemplate/edge_file.txt";
    QString originalIndexPath = "C:/Users/Dumex/Desktop/SocialNetworksCollection/SNAP_DumexTemplate/vertex_file.txt";
    QFile file(edgePath);
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    //edge are Source Target Seperated by tab \t
    //begin writing edge
    for (int i = 0; i < myEdgeList.size(); i++)
    {
        Edge * e = myEdgeList.at(i);
        Vertex * from = e->fromVertex();
        Vertex * to = e->toVertex();
        int dumex_v = myVertexList.indexOf(from);
        int dumex_u = myVertexList.indexOf(to);
        if (dumex_v == dumex_u)
        {
            qDebug() << "Error: Self Loop Edge";
            return;
        }
        out << dumex_v << "\t" << dumex_u << endl;
    }
    file.close();
    qDebug() << "Now Writing Vertex!";
    //write the node
    //Dumex Index - Original Index seperate by tab \t
    QFile vFile(originalIndexPath);
    vFile.open(QFile::WriteOnly | QFile::Text);
    QTextStream out2(&vFile);
    for (int i = 0; i < myVertexList.size(); i++)
    {
        Vertex * v = myVertexList.at(i);
        int origin_index = v->getIndex();
        out2 << i << "\t" << origin_index << endl;
    }
    vFile.close();
    qDebug() << "FINISHED! Check Files and Move to Correct Location";
}

/** REINDEXING THE GROUND TRUTH
 * @brief MainWindow::reindexing_ground_truth
 */
void MainWindow::reindexing_ground_truth()
{
    read_large_ground_truth_communities();
    // REINDEXING THE GROUND TRUTH FILE
    QString truthPath = "C:/Users/Dumex/Desktop/SocialNetworksCollection/SNAP_DumexTemplate/truth_file.txt";
    QFile truthfile(truthPath);
    truthfile.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&truthfile);
    qDebug() << "Reindexing the SNAP Ground Truth";
    //community are seperate by \n
    //memebrs of community are seperated by \t
    QMap<int, int> map;
    for (int i = 0; i < myVertexList.size(); i++)
        map.insert(myVertexList[i]->getIndex(), i); //snap_id  -> dumex_id
    //precheck ground_truth

    for (int i = 0; i < ground_truth_communities.size(); i++)
    {
        QList<int> c = ground_truth_communities[i]; //SNAP index
        for (int j = 0; j < c.size(); j++)
        {
            //find the appropriate dumex index
            int snap_id = c[j];
            if (map.contains(snap_id))
            {
                out << map.value(snap_id) << "\t";
            }
            else
            {
                qDebug() << "KEY NOT FOUND; TERMINATING!";
                return;
            }
        }
        out << endl;
    }
    truthfile.close();
    qDebug() << "DONE!";
}


/** READ SNAP GROUND TRUTH
 *
 * @brief MainWindow::read_large_ground_truth_communities
 */

void MainWindow::read_large_ground_truth_communities()
{
    qDebug() << "PARSING GROUND TRUTH COMMUNITIES";
    ground_truth_communities.clear();
    QString filePath = QFileDialog::getOpenFileName(this,
                                                     ("Select the GROUND TRUTH FILE")
                                                     , "C:/Users/Dumex/Desktop"
                                                     , ("Text Files (*.txt *csv)"));
    QFile file(filePath);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);
    int max = -1;
    while (!in.atEnd())
    {
        QStringList str = in.readLine().split('\n'); //a community
        QStringList sub_str = str[0].split('\t'); //community split by '\t'
        QList<int> community;
        for (int i = 0; i < sub_str.size(); i++)
        {
            bool ok;
            int indx = sub_str[i].toInt(&ok); //get vertex index
            if (ok)
            {
                community.append(indx);
                if (indx > max)
                    max = indx;
            }
        }
        ground_truth_communities.append(community);
    }
    //readjust if index does not start from 0
    if (max > myVertexList.size())
        qDebug() << "Did you load the SNAP or DUMEX file?";
    //check sum
    int n = 0;
    for (int i = 0; i < ground_truth_communities.size(); i++)
        n+=ground_truth_communities[i].size();

    if (n == myVertexList.size())
    {
        qDebug() << "OK!";
    }
    qDebug() << "FINISHED! Number of Comm: " << ground_truth_communities.size();
    file.close();
}


/** Require Clarification
 * We Process The OverLap vertices
 * for now, for overlap vertices, retain each vertex in the largest community
 * @brief MainWindow::large_process_overlap
 */
void MainWindow::large_process_overlap()
{
    QMap<int,int> map;
    for (int i = 0; i < ground_truth_communities.size(); i++)
    {
        QList<int> c = ground_truth_communities[i];
        for(int j = 0; j < c.size(); j++)
        {
            int id = c[j];
            if (!map.contains(id))
                map.insert(id,i);
            else
                map.insertMulti(id, i);
        }
    }
    QList<int> keys = map.keys();
    for (int k = 0; k < keys.size(); k++)
    {
        int id = keys[k];
        QList<int> comms = map.values(id);
        if (comms.size() > 1) //belong to more than 1 community
        {
            int largest_comm_size = -1, chosen_comm = -1;
            for (int i = 0; i < comms.size(); i++)
            {
                int cid = comms[i];
                int c_size = ground_truth_communities[cid].size();
                if (c_size > largest_comm_size)
                {
                    largest_comm_size = c_size;
                    chosen_comm = i;
                }
            }

            for(int i = 0; i < comms.size(); i++)
            {
                if (i != chosen_comm)
                {
                    ground_truth_communities[comms[i]].removeOne(id);
                }
            }
        }
    }
    //final check
    int n = 0;
    for (int i = 0; i < ground_truth_communities.size(); i++)
        n += ground_truth_communities[i].size();
    qDebug() << n;
}

// ----------------------------------------------------- PREDEFINED EXPERIMENTS ------------------------------
// -----------------------------------------------------------------------------------------------------------
/** Girvan and Newman Experiments with Gnp
 * fix n = 128, m = 4
 * pin is the probability for intra edge
 * pout is the probability for inter edge
 * p_in geq p_out
 * zout = pin/pout
 * fix average_average degree at 16;
 * @brief MainWindow::GN_experiment
 */
void MainWindow::GN_experiment()
{
    GRAPHICS = false;
    int n = 128, m = 4, v_per_c = n/m, average_d = 16;
    //fix z_out, then determine p_in and p_out
    //create new file
    QString filePath("C:/Users/Dumex/Desktop/GN_experiment.txt");
    QFile file(filePath);
    file.open(QIODevice::Append | QIODevice::Text);
    QTextStream out(&file);
    out.setRealNumberPrecision(5);
    out << " ***************************************** \n";
    out << " ------* A NEW RUN STARTS FROM BELOW *---- \n";
    out << "RAND\tJACCARD\tARI\tQ\tZ_out\n";
    file.close();
    for (int k = 19; k <= 19; k++)
    {
        for (int z_out = 0 ; z_out <= 12; z_out++)
        {
            //generate GN graph
            //run algorithm here
            double p_in = 0.0, p_out = 0.0;
            p_out = (double) z_out/96;
            p_in = (double) (16-z_out)/31;

            int times = 1;
            double RAND = 0.0 , JACCARD = 0.0, ARI = 0.0, Q = 0.0;
            QString mess;
            if (k == 0) continue;
            else if (k == 1) {mess.append(QString("********** I.a ************ \n"));}
            else if (k == 2) {mess.append(QString( "********** I.b ************ \n"));}
            else if (k == 3) {mess.append(QString( "********** I.c ************ \n"));}
            else if (k == 4) {mess.append(QString( "********** II.a ************ \n"));}
            else if (k == 5){mess.append(QString( "********** II.a(i) ************ \n"));}
            else if (k == 6) {mess.append(QString( "********** II.b ************ \n"));}
            else if (k == 7){mess.append(QString( "********** II.b(i) ************ \n"));}
            else if (k == 8) {mess.append(QString( "********** II.c ************ \n"));}
            else if (k == 9) {mess.append(QString( "********** II.d ************ \n"));}
            else if (k == 10) {mess.append(QString( "********** II.e ************ \n"));}
            else if (k == 11) {mess.append(QString( "********** II.f ************ \n"));}
            else if (k == 12) {mess.append(QString( "********** II.h ************ \n"));}
            else if (k == 13){mess.append(QString( "********** II.g ************ \n"));}
            else if (k == 14){mess.append(QString( "********** III.a ************ \n"));}
            else if (k == 15){mess.append(QString( "********** III.b ************ \n"));}
            else if (k == 16){mess.append(QString( "********** III.c ************ \n"));}
            else if (k == 17){mess.append(QString( "********** III.d ************ \n"));}
            else if (k == 18){mess.append(QString( "********** III.e ************ \n"));}
            else if (k == 19){mess.append(QString( "********** Betweenness Centrality Clustering ************ \n")); times = 1;}
            else{}
            for (int i = 0 ; i < times; i++)
            {
                create_Gnp_planted(p_in, p_out);
                if (k == 0) continue;
                else if (k == 1) {random_aggregate();}
                else if (k == 2) {random_aggregate_with_degree_comparison(); }
                else if (k == 3) {random_aggregate_with_weight_comparison(); }
                else if (k == 4) {random_aggregate_with_neighbour_initial_degree_bias(); }
                else if (k == 5){random_aggregate_with_neighbour_initial_degree_bias_with_comparison();}
                else if (k == 6) {random_aggregate_with_neighbour_CURRENT_degree_bias(); }
                else if (k == 7){random_aggregate_with_neighbour_CURRENT_degree_bias_with_comparison();}
                else if (k == 8) {random_aggregate_highest_CURRENT_degree_neighbour();}
                else if (k == 9) {random_aggregate_with_minimum_weight_neighbour();}
                else if (k == 10) {random_aggregate_probabilistic_lowest_degree_neighbour_destructive();}
                else if (k == 11) {random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();}
                else if (k == 12) {random_aggregate_greedy_max_weight();}
                else if (k == 13){random_aggregate_greedy_max_degree();}
                else if (k == 14){random_aggregate_retain_vertex_using_triangulation();}
                else if (k == 15){random_aggregate_retain_vertex_using_probabilistic_triangulation();}
                else if (k == 16){random_aggregate_with_highest_triangulated_vertex();}
                else if (k == 17){random_aggregate_retain_vertex_using_triangulation_times_weight();}
                else if (k == 18){random_aggregate_retain_vertex_using_triangulation_of_cluster();}
                else if (k == 19){betweenness_centrality_clutering();}
                QList<double> id = compute_pairwise_matching_efficient();
                RAND+=id[0];
                JACCARD+=id[1];
                ARI+=id[2];
                Q+=calculate_modularity_for_clusters().first;
                myEdgeList.clear();
                for (int i = 0; i < myVertexList.size(); i++)
                    delete myVertexList[i];
                reset();
            }
            file.open(QIODevice::Append | QIODevice::Text);
            QTextStream out(&file);
            mess.append(QString("%1\t%2\t%3\t%4\t%5\n").arg(RAND/times).arg(JACCARD/times).arg(ARI/times).arg(Q/times).arg(z_out));
            out << mess;
            file.close();
            //adjust p_in and p_out
            graphIsReady = false;
        }
    }
}

/** Generate the graph with
 * pin is the probability for intra edge
 * pout is the probability for inter edge
 * zout = pin/pout
 * fix average_average degree at 16;
 * @brief MainWindow::create_GN_graph
 * @param pin
 * @param pout
 * @param average_degree
 */
void MainWindow::create_GN_graph(int z_in, int z_out)
{
    myVertexList.clear();
    myEdgeList.clear();
    int n = 128, m = 4, v_per_c = n/m;
    QList<QPair<int,int> > edges;
    QList<int> intra, inter;
    for (int i = 0; i < n; i++)
        intra << 0;
    for (int i = 0; i < n; i++)
        inter << 0;
    //doing actual work
    for (int i = 0; i < n; i++)
    {
        int clus = i/v_per_c;
        QList<int> in,out;
        for (int j = 0; j < n; j++)
        {
            QPair<int,int> main_e = qMakePair(i,j);
            QPair<int,int> reverse_e = qMakePair(j,i);
            bool exists = (edges.contains(main_e) || edges.contains(reverse_e));
            if (j == i)
                continue;
            else if (j/v_per_c == clus)
            {
                if (!exists && intra[j] < z_in)
                    in << j;
            }
            else
            {
                if (!exists && inter[j] < z_out)
                    out << j;
            }
        }
        int intra_TBA = z_in - intra[i],
            inter_TBA = z_out - inter[i];
        if (in.size() < intra_TBA || out.size()< inter_TBA)
        {
            qDebug() << "- Error While Generating Clusters";
            return;
        }
        QList<int> intra_e = generate_random_edge(intra_TBA, in);
        QList<int> inter_e = generate_random_edge(inter_TBA, out);
        //create intra edges
        for (int j = 0; j < intra_e.size(); j++)
        {
            int neigh = intra_e[j];
            intra[neigh]++;
            intra[i]++;
            QPair<int,int> edge = qMakePair(i,neigh);
            edges.append(edge);
        }
        //create inter edges
        for (int j = 0; j < inter_e.size(); j++)
        {
            int neigh = inter_e[j];
            inter[neigh]++;
            inter[i]++;
            QPair<int,int> edge = qMakePair(i,neigh);
            edges.append(edge);
        }
    }
    NormalGraph g;
    for (int i = 0; i < n; i++)
        boost::add_vertex(g);
    for (int i = 0; i < edges.size(); i++)
        boost::add_edge(edges.at(i).first, edges.at(i).second, g);

    myGraph = g;
    typedef boost::graph_traits<NormalGraph>::vertex_iterator ver_iter;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    //get index map
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);

    for (std::pair<ver_iter, ver_iter> vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(index[*vp.first]);
        v->setPos(QPointF(0.0,0.0));
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
    }

    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = boost::edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        myEdgeList.append(edge);
        num_edge++;
    }
    //add ground truth
    for (int i = 0 ; i < m; i++)
    {
        QList<int> c;
        ground_truth_communities.append(c);
    }

    int c_size = n/m;
    for (int i = 0; i < n; i++)
    {
        int c = i/c_size;
        ground_truth_communities[c].append(myVertexList.at(i)->getIndex());
    }

    graphIsReady = true;
}


/**
 * @brief MainWindow::generate_random_edge
 * @param no_edge: number of edge required
 * @param index_list: list of indices
 * @return
 */

QList<int> MainWindow::generate_random_edge(int no_edge, QList<int> index_list)
{
    QList<int> result;
    for(int i = 0; i < no_edge; i++)
    {
        int size = index_list.size();
        std::uniform_int_distribution<int> distribution(0,size-1);
        int ran = distribution(generator);
        result << index_list[ran];
        index_list.removeAt(ran);
    }
    return result;
}


/** As in Fortunato's
 * @brief MainWindow::create_Gnp_planted
 * @param z_in
 * @param z_out
 */
void MainWindow::create_Gnp_planted(double p_in, double p_out)
{
    int n = 128, m = 4, v_per_c = n/m;
    QList<QPair<int,int> > e;
    std::uniform_real_distribution<double> dis(0,1);
    for (int i = 0; i < n; i++)
    {
        for (int j = i+1; j < n; j++)
        {
            double ran = dis(generator);
            bool edge = false;
            if (i == j)
                continue;
            else if (i/v_per_c == j/v_per_c)
            {
                if (ran <= p_in)
                {
                    edge = true;
                }
            }
            else
            {
                if (ran <= p_out)
                {
                    edge = true;
                }
            }

            if (edge)
            {
                QPair<int,int> edge = qMakePair(i,j);
                e.append(edge);
            }
        }
    }

    NormalGraph g;
    for (int i = 0; i < n; i++)
        boost::add_vertex(g);
    for (int i = 0; i < e.size(); i++)
        boost::add_edge(e.at(i).first, e.at(i).second, g);
    e.clear();
    myGraph = g;
    typedef boost::graph_traits<NormalGraph>::vertex_iterator ver_iter;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    //get index map
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);

    for (std::pair<ver_iter, ver_iter> vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(index[*vp.first]);
        v->setPos(QPointF(0.0,0.0));
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
    }

    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = boost::edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        myEdgeList.append(edge);
        num_edge++;
    }
    //add ground truth
    for (int i = 0 ; i < m; i++)
    {
        QList<int> c;
        ground_truth_communities.append(c);
    }

    for (int i = 0; i < n; i++)
    {
        int c = i/v_per_c;
        ground_truth_communities[c].append(myVertexList.at(i)->getIndex());
    }
    graphIsReady = true;
}



/** Generate A Nested Clustering
 * @brief MainWindow::create_Gnp_nest
 * @param levels: total number of hierarchy levels
 * @param p: a list of p for each level
 * the lowest is level is the most dense cluster, p0 > p1 > p2 ... pk
 */

void MainWindow::create_Gnp_nest(int levels, QList<double> p)
{
    if (levels != p.size())
    {
        qDebug() << "- Error While Generating Nested Gnp: Levels and number of p does not match!\n Terminating ...";
        return;
    }
    int n = 128, m = 4, v_per_c = n/m;
    QList<QPair<int,int> > e;
    std::uniform_real_distribution<double> dis(0,1);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            double ran = dis(generator);
            bool edge = false;
            for (int k = 1; k < levels; k++)
            {
                if ( (j/v_per_c*k) == (i/v_per_c*k) ) //intra
                {
                    if (ran <= p[k-1])
                    {
                        edge = true;
                        break;
                    }
                }
                if (edge)
                {
                    QPair<int,int> edge = qMakePair(i,j);
                    QPair<int,int> reverse_edge = qMakePair(j,i);
                    if (e.contains(edge) || e.contains(reverse_edge)) //to avoid duplicate edge as boost does not handle this
                    {
                        // this edge already exists;
                    }
                    else
                    {
                        e.append(edge);
                        e.append(reverse_edge);
                    }
                }
            }
        }
    }

    NormalGraph g;
    for (int i = 0; i < n; i++)
        boost::add_vertex(g);
    for (int i = 0; i < e.size(); i+=2)
        boost::add_edge(e.at(i).first, e.at(i).second, g);
    myGraph = g;
    typedef boost::graph_traits<NormalGraph>::vertex_iterator ver_iter;
    typedef boost::graph_traits<NormalGraph>::edge_iterator edge_iter;
    //get index map
    typedef boost::property_map<NormalGraph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);

    for (std::pair<ver_iter, ver_iter> vp = vertices(g); vp.first != vp.second; ++vp.first)
    {
        Vertex *v = new Vertex;
        v->reSize(VERTEX_GLOBAL_SIZE);
        v->setIndex(index[*vp.first]);
        v->setPos(QPointF(0.0,0.0));
        v->setOriginPos(v->pos());
        originVertexPos.append(v->pos());
        myVertexList.push_back(v);
    }

    int num_edge = 0;
    for (std::pair<edge_iter, edge_iter> ep = boost::edges(g); ep.first != ep.second; ++ep.first)
    {
        boost::graph_traits<NormalGraph>::edge_descriptor e_desc;
        e_desc = *ep.first;
        boost::graph_traits<NormalGraph>::vertex_descriptor u, v;
        u = source(e_desc, g);
        v = target(e_desc, g);
        Vertex *v1 = myVertexList.at(index[u]);
        Vertex *v2 = myVertexList.at(index[v]);
        Edge *edge = new Edge(v1, v2, num_edge);
        myEdgeList.append(edge);
        num_edge++;
    }
    //add ground truth
    for (int i = 0 ; i < m; i++)
    {
        QList<int> c;
        ground_truth_communities.append(c);
    }

    int c_size = n/m;
    for (int i = 0; i < n; i++)
    {
        int c = i/c_size;
        ground_truth_communities[c].append(myVertexList.at(i)->getIndex());
    }
}
