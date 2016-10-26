#ifndef GRAPHGENERATOR_H
#define GRAPHGENERATOR_H

#include <QDebug>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/topology.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/grid_graph.hpp>

typedef boost::adjacency_list<> Graph;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> NormalGraph;
typedef boost::sorted_erdos_renyi_iterator<boost::minstd_rand, NormalGraph> ERGen;

class GraphGenerator
{
public:
    GraphGenerator(){}

    NormalGraph erdos_reyni_generator(int noofVertex, double p)
    {
        boost::minstd_rand gen;
        // Create graph with noofVertex nodes and edges with probability p
        NormalGraph g(ERGen(gen, noofVertex, p), ERGen(), noofVertex);
        return g;
    }

};
#endif // GRAPHGENERATOR_H
