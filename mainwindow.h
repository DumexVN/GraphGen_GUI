#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QCoreApplication>
#include <QGraphicsView>
#include <QPainter>
#include <QInputDialog>
#include <QtCore>
#include <QTimer>

#include "vertex.h"
#include "edge.h"
#include "graphgenerator.h"
#include "lineanimator.h"

class InfoWidget;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();

private slots:// private slot functions
    void generateLayerGnp();
    void generateSimpleClusterGnp();
    void generateCycle();
    void generateArtificialSocialNetwork();
    void generateGNArtificialNetwork();
    void generateErdosReyni();
    void generateKarateClub();
    void generateWeightedKarateClub();
    void generateSimpleSample();
    void read_seperated_graph_input();
    void load_ground_truth_communities();
    void read_large_graph_with_ground_truth_communities();
    void parseGMLfile();
    void draw_aggregation_result();
    void draw_aggregation_result_by_shape();
    void draw_dense_graph_aggregation_result();
    void draw_aggregation_retain_vertex_result();
    void draw_hierarchy_tree();
    void draw_walked_directed_graph();
    //investigate bridges
    void get_bridge_stats();

    //Betweeness Centrality Clustering
    void betweenness_centrality_clutering();
    //aggregation
    void random_aggregate();
    void random_aggregate_n_steps_walk();
    void random_walking_constant_restart_remove_edges();
    void random_walking_normal_retain_edges();
    void random_walking_testing();
    //
    void random_aggregate_with_degree_comparison();
    void random_aggregate_with_weight_comparison();
    void random_aggregate_with_neighbour_initial_degree_bias();
    void random_aggregate_with_neighbour_initial_degree_bias_with_comparison();
    void random_aggregate_with_neighbour_CURRENT_degree_bias();
    void random_aggregate_with_neighbour_CURRENT_degree_bias_with_comparison();
    void random_aggregate_highest_CURRENT_degree_neighbour();
    void random_aggregate_with_minimum_weight_neighbour();
    void random_aggregate_probabilistic_lowest_degree_neighbour_destructive();
    void random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();
    void random_aggregate_with_highest_edge_weight_and_weight_comparison();
    void random_aggregate_with_edge_weight_bias_and_weight_comparison();
    void random_aggregate_with_highest_triangulated_vertex();
    void random_aggregate_greedy_max_degree();
    void random_aggregate_greedy_max_weight();
    //agg without 'removing' vertices
    void random_aggregate_retain_vertex_using_triangulation();
    void random_aggregate_retain_vertex_using_probabilistic_triangulation();
    void random_aggregate_retain_vertex_using_triangulation_times_weight();
    void random_aggregate_retain_vertex_using_triangulation_of_cluster();
    //random edge removal
    void random_edge_removal();
    //random aggregation with simulate vertex removal
    void typeI_retaining_vertex();
    void typeII_retaining_vertex();
    void typeIIIc_retaining_vertex();
    void typeIIId_removing_vertex();
    void typeIV_retaining_vertex();
    //
    void compareResultWithKarateClubData();

    //test with multiple runs
    void multirun_aggregate();
    void GN_experiment();
    void multirun_aggregate_edge_weight_graph();
    void reset();
    void exportGML();
    //indices
    void calculate_indices();
    void calculate_Q();

private:
    void create_GN_graph(int z_in, int z_out);
    void create_Gnp_planted(double z_in, double z_out);
    void create_Gnp_nest(int levels, QList<double> p);
    QList<int> generate_random_edge(int no_edge, QList<int> index_list);
    void read_ground_truth_communities();
    void colouring_ground_truth_communities();
    void read_string_ground_truth();
    void export_string_ground_truth_to_id_ground_truth(const QString &fileDir);
    void generateErdosReyni(double p);
    void debugOriginalGraph(NormalGraph g);
    //random aggregation with weight update
    void typeII_with_weight_update();
    void typeII_probabilistic_current_degree_uniform_neighbour();
    void typeII_probabilistic_degree_smallest_neighbour();
    void typeIII_with_weight_update();
    void typeIII_probabilistic_with_weight_update();
    void typeIV_with_weight_update();
    void typeV_with_weight_update();
    void typeVI_with_weight_update();
    void typeVII_with_weight_update();
    void typeVIII_with_weight_update();
    QPair<double,int> calculate_modularity_for_clusters();
    double calculate_modularity_for_ground_truth();
    void createSnapShot(int t);
    void readjust_graph(Edge * mainE, QList<int> subE, Vertex * delete_vertex,int step);
    //parse result without graphic factors
    void parse_aggregation_result();
    void parse_retain_result();

    //create Actions and Menu
    void setVertexDetails();
    void setEdgeWeight();
    void writeDetailRun(int type);
    void writeModularityReport(QString fileName, QList<QPair<double,int> > Q);
    void writeCommunityMatchingReport(QString fileName, QList<QList<double> > Indices);
    bool checkGraphCondition();
    void createActions();
    void createMenus();
    void setUpGraphInfoWidget();
    void createSeperateScene();
    void reConnectGraph();
    void resetGraphics();
    void saveAggregationResultInSeperateWindow();
    QPair<int,int> compareResultWithKarateClubDataStats();
    QPair<int,int> persistentGraphCompareWithKarate();
    //comput measures
    //pairwise
    int CC();
    QList<double> compute_pairwise_matching_efficient();
    void compute_cluster_matching(QList<QList<int> > clusters);
    double compute_RAND_index(QList<QList<int> > result_clusters);
    double compute_Jaccard_index(QList<QList<int> > result_clusters);
    quint64 calA(QList<quint64> param);
    quint64 calB(QList<quint64> param);
    quint64 calC(QList<quint64> param);
    quint64 calD(QList<quint64> param);
    double calAdRand(QList<quint64> param);
    //newman fraction
    double compute_Newman_fraction_of_classified(QList<QList<int> > result_clusters);
    bool check_merge(QList<int> result_c, QList<QList<int> > permutation);
    QList<QList<int> > generateP(int size);
    void subset_sum(QList<int> numbers, int target);
    void sum_up_recursive(QList<int> numbers, int target, QList<int> partial);
    QList<int> community_matching(QList<int> truth_c, QList<int> result_c);
    QList<QPair<int,int> > create_same_set(QList<QList<int> > C);
    QList<QPair<int, int> > create_diff_set(QList<QList<int> > C);
    // for large graph
    void reindexing();
    void reindexing_ground_truth();
    void read_large_ground_truth_communities();
    void large_process_overlap();
    void examining_bridges();

    QSequentialAnimationGroup * prepare_animation_for_one_vertex(QList<int> edges, Vertex * abs_v, Vertex * absorber);
    void draw_aggregation_animation();

    void add_components_without_graphic();
    void loadGephiLayout(QString dir);
    void layout_KarateClub();
    void layout_graph(NormalGraph &G);
    void layout_artificial_network(NormalGraph &G, int no_cluster);
    void preserveVertexPointer();

    void exportGML(int take, double p);
    void exportGML(int take, double p, double Q);


    //private variables
    NormalGraph myGraph;
    QList<QPointF> originVertexPos;
    QList<Vertex*> myVertexList;
    QList<Edge*> myEdgeList;
    QList<LineAnimator*> myLine;
    QList<Vertex*> centroids;
    QList<Vertex*> cluster_vertex;
    QList<QString> vertex_label;
    QList<int> vertex_weight;
    QMap<QString, QString> vertex_value;
    QList<int> edge_weight;
    //
    QList<QList<int> > ground_truth_communities;
    QList<QPair<int,int> > hierarchy;
    QList<QList<int> > result;
    QSet<int> large_excluded;
    //
    QGraphicsScene *myScene;
    QGraphicsView *myView;
    QMainWindow * originGraph;
    bool graphIsReady;

    //all menus
    QMenu *generateGraphMenu;
    QMenu *aggregateMenu;
    InfoWidget * myInfoWidget;

    //all actions
    QAction *generateLayerGnpAction;
    QAction *generateCycleAction;
    QAction *generateSimpleClusterGnpAction;
    QAction *generateArtificialSocialNetworkAction;
    QAction *generateGNArtificialNetworkAction;
    QAction *generateErdosReyniAction;
    QAction *generateKarateClubAction;
    QAction *generateWeightedKarateClubAction;
    QAction *generateSmallSampleAction;
    QAction *readInputAction;
    QAction *readSNAPAction;
    QAction *loadGroundTruthCommAction;
    QAction *readGMLFileAction;
    QAction *comparesKarateClubResultAction;
    QAction *BCAction;
    QAction *randomAggregateAction;
    QAction *randomAggregationNStepsWalkAction;
    QAction *randomAggregateWithInitialDegreeBiasAction;
    QAction *randomAggregationWithDegreeComparisonAction;
    QAction *randomAggregationWithWeightComparisonAction;
    QAction *randomAggregationSelectHighestEdgeWeightWithVertexWeightComparisonAction;
    QAction *randomAggregationWithEdgeWeightBiasAndVertexWeightComparisonAction;
    QAction *randomAggregationWithHighestTriangulatedVertexAndWeightComparisonAction;
    QAction *randomAggregationWithMinimumDegreeNeighbourAction;
    QAction *randomAggregationProbabilisticWeightSelectionWithMinimumDegreeNeighbourAction;
    QAction *randomAggregationII_a_i_Action;
    QAction *randomAggregationII_b_i_Action;

    QAction *randomAggregationProbabilisticTriangulationAction;
    QAction *randomAggregationRetainingVertexWithTriangulationCoeffAndWeightComparisonAction;
    QAction *randomAggregationRetainingVertexWithTriangulationxWeightIndexAction;
    QAction *randomAggregationRetainingVertexWithHighestTrianglesClusterAction;
    QAction *randomAggregationWithNeighbourCurrentDegreeBiasedAction;
    QAction *randomAggregateHighestDegreeNeighbourAction;
    QAction *randomAggregateLowestDegreeNeighbourDestructiveAction;

    QAction *randomAggregateMaxDegreeAction;
    QAction *randomAggregateMaxWeightAction;

    QAction *randomWalkingRemovingEdgesAction;
    QAction *randomWalkingNormalAction;
    QAction *randomWalkingTestAction;

    QAction *randomEdgeRemovalAction;

    QAction *typeIretainAction;
    QAction *typeIIretainAction;
    QAction *typeIIIretainAction;
    QAction *typeIIIdremoveAction;
    QAction *typeIVretainAction;

    QAction *GNGraphExperimentAction;

    QAction *resetAction;

    QAction *calculateIndicesAction;
    QAction *calculateQAction;
    //multiple test
    QAction *multipleAggregateAction;
    QAction *multipleAggregateEdgeWeightedGraphAction;
    //empirical
    QAction *examineBridgeAction;
    //for large graph
    QAction *reindexingSNAPAction;
    QAction *readDumexInputAction;
    QAction *LARGERerun;
    QAction *saveGMLAction;
};

#endif // MAINWINDOW_H
