#ifndef INFOWIDGET_H
#define INFOWIDGET_H

#include <QWidget>

class QLCDNumber;

namespace Ui {
class InfoWidget;
}

class InfoWidget : public QWidget
{
    Q_OBJECT

public:
    explicit InfoWidget(QWidget *parent = 0,
                        int G_num_vertex = 0, int G_num_edge = 0);
    ~InfoWidget();

public slots:
    void display_vertex_selected_info(QList<int> info);
    void update_vertex_index(const int &i);
    void update_vertex_degree_origin(const int &i);
    void update_vertex_weight_origin(const int &i);

    void aggregation_selected_vertex(QList<int> info);
    void aggregation_absorbed_vertex(QList<int> info);
private:
    void setLCDdisplay(QLCDNumber * lcd);
    Ui::InfoWidget *ui;
};

#endif // INFOWIDGET_H
