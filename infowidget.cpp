#include "infowidget.h"
#include "ui_infowidget.h"

#include <QDebug>
#include <QPalette>

InfoWidget::InfoWidget(QWidget *parent,
                       int G_num_vertex, int G_num_edge) :
    QWidget(parent),
    ui(new Ui::InfoWidget)
{
    ui->setupUi(this);
    //reset styleq
    for (int i = 0; i < ui->graph_info->count(); i++)
    {
        QLayoutItem * item = ui->graph_info->itemAt(i);
        QWidget * widget = item->widget();
        QLCDNumber * lcd = qobject_cast<QLCDNumber*>(widget);
        if (lcd)
        {
            setLCDdisplay(lcd);
            lcd->setSegmentStyle(QLCDNumber::Flat);
        }
    }
    for (int i = 0; i < ui->vertex_info->count(); i++)
    {
        QLayoutItem * item = ui->vertex_info->itemAt(i);
        QWidget * widget = item->widget();
        QLCDNumber * lcd = qobject_cast<QLCDNumber*>(widget);
        if (lcd)
        {
            setLCDdisplay(lcd);
            lcd->setSegmentStyle(QLCDNumber::Flat);
        }
    }
    for (int i = 0; i < ui->agg_info->count(); i++)
    {
        QLayoutItem * item = ui->vertex_info->itemAt(i);
        QWidget * widget = item->widget();
        QLCDNumber * lcd = qobject_cast<QLCDNumber*>(widget);
        if (lcd)
        {
            setLCDdisplay(lcd);
            lcd->setSegmentStyle(QLCDNumber::Flat);
        }
    }
    setLCDdisplay(ui->v_index);
    ui->v_index->setSegmentStyle(QLCDNumber::Flat);
    //set value
    ui->graphInfo_numVertex->display(G_num_vertex);
    ui->graphInfo_numEdge->display(G_num_edge);
}


InfoWidget::~InfoWidget()
{
    delete ui;
}

void InfoWidget::display_vertex_selected_info(QList<int> info)
{ //info << index << degree << weight;
    ui->v_index->display(info[0]);
    ui->v_degree_origin->display(info[1]);
    ui->v_w_origin->display(info[2]);
}

void InfoWidget::update_vertex_index(const int &i)
{
    ui->v_index->display(i);
}

void InfoWidget::update_vertex_degree_origin(const int &i)
{
    ui->v_degree_origin->display(i);
}

void InfoWidget::update_vertex_weight_origin(const int &i)
{
    ui->v_w_origin->display(i);
}


void InfoWidget::aggregation_selected_vertex(QList<int> info)
{
    ui->agg_selected_degree->display(info[0]);
    ui->agg_selected_weight->display(info[1]);
}

void InfoWidget::aggregation_absorbed_vertex(QList<int> info)
{
    ui->agg_absorbed_degree->display(info[0]);
    ui->agg_absorbed_weight->display(info[1]);
}

void InfoWidget::setLCDdisplay(QLCDNumber *lcd)
{
    lcd->setAutoFillBackground(true);
    QPalette Pal = lcd->palette();
    Pal.setColor(QPalette::Normal, QPalette::WindowText, Qt::black);
    Pal.setColor(QPalette::Normal, QPalette::Window, Qt::lightGray);
    lcd->setPalette(Pal);

}

