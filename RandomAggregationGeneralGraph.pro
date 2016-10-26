#-------------------------------------------------
#
# Project created by QtCreator 2015-10-23T12:14:33
#
#-------------------------------------------------

QT       += core gui
QT       += xml

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets


TARGET = RandomAggregationGeneralGraph
TEMPLATE = app

CONFIG += c++11

SOURCES += main.cpp\
        mainwindow.cpp \
    vertex.cpp \
    graphview.cpp \
    edge.cpp \
    arrow.cpp \
    clustercentroid.cpp \
    clustervertex.cpp \
    infowidget.cpp

HEADERS  += mainwindow.h \
    vertex.h \
    graphview.h \
    edge.h \
    GraphGenerator.h \
    arrow.h \
    clustercentroid.h \
    clustervertex.h \
    lineanimator.h \
    infowidget.h

FORMS    += mainwindow.ui \
    infowidget.ui


INCLUDEPATH += "C:\Boost\boost_1_56_0"

