/****************************************************************************
**
** Copyright (C) 2010 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial Usage
** Licensees holding valid Qt Commercial licenses may use this file in
** accordance with the Qt Commercial License Agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and Nokia.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 2.1 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU Lesser General Public License version 2.1 requirements
** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** In addition, as a special exception, Nokia gives you certain additional
** rights.  These rights are described in the Nokia Qt LGPL Exception
** version 1.1, included in the file LGPL_EXCEPTION.txt in this package.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3.0 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU General Public License version 3.0 requirements will be
** met: http://www.gnu.org/copyleft/gpl.html.
**
** If you have questions regarding the use of this file, please contact
** Nokia at qt-info@nokia.com.
** $QT_END_LICENSE$
**
****************************************************************************/

#include <QtGui>

#include "arrow.h"
#include <math.h>

const qreal Pi = 3.14;

Arrow::Arrow(Vertex *startV, Vertex *endV)  :
    Edge(startV, endV, 0)
{
    myFromVertex = startV;
    myToVertex = endV;

    myFromVertex->addEdge(this);
    myToVertex->addEdge(this);

    setFlag(QGraphicsItem::ItemIsSelectable, false);
    myColor = Qt::black;
    setPen(QPen(myColor, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

    setVisible(true);
}

Arrow::~Arrow()
{
    myFromVertex->removeEdge(this);
    myToVertex->removeEdge(this);
}

QRectF Arrow::boundingRect() const
{
  //  qreal extra = (pen().width() + 1) / 2.0;
    qreal extra = 0.0;
    return QRectF(line().p1(), QSizeF(line().p2().x() - line().p1().x(),
                                      line().p2().y() - line().p1().y()))
        .normalized()
        .adjusted(-extra, -extra, extra, extra);
}

QPainterPath Arrow::shape() const
{
    QPainterPath path = QGraphicsLineItem::shape();
    path.addPolygon(arrowHead);
    return path;
}

void Arrow::setColour(QColor cor)
{
    myColor = cor;
    update();
}

void Arrow::changePen(QColor cor, int size)
{
    myColor = cor;
    setPen(QPen(myColor, size, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    update();
}

void Arrow::updatePosition()
{
    QLineF line(mapFromItem(myFromVertex, 0, 0), mapFromItem(myToVertex, 0, 0));
    setLine(line);
}

void Arrow::paint(QPainter *painter, const QStyleOptionGraphicsItem *,
          QWidget *)
{
    if (myFromVertex->collidesWithItem(myToVertex))
        return;

    QPen myPen = pen();
    myPen.setColor(myColor);
    qreal arrowSize = 15;
    painter->setPen(myPen);
    painter->setBrush(myColor);

    int size = myToVertex->getSize();
    size *= 1.25; //readjust the bounding rectangle
    QPointF startPos(myFromVertex->x(), myFromVertex->y());
    QPointF endPos(myToVertex->x(), myToVertex->y());
    QLineF centerLine(startPos,endPos);
    QRectF endRectF(myToVertex->x() - size/2, myToVertex->y()-size/2, size, size);
    QPolygonF endPolygon(endRectF);
  //  endPolygon.removeLast();
    QPointF p1 = endPolygon.first();// + myToVertex->pos();
    QPointF p2;
    QPointF intersectPoint;
    QLineF polyLine;
    for (int i = 1; i < endPolygon.count(); ++i)
    {
        p2 = endPolygon.at(i) ; // + myToVertex->pos();
        polyLine = QLineF(p1, p2);
        QLineF::IntersectType intersectType =
            polyLine.intersect(centerLine, &intersectPoint);
        if (intersectType == QLineF::BoundedIntersection)
            break;
            p1 = p2;
    }

    setLine(QLineF(intersectPoint, startPos));

    double angle = ::acos(line().dx() / line().length());
    if (line().dy() >= 0)
        angle = (Pi * 2) - angle;

        QPointF arrowP1 = line().p1() + QPointF(sin(angle + Pi / 3) * arrowSize,
                                        cos(angle + Pi / 3) * arrowSize);
        QPointF arrowP2 = line().p1() + QPointF(sin(angle + Pi - Pi / 3) * arrowSize,
                                        cos(angle + Pi - Pi / 3) * arrowSize);

        arrowHead.clear();
        arrowHead << line().p1() << arrowP1 << arrowP2;
        painter->drawLine(line());
        painter->drawPolygon(arrowHead);
        if (isSelected()) {
            painter->setPen(QPen(myColor, 1, Qt::DashLine));
        QLineF myLine = line();
        myLine.translate(0, 4.0);
        painter->drawLine(myLine);
        myLine.translate(0,-8.0);
        painter->drawLine(myLine);
    }
}
