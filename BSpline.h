#ifndef BSPLINE_H
#define BSPLINE_H

#include <QVector>

class BSpline
{
public:
    BSpline();
    BSpline(const QVector<qreal> &knots, int order);

    void SetOrder(int newOrder);
//    qreal GetOrder() const {return order;}
    QVector<qreal> BasisFunc(qreal u, int span);
    QVector<QVector<qreal> > DerivativeBasis(qreal u, int span, int diffOrder);
    int FindSpan(qreal u);

//    static BSpline FromInterpolation(int order, const QVector<QVector3D> &points);

private:
    QVector<qreal> knots;
    int order;

    friend class BSplineCurve;
    friend class BSplineSurface;
};

#endif // BSPLINE_H
