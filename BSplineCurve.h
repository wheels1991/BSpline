#ifndef BSPLINECURVE_H
#define BSPLINECURVE_H
#include <QVector>
#include <QVector3D>
#include <qmath.h>
#include <BSpline.h>

/*!
 * \brief B样条曲线
 */
class BSplineCurve
{
public:
    BSplineCurve() {}
    BSplineCurve(int _order, QVector<qreal> knotVector, QVector<QVector3D > controlPoints);
    BSplineCurve(BSpline bs, QVector<QVector3D> controlPoints);
    qreal PointOnCurve(qreal u, Qt::Axis axis);
    QVector<qreal> CurveDerivs(qreal u, int diffOrder, int p, Qt::Axis axis);
    static BSplineCurve FromInterpolation(int order, const QVector<QVector3D> &interpolationPoints);
    static BSplineCurve FromInterpolation(int order, QVector<QVector3D> interpolationPoints, QVector3D D0, QVector3D Dn);
    static BSplineCurve FromInterpolation(int order, QVector<QVector3D> interpolationPoints, QVector3D D0, QVector3D Dn, QVector3D A0, QVector3D An);
    static QVector<QVector3D> CurveSubdivide(const QVector<QVector3D> &points, const int numOfSubdivide);

private:
    QVector<QVector3D> ctrlPoints;
    BSpline bspline;
};

#endif // BSPLINECURVE_H
