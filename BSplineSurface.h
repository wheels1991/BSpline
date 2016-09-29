#ifndef BSPLINESURFACE_H
#define BSPLINESURFACE_H
#include <QVector>
#include <QVector3D>
#include <qmath.h>
#include <BSpline.h>
/*!
 * \brief B样条曲面
 */
class BSplineSurface
{
public:
    BSplineSurface();
    BSplineSurface(int orderU, int orderV, QVector<qreal> knotVectorU, QVector<qreal> knotVectorV, QVector<QVector<QVector3D > > controlPoints);
    BSplineSurface(const BSpline &u, const BSpline &v, QVector<QVector<QVector3D> > contrlPoints);
    qreal PointOnSurface(qreal u, qreal v, Qt::Axis axis);
    QVector< QVector<qreal> > SurfaceDerivative(qreal u, qreal v, int diffOrder, Qt::Axis axis);
    QVector3D NormVector(qreal u, qreal v);

    static BSplineSurface FromInterpolation(int orderU, int orderV, const QVector<QVector<QVector3D> > &points);

private:
    bool isValid;
    QVector<QVector<QVector3D> > ctrlPoints;          //外层为u方向，内层为v方向, 最内层为点的维数
    BSpline bsplineU;
    BSpline bsplineV;
};

#endif // BSPLINESURFACE_H
