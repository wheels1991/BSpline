#include <QCoreApplication>
#include <QDebug>
#include "BSplineCurve.h"
#include "BSplineSurface.h"
#include "PathPlanning.h"
#include "MyMath.h"

#include <QTime>


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    //定义插值点
    int nu = 10;
    int nv = 10;
    int p = 3;
    int q = 3;
//    int dimension = 2;
    QVector<QVector<QVector3D> > interpolationPoints(nu + 1,
            QVector<QVector3D>(nv + 1, QVector3D(0.0, 0.0, 0.0)));
    QVector<qreal> x(nu+1, 0.0);
    QVector<qreal> y(nv+1, 0.0);
    for (int i = 0; i <= nu; i++) {
        x[i] = 200.0 / nu * i;
    }
    for (int i = 0; i <= nv; i++) {
        y[i] = 200.0 / nv * i;
    }
    for (int i = 0; i <= nu; i++) {
        for (int j =0; j <= nv; j++) {
            interpolationPoints[i][j][0] = x[i];
            interpolationPoints[i][j][1] = y[j];
            interpolationPoints[i][j][2] = qSin( x[i] * M_PI / 100.0 )
                                           * qExp( -3.0 * x[i] /200.0 )
                                           * qSin( y[j] * M_PI / 200.0 )
                                           * 100.0;
        }
    }

    BSplineSurface surface = BSplineSurface::FromInterpolation(p, q, interpolationPoints);
    qreal u = 0.5;
    qreal v = 0.3;
    qDebug() << "the point on the surface is: " << surface.PointOnSurface(u, v, Qt::ZAxis);
    qDebug() << "the derivs on the surface point is: " << surface.SurfaceDerivative(u, v, 2, Qt::ZAxis);
    qDebug() << "the norm vector is: " << surface.NormVector(u,v);
    qreal feedInterval = 20;
    qreal scanInterval = 5;
    qreal roughScanInterval = 10;
    qreal incidence = 0;
    qreal azimuth = 0;

    PathPlanning pathPlanning(surface, feedInterval, scanInterval, incidence,azimuth);
    QVector<QVector<QVector3D> > path = pathPlanning.GetPath();
    QVector<QVector<QVector3D> > pathNorm = pathPlanning.GetPathDirection();
    qDebug() << "计算完成";
    return a.exec();
}
