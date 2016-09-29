#ifndef PATHPLANNING_H
#define PATHPLANNING_H
#include "BSplineCurve.h"
#include "BSplineSurface.h"
/*!
 * \brief 对一个给定的B样条曲面，做路径规划
 *        路径规划的结果是给出TCP的轨迹，以csv的形式输出
 */
class PathPlanning
{
public:
    PathPlanning();
    //等弧长进给、等弧长细分
    PathPlanning(BSplineSurface _surface, qreal feedInterval, qreal scanInterval,
                 qreal roughScanInterval, qreal incidence, qreal azimuth);
    //等参数进给、等弧长细分
    PathPlanning(BSplineSurface _surface, qreal feedInterval, qreal scanInterval,
                 qreal incidence, qreal azimuth);

    QVector<QVector<QVector3D> > GetPath();
    QVector<QVector<QVector3D> > GetPathDirection();
private:
    QVector<QVector<qreal> > feed(QVector<QVector<qreal> > pointsList, qreal feed);
    QVector<QVector<qreal> > CurveSubdivide(QVector<QVector<qreal> > pointsList, qreal interval);


private:
    BSplineSurface surface;
    QVector<QVector<QVector3D> > path;      //外层为扫描线脚标，内层为扫描线对应的点序列脚标
    QVector<QVector<QVector3D> > pathDirection;
};

#endif // PATHPLANNING_H
