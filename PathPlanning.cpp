#include "PathPlanning.h"
#include "qmath.h"
#include <qdebug.h>

PathPlanning::PathPlanning()
{

}
/*!
 * \brief PathPlanning::PathPlanning    等弧长进给、等弧长细分
 * \param _surface                      b样条曲面
 * \param _feedInterval                 进给量
 * \param _scanInterval                 单条扫描线细分量
 * \param _roughScanInterval            单条扫描线的粗分量
 * \param _incidence                    入射角，声束与法向量的夹角
 * \param _azimuth                      方向角，声束在法向量垂直平面上的投影与切线方向的夹角
 */
PathPlanning::PathPlanning(BSplineSurface _surface, qreal feedInterval, qreal scanInterval, qreal roughScanInterval, qreal incidence, qreal azimuth)
{
    surface = _surface;
//    feedInterval = _feedInterval;
//    scanInterval = _scanInterval;
//    incidence = _incidence;
//    azimuth = _azimuth;
    //以下开始进行路径规划，沿u方向扫描，v方向进给
    //设计第一条扫描线
    QVector<QVector<qreal> > roughScanLine(2, QVector<qreal>(21, 0.0));
    for (int i = 0; i <= 20; i++) {          //第一条扫描线为v = 0
        roughScanLine[0][i] = i / 20.0;
    }

    roughScanLine = CurveSubdivide(roughScanLine, roughScanInterval);                       //对第一条扫描线进行粗分，确定进给参考点

    int ScanLineNo = 0; //测试用数据，用来记录扫描线
    /*不断进给，计算各扫描线*/
    while (1) {
        /*设置退出规划完成的条件*/
        //找到scanLine中v的最小值vmin
        qreal vmin = roughScanLine[1][0];
        for (int i = 1; i < roughScanLine[0].size(); i++) {
            vmin = vmin > roughScanLine[1][i] ? roughScanLine[1][i] : vmin;
        }
        if (vmin >= 1) {     //当粗扫描线的v的最小值不小于1时，说明扫描覆盖了所有曲面，退出规划
             break;
        }
        /*开始进行规划 */
        //扫描细分
        QVector<QVector3D> lineTemp;                                                        //临时存放一条扫描线(x, y, z)
        QVector<QVector3D> normTemp;                                                        //临时存放一条扫描线的法向量
        QVector<QVector3D> directionTemp;                                                   //临时存放一条扫描线的声束指向
        QVector<QVector<qreal> > scanLine = CurveSubdivide(roughScanLine, scanInterval);    /*!< 对当前扫描线(u, v)进行细分，得到输出扫描线scanLine：(u, v)的点序列*/
        for (int i = 0; i <= (scanLine[0].size() - 1); i++) {
            const QVector3D point(surface.PointOnSurface(scanLine[0][i], scanLine[1][i], Qt::XAxis),
                                  surface.PointOnSurface(scanLine[0][i], scanLine[1][i], Qt::YAxis),
                                  surface.PointOnSurface(scanLine[0][i], scanLine[1][i], Qt::ZAxis));  //将点(x, y, z)保存下来
            lineTemp.append(point);                                                         //将点坐标保存到扫描线的点列表中
            normTemp.append(surface.NormVector(scanLine[0][i], scanLine[1][i]));            //将点的法向量保存到扫描线的法向量列表中
            if (i >= 1) {
                QVector3D Ltemp = point - lineTemp[i - 1];                                  //切向量：上一个扫描点到下一个扫描点的向量，采用的差分法计算切向量
                Ltemp.normalize();                                                          //将切向量归一化
                const QVector3D Mtemp = QVector3D::crossProduct(Ltemp, normTemp[i - 1]);    //切向量与法向量的公法向量，这三个向量及对点扫描点构成随动坐标系
                QVector3D direction = qSin(incidence) * (qSin(azimuth) * Mtemp + qCos(azimuth) * Ltemp)
                                 + qCos(incidence) * normTemp[i - 1];                       /*!< 采用入射角和方向角确定声束在随动坐标系中的方向*/
                directionTemp.append(direction);                                            //将扫描点的声束方向添加到扫描线的声束方向列表中
            }
        }
        lineTemp.removeLast();                                                              //移除扫描线中最后一个扫描点，因为该扫描点无法利用差分计算切向量
        path.append(lineTemp);                                                              //将一条新的扫描线添加到路径中
        pathDirection.append(directionTemp);                                                //将一条新的扫描线的声束方向添加到路径声束方向中
        //扫描进给
        roughScanLine = feed(roughScanLine, feedInterval);                                                /*!< 扫描线进给*/
        ScanLineNo++;
        qDebug() << "第" << ScanLineNo << "条线点数为：" << lineTemp.size();
    }
    ScanLineNo++;
    /*最后一条扫描线，即所有点的v=0*/
    QVector<QVector3D> lineTemp;                                                            //临时存放一条扫描线(x, y, z)
    QVector<QVector3D> normTemp;                                                            //临时存放一条扫描线的法向量
    QVector<QVector3D> directionTemp;                                                       //临时存放一条扫描线的声束指向
    QVector<QVector<qreal> > scanLine = CurveSubdivide(roughScanLine, scanInterval);        //对当前扫描线(u, v)进行细分，得到输出扫描线scanLine：(u, v)的点序列
    for (int i = 0; i <= (scanLine[0].size() - 1); i++) {
        const QVector3D point(surface.PointOnSurface(scanLine[0][i], scanLine[1][i], Qt::XAxis),
                              surface.PointOnSurface(scanLine[0][i], scanLine[1][i], Qt::YAxis),
                              surface.PointOnSurface(scanLine[0][i], scanLine[1][i], Qt::ZAxis));  //将点(x, y, z)保存下来
        lineTemp.append(point);                                                             //将点坐标保存到扫描线的点列表中
        normTemp.append(surface.NormVector(scanLine[0][i], scanLine[1][i]));                //将点的法向量保存到扫描线的法向量列表中
        if (i >= 1) {
            QVector3D Ltemp = point - lineTemp[i - 1];                                      //切向量：上一个扫描点到下一个扫描点的向量，采用的差分法计算切向量
            Ltemp.normalize();                                                              //将切向量归一化
            const QVector3D Mtemp = QVector3D::crossProduct(Ltemp, normTemp[i - 1]);        //切向量与法向量的公法向量，这三个向量及对点扫描点构成随动坐标系
            QVector3D direction = qSin(incidence) * (qSin(azimuth) * Mtemp + qCos(azimuth) * Ltemp)
                             + qCos(incidence) * normTemp[i - 1];                           //采用入射角和方向角确定声束在随动坐标系中的方向
                directionTemp.append(direction);                                            //将扫描点的声束方向添加到扫描线的声束方向列表中
            }
        }
        lineTemp.removeLast();                                                              //移除扫描线中最后一个扫描点，因为该扫描点无法利用差分计算切向量
        path.append(lineTemp);                                                              //将一条新的扫描线添加到路径中
        pathDirection.append(directionTemp);                                                //将一条新的扫描线的声束方向添加到路径声束方向中
    qDebug() << "第" << ScanLineNo << "条线点数为：" << lineTemp.size();


}
/*!
 * \brief PathPlanning::PathPlanning    等参数进给、等弧长细分
 * \param _surface                      b样条曲面
 * \param _feedInterval                 进给量
 * \param _scanInterval                 单条扫描线细分量
 * \param _incidence                    入射角，声束与法向量的夹角
 * \param _azimuth                      方向角，声束在法向量垂直平面上的投影与切线方向的夹角
 */
PathPlanning::PathPlanning(BSplineSurface _surface, qreal feedInterval, qreal scanInterval, qreal incidence, qreal azimuth)
{
    surface = _surface;
    /*计算参数曲面沿V方向的最大曲线长度Lvmax，决定沿V方向的细分个数n = Lvmax / JogInterval + 1*/
    const int numTemp = 10;                                                     //沿v方向的曲线的条数，用来计算寻找最大长度Lvmax
    qreal Lvmax = 0;
    for (int i = 0; i <= numTemp; i++) {
        const qreal u = 1.0 * i / numTemp;
        qreal dL = 0;   //临时弧长
        qreal v = 0;
        qreal xtemp1 = surface.PointOnSurface(u, v, Qt::XAxis);
        qreal ytemp1 = surface.PointOnSurface(u, v, Qt::YAxis);
        qreal ztemp1 = surface.PointOnSurface(u, v, Qt::ZAxis);
        //计算一条v方向曲线的长度
        while (v < 1) {
            v = v + 0.0001;
            const qreal xtemp2 = surface.PointOnSurface(u, v, Qt::XAxis);
            const qreal ytemp2 = surface.PointOnSurface(u, v, Qt::YAxis);
            const qreal ztemp2 = surface.PointOnSurface(u, v, Qt::ZAxis);
            dL += qSqrt( qPow(xtemp2 - xtemp1, 2) +
                         qPow(ytemp2 - ytemp1, 2) +
                         qPow(ztemp2 - ztemp1, 2) );
            xtemp1 = xtemp2;
            ytemp1 = ytemp2;
            ztemp1 = ztemp2;
        }
        Lvmax = Lvmax > dL ? Lvmax : dL;                                                        //取最大曲线长度
    }
    //沿v方向进给，计算各扫描线的扫描点，按等弧长细分
    const int numOfScanLine = int(Lvmax / feedInterval);                                        //扫描线条数为numOfScanLine + 1
    for (int i = 0; i <= numOfScanLine; i++) {
        QVector<QVector3D> lineTemp;                                                            //临时存放一条扫描线(x, y, z)
        QVector<QVector3D> normTemp;                                                            //临时存放一条扫描线的法向量
        QVector<QVector3D> directionTemp;                                                       //临时存放一条扫描线的声束方向
        int numOfScanPoints = 0;                                                                //临时存放一条扫描线点的个数，随着细分进行，逐渐递增

        /*对第一条扫描线进行等弧长细分*/
        qreal u = 0;
        const qreal v = 1.0 * i / numOfScanLine;                                                //固定v，沿u方向扫描
        while (u < 1) {
            qreal xtemp1 = surface.PointOnSurface(u, v, Qt::XAxis);
            qreal ytemp1 = surface.PointOnSurface(u, v, Qt::YAxis);
            qreal ztemp1 = surface.PointOnSurface(u, v, Qt::ZAxis);
            qreal dL = 0;
            while (dL <= scanInterval && u < 1) {                                                //固定v，沿u方向计算曲线弧长，每dL的弧长，记录一个点
                u += 0.0001;
                const qreal xtemp2 = surface.PointOnSurface(u, v, Qt::XAxis);
                const qreal ytemp2 = surface.PointOnSurface(u, v, Qt::YAxis);
                const qreal ztemp2 = surface.PointOnSurface(u, v, Qt::ZAxis);
                dL += qSqrt( qPow(xtemp2 - xtemp1, 2) +
                             qPow(ytemp2 - ytemp1, 2) +
                             qPow(ztemp2 - ztemp1, 2) );
                xtemp1 = xtemp2;
                ytemp1 = ytemp2;
                ztemp1 = ztemp2;
            }
            const QVector3D point(xtemp1, ytemp1, ztemp1);                                      //记录下扫描点
            lineTemp.append(point);                                                             //将扫描点添加到扫描点列表中
            normTemp.append(surface.NormVector(u, v));                                           //将法向量添加到扫描点法向量列表中
            if (numOfScanPoints >= 1) {
                QVector3D Ltemp = point - lineTemp[numOfScanPoints - 1];                        //切向量，上一个扫描点指向下一个扫描点的向量
                Ltemp.normalize();                                                              //将切向量归一化
                const QVector3D Mtemp = QVector3D::crossProduct(Ltemp, normTemp[numOfScanPoints - 1]);//切向量与法向量的公法向量，这三个向量构成随动坐标系
                const QVector3D direction = qSin(incidence) * (qSin(azimuth) * Mtemp + qCos(azimuth) * Ltemp)
                                 + qCos(incidence) * normTemp[numOfScanPoints - 1];             /*!< 采用入射角和方向角确定声束在随动坐标系中的方向*/
                directionTemp.append(direction);                                                //将声束方向添加到扫描线声束方向列表中
            }
            numOfScanPoints++;                                                                  //扫描线上的扫描点个数加1
        }
        lineTemp.removeLast();                                                                  //移除扫描线中最后一个扫描点，因为该扫描点无法利用差分计算切向量
        path.append(lineTemp);                                                                  //将生成的单条扫描线添加到路径中
        pathDirection.append(directionTemp);                                                    //将生成的单条扫描线的声束方向添加到路径声束方向中
        qDebug() << "第" << i << "条线点数为：" << numOfScanPoints - 1;
    }
}
/*!
 * \brief PathPlanning::GetPath             获得规划路径，即扫描点序列
 * \return
 */
QVector<QVector<QVector3D> > PathPlanning::GetPath()
{
    return path;
}
/*!
 * \brief PathPlanning::GetPathNormVector   获得规划路径对应的扫描点声束方向
 * \return
 */
QVector<QVector<QVector3D> > PathPlanning::GetPathDirection()
{
    return pathDirection;
}
/*!
 * \brief PathPlanning::feed                给定一条扫描点序列，预测下一条扫描点序列
 * \param pointsList                        一个点序列(u, v)，需要借助曲面将其转换成坐标值
 * \return                                  返回下一条扫描点序列(u, v)
 */
QVector<QVector<qreal> > PathPlanning::feed(QVector<QVector<qreal> > pointsList, qreal feed)
{
    QVector<QVector<qreal> > nextList(2, QVector<qreal>(pointsList[0].size(), 0.0));    //声明下一条扫描线, 与当前扫描线点数相同
    /*对当前扫描线上的点沿(u, v) = (0, 1)方向进给feedInterval */
    for (int i = 0; i < pointsList[0].size(); i++) {                                    //共pointsList[0].size()个点，分别进行沿v方向的进给
        qreal utemp = pointsList[0][i];                                                 //临时变量，用来对弧长进行积分
        qreal vtemp = pointsList[1][i];
        qreal dl = 0.0;
//      qreal du = 0.0;                                                                 //(du, dv)即进给的方向上的步长
        qreal dv = 0.0001;                                                              //即沿v方向进给

        qreal xtemp1 = surface.PointOnSurface(utemp, vtemp, Qt::XAxis);                 //临时变量，用来计算弧长
        qreal ytemp1 = surface.PointOnSurface(utemp, vtemp, Qt::YAxis);
        qreal ztemp1 = surface.PointOnSurface(utemp, vtemp, Qt::ZAxis);
        while (dl < feed && vtemp < 1) {
//          utemp += du;
            vtemp += dv;
            const qreal xtemp2 = surface.PointOnSurface(utemp, vtemp, Qt::XAxis);       //临时变量，用来计算弧长
            const qreal ytemp2 = surface.PointOnSurface(utemp, vtemp, Qt::YAxis);
            const qreal ztemp2 = surface.PointOnSurface(utemp, vtemp, Qt::ZAxis);
            dl += qSqrt( qPow(xtemp2 - xtemp1, 2) +                                     //计算弧长
                         qPow(ytemp2 - ytemp1, 2) +
                         qPow(ztemp2 - ztemp1, 2) );
            vtemp = vtemp > 1 ? 1 : vtemp;                                              //对vtemp进行控制，防止越界[0, 1]
            vtemp = vtemp < 0 ? 0 : vtemp;
            xtemp1 = xtemp2;
            ytemp1 = ytemp2;
            ztemp1 = ztemp2;
        }
        nextList[0][i] = utemp;                                                         //将进给后的点保存下来
        nextList[1][i] = vtemp;

    }
    nextList[0][0] = pointsList[0][0];                                                  //修正起点参数u，防止进给过程中出现收缩
    nextList[0][pointsList[0].size() - 1] = pointsList[0][pointsList[0].size() - 1];

    return nextList;
}
/*!
 * \brief PathPlanning::CurveSubdivide      给定一条点序列，对其进行等弧长细分，细分间隔为L
 * \param pointsList                        一个点序列(u, v)，需要借助曲面将其转换成坐标值
 * \param interval                          细分弧长间隔
 * \return                                  返回细分后的点序列(u, v)
 */
QVector<QVector<qreal> > PathPlanning::CurveSubdivide(QVector<QVector<qreal> > pointsList, qreal interval)
{
    /*对（u, v)做曲线拟合*/
    //BSplineCurveInterpolation类只支持3维，因此uvCurve被修正成3D曲线，本质还是二维曲线
    QVector<QVector3D> uvCurve(pointsList[0].size(), QVector3D(0.0, 0.0, 0.0));
    for (int i = 0; i < pointsList[0].size(); i++) {
        uvCurve[i] = QVector3D(pointsList[0][i], pointsList[1][i], 0.0);
    }
    int p = 2;
    BSplineCurve curve = BSplineCurve::FromInterpolation(p, uvCurve);    //将(u, v)拟合成一条p次曲线curve

    /* 对曲线进行细分, 细分间隔为dt */

    QVector<QVector<qreal> > uvOut(2, QVector<qreal>());    //声明返回值变量

    qreal t = 0.0;
    qreal dt = 0.0002;
    qreal utemp = curve.PointOnCurve(t, Qt::XAxis);
    qreal vtemp = curve.PointOnCurve(t, Qt::YAxis);
    qreal xtemp1 = surface.PointOnSurface(utemp, vtemp, Qt::XAxis);
    qreal ytemp1 = surface.PointOnSurface(utemp, vtemp, Qt::YAxis);
    qreal ztemp1 = surface.PointOnSurface(utemp, vtemp, Qt::ZAxis);
    /*计算弧长，进行细分*/
    while (t < 1) {
        qreal dl = 0.0;
        utemp = utemp > 1 ? 1 : utemp;                                              //对utemp, vtemp进行控制，防止越界[0, 1]
        utemp = utemp < 0 ? 0 : utemp;
        vtemp = vtemp > 1 ? 1 : vtemp;
        vtemp = vtemp < 0 ? 0 : vtemp;
        uvOut[0].append(utemp);                                                     //保存u值
        uvOut[1].append(vtemp);                                                     //保存v值
        while (dl < interval && t < 1) {                                            //计算每一段弧长interval
            t += dt;
            t = t > 1 ? 1 : t;
            utemp = curve.PointOnCurve(t, Qt::XAxis);
            vtemp = curve.PointOnCurve(t, Qt::YAxis);
            const qreal xtemp2 = surface.PointOnSurface(utemp, vtemp, Qt::XAxis);
            const qreal ytemp2 = surface.PointOnSurface(utemp, vtemp, Qt::YAxis);
            const qreal ztemp2 = surface.PointOnSurface(utemp, vtemp, Qt::ZAxis);
            dl += qSqrt( qPow(xtemp2 - xtemp1, 2) +                                 //计算弧长
                         qPow(ytemp2 - ytemp1, 2) +
                         qPow(ztemp2 - ztemp1, 2) );
            xtemp1 = xtemp2;
            ytemp1 = ytemp2;
            ztemp1 = ztemp2;
        }
    }
    /*最后一个点, 不能保证最后一个点与倒数第二个点间隔为interval*/
    t = 1.0;                                                                        //补上最后一个点t = 1
    utemp = curve.PointOnCurve(t, Qt::XAxis);
    vtemp = curve.PointOnCurve(t, Qt::YAxis);
    uvOut[0].append(utemp);                                                         //保存u值，最后一个点
    uvOut[1].append(vtemp);                                                         //保存v值

    return uvOut;
}

