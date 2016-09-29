#include "BSplineCurve.h"
#include "MyMath.h"
#include <QDebug>

BSplineCurve::BSplineCurve(int _order, QVector<qreal> knotVector, QVector<QVector3D> controlPoints) :
    bspline(knotVector, _order)
{
    ctrlPoints = controlPoints;
}

BSplineCurve::BSplineCurve(BSpline bs, QVector<QVector3D> controlPoints)
{
    bspline = bs;
    ctrlPoints = controlPoints;
}

qreal BSplineCurve::PointOnCurve(qreal u, Qt::Axis axis)
{
    u = qBound<qreal>(0, u, 1);
    int span = bspline.FindSpan(u);
    QVector<qreal> basis = bspline.BasisFunc(u, span);
    qreal CP = 0;
    for (int i = 0; i <= bspline.order; i++) {
        CP += basis[i] * ctrlPoints[span - bspline.order + i][axis];
    }
    return CP;
}

/*!
 * \brief CurveDerivs::BasisFunc            计算指定参数的axis轴坐标的多阶导数
 * \param u                                 指定参数
 * \param diffOrder                         导数阶次
 * \param axis                              轴坐标标志
 * \param p                                 基函数阶次
 * \return                                  指定参数的axis轴坐标的多阶导数，深度为diffOrder + 1
 */
QVector<qreal> BSplineCurve::CurveDerivs(qreal u, int diffOrder, int p, Qt::Axis axis)
{
    u = qBound<qreal>(0, u, 1);
    QVector<qreal> der(diffOrder + 1, 0);       //长度为：dirrOrder + 1
    int du = p > diffOrder ? diffOrder : p;     // du = min(p, diffOrder)
    for (int k = p + 1; k <= diffOrder; k++) {
        //将高于p将的导数清零
        der[k] = 0;
    }
    int span = bspline.FindSpan(u);
    QVector<QVector<qreal> > nders = bspline.DerivativeBasis(u, span, diffOrder);
    for (int k = 0; k <= du; k++) {
        der[k] = 0;
        for (int j = 0; j <= p; j++) {
            der[k] = der[k] + nders[k][j] * ctrlPoints[span - p + j][axis];
        }
    }
    return der;
}

/*!
 * \brief BSplineCurve::FromInterpolation 不指定端点导矢的曲线拟合
 * \param _order
 * \param _interpolationPoints
 * \return
 */
BSplineCurve BSplineCurve::FromInterpolation(int order, const QVector<QVector3D> &interpolationPoints)
{
    BSplineCurve curve;
    BSpline &bspline = curve.bspline;
    bspline.order = order;
    int lengthOfInterPoints = interpolationPoints.size() - 1;                      //插值点的个数
    int lengthOfCtrlPoints = lengthOfInterPoints;                                   //不指定端点导矢的情况下，控制点与插值点个数相同
    int dimension = 2;//interpolationPoints[0].size() - 1;                       //插值点的维度, 在QVector3D下dimension = 2
/*在这里就完成曲面的拟合，计算出插值点对应的节点、节点矢量、控制点*/

    /*下面计算插值点对应的参数值knotsInterpolation*/
    QVector<qreal> knotsInterpolation(lengthOfInterPoints + 1, 0);              //定义插值点对应的参数值
    QVector<qreal> temp(lengthOfInterPoints, 0);                              //保存各相邻插值点间的距离
    qreal d = 0;                                                                //相邻插值点间的距离之和
    for (int i = 0; i < lengthOfInterPoints; i++) {
        temp[i] = 0;
        for (int j = 0; j <= dimension; j++) {
            temp[i] += qPow( interpolationPoints[i][j] - interpolationPoints[i + 1][j], 2);
        }
        temp[i] = qSqrt(temp[i]);                                            //计算相邻两点的距离
        d += temp[i];                                                       //总距离
    }
    /*利用9.5式计算插值点对应的参数knotsInterpolation      */
    for (int i = 1; i <= lengthOfInterPoints; i++) {
        knotsInterpolation[i] = knotsInterpolation[i - 1] + temp[i - 1] / d;    //knotsInterpolation[0] = 0
    }
    /*配置节点矢量*/
    int lengthOfKnots = lengthOfCtrlPoints + bspline.order + 1;                    //节点矢量的长度为lengthOfKnots+1
    bspline.knots.resize(lengthOfKnots + 1);
    for (int i = 0; i <= lengthOfKnots; i++ ) {
        if (i <= bspline.order) {
            bspline.knots[i] = 0;
        } else if (i >= lengthOfKnots - bspline.order) {
            bspline.knots[i] = 1;
        } else {
            qreal sum_temp = 0;
            for (int k = 1; k <= bspline.order; k++) {
                sum_temp += knotsInterpolation[i - k];
            }
            bspline.knots[i] = 1.0 / bspline.order * sum_temp;
        }
    }
    /*建立系数矩阵A */
    QVector<QVector<qreal> > A(lengthOfCtrlPoints + 1, QVector<qreal>(lengthOfCtrlPoints + 1, 0) );      //系数矩阵为(lengthOfCtrlPoints+1)*(lengthOfCtrlPoints+1)
    for (int i = 0; i <= lengthOfInterPoints; i++) {
        int span = bspline.FindSpan(knotsInterpolation[i]);
        QVector<qreal> basis = bspline.BasisFunc(knotsInterpolation[i], span);
        for (int j = 0; j <= bspline.order; j++) {
            A[i][span-bspline.order+j] = basis[j];
        }
    }
    /*利用高斯消元法求解控制点: 对各个Axis分别进行求解  */
    curve.ctrlPoints = QVector<QVector3D >(lengthOfCtrlPoints + 1, QVector3D(0, 0, 0));
    for (int i = 0; i <= dimension; i++) {
        QVector<qreal> b(lengthOfCtrlPoints+1, 0);
        for (int j = 0; j <= lengthOfCtrlPoints; j++) {  //提取某个axis的坐标
            b[j] = interpolationPoints[j][i];
        }
        QVector<qreal> x = Gauss(A, b);                  //高斯求解
        for (int j = 0; j <= lengthOfCtrlPoints; j++) {  //复合某个axis的坐标
            curve.ctrlPoints[j][i] = x[j];
        }
    }
    return curve;
}
/*!
 * \brief BSplineCurve::FromInterpolation       指定端点一阶导矢的曲线拟合
 * \param _order
 * \param _interpolationPoints
 * \param D0
 * \param Dn
 * \return
 */
BSplineCurve BSplineCurve::FromInterpolation(int order, QVector<QVector3D> interpolationPoints, QVector3D D0, QVector3D Dn)
{
    BSplineCurve curve;
    BSpline &bspline = curve.bspline;
    order = qMax(order, 2);                                                     //指定端点一阶导矢的曲线，阶次不能少于2阶，否则不存在这样的曲线
    bspline.order = order;
    int lengthOfInterPoints = interpolationPoints.size() - 1;                   //插值点的个数
    int lengthOfCtrlPoints = lengthOfInterPoints + 2;                           //指定端点一阶导矢的情况下，控制点比插值点个数多2个
    int dimension = 2;//interpolationPoints[0].size() - 1;                      //插值点的维度, 在QVector3D下dimension = 2

    /*在这里就完成曲面的拟合，计算出插值点对应的节点、节点矢量、控制点*/

    /*下面计算插值点对应的参数值knotsInterpolation*/
    QVector<qreal> knotsInterpolation(lengthOfInterPoints + 1, 0);              //定义插值点对应的参数值的长度
    QVector<qreal> temp(lengthOfInterPoints, 0);                                 //保存各相邻插值点间的距离
    qreal d = 0;                                                                //相邻插值点间的距离之和
    for (int i = 0; i < lengthOfInterPoints; i++) {
        temp[i] = 0;
        for (int j = 0; j <= dimension; j++) {
            temp[i] += qPow( interpolationPoints[i][j] - interpolationPoints[i+1][j], 2);
        }
        temp[i] = qSqrt(temp[i]);                                               //计算相邻两点的距离
        d += temp[i];                                                           //总距离
    }
    /*利用9.5式计算插值点对应的参数knotsInterpolation      */
    for (int i = 1; i <= lengthOfInterPoints; i++ ) {
        knotsInterpolation[i] = knotsInterpolation[i-1] + temp[i-1] / d;        //knotsInterpolation[0] = 0
    }
    /*配置节点矢量 */
    int lengthOfKnots = lengthOfCtrlPoints + bspline.order + 1;                //节点矢量的长度为lengthOfKnots+1
    bspline.knots.resize(lengthOfKnots + 1);
    for (int i = 0; i <= lengthOfKnots; i++ ) {
        if (i <= order) {
             bspline.knots[i] = 0;
        }
        else if (i >= lengthOfKnots - order) {
             bspline.knots[i] = 1;
        }
        else{
            qreal sum_temp = 0;
            for (int k = 1; k <= order; k++ ) {
                sum_temp += knotsInterpolation[i - k - 1];
            }
             bspline.knots[i] = 1.0 / order * sum_temp;
        }
    }
    /*建立系数矩阵A */
    QVector<QVector<qreal> > A( lengthOfCtrlPoints+1, QVector<qreal>(lengthOfCtrlPoints+1, 0));      //系数矩阵为(lengthOfCtrlPoints+1)*(lengthOfCtrlPoints+1)
    //定系数矩阵第一行
    int span = bspline.FindSpan(knotsInterpolation[0]);
    QVector<qreal> basis = bspline.BasisFunc(knotsInterpolation[0], span);
    for (int i = 0; i <= order; i++) {
        A[0][span-order+i] = basis[i];
    }
    //定系数矩阵第二行
    A[1][0] = order / bspline.knots[order + 1] *(-1);
    A[1][1] = order / bspline.knots[order + 1];
    //定系数矩阵中间部分
    for (int i = 2; i <= lengthOfCtrlPoints-2; i++) {
        span = bspline.FindSpan(knotsInterpolation[i - 1]);
        basis = bspline.BasisFunc(knotsInterpolation[i - 1], span);
        for (int j = 0; j <= order; j++) {
            A[i][span-order+j] = basis[j];
        }
    }
    //定系数矩阵倒数第二行
    A[lengthOfCtrlPoints-1][lengthOfCtrlPoints-1] = (-1.0) * order / (1 - bspline.knots[lengthOfKnots - order - 1]);
    A[lengthOfCtrlPoints-1][lengthOfCtrlPoints] = 1.0 * order / (1 - bspline.knots[lengthOfKnots - order - 1]);
    //定系数矩阵最后一行
    span = bspline.FindSpan(knotsInterpolation[lengthOfInterPoints]);
    basis = bspline.BasisFunc(knotsInterpolation[lengthOfInterPoints], span);
    for (int i = 0; i <= order; i++) {
        A[lengthOfCtrlPoints][span - order + i] = basis[i];
    }

    //重新构造插值点
    QVector<QVector3D> interPointsTemp(lengthOfCtrlPoints+1, QVector3D(0, 0, 0));
    interPointsTemp[0] = interpolationPoints[0];
    interPointsTemp[1] = D0;
    for (int i = 1; i <= lengthOfInterPoints - 1; i++) {
        interPointsTemp[i+1] = interpolationPoints[i];
    }
    interPointsTemp[lengthOfCtrlPoints-1] = Dn;
    interPointsTemp[lengthOfCtrlPoints] = interpolationPoints[lengthOfInterPoints];

    /*利用高斯消元法求解控制点: 对各个Axis分别进行求解  */
    curve.ctrlPoints = QVector<QVector3D >(lengthOfCtrlPoints + 1, QVector3D(0, 0, 0));
    for (int i = 0; i <= dimension; i++) {
        QVector<qreal> b(lengthOfCtrlPoints+1, 0);
        for (int j = 0; j <= lengthOfCtrlPoints; j++) {  //提取某个axis的坐标
            b[j] = interPointsTemp[j][i];
        }
        QVector<qreal> x = Gauss(A, b);                  //高斯求解
        for (int j = 0; j <= lengthOfCtrlPoints; j++) {  //复合某个axis的坐标
            curve.ctrlPoints[j][i] = x[j];
        }
    }
    return curve;
}

BSplineCurve BSplineCurve::FromInterpolation(int order, QVector<QVector3D> interpolationPoints, QVector3D D0, QVector3D Dn, QVector3D A0, QVector3D An)
{
    BSplineCurve curve;
    BSpline &bspline = curve.bspline;
    order = qMax(order, 3);
    bspline.order = order;                                                      //指定端点二阶导矢的曲线，拟合阶次不能低于3次，否则不存在这样的曲线
    int lengthOfInterPoints = interpolationPoints.size() - 1;                   //插值点的个数
    int lengthOfCtrlPoints = lengthOfInterPoints + 4;                           //指定端点二阶导矢的情况下，控制点比插值点个数多4个
    int dimension = 2;//interpolationPoints[0].size() - 1;                      //插值点的维度

    /*在这里就完成曲面的拟合，计算出插值点对应的节点、节点矢量、控制点*/

    /*下面计算插值点对应的参数值knotsInterpolation
      利用9.5式计算插值点对应的参数knotsInterpolation*/
    QVector<qreal> knotsInterpolation(lengthOfInterPoints + 1, 0);              //定义插值点对应的参数值的长度
    QVector<qreal> temp(lengthOfInterPoints, 0);                                 //保存各相邻插值点间的距离
    qreal d = 0;                                                                //相邻插值点间的距离之和
    for (int i = 0; i < lengthOfInterPoints; i++) {
        temp[i] = 0;
        for (int j = 0; j <= dimension; j++) {
            temp[i] += qPow( interpolationPoints[i][j] - interpolationPoints[i+1][j], 2);
        }
        temp[i] = qSqrt(temp[i]);                                               //计算相邻两点的距离
        d += temp[i];                                                           //总距离
    }
    for (int i = 1; i <= lengthOfInterPoints; i++ ) {
        knotsInterpolation[i] = knotsInterpolation[i-1] + temp[i-1] / d;        //knotsInterpolation[0] = 0
    }

    /*构造一个临时插值点对应的节点，方法是将原节点向两头拓展一个，即左侧加0，右侧加1*/
    QVector<qreal> knotInterTemp(lengthOfInterPoints + 3, 0);
    knotInterTemp[0] = 0;
    for (int i = 1; i <= lengthOfInterPoints + 1; i++) {
        knotInterTemp[i] = knotsInterpolation[i - 1];
    }
    knotInterTemp[lengthOfInterPoints + 2] = 1.0;

    /*配置节点矢量 */
    int lengthOfKnots = lengthOfCtrlPoints + order + 1;                         //节点矢量的长度为lengthOfKnots+1
    bspline.knots.resize(lengthOfKnots + 1);
    for (int i = 0; i <= lengthOfKnots; i++ ) {
        if (i <= order) {
            bspline.knots[i] = 0;
        }
        else if (i >= lengthOfKnots-order) {
            bspline.knots[i] = 1;
        }
        else{
            qreal sum_temp = 0;
            for (int k = 1; k <= order; k++ ) {
                sum_temp += knotInterTemp[i - k - 1];
            }
            bspline.knots[i] = 1.0 / order * sum_temp;
        }
    }
    /*建立系数矩阵A */
    QVector<QVector<qreal> > A( lengthOfCtrlPoints+1, QVector<qreal>(lengthOfCtrlPoints+1, 0) );      //系数矩阵为(lengthOfCtrlPoints+1)*(lengthOfCtrlPoints+1)
    //定系数矩阵第一行
    int span = bspline.FindSpan(knotsInterpolation[0]);
    QVector<qreal> basis = bspline.BasisFunc(knotsInterpolation[0], span);
    for (int i = 0; i <= order; i++) {
        A[0][span-order+i] = basis[i];
    }
    //定系数矩阵第二行
    A[1][0] = order / bspline.knots[order + 1] *(-1);
    A[1][1] = order / bspline.knots[order + 1];
    //定系数矩阵第三行
    A[2][0] = order * (order - 1) / bspline.knots[order + 1] / bspline.knots[order + 1];
    A[2][1] = order * (order - 1) / bspline.knots[order + 1] * ( -1 * bspline.knots[order + 1] - bspline.knots[order + 2] ) / ( bspline.knots[order + 1] * bspline.knots[order + 2]);
    A[2][2] = order * (order - 1) / bspline.knots[order + 1] / bspline.knots[order + 2];
    //定系数矩阵中间部分
    for (int i = 3; i <= lengthOfCtrlPoints - 3; i++) {
        span = bspline.FindSpan(knotsInterpolation[i - 2]);
        basis = bspline.BasisFunc(knotsInterpolation[i - 2], span);
        for (int j = 0; j <= order; j++) {
            A[i][span - (order - j)] = basis[j];
        }
    }
    //定系数矩阵倒数第三行
    A[lengthOfCtrlPoints-2][lengthOfCtrlPoints - 2] =
            order * (order - 1)
            / (1 - bspline.knots[lengthOfKnots - order - 1])
            / (1 - bspline.knots[lengthOfKnots - order - 2]);

    A[lengthOfCtrlPoints-2][lengthOfCtrlPoints-1] =
            order
            * (order - 1)
            / (1 - bspline.knots[lengthOfKnots - order - 1])
            * (bspline.knots[lengthOfKnots - order -1] + bspline.knots[lengthOfKnots - order - 2] -2)
            / (1 - bspline.knots[lengthOfKnots - order - 1])
            / (1 - bspline.knots[lengthOfKnots - order - 2]);

    A[lengthOfCtrlPoints-2][lengthOfCtrlPoints] =
            order
            * (order - 1)
            / (1 - bspline.knots[lengthOfKnots - order - 1])
            / (1 - bspline.knots[lengthOfKnots - order -1]);

    //定系数矩阵倒数第二行
    A[lengthOfCtrlPoints - 1][lengthOfCtrlPoints - 1] =
            (-1.0) * order
            / (1 - bspline.knots[lengthOfKnots - order - 1]);
    A[lengthOfCtrlPoints-1][lengthOfCtrlPoints] =
            1.0 * order
            / (1 - bspline.knots[lengthOfKnots - order - 1]);
    //定系数矩阵最后一行
    span = bspline.FindSpan(knotsInterpolation[lengthOfInterPoints]);
    basis = bspline.BasisFunc(knotsInterpolation[lengthOfInterPoints], span);
    for (int i = 0; i <= order; i++) {
        A[lengthOfCtrlPoints][span - order + i] = basis[i];
    }

    //重新构造插值点
    QVector<QVector3D > interPointsTemp(lengthOfCtrlPoints+1, QVector3D(0, 0, 0));
    interPointsTemp[0] = interpolationPoints[0];
    interPointsTemp[1] = D0;
    interPointsTemp[2] = A0;
    for (int i = 3; i <= lengthOfCtrlPoints - 3; i++) {
        interPointsTemp[i] = interpolationPoints[i-2];
    }
    interPointsTemp[lengthOfCtrlPoints-2] = An;
    interPointsTemp[lengthOfCtrlPoints-1] = Dn;
    interPointsTemp[lengthOfCtrlPoints] = interpolationPoints[lengthOfInterPoints];

    /*利用高斯消元法求解控制点: 对各个Axis分别进行求解  */
    curve.ctrlPoints = QVector<QVector3D >(lengthOfCtrlPoints+1, QVector3D(0, 0, 0));
    for (int i = 0; i <= dimension; i++) {
        QVector<qreal> b(lengthOfCtrlPoints+1, 0);
        for (int j = 0; j <= lengthOfCtrlPoints; j++) {         //提取某个axis的坐标
            b[j] = interPointsTemp[j][i];
        }
        QVector<qreal> x = Gauss(A, b);                          //高斯求解
        for (int j = 0; j <= lengthOfCtrlPoints; j++) {         //复合某个axis的坐标
            curve.ctrlPoints[j][i] = x[j];
        }
    }
    return curve;
}
/*! 对点序列进行重新拟合并细分成numOfSubdivide个新点
*/
QVector<QVector3D> BSplineCurve::CurveSubdivide(const QVector<QVector3D> &points, const int numOfSubdivide)
{
    QVector<QVector3D> pointsOut;// = QVector(numOfSubdivide);
    pointsOut.resize(numOfSubdivide);
    BSplineCurve curve = BSplineCurve::FromInterpolation(3, points);
    //求曲线总长度
    qreal u = 0;
    qreal lengthOfCurve = 0;
    QVector<qreal> tempX = curve.CurveDerivs(u, 1, 3, Qt::XAxis);
    QVector<qreal> tempY = curve.CurveDerivs(u, 1, 3, Qt::YAxis);
    QVector<qreal> tempZ = curve.CurveDerivs(u, 1, 3, Qt::ZAxis);
    while (u < 1) {
        u += 0.0002;
        tempX = curve.CurveDerivs(u, 1, 3, Qt::XAxis);
        tempY = curve.CurveDerivs(u, 1, 3, Qt::YAxis);
        tempZ = curve.CurveDerivs(u, 1, 3, Qt::ZAxis);

        lengthOfCurve +=  0.0002 * qSqrt(qPow(tempX[1], 2) + qPow(tempY[1], 2) + qPow(tempZ[1], 2));
    }
    qreal interval = lengthOfCurve / (numOfSubdivide - 1); /*点间隔 */
    /* 计算各细分点 */
    u = 0;
    int numTemp = 1;
    qreal dL = 0;
    pointsOut[0] = points[0];       /*指定第一个点 */
    u = 0;
    while (u < 1 && numTemp < numOfSubdivide) {
        while (u < 1 && dL < numTemp * interval) {
            u += 0.0002;
            tempX = curve.CurveDerivs(u, 1, 3, Qt::XAxis);
            tempY = curve.CurveDerivs(u, 1, 3, Qt::YAxis);
            tempZ = curve.CurveDerivs(u, 1, 3, Qt::ZAxis);
            dL += 0.0002 * qSqrt(qPow(tempX[1], 2) + qPow(tempY[1], 2) + qPow(tempZ[1], 2));
        }
        pointsOut[numTemp] = QVector3D(tempX[0], tempY[0], tempZ[0]);
        numTemp++;
    }
    pointsOut[numTemp - 1] = points.last(); /*指定最后一个点 */
    return pointsOut;
}
