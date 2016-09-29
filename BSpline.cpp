#include "BSpline.h"
#include <QtMath>

BSpline::BSpline()
{
    SetOrder(2);
}

BSpline::BSpline(const QVector<qreal> &knots, int order) :
    knots(knots)
{
    qSort(this->knots);
    SetOrder(order);
}

void BSpline::SetOrder(int newOrder)
{
    order = qMax(1, newOrder);
    if (knots.size() < 2 * order) {
        knots.resize(2 * order);
    }
}

QVector<qreal> BSpline::BasisFunc(qreal u, int span)
{
    u = qBound<qreal>(0, u, 1);
    QVector<qreal> basis(order + 1, 0);        //长度为：p+1
    QVector<qreal> left(order + 1, 0);
    QVector<qreal> right(order + 1, 0);
    basis[0] = 1;
    for (int j = 1; j <= order; j++) {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        qreal saved = 0;
        for (int r = 0; r < j; r++) {
            qreal temp = basis[r] / (right[r + 1] + left[j - r]);
            basis[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        basis[j] = saved;
    }
    return basis;
}

/*! 基函数的0-n阶导数 */
QVector<QVector<qreal> > BSpline::DerivativeBasis(qreal u, int span, int diffOrder)
{
    u = qBound<qreal>(0, u, 1);
    QVector<QVector<qreal> > ders(diffOrder + 1, QVector<qreal>(order + 1, 0));      //定义一个[dirrOrder + 1]*[p + 1] 的数组
    QVector<QVector<qreal> > ndu(order + 1, QVector<qreal>(order + 1, 0));       //临时变量，用来存储基函数和节点之差
    QVector<qreal> left(order + 1, 0);
    QVector<qreal> right(order + 1, 0);
    ndu[0][0] = 1;
    for (int j = 1; j <= order; ++j) {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        qreal saved = 0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            qreal temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    for (int j = 0; j <= order; ++j) {                                          //载入基函数的值
        ders[0][j] = ndu[j][order];
    }
    /* 计算导数 */
    QVector<QVector<qreal> > a(2, QVector<qreal>(diffOrder + 1, 0));
    for (int r = 0; r <= order; r++) {
        int s1 = 0;
        int s2 = 1;
        a[0][0] = 1;
        /*循环计算k阶导数, k = 1, 2, ..., p */
        for (int k = 1; k <= diffOrder; ++k) {
            qreal d = 0;
            int rk = r - k;
            int pk = order - k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r - 1 <= pk) ? k -1 : order - r;
            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            qSwap(s1, s2);
        }

    }
    /* 对结果乘以正确的因子 */
    int r = order;
    for (int k = 1; k <= diffOrder; ++k) {
        for (int j = 0; j <= order; ++j) {
            ders[k][j] *= r;
        }
        r *= (order - k);
    }
    return ders;
}

///*! 计算参数u在节点矢量中的位置 */
//int BSpline::FindSpan(qreal u)
//{
//    const int maxIndex = knots.size() - 1 - order;                               /*!< 控制点的个数 */
//    Q_ASSERT_X(maxIndex >= 0 && maxIndex < knots.size() - 1, "BSpline::FindSpan", "Index of control points out of bound");

//    u = qBound<qreal>(0, u, 1);
//    if (qFuzzyCompare(u, knots[maxIndex])) {
//        return maxIndex - 1;
//    }
//    int low = order;
//    int high = maxIndex;
//    if (low >= high) {
//        return high;
//    }

//    int mid;
//    do {
//        mid = (low + high) / 2;
//        if (u < knots[mid]) {
//            high = mid;
//        } else {
//            low = mid;
//        }
//    } while ((high - low) > 1 && !(u >= knots[mid] && u < knots[mid + 1]));
//    return mid;
//}
int BSpline::FindSpan(qreal u)
{
    int n = knots.size() - order - 2;
    u = u > 1 ? 1 : u;                      //保证u介于[0, 1]
    u = u < 0 ? 0 : u;
    if( u >= knots[ n + 1 ]){
        return n;
    }
    int low = order;
    int high = n + 1;
    int mid = ( low + high ) / 2;
    while( u < knots[mid] || u >= knots[mid + 1] )
    {
        if( u < knots[mid]){
            high = mid;
        }
        else{
            low = mid;
        }
        mid = ( low + high ) / 2;
    }
    return mid;
}

//BSpline BSpline::FromInterpolation(int order, const QVector<QVector3D> &points)
//{
//    BSpline spline;
//    /*配置节点矢量knotsInU*/
//    int lengthOfKnotsInU = lengthOfInterPointsInU + spline.order + 1;                       //节点矢量的长度为lengthOfKnotsInU+1
//    QVector<qreal> knotsInU(lengthOfKnotsInU + 1, 0);
//    for (int i = 0; i <= lengthOfKnotsInU; i++) {
//        if (i <= spline.order) {                                                             //前orderInU+1个节点是0
//            knotsInU[i] = 0;
//        } else if (i >= lengthOfKnotsInU - spline.order) {                                  //后orderInU+1个节点是1
//            knotsInU[i] = 1.0;
//        } else {
//            qreal sum_temp = 0;
//            for (int k = 1; k <= spline.order; k++) {
//                sum_temp += knotsInterInU[i - k];
//            }
//            knotsInU[i] = 1.0 / spline.order * sum_temp;
//        }
//    }
//}
