#include "MyMath.h"

/*!
 * \brief Gauss             利用高斯消元法求解线性方程组Ax=b（不是列主元法，因此对于对角线元素为零的情况，可能无法求解）
 * \param A                 系数矩阵
 * \param b
 * \return
 */
QVector<qreal> Gauss(QVector<QVector<qreal> > A, QVector<qreal> b)
{
    /*判断A与b的维度是否一致*/
    int row0 = A.size() - 1;
    int col0 = A[0].size() - 1;
    int n = b.size() - 1;
    if (row0 != col0 || row0 != n) {
        return QVector<qreal>();            //返回一个空的值，在外面判断是否求解成功
    }
    QVector<qreal> x(n+1, 0.0);
    /*将系数矩阵上三角化*/
    for (int row1 = 0; row1 <= n-1; row1++) {              //从第二行开始消元
        for (int row2 = row1+1; row2 <= n; row2++) {
            qreal m = A[row2][row1]/A[row1][row1];
            for (int col = row1; col <= n; col++) {
                A[row2][col] = A[row2][col] - m * A[row1][col];
            }
            b[row2] = b[row2] - m * b[row1];
        }
    }
    x[n] = b[n] / A[n][n];
    for (int i = n - 1; i >= 0; i--) {
        if (qFuzzyIsNull(A[i][i])) {             //对角化后主元为0的情况，无法求解
            return QVector<qreal>();
        }
        qreal sum = 0.0;
        for (int j = i + 1; j <= n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
    return x;

}

