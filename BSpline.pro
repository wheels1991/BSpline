QT += core
QT += gui widgets

#QT -= gui

CONFIG += C++11
TARGET = BSpline
#CONFIG += console
#CONFIG -= app_bundle

TEMPLATE = app

INCLUDEPATH += ../../../Source

SOURCES += main.cpp \
    BSplineCurve.cpp \
    BSplineSurface.cpp \
    MyMath.cpp \
    BSpline.cpp \
    PathPlanning.cpp \

HEADERS += \
    BSplineCurve.h \
    BSplineSurface.h \
    MyMath.h \
    BSpline.h \
    PathPlanning.h \

