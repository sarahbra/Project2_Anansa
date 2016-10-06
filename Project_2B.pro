QT += core
QT -= gui

CONFIG += c++11

TARGET = Project_2B
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    jakobi_method.cpp \
    A.cpp \
    2b.cpp

HEADERS += \
    A.h
