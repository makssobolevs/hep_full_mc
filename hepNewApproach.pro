TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

VPATH += ./src/
SOURCES += main.cpp \
        definitions.cpp \
        matrixElement.cpp \
        msimplified.cpp \

HEADERS += definitions.h \
        matrixElement.h \
        msimplified.h \
        
OBJECTS_DIR = ../obj/

QMAKE_CXXFLAGS += -O3
