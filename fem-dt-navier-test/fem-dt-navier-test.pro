TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp

LIBS += -ldolfin

OTHER_FILES += \
    navier.ufl \
    navier_s5.ufl \
    navier_s4.ufl \
    navier_s3.ufl \
    navier_s2.ufl \
    navier_s1.ufl \
    navier_n3.ufl \
    navier_n2.ufl \
    navier_n1.ufl \
    navier_d3.ufl \
    navier_d2.ufl \
    navier_d1.ufl \
    navier_velocity.ufl \
    navier_pressure.ufl \
    navier_psi.ufl
