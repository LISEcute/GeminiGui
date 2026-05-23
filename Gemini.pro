TARGET = Gemini
TEMPLATE = app
CONFIG += c++17
QT       += widgets core gui printsupport sql

UI_DIR = ui
MOC_DIR = moc
RCC_DIR = rcc
OBJECTS_DIR = obj


win32-g++ {
DESTDIR = c:/Gemini/_install
}
win32-msvc {
DESTDIR = c:/Gemini/_install_MSVC
}

win32:VERSION = 3.1.3.1 # major.minor.patch.build
else:VERSION  = 3.1.3   # major.minor.patch

win32 {
        QMAKE_TARGET_COPYRIGHT = "LISE group at FRIB/MSU"
        QMAKE_TARGET_COMPANY   = "LISE group at FRIB/MSU"
        }

#==================================================
# The following define makes your compiler warn you if you use any
# feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Input


SOURCES += \
    g_Gemini/source/Gdr.cpp \
    g_Gemini/test/testDecay.cpp \
    g_Gemini/test/testWidth.cpp \
    g_Gemini/test/testFusion.cpp \
    g_Gemini/source/gemini.cpp \
    g_Gemini/source/Angle.cpp \
    g_Gemini/source/AngleDist.cpp \
    g_Gemini/source/Chart.cpp \
    g_Gemini/source/Evap.cpp \
    g_Gemini/source/Fus.cpp \
    g_Gemini/source/History.cpp \
    g_Gemini/source/LevelDensity.cpp \
    g_Gemini/source/LightP.cpp \
    g_Gemini/source/Mass.cpp \
    g_Gemini/source/Nucleus.cpp \
    g_Gemini/source/Nuclide.cpp \
    g_Gemini/source/Random.cpp \
    g_Gemini/source/Scission.cpp \
    g_Gemini/source/SigBarDist.cpp \
    g_Gemini/source/SigCharged.cpp \
    g_Gemini/source/TlArray.cpp \
    g_Gemini/source/TlBarDist.cpp \
    g_Gemini/source/Weight.cpp \
    g_Gemini/source/Yrast.cpp \
    g_Gemini/windows/gm_about.cpp \
    g_Gemini/windows/gm_angular_distribution.cpp \
    g_Gemini/windows/gm_cs_file.cpp \
    g_Gemini/windows/gm_fileReadWrite.cpp \
    g_Gemini/windows/gm_main.cpp \
    g_Gemini/windows/gm_mainwindow.cpp \
    g_Gemini/windows/gm_mwExecCompound.cpp \
    g_Gemini/windows/gm_mwExecFusion.cpp \
    g_Gemini/windows/gm_results.cpp \
    w_Stuff/w_Label_clickable.cpp

HEADERS  += \
    g_Gemini/source/CRun.h \
    g_Gemini/source/CAngle.h \
    g_Gemini/source/CAngleDist.h \
    g_Gemini/source/CChart.h \
    g_Gemini/source/CEvap.h \
    g_Gemini/source/CFus.h \
    g_Gemini/source/CGdr.h \
    g_Gemini/source/CHistory.h \
    g_Gemini/source/CLevelDensity.h \
    g_Gemini/source/CLightP.h \
    g_Gemini/source/CMass.h \
    g_Gemini/source/CNucleus.h \
    g_Gemini/source/CNuclide.h \
    g_Gemini/source/CRandom.h \
    g_Gemini/source/CRunThick.h \
    g_Gemini/source/CScission.h \
    g_Gemini/source/CSigBarDist.h \
    g_Gemini/source/CSigCharged.h \
    g_Gemini/source/CTlArray.h \
    g_Gemini/source/CTlBarDist.h \
    g_Gemini/source/CWeight.h \
    g_Gemini/source/CYrast.h \
    g_Gemini/source/SStoreEvap.h \
    g_Gemini/windows/gm_about.h \
    g_Gemini/windows/gm_angular_distribution.h \
    g_Gemini/windows/gm_ftype.h \
    g_Gemini/windows/gm_mainwindow.h \
    g_Gemini/windows/gm_results.h \
    w_Stuff/w_Label_clickable.h

FORMS    += \
    g_Gemini/windows/gm_about.ui \
    g_Gemini/windows/gm_mainwindow.ui

OTHER_FILES += \
    g_Gemini/tbl/evap.inp \
    g_Gemini/tbl/gemini.inp \
    g_Gemini/tbl/GDR.tbl \
    g_Gemini/tbl/sad.tbl \
    g_Gemini/tl/alpha.inv \
    g_Gemini/tl/alpha.tl \
    g_Gemini/tl/alphaM.inv \
    g_Gemini/tl/alphaM.tl \
    g_Gemini/tl/alphaP.inv \
    g_Gemini/tl/alphaP.tl \
    g_Gemini/tl/be7.tl \
    g_Gemini/tl/be7M.tl \
    g_Gemini/tl/be7P.tl \
    g_Gemini/tl/be8.tl \
    g_Gemini/tl/be8M.tl \
    g_Gemini/tl/be8P.tl \
    g_Gemini/tl/be9.tl \
    g_Gemini/tl/be9M.tl \
    g_Gemini/tl/be9P.tl \
    g_Gemini/tl/be10.tl \
    g_Gemini/tl/be10M.tl \
    g_Gemini/tl/be10P.tl \
    g_Gemini/tl/deuteron.inv \
    g_Gemini/tl/deuteron.tl \
    g_Gemini/tl/deuteronM.inv \
    g_Gemini/tl/deuteronM.tl \
    g_Gemini/tl/deuteronP.inv \
    g_Gemini/tl/deuteronP.tl \
    g_Gemini/tl/he3.inv \
    g_Gemini/tl/he3.tl \
    g_Gemini/tl/he3M.inv \
    g_Gemini/tl/he3M.tl \
    g_Gemini/tl/he3P.inv \
    g_Gemini/tl/he3P.tl \
    g_Gemini/tl/he6.inv \
    g_Gemini/tl/he6.tl \
    g_Gemini/tl/he6M.inv \
    g_Gemini/tl/he6M.tl \
    g_Gemini/tl/he6P.inv \
    g_Gemini/tl/he6P.tl \
    g_Gemini/tl/li6.inv \
    g_Gemini/tl/li6.tl \
    g_Gemini/tl/li6M.inv \
    g_Gemini/tl/li6M.tl \
    g_Gemini/tl/li6P.inv \
    g_Gemini/tl/li6P.tl \
    g_Gemini/tl/li7.inv \
    g_Gemini/tl/li7.tl \
    g_Gemini/tl/li7M.inv \
    g_Gemini/tl/li7M.tl \
    g_Gemini/tl/li7P.inv \
    g_Gemini/tl/li7P.tl \
    g_Gemini/tl/li8.inv \
    g_Gemini/tl/li8.tl \
    g_Gemini/tl/li8M.inv \
    g_Gemini/tl/li8M.tl \
    g_Gemini/tl/li8P.inv \
    g_Gemini/tl/li8P.tl \
    g_Gemini/tl/neutron.tl \
    g_Gemini/tl/neutronM.tl \
    g_Gemini/tl/neutronP.tl \
    g_Gemini/tl/proton.inv \
    g_Gemini/tl/proton.tl \
    g_Gemini/tl/protonM.inv \
    g_Gemini/tl/protonM.tl \
    g_Gemini/tl/protonP.inv \
    g_Gemini/tl/protonP.tl \
    g_Gemini/tl/triton.inv \
    g_Gemini/tl/triton.tl \
    g_Gemini/tl/tritonM.inv \
    g_Gemini/tl/tritonM.tl \
    g_Gemini/tl/tritonP.inv \
    g_Gemini/tl/tritonP.tl

RESOURCES += \
    g_Gemini/gm_resources.qrc \
    g_Gemini/icons/gm_icons.qrc \
    g_Gemini/icons/icons.qrc \
    lise.qrc


#DISTFILES +=   g_Gemini/icons/Gemini_logo.ico


RC_ICONS += \
       g_Gemini/icons/Gemini_logo.ico \
	g_Gemini/icons/Gemini_logo16.ico

DEFINES += __APPLE_
ICON = ../Icons_macos/gemini.icns
QMAKE_INFO_PLIST = ../Info.plist
QT += core concurrent

