#ifndef FTYPE_H
#define FTYPE_H

#define Gemini_version "Version 3.1.3"
#define Gemini_date    "16-AUG-2023"


#endif // FTYPE_H

// 2.3.2  12/21/2020  updated
// 2.4.1  02/28/2021  About dialog
// 2.5.1  03/03/2021
// Global revision, files location, AME masses
// 2.6.1    03/05/2021  Final
// 2.6.2  03/19/2021 modifications with lise2016 path
// 2.6.3  04/18/2021 additional basePATH2

// 2.6.4  04/25/2021
// revision of paths, creation of arg-file reading utiltiy

// 2.6.5  04/25/2021
// revision of open and save dialogs

// 2.6.6  04/26/2021
// LISErootPATH = QCoreApplication::applicationDirPath();

// 2.6.7  05/25/2021
// typos correction

// 2.6.8  06/03/2021
// corrections in HighDpiScaling attribute at start

// 2.7.1  06/09/2021
// filename corrections in the case of net names //****.***.**/**

// 2.7.2  06/10/2021
// adaptation to Qt6 by Oleg

// 2.7.3  07/12/2021  corrrection for lise.ini path
// 2.7.4  12/03/21  modified for non-latin ame2016 path
// 2.7.5  12/04/21  read/write modification for const QString&
// 2.7.6  12/29/21  compatability with linux

// 2.7.7  01/05/22  test functions
// 2.7.8  01/05/22  print GEMINI properties

// 3.0.0  01/06/22  update with original current GEMINI version

// 3.0.2  07/13/22  migrating to Qt63
// 3.0.3  06/18/23  Changed the following files to .accdb file types- chart.tbl, mass.tbl, mass_tf.tbl, and lise2016.dbf
// 3.0.4  06/26/23  DAK - database updates, pallete changed
// 3.0.5  06/26/23  DAK - bug fixed
// 3.0.6  06/26/23  commented debug lines
// 3.0.7  06/26/23  no DBF lines in the code
// --------------------------------------  06/26/23  difference in IMF oberved!!!
// 3.0.8  06/26/23  correction in FRDMFinder for "return" location
// 3.0.9  06/26/23  optimization of "Mass" class
// 3.0.10 07/01/23  compiled with MSVC
// --------------------------------------
// 3.0.11 07/02/23  added two different database connections
// 3.0.12 07/02/23  removed dak_chart and chart.cpp now uses sql
// 3.0.13 07/04/23  issue with tabs: mode -> tab_mode. Modifications in Mass.cpp fordatabase closing
//                  still difference between versions
// 3.0.14 07/04/23  files chart.tbl, mass.tbl, mass_tf.tbl -- erased from package.
// 3.0.15 07/04/23  temporary tab_mode=1. It looks like it'some Qt 6.5.1 bug.
// 3.1.01 07/04/23  middle version has been changed
// 3.1.02 07/12/23  SQLite compatibaly established
// 3.1.08 08/16/23  Corrections with database path
