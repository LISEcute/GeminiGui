#include <QFile>
#include <QTextStream>
//#include <QDebug>

#include "CChart.h"
#include "qdebug.h"

int const CChart::iZmax = 136;
CChart* CChart::fInstance = 0;

#include <QSqlQuery>
#include <QSqlRecord>

extern QSqlDatabase Geminidb;
extern std::vector<std::vector<QVariant>> dakChart;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//****************************************************
  /**
   * Constructor reads in files with neutron and 
   * proton rick limits
   */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CChart::CChart()
{    
   if (Geminidb.isValid())
    {
        dakChart.clear();
        QSqlQuery query(Geminidb);
        query.prepare("SELECT * FROM chart ORDER BY Z ASC");
        query.exec();
        //====================================================================================
        if (!query.exec())
        {
            qDebug() << "Search Didn't work";
        }
        //====================================================================================
        while (query.next())
        {
            std::vector<QVariant> row; // this creates a 1d array
            for (int i= 0; i < query.record().count();i++)
            {
                row.push_back(query.value(i)); // here we add the values into the 1d vecotr
            }
            dakChart.push_back(row); // we append the 1d vector into the 2d vector
        }
        query.finish();
    }
   else
    {
        qDebug() << "Gemini DB connection is not valid.";
    }

  isotope = new SIsotope[iZmax+1];
  iZindex = new int [iZmax+1];

  int iZ,iAmin,iAmax,i = 0;
  for (;;)
      {
      iZ = dakChart[i][0].toInt();
      iAmin = dakChart[i][1].toInt();
      iAmax = dakChart[i][2].toInt();
      if (iZ >= iZmax) break;
      isotope[iZ].iAmin = iAmin;
      isotope[iZ].iAmax = iAmax;
      i++;
      }

  //construct index file
  iMassDim = 0;
  for ( int i=0;i<iZmax+1;i++)
      {
      iZindex[i] = iMassDim;
      iMassDim += isotope[i].iAmax - isotope[i].iAmin + 1;
      }

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CChart* CChart::instance()
{
    if (fInstance == 0) {fInstance = new CChart;}
    return fInstance;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//***********************************************************
  /**
   *descructor
   */
CChart::~CChart()
{
  //descructor
  delete [] isotope;
  delete [] iZindex;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//**************************************************
  /**
   *returns the maxium iA for a given element that will be considered 
   *in the decay 
   \param iZ is the proton number
   */
int CChart::getAmin(int iZ)
{
  //returns the minimum A value for a given Z
  if (iZ > iZmax) 
    {
      cout << "Chart above its limits" << endl;
      cout << "iZ=" << iZ << endl;
      abort();
    }
  return isotope[iZ].iAmin;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//**************************************************
  /**
   * Returns the minimum iA for a given element that will be considered 
   *in the decay 
  \param iZ is the proton number
  */

int CChart::getAmax(int iZ)
{
  //returns the maximum A value for a given Z
  if (iZ > iZmax) 
    {
      cout << "CChart above its limits" << endl;
      cout << "iZ=" << iZ << endl;
      abort();
    }
  return isotope[iZ].iAmax;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*******************************************************
/**
 * Returns the index number of a particular nuclide
\param iZ is the proton number
\param iA is the mass number
 */
int CChart::getIndex(int iZ, int iA)
{
  //returns the index in the mass table for a specified Z,A
  if (iZ < 0) 
      {
      cout << " Z < 0 in chart" <<endl;
      return -1;
      }
  else if (iZ > iZmax)
      {
      //cout << " Z > iZmax in CChart " << endl;
      return -1;
      }
  if (iA < isotope[iZ].iAmin || iA > isotope[iZ].iAmax)
      {
      //cout << " outside chart of nuclides defined in CChart" << endl;
      return -1;
      }
 //qDebug() << "iZmax = " << iZmax;
 //qDebug() << "iZindex[" << iZ << "] = " << iZindex[iZ];// + iA - isotope[iZ].iAmin;
//      qDebug() << iZindex[iZ] << iA << isotope[iZ].iAmin << "get index";
  return iZindex[iZ] + iA - isotope[iZ].iAmin; 
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
