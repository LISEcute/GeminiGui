#include "gm_mainwindow.h"

extern FILE *mfopen(const QString& filename, const char* operand);

void MainWindow::write_cs_file(const QString & filename_cs)
{
   // int II = 1;
    float _SIGMA = 100;
    FILE *file_cs;
    file_cs=mfopen(filename_cs, "wt");
    int _INPUT = 2;
    /*
    fJ
*/
    _INPUT = 2;
    if(file_cs) {
       fprintf(file_cs,"!-title- Zcomp  Ncomp  Energy  CSfus  Input _AJNUC\n");
       fprintf(file_cs,"!Gemini %d %d %.3f %.3g %d %.3g",iZCN,iACN-iZCN,fEx,_SIGMA, _INPUT, (_INPUT==1 ? 0: l0));
    }

          fclose(file_cs);
}

