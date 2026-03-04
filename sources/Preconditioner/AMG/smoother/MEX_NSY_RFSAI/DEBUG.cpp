#include "DEBUG.h"

FILE *dbfile;

void Open_DebugLog(){

   if (DEBUG){
      dbfile = fopen("DEBUG","w");
      if (dbfile == nullptr) cout << "ERROR OPENING DEBUG LOGFILE" << endl;
      fprintf(dbfile,"********** INIZIO DEBUG ***********\n\n");
   }

}

void Close_DebugLog(){

   if (DEBUG) fclose(dbfile);

}
