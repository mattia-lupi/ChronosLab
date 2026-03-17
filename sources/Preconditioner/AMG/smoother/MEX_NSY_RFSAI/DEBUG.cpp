#include "DEBUG.h"
#include <iostream>

FILE *dbfile;

void Open_DebugLog(){

   if (DEBUG){
      dbfile = fopen("DEBUG","w");
      if (dbfile == nullptr) std::cout << "ERROR OPENING DEBUG LOGFILE" << std::endl;
      fprintf(dbfile,"********** INIZIO DEBUG ***********\n\n");
   }

}

void Close_DebugLog(){

   if (DEBUG) fclose(dbfile);

}
