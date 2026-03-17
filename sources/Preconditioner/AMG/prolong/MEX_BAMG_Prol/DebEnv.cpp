#include "DebEnv.h"

//----------------------------------------------------------------------------------------

// Creates the object
DebugEnvironment::DebugEnvironment(){

}

//----------------------------------------------------------------------------------------

// Deletes the object
DebugEnvironment::~DebugEnvironment(){

   // Close log files
   CloseDebugLog();

}

//----------------------------------------------------------------------------------------

// Set-up of the object.
void DebugEnvironment::SetDebEnv(const int nthreads_in, const char *mode){

   nthreads = nthreads_in;

   if (DEBUG){

      // Open log files
      OpenDebugLog(mode);

      // Set AMG levels to be printed (aAMG level indexes start from 0)
      iLevPrint.resize(3);
      iLevPrint[0] = 0;
      iLevPrint[1] = 1111;
      iLevPrint[2] = 1112;

   }

}

//----------------------------------------------------------------------------------------

// Check if the level has to be printed
bool DebugEnvironment::ChkPrtLevel(iReg ilev){

   if (DEBUG){

      iExt size = iLevPrint.size();
      for ( iExt i = 0; i < size; i++ ) {
         if ( iLevPrint[i] == ilev ) return true;
      }

   }

   return false;

}

//----------------------------------------------------------------------------------------

// Open log files
void DebugEnvironment::OpenDebugLog(const char *mode){

   if (DEBUG && !OPEN_LOG){

      OPEN_LOG = true;

      // Open the log unit associated to ranks
      std::string logfile_name = "M3ELib_Logfile";
      r_logfile = fopen(logfile_name.data(),mode);

      // Resize the log unit for threads
      t_logfile.resize(nthreads);

      // Open the log unit associated to threads and ranks
      for (iReg i = 0; i < nthreads; i++){
         std::stringstream ss;
         ss << std::setw(2) << std::setfill('0') << i;
         std::string myid_label = ss.str();
         std::string logfile_name = "t_M3ELib_Logfile." + myid_label;
         t_logfile[i] = fopen(logfile_name.data(),mode);
      }

   }

}

//----------------------------------------------------------------------------------------

// Close log files
void DebugEnvironment::CloseDebugLog(){

   if (DEBUG && OPEN_LOG){
      OPEN_LOG = false;
      fclose(r_logfile);
      for (iReg i = 0; i < nthreads; i++) fclose(t_logfile[i]);
   }

}

//----------------------------------------------------------------------------------------

DebugEnvironment DebEnv;
