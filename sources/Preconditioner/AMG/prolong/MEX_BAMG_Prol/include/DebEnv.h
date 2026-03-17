
#pragma once

#include <sstream>  // std::stringstream
#include <string>   // std::string
#include <iomanip>  // std::setw
#include <vector>   // std::vector

#include "precision.h"

//----------------------------------------------------------------------------------------

//#define DEBUG true
#define DEBUG false

// Flag to activate debug locally in BAMG
#define BAMG_DEBUG true

#define VERB_LEV 2

//----------------------------------------------------------------------------------------

/**
 * class DebugEnvironment.
 * @brief This class is used to manage the Debug Environment.
 */
class DebugEnvironment {

   //-------------------------------------------------------------------------------------

   // preivate members
   private:
   bool OPEN_LOG = false;

   // public object
   public:

   type_OMP_iReg nthreads = 1;

   /**
    * @brief MPI rank log files.
    */
   FILE *r_logfile;

   /**
    * @brief OMP threads log files.
    */
   std::vector<FILE*> t_logfile;

   /**
    * @brief List of AMG levels to be printed.
    */
   std::vector<iReg> iLevPrint;

   //-------------------------------------------------------------------------------------

   // public functions
   public:

   /**
    * @brief Open log files.
    */
   void OpenDebugLog(const char *mode);

   /**
    * @brief Close log files.
    */
   void CloseDebugLog();

   /**
    * @brief Constructs the object.
    */
   DebugEnvironment();

   /**
    * @brief Deletes the object.
    */
   ~DebugEnvironment();

   /**
    * @brief Set-up of the object.
    */
   void SetDebEnv(const int ntreads_in, const char *mode);

   /**
    * @brief Check if the level has to be printed
    */
   bool ChkPrtLevel(iReg ilev);

};

// Set DebEnv as global
extern DebugEnvironment DebEnv;

//----------------------------------------------------------------------------------------
