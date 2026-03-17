#include <iomanip>
#include <omp.h>

#include "precision.h"
#include "linsol_error.h"

//----------------------------------------------------------------------------------------

/**
 * @brief Print error message.
 */
linsol_error::linsol_error(const std::string & function_name,
                           const std::string & error_message) {

   std::cout << "FUNCTION: " << function_name << " - ERROR: " << error_message << std::endl;

}

//----------------------------------------------------------------------------------------

/**
 * @brief Print memory error message.
 */
linsol_error::linsol_error(const std::string & function_name,
                           const std::string & error_message,
                           unsigned long int nbytes_requested) {

   #pragma omp critical (error_print)
   {

      // Retrieve thread id
      type_OMP_iReg mythid = omp_get_thread_num();

      // Dump error
      std::cout << "THREAD ID: " << mythid << std::endl;
      std::cout << "FUNCTION: " << function_name << " - ERROR: "<< error_message << std::endl;
      rExt req = static_cast<rExt>(nbytes_requested) / static_cast<rExt>(MemUM);
      std::cout << std::fixed << std::setprecision(5);
      std::cout << "REQUESTED:           " << req << " " << UM_label << std::endl;
      std::cout.flush();

   }

}
