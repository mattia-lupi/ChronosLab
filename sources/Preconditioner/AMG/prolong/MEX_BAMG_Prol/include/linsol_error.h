/**
 * @file linsol_error.h
 * @brief This header is used to manage the error messages.
 * @date July 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent
 */

#pragma once

#include<string>
#include<iostream>

#define kB 1024L
#define MB 1024L*kB
#define GB 1024L*MB

#define INFO_SZ 6

#define MemUM MB

#if MemUM == kB
#define UM_label "kb"
#elif MemUM == MB
#define UM_label "Mb"
#elif MemUM == GB
#define UM_label "Gb"
#endif

/**
 * class linsol_error.
 * @brief This class is used to print error messages.
 *
 */
class linsol_error{

   public:

   /**
    * @brief Print error message.
    */
   linsol_error(const std::string & function_name,
                const std::string & error_message);

   /**
    * @brief Print memory error message.
    */
   linsol_error(const std::string & function_name,
                const std::string & error_message,
                unsigned long int nbytes_requested);

};
