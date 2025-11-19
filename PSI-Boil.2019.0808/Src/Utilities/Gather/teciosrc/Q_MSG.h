 #ifndef Q_MSG_H
 #define Q_MSG_H
 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___3257
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
 #define MAX_STATUS_LINE_MSG_LEN 255
 #define MAX_RUNNING_COORDS_TEXT_LEN 80
#include "TranslatedString.h"
EXTERN ___372 ___4476(const char  *___2871, char       **___2698); EXTERN void Information(tecplot::___4216 format, ...); EXTERN void ___4445(tecplot::___4216 format, ...); EXTERN void ___1175(tecplot::___4216 format, ...);
 #endif 
