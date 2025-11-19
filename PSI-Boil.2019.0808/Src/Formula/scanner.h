#ifndef SCANNER_H
#define SCANNER_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cstring>

enum {LEXERR,
      SYMBOL,
      INTLIT,
      REALLIT,
      EXPLIT,
      STRLIT,
      LPAREN,
      RPAREN,
      SEMIC,
      COLON,
      COMMA,
      PERIOD,
      APOST,
      PLUSOP,
      MINUSOP,
      MUXOP,
      DIVOP,
      SINF,
      COSF,
      TANF,
      ASINF,
      ACOSF,
      ATANF,
      SQRT,
      POWOP,
      ASSIGNOP,
      EOL,
      UNDEF};

void scanner(char **text, char *token, int &ttype);

#endif
