#include "scanner.h"

/******************************************************************************/
void scanner(char **text, char *token, int &ttype)
/*----------------------------------------------------------------------------*/
{
 char *token0;

 token0=token;

 for( ; **text == ' ' || **text == '\t' || **text == '\n'; (*text)++);

 if(**text == '\0')
  {
   ttype=EOL;
   return;
  }

 if( (**text >= 'A' && **text <= 'Z') || (**text >= 'a' && **text <= 'z') )
  {
   ttype=SYMBOL;
   while( (**text >= 'A' && **text <= 'Z') ||
	  (**text >= 'a' && **text <= 'z') ||
	  (**text >= '0' && **text <= '9') )
    {
     *token++=*(*text)++;
    }
   *token='\0';
   if(0==strcmp(token0,"sin"))  ttype=SINF;
   if(0==strcmp(token0,"cos"))  ttype=COSF;
   if(0==strcmp(token0,"tan"))  ttype=TANF;
   if(0==strcmp(token0,"asin")) ttype=ASINF;
   if(0==strcmp(token0,"acos")) ttype=ACOSF;
   if(0==strcmp(token0,"atan")) ttype=ATANF;
   if(0==strcmp(token0,"sqrt")) ttype=SQRT;
   if(0==strcmp(token0,"tanh")) ttype=TANHF;
   if(0==strcmp(token0,"abs")) ttype=ABSF;
   return;
  }

 if(**text == '"')
  {
   ttype=STRLIT;
   (*text)++;     /* preskoci navodnik */
   while( **text != '"' && **text )
    {
     *token++=*(*text)++;
    }
   (*text)++;     /* preskoci navodnik */
   *token='\0';
   return;
  }

 if(**text >= '0' && **text <= '9')
  {
   ttype=INTLIT;
   while(**text >= '0' && **text <= '9')
    {
     *token++=*(*text)++;
     if(**text=='.')
      {
       ttype=REALLIT;
       *token++=*(*text)++;
      }
     //if(ttype==REALLIT && (**text=='e' || **text=='E'))
     if(**text=='e' || **text=='E')
      {
       ttype=EXPLIT;
       *token++=*(*text)++;
      }
     if(ttype==EXPLIT && (**text=='+' || **text=='-'))
      {
       *token++=*(*text)++;
      }
    }
   *token='\0';
   return;
  }

 if(**text == '(')
  {
   ttype=LPAREN;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == ')')
  {
   ttype=RPAREN;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == ';')
  {
   ttype=SEMIC;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == ':')
  {
   ttype=COLON;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == ',')
  {
   ttype=COMMA;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == '.')
  {
   ttype=PERIOD;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == '\'')
  {
   ttype=APOST;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == '+')
  {
   ttype=PLUSOP;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == '-')
  {
   ttype=MINUSOP;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == '*')
  {
   ttype=MUXOP;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == '/')
  {
   ttype=DIVOP;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == '^')
  {
   ttype=POWOP;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 if(**text == '=')
  {
   ttype=ASSIGNOP;
   *token++=*(*text)++;
   *token='\0';
   return;
  }

 ttype=LEXERR;
 return;
}
