#include "stdafx.h"
#include "MASTER.h"
 #define ___3858
#include "GLOBAL.h"
#include "CodeContract.h"
#include "STRUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "STRLIST.h"
using tecplot::___1095; static ___372 StringListItemDestructor(void*      ___2096, ___90 ___492) { REQUIRE(VALID_REF(___2096)); REQUIRE(VALID_REF(*static_cast<char**>(___2096)) || *static_cast<char**>(___2096) == NULL); ___4276(___492); char** StringRef = static_cast<char**>(___2096); if (*StringRef != NULL) { ___1528(*StringRef, "string"); *StringRef = NULL; } ENSURE(*StringRef == NULL); return ___4224; } static ___372 StringListItemDuplicator(void*      ___3947, void*      ___3643, ___90 ___492) { REQUIRE(VALID_REF(___3947)); REQUIRE(VALID_REF(___3643)); REQUIRE(VALID_REF(*static_cast<char**>(___3643)) || *static_cast<char**>(___3643) == NULL); ___4276(___492); ___372 ___2038 = ___4224; char** TargetStringRef = static_cast<char**>(___3947); char** SourceStringRef = static_cast<char**>(___3643); if (*SourceStringRef != NULL) ___2038 = ((*TargetStringRef = ___1133(___1095(*SourceStringRef))) != NULL); else *TargetStringRef = NULL; ENSURE(VALID_REF(*TargetStringRef) || *TargetStringRef == NULL); ENSURE(VALID_BOOLEAN(___2038)); return ___2038; } ___372 ___3846(___3837 ___3817) { ___372 ___2065 = ArrayListIsValid(reinterpret_cast<___134>(___3817)); if (___2065) { ___2225 stringCount = ___101(reinterpret_cast<___134>(___3817));
 #if defined PERFORM_EXPENSIVE_STRLIST_TESTS
{ for (___2225 index = 0; index < stringCount; index++) { char* string = ___100(reinterpret_cast<___134>(___3817), index); if (string != NULL && !VALID_REF(string)) { ___2065 = ___1303; break; } } }
 #else
{ if (stringCount > 0) { char* string = ___100(reinterpret_cast<___134>(___3817), 0); if (string != NULL && !VALID_REF(string)) { ___2065 = ___1303; } } if (___2065 && stringCount > 1) { char* string = ___100(reinterpret_cast<___134>(___3817), stringCount - 1); if (string != NULL && !VALID_REF(string)) { ___2065 = ___1303; } } }
 #endif 
} ENSURE(VALID_BOOLEAN(___2065)); return ___2065; } void ___3822(___3837 ___3818) { REQUIRE(___3846(___3818)); ArrayListDeleteAllItems(reinterpret_cast<___134>(___3818), StringListItemDestructor, 0); ENSURE(___3846(___3818) && ___3824(___3818) == 0); } void ___3839(___3837 ___3818, ___2225     ___3852, ___2225     ___682) { REQUIRE(___3846(___3818)); REQUIRE(0 <= ___3852 && ___3852 <= ___3824(___3818) - 1); REQUIRE(1 <= ___682 && ___3852 + ___682 <= ___3824(___3818)); ArrayListDeleteItems(reinterpret_cast<___134>(___3818), ___3852, ___682, StringListItemDestructor, 0); ENSURE(___3846(___3818)); } void ___3838(___3837 ___3818, ___2225     ___3852) { REQUIRE(___3846(___3818)); REQUIRE(0 <= ___3852 && ___3852 <= ___3824(___3818) - 1); ArrayListDeleteItems(reinterpret_cast<___134>(___3818), ___3852, 1, StringListItemDestructor, 0); ENSURE(___3846(___3818)); } void ___3826(___3837* ___3818) { REQUIRE(VALID_REF(___3818)); REQUIRE(*___3818 == NULL || ___3846(*___3818)); if (*___3818 != NULL) ArrayListDealloc(reinterpret_cast<___134*>(___3818), StringListItemDestructor, 0); ENSURE(*___3818 == NULL); } ___3837 ___3819() { ___3837 ___3357 = reinterpret_cast<___3837>(ArrayListAlloc(0, ArrayListType_CharPtr, NULL, 0)); ENSURE(___3357 == NULL || ___3846(___3357)); return ___3357; } ___372 ___3821(___3837 ___3818, char const*   ___3811) { REQUIRE(___3846(___3818)); REQUIRE(___3811 == NULL || VALID_REF(___3811)); ___372 ___2038 = ___3841(___3818, ___3824(___3818), ___3811); ENSURE(___3846(___3818)); ENSURE(VALID_BOOLEAN(___2038)); return ___2038; } ___2225 ___3824(___3837 ___3818) { REQUIRE(___3846(___3818)); ___2225 ___3357 = ___101(reinterpret_cast<___134>(___3818)); ENSURE(___3357 >= 0); return ___3357; } char* ___3832(___3837 ___3818, ___2225     ___3852) { REQUIRE(___3846(___3818)); REQUIRE(0 <= ___3852 && ___3852 <= ___3824(___3818) - 1); char* ___3357; char const* StringRef = ___3833(___3818, ___3852); if (StringRef == NULL) ___3357 = NULL; else ___3357 = ___1133(___1095(StringRef)); ENSURE(___3357 == NULL || VALID_REF(___3357)); return ___3357; }
 #if !defined USE_MACROS_FOR_FUNCTIONS
char const* ___3834(___3837 ___3818, ___2225     ___3852) { REQUIRE(___3846(___3818)); REQUIRE(0 <= ___3852 && ___3852 <= ___3824(___3818) - 1); char const* ___3357 = ___3835(___3818, ___3852); ENSURE(___3357 == NULL || VALID_REF(___3357)); return ___3357; }
 #endif
___372 ___3841(___3837 ___3818, ___2225     ___3852, char const*   ___3811) { REQUIRE(___3846(___3818)); REQUIRE(___3852 >= 0); REQUIRE(___3811 == NULL || VALID_REF(___3811)); ___372       ___2038; ArrayListItem_u ItemCopy; if (___3811 != NULL) { ItemCopy.___472 = ___1133(___1095(___3811)); ___2038 = (ItemCopy.___472 != NULL); } else { ItemCopy.___472 = NULL; ___2038 = ___4224; } if (___2038) ___2038 = ArrayListSetItem(reinterpret_cast<___134>(___3818), ___3852, ItemCopy, StringListItemDestructor, 0); ENSURE(___3846(___3818)); ENSURE(VALID_BOOLEAN(___2038)); return ___2038; } ___372 ___3836(___3837 ___3818, ___2225     ___3852, char const*   ___3811) { REQUIRE(___3846(___3818)); REQUIRE(___3852 >= 0); REQUIRE(___3811 == NULL || VALID_REF(___3811)); ___372       ___2038; ArrayListItem_u ItemCopy; if (___3811 != NULL) { ItemCopy.___472 = ___1133(___1095(___3811)); ___2038 = (ItemCopy.___472 != NULL); } else { ItemCopy.___472 = NULL; ___2038 = ___4224; } if (___2038) ___2038 = ArrayListInsertItem(reinterpret_cast<___134>(___3818), ___3852, ItemCopy); ENSURE(___3846(___3818)); ENSURE(VALID_BOOLEAN(___2038)); return ___2038; } ___3837 ___3823(___3837 ___3818) { REQUIRE(___3846(___3818)); ___3837 ___3357 = reinterpret_cast<___3837>(ArrayListCopy(reinterpret_cast<___134>(___3818), StringListItemDuplicator, 0)); ENSURE(___3357 == NULL || (___3846(___3357) && ___3824(___3357) == ___3824(___3818))); return ___3357; } ___372 ___3820(___3837 ___3944, ___3837 ___3640) { REQUIRE(___3846(___3944)); REQUIRE(___3846(___3640)); ___3837 SourceCopy = ___3823(___3640); ___372 ___2038 = (SourceCopy != NULL); if (___2038) { ArrayListAppend(reinterpret_cast<___134>(___3944), reinterpret_cast<___134>(SourceCopy)); ArrayListDealloc(static_cast<___134*>(static_cast<void*>(&SourceCopy)), NULL, 0); } ENSURE(___3846(___3944)); ENSURE(VALID_BOOLEAN(___2038)); return ___2038; } char* ___3845(___3837 ___3818) { REQUIRE(___3846(___3818)); size_t ___2222 = 0; ___2225 ___682 = ___3824(___3818); if (___682 >= 1) { ___2225 ___1924; for (___1924 = 0, ___2222 = strlen("\n") * (___682 - 1); ___1924 < ___682; ___1924++) { char* ___3811 = ___100(reinterpret_cast<___134>(___3818), ___1924); if (___3811 != NULL) ___2222 += strlen(___3811); } } char* ___3357 = ___23(___2222 + 1, char, "new line separated string"); if (___3357 != NULL) { ___2225 ___1924; for (___1924 = 0, strcpy(___3357, ""); ___1924 < ___682; ___1924++) { char* ___3811 = ___100(reinterpret_cast<___134>(___3818), ___1924); if (___1924 != 0) strcat(___3357, "\n"); if (___3811 != NULL) strcat(___3357, ___3811); } } ENSURE(___3357 == NULL || VALID_REF(___3357)); return ___3357; } ___3837 ___3829(char const* ___3811)
{ REQUIRE(VALID_REF(___3811)); ___3837 ___3357 = ___3819(); int ___3683; int EndIndex; for (___3683 = EndIndex = 0; ___3357 != NULL; EndIndex++) { if (___3811[EndIndex] == '\n' || ___3811[EndIndex] == '\0') { int ___2222 = EndIndex - ___3683; char* SubString = ___23(___2222 + 1, char, "sub string"); if (SubString != NULL) { ___676(SubString, ___3811, ___3683, ___2222); ___3821(___3357, SubString); ___1528(SubString, "sub string"); if (___3811[EndIndex] != '\0') ___3683 = EndIndex + 1; else break; } else { ___3826(&___3357); ___3357 = NULL; break; } } } ENSURE(___3357 == NULL || ___3846(___3357)); return ___3357; } char** ___3843(___3837 ___3818) { REQUIRE(___3846(___3818)); char** ___3357 = static_cast<char**>(ArrayListToArray(reinterpret_cast<___134>(___3818), StringListItemDuplicator, 0)); ENSURE(___3357 == NULL || VALID_REF(___3357)); return ___3357; } ___3837 ___3827(char const** ___3814, ___2225    ___682) { REQUIRE((___682 == 0 && ___3814 == NULL) || (___682 >= 1 && VALID_REF(___3814))); ___3837 ___3357 = reinterpret_cast<___3837>(ArrayListFromArray(static_cast<void*>(___3814), ___682, ArrayListType_CharPtr, StringListItemDuplicator, 0)); ENSURE(___3357 == NULL || ___3846(___3357)); return ___3357; }
 #define ISJOINCHAR(c) ((c == ';') || (c == '+'))
static void SkipWhiteSpaceOrComma(char const** ___683) { REQUIRE(VALID_REF(___683) && VALID_REF(*___683)); while (___2080(**___683) || (**___683 == ',')) (*___683)++; } static ___372 GetNextSubString(char const** OriginalCPtr, char**       NextSubString) { REQUIRE(VALID_REF(OriginalCPtr) && (VALID_REF(*OriginalCPtr))); REQUIRE(VALID_REF(NextSubString)); ___372 ___2038 = ___4224; *NextSubString = NULL; char const* ___683 = *OriginalCPtr; SkipWhiteSpaceOrComma(&___683); char InsideDelimiter = '\0'; if (*___683 == '"'|| *___683 == '\'') { InsideDelimiter = *___683; ___683++; } char const* CStart = ___683; while (*___683 && ((InsideDelimiter && (*___683 != InsideDelimiter)) || (!InsideDelimiter && (*___683 != ',')       && !ISJOINCHAR(*___683)  && !___2080(*___683)))) { if (InsideDelimiter  && (*___683 == '\\')  && (*(___683 + 1) == InsideDelimiter || *(___683 + 1) == '\\')) ___683 += 2; else ___683++; } if (InsideDelimiter && (*___683 != InsideDelimiter)) ___2038 = ___1303; if (___2038 && CStart < ___683) { size_t StrLen = static_cast<size_t>(___683 - CStart); *NextSubString = ___23(StrLen + 1, char, "GetNextSubString: NextSubString"); if (*NextSubString) { char* NPtr = *NextSubString; while (CStart < ___683) { if (*CStart == '\\') { if (*(CStart + 1) == InsideDelimiter) CStart++; else *NPtr++ = *CStart++; } *NPtr++ = *CStart++; } *NPtr = '\0'; } else ___2038 = ___1303; } if (___2038) { if (InsideDelimiter) ___683++; SkipWhiteSpaceOrComma(&___683); *OriginalCPtr = ___683; } ENSURE(VALID_BOOLEAN(___2038)); return ___2038; } ___3837 ___3828(char const* ___3811) { REQUIRE(VALID_REF(___3811)); SkipWhiteSpaceOrComma(&___3811); REQUIRE(!ISJOINCHAR(*___3811)); ___372 ___2038 = ___4224; ___3837 ___3357 = ___3819(); char const*   ___683   = ___3811; char* CurString = NULL; while (___2038 && *___683 != '\0') { char*     NextSubString = NULL; ___372 WantsToJoin   = ___1303; if (ISJOINCHAR(*___683)) { WantsToJoin = ___4224; ___683++; SkipWhiteSpaceOrComma(&___683); } ___2038 = GetNextSubString(&___683, &NextSubString); if (___2038) { if (WantsToJoin) ___3937(&CurString, '\n'); if (NextSubString != NULL && strlen(NextSubString) != 0) ___2038 = ___3939(&CurString, NextSubString, ___1303, ___1303); else if (CurString == NULL) CurString = ___1133(___1095("")); } if (NextSubString != NULL) ___1528(NextSubString, "StringListFromCompound: NextSubString"); if (___2038 && !ISJOINCHAR(*___683)) { ___3821(___3357, CurString); if (CurString != NULL) ___1528(CurString, "current string"); CurString = NULL; } } if (CurString != NULL) ___1528(CurString, "current string"); if (!___2038) ___3826(&___3357); ENSURE(___3357 == NULL || ___3846(___3357)); return ___3357; } char *___3844(___3837 ___3818, char          ___1816, char const*   ___473) { REQUIRE(___3846(___3818)); REQUIRE(___3824(___3818) >= 1); REQUIRE(ISJOINCHAR(___1816)); REQUIRE(VALID_REF(___473)); char* ___3357 = NULL; ___372 ___2038 = ___4224; ___2225 ___1924; ___2225 ___682;
for (___1924 = 0, ___682 = ___3824(___3818), ___2038 = ___4224; ___1924 < ___682 && ___2038; ___1924++) { char* ___3811 = ___3832(___3818, ___1924); if (___3811 != NULL && strlen(___3811) != 0) { char*       CStart = NULL; char*       CEnd = NULL; char*       EscapedString = NULL; char const* EscChar = NULL; char*       StrChar = NULL; for (StrChar = ___3811; *StrChar != '\0'; StrChar++) { for (EscChar = ___473; *EscChar != '\0'; EscChar++) if (*StrChar == *EscChar) { ___2038 = ___3937(&EscapedString, '\\'); ___2038 = ___3937(&EscapedString, '\\'); break; } ___2038 = ___3937(&EscapedString, *StrChar); } CEnd = EscapedString; while (___2038 && CEnd && *CEnd != '\0') { int ___2221 = 0; CStart = CEnd; while (*CEnd != '\0' && *CEnd != '\n') { ___2221++; if (*CEnd == '"') ___2221++; CEnd++; } char* TString = ___23(___2221 + 4, char, "temp compound sub-string"); if (TString != NULL) { if (CStart == EscapedString) { if (___1924 != 0) ___2038 = ___3937(&___3357, ' '); } else { ___2038 = ___3937(&___3357, ___1816); } char* TStr = TString; *TStr++ ='"'; while (CStart && CStart != CEnd) { if (*CStart == '"') *TStr++ = '\\'; *TStr++ = *CStart++; } *TStr++ = '"'; *TStr = '\0'; ___3939(&___3357, TString, ___1303, ___1303); ___1528(TString, "___3844"); TString = NULL; if (*CEnd) CEnd++; } else { ___2038 = ___1303; } } if (EscapedString != NULL) ___1528(EscapedString, "escaped string"); } else { if (___1924 == 0) ___3939(&___3357, "\"\"", ___1303, ___1303); else ___3939(&___3357, " \"\"", ___1303, ___1303); } if (___3811 != NULL) ___1528(___3811, "string list ___2083"); } if (!___2038) { if (___3357 != NULL) { ___1528(___3357, "___3844"); ___3357 = NULL; } } ENSURE(___3357 == NULL || VALID_REF(___3357)); return ___3357; }