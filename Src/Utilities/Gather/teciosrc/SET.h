 #pragma once
 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___3498
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
#include "ThirdPartyHeadersBegin.h"
#include <algorithm>
#include <set>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "RawArray.h"
#include "CodeContract.h"
 #define ___2894(X,Y)           ((int)(((X)-1)/(Y)+1)*(Y))
static size_t const ___3477 = 8*sizeof(___3481); static size_t const ___3496 = static_cast<___3478>(1)<<(___3477-1);
 #if defined _DEBUG
 #  define USE_FUNCTIONS_FOR_SETS
 #endif
struct ___3500 { ___3491 size; ___3479 data; };
 #define ___2054(___3474) ((___3474)==NULL)
inline size_t ___3480(___3499 ___3474) { REQUIRE(VALID_REF(___3474)); return ___3474->size / ___3477 * sizeof(___3481); } EXTERN ___3499 ___29(___372 ___3572); EXTERN void ___937(___3499 *___3474); EXTERN ___372 ___3494(void       *___2096, ___90  ___492); EXTERN ___372 ___1199(___3499     ___3474, ___3491 max_val, ___372  ___3572); EXTERN ___372 AllocAndCopySet(___3499&  ___1119, ___3499   ___3654); EXTERN ___372 ___674(___3499    ___1119, ___3499    ___3654, ___372 ___3572); EXTERN ___372 ___83(___3499 ___1119, ___3499 ___3654, ___372 ___3572); EXTERN void ___491(___3499 ___3474);
 #if defined USE_FUNCTIONS_FOR_SETS
EXTERN ___372 ___17(___3499     ___3474, ___3491 ___2399, ___372  ___3572);
 #else
 #if defined __cplusplus
inline ___372 ___17(___3499     ___3474, ___3491 ___2399, ___372  ___3572) { if (___3474 && (___2399 + 1 <= ___3474->size || ___1199(___3474, ___2399 + 1, ___3572))) { ___3491 word = ___2399 / ___3477; ___3478 bit = (___3478)1 << (___2399 % ___3477); ___3474->data[word] |= bit; return ___4224; } else return ___1303; }
 #elif defined TECPLOTKERNEL
 #define ___17(___3474,___2399,___3572) \
 (((___3474) && \
 ((___2399)+1 <= (___3474)->size || \
 ___1199((___3474), (___2399)+1, (___3572)))) \
 ? (((___3474)->data[(___2399) / ___3477].___1344((___3478)1 << ((___2399) % ___3477))), ___4224) \
 : ___1303)
 #else
 #define ___17(___3474,___2399,___3572) \
 (((___3474) && \
 ((___2399)+1 <= (___3474)->size || \
 ___1199((___3474), (___2399)+1, (___3572)))) \
 ? (((___3474)->data[(___2399) / ___3477] |= (___3478)1 << ((___2399) % ___3477)), ___4224) \
 : ___1303)
 #endif
 #endif
EXTERN ___372 AddRangeToSet(___3499 ___3474, ___3491 rangeStart, ___3491 rangeEnd); EXTERN void ___3332(___3499     ___3474, ___3491 ___2399); EXTERN void ___955(___3499     ___3474, ___3491 ___2400); EXTERN ___372 ___1953(___3499     ___3474, ___3491 ___2400, ___372  ___3569);
 #if defined USE_FUNCTIONS_FOR_SETS
EXTERN ___372 ___1954(___3499     ___3474, ___3491 ___2399);
 #else
 #if defined __cplusplus
inline ___372 ___1954(___3499     ___3474, ___3491 ___2399) { if (___3474 && (0 <= ___2399 && ___2399 < ___3474->size)) { ___3491 word = ___2399 / ___3477; ___3478 bit = (___3478)1 << (___2399 % ___3477); return (___3474->data[word]&bit) != 0; } else return ___1303; }
 #elif defined TECPLOTKERNEL
 #define ___1954(___3474,___2399)  ((___3474 && (0<=(___2399) && (___2399)<(___3474)->size)) \
 ? ((___3474)->data[(___2399)/___3477].load()&((___3478)1<<((___2399)%___3477)))!=0 \
 : ___1303)
 #else
 #define ___1954(___3474,___2399)  ((___3474 && (0<=(___2399) && (___2399)<(___3474)->size)) \
 ? ((___3474)->data[(___2399)/___3477]&((___3478)1<<((___2399)%___3477)))!=0 \
 : ___1303)
 #endif
 #endif
EXTERN ___372 ___2013(___3499 ___3474); EXTERN ___372 ___1820(___3499 ___3474); EXTERN ___3491 ___2401(___3499 ___3474); EXTERN ___372 ___2031(___3499 ___3474); EXTERN ___3491 ___1759(___3499     ___3474, ___3491 ___3680); EXTERN ___3491 ___1767(___3499     ___3474, ___3491 ___3680); EXTERN ___372 ___1173(___3499  ___3475, ___3499  ___3476); ___3499 intersection( ___3499 ___3475, ___3499 ___3476); EXTERN ___372 ___2060(___3499 ___484, ___3499 ___2971); EXTERN ___3491 ___2402(___3499 ___3474, ___3491    ___2400); EXTERN ___3491 ___2865(___3499 ___3474, ___3491    ___2864); EXTERN ___372 ___675(___3499     ___1124, ___3491 ___1123, ___3499     ___3661, ___3491 ___3660); EXTERN void ___3558(___3499     ___3474, ___3491 ___3556, ___3491 ___3557, ___3491 ___3554);
 #define ___1744(___3474) (___1759((___3474), ___333))
 #define ___1749(___3474)  (___1767((___3474), ___333))
 #define ___1470(___2400, ___3474) \
 for (___2400 = ___1744((___3474)); \
 ___2400 != ___333; \
 ___2400 = ___1759((___3474), (___2400)))
 #define ForAllMembersInEntIndexSet(___2400, ___3474) \
 for (___2400 = static_cast<___1170>(___1744((___3474))); \
 ___2400 != static_cast<___1170>(___333); \
 ___2400 = static_cast<___1170>(___1759((___3474), (___2400))))
 #define ___1469(___2400, ___3474) \
 for (___2400 = ___1749((___3474)); \
 ___2400 != ___333; \
 ___2400 = ___1767((___3474), (___2400)))
namespace tecplot { template <typename T> std::vector<T> ___4192(___3499 ___2098) { REQUIRE(VALID_REF(___2098) || ___2098 == 0); std::vector<T> ___3356; size_t const count = ___2401(___2098); if (count != 0) { ___3356.reserve(count); ___3491 ___2083; ___1470(___2083,___2098) ___3356.push_back(static_cast<T>(___2083)); } return ___3356; } template <typename T> inline std::set<T> ___4184(___3499 const set) { REQUIRE(VALID_REF_OR_NULL(set)); ___1170 ___4312; std::set<T> ___3356; if (set != NULL) { ___1470(___4312, set) { ___3356.insert(static_cast<T>(___4312)); } } return ___3356; } template <typename CONTAINER> ___3499 ___4184( CONTAINER const& ___2097, bool             isSorted = true) { REQUIRE(IMPLICATION(isSorted && !___2097.empty(), ___2097[___2097.size()-1] == ___333 || ___2097[___2097.size()-1] == *std::max_element(&___2097[0], &___2097[0]+___2097.size()))); ___3499 ___3356 = ___29(___1303); if (___3356 == NULL) throw std::bad_alloc(); if (!___2097.empty()) { typename CONTAINER::value_type largestMember = static_cast<typename CONTAINER::value_type>(___333); if (isSorted) { for (typename CONTAINER::value_type const* iter = &___2097[___2097.size()-1]; iter >= &___2097[0]; --iter) if ((largestMember = *iter) != static_cast<typename CONTAINER::value_type>(___333)) break; } else { largestMember = *std::max_element(&___2097[0], &___2097[0]+___2097.size()); } if (largestMember != static_cast<typename CONTAINER::value_type>(___333)) { if (!___1199(___3356, static_cast<___3491>(largestMember + 1), ___1303)) { ___937(&___3356); throw std::bad_alloc(); } typename CONTAINER::value_type const* itemsArray = &___2097[0]; size_t const ___2810 = ___2097.size(); for (size_t ___1990 = 0; ___1990 < ___2810; ++___1990) if (itemsArray[___1990] != static_cast<typename CONTAINER::value_type>(___333)) (void)___17(___3356,static_cast<___3491>(itemsArray[___1990]),___1303); } } ENSURE(VALID_REF(___3356)); return ___3356; } template <typename T> void ___4183( ___3499       ___2098, ___3267<T>& ___3356) { REQUIRE(VALID_REF(___2098) || ___2098 == 0); size_t const count = ___2401(___2098); if (count != 0) { ___3356.reserve(count); ___3356.___3501(count); T* ___3358 = &___3356[0]; size_t ___2863 = 0; ___3491 ___2083; ___1470(___2083,___2098) ___3358[___2863++] = static_cast<T>(___2083); } else { ___3356.___3501(0); } } }
