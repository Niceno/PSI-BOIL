 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <string>
#include <stack>
#include "ThirdPartyHeadersEnd.h"
#include "StandardIntegralTypes.h"
 #if defined MAKEARCHIVE
 #define common_strutil_API 
 #else
#include "common_strutil_Exports.h" 
 #endif
namespace tecplot { class common_strutil_API ___3439 { public: ___3439(); ___3439(___3439 const& ___2886); ___3439(std::string const& source); ___3439( std::string const& source, char               delimiter); ___3439( std::string const& source, std::string const& delimiters); virtual ~___3439(); ___3439& operator=(___3439 const& ___2886); virtual void append(std::string const& source); virtual void reset(); virtual void reset(std::string const& source); virtual void reset( std::string const& source, std::string const& delimiters); virtual std::string const& source() const; virtual std::string::size_type position() const; virtual void savePosition(); virtual void forgetPosition(); virtual void eraseToPosition(); virtual void restorePosition(); virtual bool hasNextToken() const; virtual bool hasNextToken(char customDelimiter) const; virtual bool hasNextToken(std::string const& customDelimiters) const; virtual bool hasNextInt64Token() const; virtual bool hasNextInt64Token(char customDelimiter) const; virtual bool hasNextInt64Token(std::string const& customDelimiters) const; virtual bool hasNextDoubleToken() const; virtual bool hasNextDoubleToken(char customDelimiter) const; virtual bool hasNextDoubleToken(std::string const& customDelimiters) const; virtual bool hasNextTokenBeginningWith(std::string const& pattern) const; virtual bool hasNextTokenBeginningWith( std::string const& pattern, char               customDelimiter) const; virtual bool hasNextTokenBeginningWith( std::string const& pattern, std::string const& customDelimiters) const; virtual bool hasNextTokenEndingWith(std::string const& pattern) const; virtual bool hasNextTokenEndingWith( std::string const& pattern, char               customDelimiter) const; virtual bool hasNextTokenEndingWith( std::string const& pattern, std::string const& customDelimiters) const; virtual bool hasNextTokenContaining(std::string const& pattern) const; virtual bool hasNextTokenContaining( std::string const& pattern, char               customDelimiter) const; virtual bool hasNextTokenContaining( std::string const& pattern, std::string const& customDelimiters) const; virtual bool hasNextChar() const; virtual bool hasNextChar(char matchChar) const; virtual bool hasNext(std::string const& pattern) const; virtual bool hasNextInt64() const; virtual bool hasNextDouble() const; virtual bool hasNextLine() const; virtual std::string nextToken(); virtual std::string nextToken(char customDelimiter); virtual std::string nextToken(std::string const& customDelimiters); virtual std::string peekNextToken() const; virtual std::string peekNextToken(char customDelimiter) const; virtual std::string peekNextToken(std::string const& customDelimiters) const; virtual int64_t nextInt64Token(); virtual int64_t nextInt64Token(char customDelimiter); virtual int64_t nextInt64Token(std::string const& customDelimiters); virtual double nextDoubleToken(); virtual double nextDoubleToken(char customDelimiter); virtual double nextDoubleToken(std::string const& customDelimiters);
virtual bool skipTokens(size_t NumTokensToSkip); virtual void skipWhitespace(); virtual void skipDelimiters(); virtual void skipDelimiters(char customDelimiter); virtual void skipDelimiters(std::string const& customDelimiters); virtual char nextChar(); virtual char peekNextChar() const; virtual std::string peekNext(std::string::size_type numChars) const; virtual std::string next(std::string::size_type numChars); virtual std::string remaining(); virtual std::string peekRemaining() const; virtual std::string::size_type numCharsRemaining() const; virtual int64_t nextInt64(); virtual double nextDouble(); virtual std::string nextLine(); virtual bool findInLine(std::string const& pattern); virtual bool find(std::string const& pattern); static std::string const noDelimiters; static std::string const whitespaceDelimiters; static std::string const newlineDelimiters; private: void ___1930( std::string const& source, std::string const& delimiters); std::string                        m_sourceStr; std::string::size_type             m_sourcePos; std::string                        m_defaultDelimiters; std::stack<std::string::size_type> m_savedPositions; }; }
