#ifndef BRIANEX_H
#define BRIANEX_H

/*
 *
 * brianException.h: Exception handling in BRIAN
 * BRIAN Software Package Version 3.0
 *
 * $Id: brianException.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.30 (20/01/13): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions and classes for exception handling in BRIAN.
*/

//! Handles a floating point exception.

struct fpException : public std::exception {
	char msg[256];							//!< holds message

	fpException(const int c);
	fpException(const fpException& e)
		{ memcpy(msg,e.msg,256); }	
	~fpException() _GLIBCXX_USE_NOEXCEPT { }
	const char *what () const _GLIBCXX_USE_NOEXCEPT			//!< shows message
		{ return msg; }
};

//! Handles an option exception.

struct optException : public std::exception {
	char msg[256];							//!< holds message

	optException(const char* m);
	optException(const char* m, const int c);
	optException(const char* m, const int c0, const int c1);
	optException(const char* m, const char* s0, const char* s1);
	optException(const char* m, const int c, const char* s);
	optException(const char* m, const char* c);
	optException(const optException& e)
		{ memcpy(msg,e.msg,256); }	
	~optException()  _GLIBCXX_USE_NOEXCEPT { }
	const char *what () const _GLIBCXX_USE_NOEXCEPT			//!< shows message
		{ return msg; }
};

//! Handles a runtime exception.

struct rtException : public std::exception {
	char msg[256];							//!< holds message

	rtException(const char* m);
	rtException(const char* m, const char* s);
	rtException(const char* m, const int c);
	rtException(const char* m, const int c0, const int c1);
	rtException(const rtException& e)
		{ memcpy(msg,e.msg,256); }	
	~rtException()  _GLIBCXX_USE_NOEXCEPT { }
	const char *what () const  _GLIBCXX_USE_NOEXCEPT		//!< shows message
		{ return msg; }
};

extern FILE* openFile(const char* fname, const char* how);
extern void closeFile(FILE *fp);
extern void clearFpException();
extern int testFpException();
extern void checkFpException();
extern bool brianException();
extern char *__progname;

#endif
