#ifndef CMDLINE_H
#define CMDLINE_H

/*
 *
 * cmdline.h: parse command line of a BRIAN program
 * BRIAN Software Package Version 3.0
 *
 * $Id: cmdline.h 509 2017-03-27 20:15:06Z kruggel $
 *
 * 0.10 (07/10/09): initial version
 * 0.20 (04/09/10): transparent compression implemented
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (10/01/14): documented, introduced class, split from file
 * 0.41 (02/03/15): init argVector std::explicitly
 * v406 (28/09/16): bumped to version 3.0
 *
 * some stuff retained from the Vista library, we thankfully acknowledge:
 *
 * Copyright 1994 University of British Columbia
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appears in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation. UBC makes no representations about the suitability of this
 * software for any purpose. It is provided "as is" without std::express or
 * implied warranty.
 *
 * Author: Arthur Pope, UBC Laboratory for Computational Intelligence
 *
 */

/*! \file
    \brief Provides an command line interface for BRIAN programs.
*/

/*! \enum repnType
    \brief Symbolic data representation.
*/

enum class repnType {
	unknown = 0,				///< unknown representation (=error)
	bit = 1,				///< 1-bit integer, [0, 1]
	byte = 2,				///< 8-bit integer, [0, 255]
	shrt = 3,				///< 16-bit integer, [-32768, 32767]
	sint = 4,				///< 32-bit integer, [-2**31, 2**31-1]
	flt = 5,				///< 32-bit IEEE floating point
	string = 6,				///< null-terminated string
	image = 7,				///< image
	graph = 8				///< graph
};

//! Collects information about a data representation

struct repnInfo {
	const char*	name;			//!< name string
	const size_t	size;			//!< size, in bytes
	const unsigned int precision;		//!< precision, in bits
	const repnType	repn;			//!< representation code
};

extern const repnInfo repnTable[];
extern repnType lookupRepn(const char *name);
extern size_t lookupRepnSize(repnType repn);

//! Collects information about a data dictionary entry

struct dictEntry {
	const char* key;			//!< keyword
	const int value;			//!< value
};

using vDict = std::vector<dictEntry>;

extern const vDict boolDict;

//! Represents command line options taking multiple values as a vector

struct argVector {
	unsigned int number;			//!< number of arguments
	void* 	vector;				//!< vector of arguments

	argVector()
		: number(0), vector(nullptr) { }
	~argVector()
		{ delete [] static_cast<unsigned char *>(vector); }
	argVector(const argVector& b) = delete;
	argVector(argVector&& b)
		: number(b.number), vector(b.vector) { b.vector = nullptr; }
	argVector& operator=(const argVector& b) = delete;
	argVector& operator=(argVector&& b)
		{ assert(this != &b); delete [] static_cast<unsigned char *>(vector);
		  vector = b.vector; b.vector = nullptr; number = b.number; return *this; }
};

//! Collects information about a command line option

class option {
	const vDict* getDict() const
		{ return (repn == repnType::bit && dict == nullptr)? &boolDict: dict; }
	int	convertDictToken(const vDict& dc, const char* s) const;
	void	convertToken(void* p, const vDict* dp, const char* s) const;
	void	printDict(const vDict& dc, const char* s) const;
	void	printDictValue(const vDict& dc, const int v) const;
	void	printValue(const vDict* dp, const void* p) const;
public:
	static bool required;			//!< provide distinguished values for the "found" field of an option
	static bool optional;			// exists just once
	const char* key;			//!< keyword signaling option
	repnType repn;				//!< type of value supplied by option
	unsigned int number;			//!< number of values supplied
	void*	value;				//!< location for storing value(s)
	bool*	found;				//!< whether optional arg found
	const vDict* dict;			//!< optional dict of value keywords
	const char* blurb;			//!< on-line help blurb

	bool	matches(const char* k) const
		{ return k && strcmp(key,k) == 0; }
	void	checkRequired(const bool seen) const					//! checks that a required option was seen.
		{ if (seen || found != &required) return;
		  throw optException("Option -%s must be specified", key); }
	void	print() const;
	void	convert(int& arg, const int argc, char **argv);
	void	cleanDuplicateOption()
		{ printf("Duplicate -%s option; using last.\n", key); if (number) return;
		  argVector *v = reinterpret_cast<argVector*>(value);			// delete its previous value
		  unsigned char *p = reinterpret_cast<unsigned char *>(v->vector);
		  delete [] p; v->vector = nullptr; v->number = 0; }
	void	getFiles(const int fd, int& argc, char **argv);
};

//! Collects information for command line parsing

class cmdline {
	option*	opts;				//!< vector of options
	unsigned int nopts;
	int	argc;				//!< number of command line tokens
	char**	argv;				//!< array of command line tokens

	unsigned int lookupOption(const char* key) const				//! returns option id for a key
		{ for (unsigned int i = 0; i < nopts; i++) {
			if (opts[i].matches(key)) return i; };
		  throw optException("Option -%s not defined in option table", key); }
public:
	cmdline(int _argc, char** _argv, const unsigned int _nopts = 0, option* _opts = nullptr)
		: opts(_opts), nopts(_nopts), argc(_argc), argv(_argv)
		{ setenv("LANG", "C", 1); }
	void	usage(const char* others = nullptr)					//! prints information about how to use a program.
		{ const char *build = "$Rev: 440 $";
		  printf("\nUsage: %s (rev %d) <options>", argv[0], atoi(build+5));
		  if (others) printf(" %s", others);
		  printf(", where <options> includes:\n    -help\n\tDisplay this help and exit.\n");
		  for (unsigned int i = 0; i < nopts; i++) opts[i].print(); }
	bool	parse();
	void	parse(FILE** inp, FILE** outp);
	int	args() const								//! returns number of remaining tokens
		{ return argc; }
};

#endif

