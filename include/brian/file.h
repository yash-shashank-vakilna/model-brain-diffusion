#ifndef FILE_H
#define FILE_H

/*
 *
 * file.h: read & write data to external media in BRIAN format
 * BRIAN Software Package Version 3.0
 *
 * $Id: file.h 439 2016-11-12 23:25:07Z frithjof $
 *
 * 0.10 (07/10/09): initial version
 * 0.20 (04/09/10): transparent compression implemented
 * 0.30 (31/12/12): released version 2.4
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
    \brief Provides an file interface for BRIAN data sets.
*/


//! Represents a key-value pair in a BRIAN file

class attribute {
protected:
	const char* strsave(const char *s) const
		{ if (s == nullptr || *s == 0) return nullptr; 
		  char *t = new char [strlen(s)+1]; strcpy(t,s); return t; }
	void	clear()
		{ delete [] key; delete [] value; } 
	void	assign(const char* _key, const char* _value)
		{ key = strsave(_key); value = strsave(_value); } 
public:
	const char* key;			//!< key string
	const char* value;			//!< value string

	attribute()									//! allocates an empty attribute.
		: key(nullptr), value(nullptr) { }
	attribute(const char* _key, const char* _value)					//! allocates an attribute for key and value.
		: key(strsave(_key)), value(strsave(_value)) { }
	attribute(const attribute& b)							//! copies an attribute.
		: key(strsave(b.key)), value(strsave(b.value)) { }
	attribute(attribute&& b)							//! moves an attribute.
		: key(b.key), value(b.value) { b.key = nullptr; b.value = nullptr; }
	virtual ~attribute() { clear(); }
	attribute& operator=(const attribute& b)					//! assign attribute from b.
		{ if (this != &b) { clear(); assign(b.key, b.value); }; return *this; }
	attribute& operator=(attribute&& b)						//! move assign attribute from b.
		{ assert(this != &b); delete [] key; key = b.key; b.key = nullptr;
		  delete [] value; value = b.value; b.value = nullptr; return *this; }
	bool	matchesKey(const char* k) const
		{ return strcmp(key, k) == 0; }
	void	updateValue(const char* val)
		{ delete [] value; value = strsave(val); }
	void	print() const
		{ printf("%s: %s\n", key, value); }
};

//! Represents a list of attributes in a BRIAN file

class attrList {
	std::list<attribute> ls;		//!< the attribute list.
public:
	attrList()									//! constructs an empty attribute list.
		: ls() { }
	attrList(const attrList& b)							//! copies from attribute list b.
		: ls(b.ls) { }
	attrList(attrList&& b)								//! moves from attribute list b.
		: ls(std::move(b.ls)) { }
	~attrList() = default;
	attrList& operator=(const attrList& b)						//! assigns from attribute list b.
		{ if (this != &b) ls = b.ls; return *this; }
	attrList& operator=(attrList&& b)						//! move assigns from attribute list b.
		{ assert(this != &b); ls = std::move(b.ls); return *this; }
	const char* lookup(const char* key) const					//! given key, looks up value.
		{ for (const auto& at : ls) if (at.matchesKey(key)) return at.value;
		  return nullptr; }
	unsigned int lookupNumber(const char* key) const
		{ for (const auto& at : ls) if (at.matchesKey(key)) return ATOU(at.value);
		  return 0; }
	void	remove(const char* key)							//! given key, removes the attribute from list.
		{ for (std::list<attribute>::iterator it = ls.begin(); it != ls.end(); it++)
			if (it->matchesKey(key)) { ls.erase(it); return; } }
	void	update(const char* key, const char* val)				//! given key, updates value with nvalue else adds an attribute.
		{ for (auto& at : ls) if (at.matchesKey(key)) { at.updateValue(val); return; };
		  ls.push_back({key, val}); }
	void	write(FILE* fp) const;
	void	read(FILE* fp);
	void	clear()
		{ ls.clear(); }
	void	add(const attribute& a)
		{ ls.push_back(a); }
	void	writeString(FILE *fp, const char *str) const;
	bool	readString(FILE *fp, int ch, char* value);
	void	print() const
		{ for (const auto& a: ls) a.print(); }
};

//! Represents a data object in a BRIAN file

struct bundle : public attribute {
	attrList at;				//!< list of attributes
	unsigned char* data;			//!< points to anon data
	size_t	length;				//!< length to anon data
	repnType repn;				//!< representation of anon data

	bundle()									//! allocates an empty bundle.
		 : attribute(), at(), data(nullptr), length(0), repn(repnType::unknown) { }
	bundle(const char* key, const char* value, const attrList& ls)			//! allocates a bundle for (key,value) and attribute list l.
		 : attribute(key, value), at(ls), data(nullptr), length(0), repn(repnType::unknown) { }
	bundle(const bundle& b)								//! copies from bundle b.
		 : attribute(b), at(b.at), data(nullptr), length(b.length), repn(b.repn)
		{ if (length) { data = new unsigned char [length]; memcpy(data, b.data, length); }; }
	bundle(bundle&& b)								//! moves from bundle b.
		 : attribute(b), at(b.at), data(b.data), length(b.length), repn(b.repn)
		{ b.data = nullptr; }
	~bundle() { delete [] data; data = nullptr; }
	bundle&	operator=(const bundle& b)						//! assigns from bundle b.
		{ if (this != &b) { clear(); delete [] data; assign(b.key,b.value); 	// clear old info
		  	at = b.at; length = b.length; repn = b.repn;	// copy new info
		  	if (length) { data = new unsigned char [length];
				memcpy(data, b.data, length); } };
		  return *this; }
	bundle&	operator=(bundle&& b)							//! move assigns from bundle b.
		{ assert(this != &b); attribute(std::move(b)); at = std::move(b.at);
		  length = b.length; repn = b.repn;
		  delete [] data; data = b.data; b.data = nullptr; return *this; }
	void	copyAttributes(const attrList& _at)					//! copies attribute list at to this bundle.
		{ at = _at; }
	void	resize(size_t l)							//! resize data container to size l bytes.
		{ if (data) { if (l == length) return; delete [] data; data = nullptr; };
		  if (l) data = new unsigned char [l];
		  length = l; }
	void	copyStream(const std::vector<unsigned char>& v)
		{ resize(v.size()); memcpy(data,&v[0],length); }
	void	read(FILE* fp, size_t& offset);
	void	write(FILE* fp, size_t& offset, const bool comp);
	void	print() const
		{ at.print(); }
};	

using bundleList = std::list<bundle>;

extern bool testHeader(FILE* fp, bundleList& list);
extern void readHeader(FILE* fp, bundleList& list);
extern void readFile(FILE* fp, bundleList& list);
extern void writeFile(FILE* fp, bundleList& list);

#endif

