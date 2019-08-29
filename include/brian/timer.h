#ifndef TIMER_H
#define TIMER_H

/*
 *
 * timer.h: measure micro-second time intervals
 * BRIAN Software Package Version 3.0
 *
 * $Id: timer.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (01/11/09): initial version
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Functions for measuring micro-second time intervals.
*/

#include <sys/time.h>

//! Timer class

class timer {
	struct timeval base;			//!< holds time record
	void	set()									//! sets current time.
		{ gettimeofday(&base, 0); }
public:
	timer()										//! allocates a timer and sets current time.
		 : base() { set(); }
	void	start()									//! sets current time.
		{ set(); }
	double	elapsed() const								//! evaluates time interval from start to current time.
		{ struct timeval now; gettimeofday(&now, 0);
		  long int s = now.tv_sec-base.tv_sec, f = now.tv_usec-base.tv_usec;
		  if (f < 0) { s--; f += 1000000.0; };
		  double t = s+double(f/1000000.0); return t; }
};

//! \relates timer
inline void print(const timer& b)
//! prints time interval from start to current time.
{ 
	printT(b.elapsed()); printf("s\n");
}

#endif
