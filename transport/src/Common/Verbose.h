#ifndef VERBOSE_H_
#define VERBOSE_H_


#define VERBOSE

#define VERBOSE_DYN

#ifdef VERBOSE
	#include<R.h>
	#include<cstdio>
	#ifdef VERBOSE_DYN
		extern bool verbose_mode;
		//#define eprintf(format, ...) if(verbose_mode) Rprintf(format, ##__VA_ARGS__);
		#define eprintf(...) if(verbose_mode) Rprintf(__VA_ARGS__);
	#else
		//#define eprintf(format, ...) Rprintf(format, ##__VA_ARGS__);
		#define eprintf(...) Rprintf(__VA_ARGS__);
	#endif
#else
	//#define eprintf(format, ...);
	#define eprintf(...);
#endif





#endif
