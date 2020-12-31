/*
 * MTB_DoubleMetaphone.cpp
 *
 *  Created on: 06.07.2017
 *      Author: Rukasz
 *
 *      Copyright 2007, Stephen Lacy slacy@slacy.com
 *
 * This code is a derivative work from an implementation by Maurice Aubrey maurice@hevanet.com,
 * and modified to use STL vector and string classes instead of bare pointers.
 *
 * Original Comments by Maurice Aubrey:
 *
 * All rights reserved.
 *
 * This code is based heavily on the C++ implementation by Lawrence Philips and incorporates several bug fixes courtesy of Kevin Atkinson kevina@users.sourceforge.net.
 *
 * This module is free software; you may redistribute it and/or modify it under the same terms as Perl itself.
 */

#include "MTB_DoubleMetaphone.h"


string MTB_DoubleMetaphone::getName() {
	return "Exact";
}

double MTB_DoubleMetaphone::getRelativeValue(string var1, string var2) {
	if (var1.compare(var2) == 0) {
		return 1.0;
	}
	return 0.0;
}

double MTB_DoubleMetaphone::getAbsoluteValue(string var1, string var2) {

	return getRelativeValue(var1, var2);
}

void CreateUpper(string &s) {
	for (unsigned int i = 0; i < s.length(); i++) {
		s[i] = toupper(s[i]);
	}
}

int IsVowel(string &s, unsigned int pos) {
	char c;

	if (pos >= s.length())
		return 0;

	c = s[pos];
	if ((c == 'A') || (c == 'E') || (c == 'I') || (c == 'O') || (c == 'U')
			|| (c == 'Y')) {
		return 1;
	}

	return 0;
}

int SlavoGermanic(string &s) {
	if ((char *) strstr(s.c_str(), "W"))
		return 1;
	else if ((char *) strstr(s.c_str(), "K"))
		return 1;
	else if ((char *) strstr(s.c_str(), "CZ"))
		return 1;
	else if ((char *) strstr(s.c_str(), "WITZ"))
		return 1;
	else
		return 0;
}

char GetAt(string &s, unsigned int pos) {
	if (pos >= s.length()) {
		return '\0';
	}

	return s[pos];
}

void SetAt(string &s, unsigned int pos, char c) {
	if (pos >= s.length()) {
		return;
	}

	s[pos] = c;
}

/*
 Caveats: the START value is 0 based
 */
int StringAt(string &s, unsigned int start, unsigned int length, ...) {
	char *test;
	const char *pos;
	va_list ap;

	if (start >= s.length()) {
		return 0;
	}

	pos = (s.c_str() + start);
	va_start(ap, length);

	do {
		test = va_arg(ap, char *);
		if (*test && (strncmp(pos, test, length) == 0)) {
			return 1;
		}
	} while (strcmp(test, ""));

	va_end(ap);

	return 0;
}

// void doubleMetaphone(const string &str, vector<string> *codes) {
//
// }

