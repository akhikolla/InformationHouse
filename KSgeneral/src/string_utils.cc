#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <climits>
#include <sstream>

#include "string_utils.h"

using namespace std;

long string_to_long(const string& s)
{
    const char* p = s.c_str();
    char* endptr = NULL;

    errno = 0;
    long value = strtol(p, &endptr, 10);

    if ((errno == ERANGE && (value == LONG_MAX || value == LONG_MIN)) || (errno != 0 && value == 0)) {
        throw runtime_error(string("Error converting string to long: ") + s + "\n" + strerror(errno));
    }

    if (endptr == p) {
        throw runtime_error(string("Error converting string to long: ") + s + "\n" + "No digits were found.");
    }

    if (*endptr != '\0') {
        throw runtime_error(string("Trailing characters during conversion of string to long:") + s + "\n");
    }

    return value;
}

double string_to_double(const string& s)
{
    const char* p = s.c_str();
    char* endptr = NULL;

    errno = 0;
    double value = strtod(p, &endptr);

    if ((value == 0.0) && (endptr == p)) {
        throw runtime_error(string("Error converting string to double '") + s + "'");
    }

    if (errno != 0) {
        throw runtime_error(string("Error converting string to double '") + s + "'\n" + strerror(errno));
    }

    //if (*endptr != '\0') {
    //    throw runtime_error(string("Trailing characters during conversion of string to double:") + s + "\n");
    //}

    return value;
}

vector<string> split(const string& s, char delimiter)
{
    vector<string> substrings;
    int current_substring_start = 0;

    for (int i = 0; i < (int)s.size(); ++i) {
        if (s[i] == delimiter) {
            substrings.push_back(s.substr(current_substring_start, i-current_substring_start));
            current_substring_start = i+1;
        }
    }

    substrings.push_back(s.substr(current_substring_start, s.size()-current_substring_start));

    return substrings;
}

vector<double> read_comma_delimited_doubles(const string& line)
{
    vector<string> substrings = split(line, ',');

    // Ignore optional trailing comma
    if ((substrings.size() >= 1) && (substrings.back() == "")) {
        substrings.resize(substrings.size()-1);
    }

    vector<double> numbers(substrings.size());
    transform(substrings.begin(), substrings.end(), numbers.begin(), string_to_double);
    return numbers;
}

string vector_to_string(const vector<double>& v)
{
    stringstream ss;
    for (int i = 0; i < (int)v.size(); ++i) {
        ss << v[i];
        if (i != int(v.size())-1) {
            ss << ", ";
        }
    }
    ss << endl;
    return ss.str();
}
