
#ifndef WK_ERROR_FORMATTER
#define WK_ERROR_FORMATTER

// https://stackoverflow.com/questions/12261915/how-to-throw-stdexceptions-with-variable-messages
#include <stdexcept>
#include <sstream>
#include <stdexcept>

class ErrorFormatter {
public:
    ErrorFormatter() {}
    ~ErrorFormatter() {}

    template <typename Type>
    ErrorFormatter & operator << (const Type & value) {
        stream_ << value;
        return *this;
    }

    std::string str() const { return stream_.str(); }
    operator std::string () const { return stream_.str(); }

    enum ConvertToString {
        to_str
    };
    std::string operator >> (ConvertToString) { return stream_.str(); }

private:
    std::stringstream stream_;

    ErrorFormatter(const ErrorFormatter &);
    ErrorFormatter & operator = (ErrorFormatter &);
};

# endif
