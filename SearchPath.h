#ifndef SEARCHPATH_H
#define SEARCHPATH_H

#include "Stringable.h"
#include "Definitions.h"

#include <fstream>
#include <sstream>

class SearchPath: public Stringable, public std::vector<unsigned int>
{
    public:
        SearchPath();
        explicit SearchPath(const std::vector<unsigned int> &path);
        virtual ~SearchPath();

        bool operator==(const SearchPath &other) const;
        std::vector<unsigned int> operator()(unsigned int start, unsigned int end) const;
        virtual std::string toString() const;

        void rewire(const std::vector<Connection> &connections, unsigned int newNode, unsigned int patternLength);
};

#endif
