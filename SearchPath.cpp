#include "SearchPath.h"

using std::string;
using std::vector;
using std::ostream;
using std::ostringstream;
using std::endl;

SearchPath::SearchPath()
{
}

SearchPath::SearchPath(const vector<unsigned int> &path)
:vector<unsigned int>(path)
{
}

SearchPath::~SearchPath()
{
}

bool SearchPath::operator==(const SearchPath &other) const
{
    if(size() != other.size())
        return false;

    for(unsigned int i = 0; i < size(); i++)
        if(at(i) != other[i])
            return false;

    return true;
}

vector<unsigned int> SearchPath::operator()(unsigned int start, unsigned int end) const
{
    assert(start <= end);
    assert(end < size());

    return vector<unsigned int>(begin() + start, begin() + end + 1);
}

string SearchPath::toString() const
{
    ostringstream sout;
    sout << "[";
    for(unsigned int i = 0; i < size() - 1; i++)
        sout << at(i) << " --> ";
    sout << back() << "]";

    return sout.str();
}

void SearchPath::rewire(const vector<Connection> &connections, unsigned int newNode, unsigned int patternLength)
{
    unsigned int count = 0;
    unsigned int lastRewired = 0;
    for(unsigned int i = 0; i < connections.size(); i++)
    {
        if(i > 0) assert(connections[i].second > connections[i-1].second);
        if((count > 0) && (connections[i].second < (connections[lastRewired].second + patternLength))) continue;
                     
        unsigned int offset = (patternLength - 1) * count;
        unsigned int start = connections[i].second - offset;
        unsigned int end = connections[i].second + patternLength - offset;

        assert((size() >= start) && (size() >= end));
        erase(begin() + start, begin() + end);
        insert(begin() + start, newNode);
        lastRewired = i;
        count++;
    }
}
