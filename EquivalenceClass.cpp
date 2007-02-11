#include "EquivalenceClass.h"

using std::vector;
using std::string;
using std::ostringstream;

EquivalenceClass::EquivalenceClass()
{
}

EquivalenceClass::EquivalenceClass(const vector<unsigned int> &units)
{
    clear();
    for(unsigned int i = 0; i < units.size(); i++)
        push_back(units[i]);
}

EquivalenceClass::~EquivalenceClass()
{
}

EquivalenceClass EquivalenceClass::computeOverlapEC(const EquivalenceClass &other) const
{
    EquivalenceClass overlap;
    for(unsigned int i = 0; i < size(); i++)
        for(unsigned int j = 0; j < other.size(); j++)
            if(at(i) == other[j])
                overlap.add(at(i));

    return overlap;
}

double EquivalenceClass::computeOverlapRatio(const EquivalenceClass &other) const
{
    return static_cast<double>(computeOverlapEC(other).size()) / other.size();
}

bool EquivalenceClass::operator==(const EquivalenceClass &other) const
{
    if(size() != other.size())
        return false;

    for(unsigned int i = 0; i < size(); i++)
        if(at(i) != other[i])
            return false;

    return true;
}

bool EquivalenceClass::has(unsigned int unit) const
{
    for(unsigned int i = 0; i < size(); i++)
        if(at(i) == unit)
            return true;

    return false;
}

bool EquivalenceClass::add(unsigned int unit)
{
    for(unsigned int i = 0; i < size(); i++)
        if(at(i) == unit)
            return false;

    push_back(unit);
    return true;
}

LexiconUnit* EquivalenceClass::makeCopy() const
{
    return new EquivalenceClass(*this);
}

string EquivalenceClass::toString() const
{
    ostringstream sout;

    sout << "E[";
    if(size() > 0)
    {
        for(unsigned int i = 0; i < size() - 1; i++)
            sout << "P" << at(i) << " | ";
        if(size() > 0) sout << "P" << back();
    }
    sout << "]";

    return sout.str();
}
