#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "NRMatrix.h"

#include <utility>
#include <vector>

typedef std::pair<unsigned int, unsigned int> Connection;
typedef NRMatrix<std::vector<Connection> > ConnectionMatrix;
typedef std::pair<double, double> SignificancePair;
typedef std::pair<unsigned int, unsigned int> Range;

#endif
