#ifndef RDSGRAPH_H
#define RDSGRAPH_H

#include "RDSNode.h"
#include "ADIOSUtils.h"
#include "global.h"
#include "MiscUtils.h"
//#include "ParseTree.h"

#include <string>
#include <sstream>

// true, if both pvalues are less than alpha
bool isPatternSignificant(const SignificancePair &pvalues, double alpha);
bool operator<(const SignificancePair &a, const SignificancePair &b);

class RDSGraph: public Stringable
{
    public:
        RDSGraph();
        explicit RDSGraph(const std::vector<std::vector<std::string> > &sequences);

        std::vector<std::string> generate() const;
        std::vector<std::string> generate(const SearchPath &searchPath) const;
        std::vector<std::string> generate(unsigned int node) const;
        void distill(const ADIOSParams &params);

        virtual std::string toString() const;

    private:
        unsigned int corpusSize;
        std::vector<SearchPath> paths;
        std::vector<RDSNode> nodes;

        void buildInitialGraph(const std::vector<std::vector<std::string> > &sequences);
        bool distill(const SearchPath &searchPath, const ADIOSParams &params);
        bool generalise(const SearchPath &searchPath, const ADIOSParams &params);
        bool bootstrapStage(SignificantPatternInfo &bestPatternInfo, EquivalenceClassInfo &bestECInfo, const SearchPath &searchPath, const ADIOSParams &params, const Range &context, SearchPathInfo &searchPathInfo);
        bool generalisationStage(SignificantPatternInfo &bestPatternInfo, const SearchPath &searchPath, const ADIOSParams &params, const ConnectionMatrix &connections);

        // generalise and bootstrap
        SearchPath computeGeneralisedSubpaths(EquivalenceClass &ec, ConnectionMatrix &connections, const SearchPath &searchPath, unsigned int prefixStart, unsigned int slotIndex, unsigned int postfixEnd);
        SearchPath bootstrap(BootstrapInfo &bootstrapInfo, const SearchPath &searchPath, double overlapThreshold) const;

        // compute matrix and pattern searching function
        void computeConnectionMatrix(ConnectionMatrix &connections, const SearchPath &searchPath) const;
        void computeDescentsMatrix(NRMatrix<double> &flows, NRMatrix<double> &descents, const ConnectionMatrix &connections) const;
        bool findSignificantPatterns(std::vector<Range> &patterns, std::vector<SignificancePair> &pvalues, const ConnectionMatrix &connections, const NRMatrix<double> &flows, const NRMatrix<double> &descents, double eta, double alpha) const;

        // rewiring and update functions
        void updateAllConnections();
        void rewire(const std::vector<Connection> &connections, unsigned int ec);
        void rewire(const std::vector<Connection> &connections, EquivalenceClass *ec);
        void rewire(const std::vector<Connection> &connections, SignificantPattern *pattern);
        std::vector<Connection> getRewirableConnections(const ConnectionMatrix &connections, const Range &bestSP, double alpha) const;

        // pattern searching auxiliary functions
        bool findRightDescentRow(unsigned int &descentRow, const NRMatrix<double> &descents, const Range &range, unsigned int column, double eta) const;
        bool findLeftDescentRow(unsigned int &descentRow, const NRMatrix<double> &descents, const Range &range, unsigned int column, double eta) const;
        double computeRightSignificance(const ConnectionMatrix &connections, const NRMatrix<double> &flows, const std::pair<unsigned int, unsigned int> &descentPoint, double eta) const;
        double computeLeftSignificance(const ConnectionMatrix &connections, const NRMatrix<double> &flows, const std::pair<unsigned int, unsigned int> &descentPoint, double eta) const;
        double findBestRightDescentColumn(unsigned int &bestColumn, NRMatrix<double> &pvalueCache, const ConnectionMatrix &connections, const NRMatrix<double> &flows, const NRMatrix<double> &descents, const Range &pattern, double eta) const;
        double findBestLeftDescentColumn(unsigned int &bestColumn, NRMatrix<double> &pvalueCache, const ConnectionMatrix &connections, const NRMatrix<double> &flows, const NRMatrix<double> &descents, const Range &pattern, double eta) const;

        std::vector<Connection> filterConnections(const std::vector<Connection> &init_cons, unsigned int start_offset, const SearchPath &search_path) const;
        std::vector<Connection> getAllNodeConnections(unsigned int nodeIndex) const;
        std::vector<Connection> getEquivalenceConnections(const EquivalenceClass &ec) const;

        // print functions
        std::string printSignificantPattern(const SignificantPattern &sp) const;
        std::string printEquivalenceClass(const EquivalenceClass &ec) const;
        std::string printNode(unsigned int node) const;
        std::string printPath(const SearchPath &path) const;
};

void printInfo(const ConnectionMatrix &connections, const NRMatrix<double> &flows, const NRMatrix<double> &descents);

#endif
