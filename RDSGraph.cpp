#include "RDSGraph.h"

using std::vector;
using std::pair;
using std::string;
using std::ostream;
using std::ostringstream;
using std::endl;
#include <iostream>
bool isPatternSignificant(const SignificancePair &pvalues, double alpha)
{
    return (pvalues.first < alpha) && (pvalues.second < alpha);
}

bool operator<(const SignificancePair &a, const SignificancePair &b)
{
    double maxA = max(a.first, a.second);
    double maxB = max(b.first, b.second);

    return maxA < maxB;
}

ADIOSParams::ADIOSParams(double eta, double alpha, unsigned int contextSize, double overlapThreshold)
{
    assert((eta >= 0.0) && (eta <= 1.0));
    assert((alpha >= 0.0) && (alpha <= 1.0));

    this->eta = eta;
    this->alpha = alpha;
    this->contextSize = contextSize;
    this->overlapThreshold = overlapThreshold;
}

RDSGraph::RDSGraph()
{
    corpusSize = 0;
}

RDSGraph::RDSGraph(const vector<vector<string> > &sequences)
{
    srand(getSeedFromTime());
    buildInitialGraph(sequences);
}

void RDSGraph::distill(const ADIOSParams &params)
{
    std::cout << "eta = " << params.eta << endl;
    std::cout << "alpha = " << params.alpha << endl;
    std::cout << "contextSize = " << params.contextSize << endl;
    std::cout << "overlapThreshold = " << params.overlapThreshold << endl;

    while(true)
    {
        bool foundPattern = false;
        for(unsigned int i = 0; i < paths.size(); i++)
        {
            std::cout << "--------------------------- working on Path " << i << " of length " << paths[i].size() << " ----------------------------------" << endl;
            std::cout << printPath(paths[i]) << endl;

            if((params.contextSize < 3) || (paths[i].size() < params.contextSize))
            {
                bool foundAnotherPattern = distill(paths[i], params);
                foundPattern = foundAnotherPattern || foundPattern;
            }
            else
            {
                bool foundAnotherPattern = generalise(paths[i], params);
                foundPattern = foundAnotherPattern || foundPattern;
            }
        }
        if(!foundPattern)
            break;
    }
}

vector<string> RDSGraph::generate() const
{
    unsigned int pathIndex = static_cast<unsigned int>(floor(uniformRand() * paths.size()));
    return generate(paths[pathIndex]);
}

vector<string> RDSGraph::generate(const SearchPath &searchPath) const
{
    vector<string> sequence;
    for(unsigned int i = 0; i < searchPath.size(); i++)
    {
        vector<string> segment = generate(searchPath[i]);
        sequence.insert(sequence.end(), segment.begin(), segment.end());
    }

    return sequence;
}

vector<string> RDSGraph::generate(unsigned int node) const
{
    vector<string> sequence;

    if(nodes[node].type == LexiconTypes::Start)
        sequence.push_back("*");
    else if(nodes[node].type == LexiconTypes::End)
        sequence.push_back("#");
    else if(nodes[node].type == LexiconTypes::Symbol)
        sequence.push_back((static_cast<BasicSymbol *>(nodes[node].lexicon))->getSymbol());
    else if(nodes[node].type == LexiconTypes::EC)
    {
        EquivalenceClass *ec = static_cast<EquivalenceClass *>(nodes[node].lexicon);
        unsigned int numberOfUnits = ec->size();
        unsigned int randomUnit = static_cast<unsigned int>(floor(numberOfUnits * uniformRand()));
        vector<string> segment = generate(ec->at(randomUnit));
        sequence.insert(sequence.end(), segment.begin(), segment.end());
    }
    else if(nodes[node].type == LexiconTypes::SP)
    {
         SignificantPattern *SP = static_cast<SignificantPattern *>(nodes[node].lexicon);
         for(unsigned int i = 0; i < SP->size(); i++)
         {
             vector<string> segment = generate((*SP)[i]);
             sequence.insert(sequence.end(), segment.begin(), segment.end());
         }
    }
    else
        assert(false);

    assert(sequence.size() > 0);

    return sequence;
}

bool RDSGraph::distill(const SearchPath &searchPath, const ADIOSParams &params)
{
    Range range(0, searchPath.size() - 1);

    // look possible significant pattern found with help of equivalence class
    ConnectionMatrix connections;
    NRMatrix<double> flows, descents;
    computeConnectionMatrix(connections, searchPath, range);
    computeDescentsMatrix(flows, descents, connections, range);

    vector<Range> patterns;
    vector<SignificancePair> pvalues;
    if(!findSignificantPatterns(patterns, pvalues, connections, flows, descents, params.eta, params.alpha))
        return false;

    SignificantPattern bestPattern(searchPath(patterns.front().first, patterns.front().second));
    vector<Connection> connectionsToRewire = getRewirableConnections(connections, patterns.front(), params.alpha);
    rewire(connectionsToRewire, new SignificantPattern(bestPattern));

    std::cout << "BEST PATTERN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    std::cout << "RANGE = [" << patterns.front().first << " " << patterns.front().second << "]" << endl;
    std::cout << bestPattern << " with " << "[" << pvalues.front().first << " " << pvalues.front().second << "]" << endl;
    std::cout << connectionsToRewire.size() << " connections rewired." << endl;
    std::cout << "END BEST PATTERN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

    return true;
}

bool RDSGraph::generalise(const SearchPath &searchPath, const ADIOSParams &params)
{
    Range range(0, searchPath.size() - 1);

    // Significant Pattern Info
    bool patternFound = false;
    BootstrapInfo bestBootstrapInfo;
    SignificantPatternInfo bestPatternInfo;
    EquivalenceClassInfo bestECInfo;

    SearchPathInfo searchPathInfo;
    searchPathInfo.searchPath = searchPath;
    searchPathInfo.alreadyTested = false;
    for(unsigned int i = 0; (i + params.contextSize) <= searchPath.size(); i++)
    {
        BootstrapInfo bootstrapInfo;
        bootstrapInfo.context = Range(i, i + params.contextSize - 1);
        SearchPath boostedPath = bootstrap(bootstrapInfo, searchPath, bootstrapInfo.context.first, bootstrapInfo.context.second, params.overlapThreshold, range);

        std::cout << "[" << bootstrapInfo.context.first << " - " << bootstrapInfo.context.second << "] ";
        //std::cout << "+++++++++++++ " << boostedPath << endl;
        //std::cout << "Context = [" << bootstrapInfo.context.first << " " << bootstrapInfo.context.second << "]" << endl;
        SignificantPatternInfo patternInfo;
        EquivalenceClassInfo ecInfo;
        if(bootstrapStage(patternInfo, ecInfo, boostedPath, params, bootstrapInfo.context, searchPathInfo))
            if((!patternFound) || (patternInfo.pvalues < bestPatternInfo.pvalues))
            {
                patternFound = true;
                bestPatternInfo = patternInfo;
                bestECInfo = ecInfo;
                bestBootstrapInfo = bootstrapInfo;
            }
    }
    std::cout << endl;

    if(patternFound)
    {
    std::cout << "BEST PATTERN###########################################################################################################" << endl;
    std::cout << "RANGE = [" << bestPatternInfo.patternRange.first << " " << bestPatternInfo.patternRange.second << "]" << endl;
    std::cout << "P[" << printSignificantPattern(bestPatternInfo.pattern) << "] with " << "[" << bestPatternInfo.pvalues.first << " " << bestPatternInfo.pvalues.second << "]" << endl;
    std::cout << bestPatternInfo.connections.size() << " connections" << endl;
    if(bestECInfo.ec.size() > 1)
        std::cout << "Slot " << bestECInfo.slot << " = E[" << printEquivalenceClass(bestECInfo.ec) << "]" << endl;
    std::cout << "Context = [" << bestBootstrapInfo.context.first << " " << bestBootstrapInfo.context.second << "]" << endl;
    for(unsigned int i = 0; i < bestBootstrapInfo.encounteredECs.size(); i++)
        std::cout << i << ":  {E[" << printEquivalenceClass(bestBootstrapInfo.encounteredECs[i]) << "] <---> "
                  << printNode(bestBootstrapInfo.overlapECs[i]) << "} ---> "
                  << bestBootstrapInfo.overlapRatios[i] << endl;
    std::cout << "END BEST PATTERN#######################################################################################################" << endl;
    }

    if(!patternFound) return false;

    // rewiring process, equivalence classes are not rewired with any connections because they are all overriden by the pattern
    std::cout << "STARTS REWIRING" << endl;
    for(unsigned int i = bestBootstrapInfo.context.first + 1; i < bestBootstrapInfo.context.second; i++)
    {
        // test if current slot is in the pattern
        if((i < bestPatternInfo.patternRange.first) || (i > bestPatternInfo.patternRange.second)) continue;

        unsigned int tempSlot = i - (bestBootstrapInfo.context.first + 1);
        if(i == bestECInfo.slot)
        {
            if(bestECInfo.ec.size() > 1)
            {
                std::cout << "NEW EC CREATED: E[" << printEquivalenceClass(bestECInfo.ec) << "]" << endl;
                rewire(vector<Connection>(), new EquivalenceClass(bestECInfo.ec));
                bestPatternInfo.pattern[i - bestPatternInfo.patternRange.first] = nodes.size() - 1;
            }
        }
        else if((nodes[bestBootstrapInfo.overlapECs[tempSlot]].type == LexiconTypes::EC) && (bestBootstrapInfo.overlapRatios[tempSlot] > 0.65))
        {
            if(bestBootstrapInfo.overlapRatios[tempSlot] < 1.0)
            {
                EquivalenceClass *otherEC = static_cast<EquivalenceClass *>(nodes[bestBootstrapInfo.overlapECs[tempSlot]].lexicon);
                EquivalenceClass overlapEC = bestBootstrapInfo.encounteredECs[tempSlot].getOverlapEC(*otherEC);
                std::cout << "NEW OVERLAP EC USED: E[" << printEquivalenceClass(overlapEC) << "]" << endl;
                rewire(vector<Connection>(), new EquivalenceClass(overlapEC));
                bestPatternInfo.pattern[i - bestPatternInfo.patternRange.first] = nodes.size() - 1;
            }
            else
            {
                std::cout << "OLD OVERLAP EC USED: E[" << printNode(bestBootstrapInfo.overlapECs[tempSlot]) << "]" << endl;
                rewire(vector<Connection>(), bestBootstrapInfo.overlapECs[tempSlot]);
                bestPatternInfo.pattern[i - bestPatternInfo.patternRange.first] = bestBootstrapInfo.overlapECs[tempSlot];
            }
        }
    }
    std::cout << bestPatternInfo.connections.size() << " occurences rewired" << endl;
    rewire(bestPatternInfo.connections, new SignificantPattern(bestPatternInfo.pattern));
    std::cout << "ENDS REWIRING" << endl;

    return patternFound;
}

bool RDSGraph::bootstrapStage(SignificantPatternInfo &bestPatternInfo, EquivalenceClassInfo &bestECInfo, const SearchPath &searchPath, const ADIOSParams &params, const Range &context, SearchPathInfo &searchPathInfo)
{
    Range range(0, searchPath.size() - 1);

    bool patternFound = false;
    vector<SearchPath> previousPaths;
    for(unsigned int i = context.first + 1; i < context.second; i++)
    {
        // look for possible equivalence class
        EquivalenceClass ec;
        ConnectionMatrix connections;
        SearchPath generalPath = computeGeneralisedSubpaths(ec, connections, searchPath, context.first, i, context.second, range);

        // see if paths already tested
        bool repeated = false;
        for(unsigned int j = 0; j < previousPaths.size(); j++)
            if(previousPaths[j] == generalPath)
            {
                repeated = true;
                break;
            }

        if(repeated) continue;

        if(generalPath == searchPathInfo.searchPath)
            if(searchPathInfo.alreadyTested)
                continue;
            else
                searchPathInfo.alreadyTested = true;

        // find best pattern and ec so far
        //std::cout << "------------- " << generalPath << endl;
        //std::cout << "Slot " << i << endl;
        previousPaths.push_back(generalPath);
        SignificantPatternInfo patternInfo;
        if(generalisationStage(patternInfo, generalPath, params, connections))
            if((!patternFound) || (patternInfo.pvalues < bestPatternInfo.pvalues))
            {
                patternFound = true;
                bestPatternInfo = patternInfo;

                if((i >= bestPatternInfo.patternRange.first) && (i <= bestPatternInfo.patternRange.second))
                {
                    bestECInfo.ec = ec;
                    bestECInfo.slot = i;
                }
                else
                    bestECInfo.ec = EquivalenceClass();

                unsigned int patternSize = bestPatternInfo.pattern.size();
                for(unsigned int j = 0; j < patternSize; j++)
                    if((ec.size() > 1) && (ec.has(bestPatternInfo.pattern[j])))
                    {
                        bestECInfo.ec = ec;
                        bestECInfo.slot = i;
                        break;
                    }
            }
    }
/*
    if(patternFound)
    {
    std::cout << "BEST PATTERN???????????????????????????????????????????????????????????????????????????????????????????????????????????" << endl;
    std::cout << "RANGE = [" << bestPatternInfo.patternRange.first << " " << bestPatternInfo.patternRange.second << "]" << endl;
    std::cout << bestPatternInfo.pattern << " with " << "[" << bestPatternInfo.pvalues.first << " " << bestPatternInfo.pvalues.second << "]" << endl;
    std::cout << bestPatternInfo.connections.size() << " connections" << endl;
    if(bestECInfo.ec.getUnits().size() > 1)
        std::cout << "Slot " << bestECInfo.slot << " = " << bestECInfo.ec << endl;
    std::cout << "END BEST PATTERN???????????????????????????????????????????????????????????????????????????????????????????????????????" << endl;
    }
*/
    return patternFound;
}

bool RDSGraph::generalisationStage(SignificantPatternInfo &bestPatternInfo, const SearchPath &searchPath, const ADIOSParams &params, const ConnectionMatrix &connections)
{
    Range range(0, searchPath.size() - 1);

    // look possible significant pattern found with help of equivalence class
    vector<Range> patterns;
    vector<SignificancePair> pvalues;
    NRMatrix<double> flows, descents;
    computeDescentsMatrix(flows, descents, connections, range);
    if(!findSignificantPatterns(patterns, pvalues, connections, flows, descents, params.eta, params.alpha))
        return false;

    // return best pattern info
    bestPatternInfo.patternRange = patterns.front();
    bestPatternInfo.pattern = SignificantPattern(searchPath(patterns.front().first, patterns.front().second));
    bestPatternInfo.pvalues = pvalues.front();
    bestPatternInfo.connections = getRewirableConnections(connections, patterns.front(), params.alpha);
/*
    std::cout << "BEST PATTERN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    std::cout << "RANGE = [" << bestPatternInfo.patternRange.first << " " << bestPatternInfo.patternRange.second << "]" << endl;
    std::cout << bestPatternInfo.pattern << " with " << "[" << bestPatternInfo.pvalues.first << " " << bestPatternInfo.pvalues.second << "]" << endl;
    std::cout << bestPatternInfo.connections.size() << " connections" << endl;
    std::cout << "END BEST PATTERN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
*/
    return true;
}

string RDSGraph::toString() const
{
    ostringstream sout;

    sout << "Search Paths" << endl;
    for(unsigned int i = 0; i < paths.size(); i++)
        sout << printPath(paths[i]) << endl;

    sout << endl << "RDS Graph Nodes " << nodes.size() << endl;
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        sout << "Lexicon " << i << ": " << printNode(i) << "   ------->  " << nodes[i].parents.size() << "  [";
        for(unsigned int j = 0; j < nodes[i].parents.size(); j++)
        {
            sout << nodes[i].parents[j].first;// << "." << nodes[i].parents[j].second;
            if(j < (nodes[i].parents.size() - 1)) sout << "   ";
        }
        sout << "]" << endl;
    }

    return sout.str();
}

void RDSGraph::buildInitialGraph(const vector<vector<string> > &sequences)
{   //pad the temporary lexicon vector with empty element to acount for special start and end state
    vector<string> lexicon;
    lexicon.push_back("");
    lexicon.push_back("");

    //insert the special symbols
    nodes.push_back(RDSNode(new StartSymbol(), LexiconTypes::Start));
    nodes.push_back(RDSNode(new EndSymbol(), LexiconTypes::End));
    for(unsigned int i = 0; i < sequences.size(); i++)
    {
        vector<unsigned int> currentPath;

        //insert start state
        currentPath.push_back(0);

        //create the main part of the graph
        for(unsigned int j = 0; j < sequences[i].size(); j++)
        {
            vector<string>::iterator foundPosition = find(lexicon.begin(), lexicon.end(), sequences[i][j]);
            if(foundPosition == lexicon.end())
            {
                lexicon.push_back(sequences[i][j]);
                nodes.push_back(RDSNode(new BasicSymbol(sequences[i][j]), LexiconTypes::Symbol));
                currentPath.push_back(lexicon.size() - 1);
            }
            else
                currentPath.push_back(foundPosition - lexicon.begin()); //+2 offset for start and end states
        }

        //insert end state
        currentPath.push_back(1);

        paths.push_back(SearchPath(currentPath));
    }

    updateAllConnections();
}

void RDSGraph::computeConnectionMatrix(ConnectionMatrix &connections, const SearchPath &searchPath, const Range &range) const
{
    assert(range.first <= range.second);
    assert(range.second < searchPath.size());

    // calculate subpath distributions, symmetrical matrix
    connections = ConnectionMatrix(searchPath.size(), searchPath.size());
    for(unsigned int i = range.first; i <= range.second; i++)
    {
        connections(i, i) = getAllNodeConnections(searchPath[i]);

        // compute the column from the diagonal
        for(unsigned int j = i + 1; j <= range.second; j++)
        {
            vector<Connection> tempConnections = connections(j - 1, i);
            filterConnections(tempConnections, searchPath, Range(j, j), i);
            connections(j, i) = tempConnections;
            connections(i, j) = tempConnections;
        }
    }
}

SearchPath RDSGraph::computeGeneralisedSubpaths(EquivalenceClass &ec, ConnectionMatrix &connections, const SearchPath &searchPath, unsigned int prefixStart, unsigned int slotIndex, unsigned int postfixEnd, const Range &range)
{
    assert(range.first <= range.second);
    assert(range.second < searchPath.size());
    assert(prefixStart < slotIndex);
    assert(slotIndex < postfixEnd);

    // get the candidate connections
    vector<Connection> equivalenceConnections = getAllNodeConnections(searchPath[prefixStart]);
    filterConnections(equivalenceConnections, searchPath, Range(prefixStart, slotIndex - 1), prefixStart);
    filterConnections(equivalenceConnections, searchPath, Range(slotIndex + 1, postfixEnd), prefixStart);

    //build equivalence class
    ec = EquivalenceClass();
    for(unsigned int i = 0; i < equivalenceConnections.size(); i++)
    {
        unsigned int currentPath = equivalenceConnections[i].first;
        unsigned int currentStart = equivalenceConnections[i].second;

        equivalenceConnections[i].second = currentStart + slotIndex - prefixStart;
        ec.add(paths[currentPath][equivalenceConnections[i].second]);
    }

    //if no generalisation possible, return early
    SearchPath generalPath = searchPath;
    if(ec.size() > 1)
    {
        bool existingECFound = false;
        for(unsigned int i = 0; i < nodes.size(); i++)
            if(nodes[i].type == LexiconTypes::EC)
                if(ec.computeOverlapRatio(*(static_cast<EquivalenceClass *>(nodes[i].lexicon))) >= 1.0)
                {
                    existingECFound = true;
                    ec = EquivalenceClass();
                    generalPath[slotIndex] = i;
                    break;
                }

         //construct temporary RDSGraph to represent rewiring of the EquivalenceClass
        if(!existingECFound)
        {
            //RDSGraph tempGraph(*this);
            //tempGraph.rewire(vector<Connection>(), new EquivalenceClass(ec));
            //generalPath[slotIndex] = tempGraph.nodes.size() - 1;
            //tempGraph.computeConnectionMatrix(connections, generalPath, range);
            rewire(vector<Connection>(), new EquivalenceClass(ec));
            generalPath[slotIndex] = nodes.size() - 1;
            computeConnectionMatrix(connections, generalPath, range);
            nodes.pop_back();
            updateAllConnections();
        }
        else
            computeConnectionMatrix(connections, generalPath, range);
    }
    else
        computeConnectionMatrix(connections, generalPath, range);

    return generalPath;
}

SearchPath RDSGraph::bootstrap(BootstrapInfo &bootstrapInfo, const SearchPath &searchPath, unsigned int prefix, unsigned int postfix, double overlapThreshold, const Range &range) const
{
    assert(range.first <= range.second);
    assert(range.second < searchPath.size());

    // find all possible connections
    vector<Connection> equivalenceConnections = getAllNodeConnections(searchPath[prefix]);
    filterConnections(equivalenceConnections, searchPath, Range(postfix, postfix), prefix);

    // find potential ECs
    bootstrapInfo.encounteredECs.clear();
    for(unsigned int i = prefix + 1; i < postfix; i++)
    {
        bootstrapInfo.encounteredECs.push_back(EquivalenceClass());
        for(unsigned int j = 0; j < equivalenceConnections.size(); j++)
        {
            unsigned int currentPath = equivalenceConnections[j].first;
            unsigned int currentStart = equivalenceConnections[j].second;

            bootstrapInfo.encounteredECs.back().add(paths[currentPath][currentStart+i-prefix]);
        }
    }

    // init bootstrap data
    bootstrapInfo.overlapECs = searchPath(prefix + 1, postfix - 1);
    bootstrapInfo.overlapRatios = vector<double>(postfix - prefix - 1, 0.0);

    // bootstrap search path
    SearchPath bootstrapPath = searchPath;
    for(unsigned int i = 0; i < bootstrapInfo.encounteredECs.size(); i++)
    {
        for(unsigned int j = 0; j < nodes.size(); j++)
            if(nodes[j].type == LexiconTypes::EC)
            {
                double overlap = bootstrapInfo.encounteredECs[i].computeOverlapRatio(*static_cast<EquivalenceClass *>(nodes[j].lexicon));
                if((overlap > bootstrapInfo.overlapRatios[i]) && (overlap > overlapThreshold))
                {
                    bootstrapInfo.overlapECs[i] = j;
                    bootstrapInfo.overlapRatios[i] = overlap;
                }
            }
        bootstrapPath[prefix + i + 1] = bootstrapInfo.overlapECs[i];
    }

    return bootstrapPath;
}

void RDSGraph::computeDescentsMatrix(NRMatrix<double> &flows, NRMatrix<double> &descents, const ConnectionMatrix &connections, const Range &range) const
{
    assert(range.first <= range.second);
    assert(range.second < connections.numberOfRows());

    // calculate P_R and P_L
    flows = NRMatrix<double>(connections.numberOfRows(), connections.numberOfColumns(), -1.0);
    for(unsigned int i = range.first; i <= range.second; i++)
        for(unsigned int j = range.first; j <= range.second; j++)
            if(i > j)
                flows(i, j) = static_cast<double>(connections(i, j).size()) / connections(i-1, j).size();
            else if(i < j)
                flows(i, j) = static_cast<double>(connections(i, j).size()) / connections(i+1, j).size();
            else
                flows(i, j) = static_cast<double>(connections(i, j).size()) / corpusSize;

    // calculate D_R and D_L
    descents = NRMatrix<double>(connections.numberOfRows(), connections.numberOfColumns(), -1.0);
    for(unsigned int i = range.first; i <= range.second; i++)
        for(unsigned int j = range.first; j <= range.second; j++)
            if(i > j)
                descents(i, j) = flows(i, j) / flows(i-1, j);
            else if(i < j)
                descents(i, j) = flows(i, j) / flows(i+1, j);
            else
                descents(i, j) = 1.0;
}

bool RDSGraph::findSignificantPatterns(vector<Range> &patterns, vector<SignificancePair> &pvalues, const ConnectionMatrix &connections, const NRMatrix<double> &flows, const NRMatrix<double> &descents, double eta, double alpha) const
{
    patterns.clear();
    pvalues.clear();

    //find candidate pattern start and ends
    unsigned int pathLength = descents.numberOfRows();
    vector<unsigned int> candidateEndRows;
    vector<unsigned int> candidateStartRows;
    for(unsigned int i = 0; i < descents.numberOfRows(); i++)
    {
        for(int j = i - 1; j >= 0; j--)
            if(descents(i, j) < eta)
            {
                candidateEndRows.push_back(i - 1);
                break;
            }

        for(unsigned int j = i + 1; j < descents.numberOfColumns(); j++)
            if(descents(i, j) < eta)
            {
                candidateStartRows.push_back(i + 1);
                break;
            }
    }

    //find candidate patterns;
    vector<Range> candidatePatterns;
    for(unsigned int i = 0; i < candidateStartRows.size(); i++)
        for(unsigned int j = 0; j < candidateEndRows.size(); j++)
            if(candidateStartRows[i] < candidateEndRows[j])
                candidatePatterns.push_back(Range(candidateStartRows[i], candidateEndRows[j]));

    //for(unsigned int i = 0; i < candidatePatterns.size(); i++)
    //    std::cout << "Candidate Pattern " << i << " = " << candidatePatterns[i].first << " " << candidatePatterns[i].second << endl;

    NRMatrix<double> pvalueCache(pathLength, pathLength, 2.0);
    for(unsigned int i = 0; i < candidatePatterns.size(); i++)
    {   //std::cout << "START+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        //std::cout << "Testing pattern at [" << candidatePatterns[i].first << " -> " << candidatePatterns[i].second << "]" << endl;

        SignificancePair tempPvalues(2.0, 2.0);
        pair<unsigned int, unsigned int> descentContexts;
        tempPvalues.second = findBestRightDescentColumn(descentContexts.second, pvalueCache, connections, flows, descents, candidatePatterns[i], eta);
        tempPvalues.first = findBestLeftDescentColumn(descentContexts.first, pvalueCache, connections, flows, descents, candidatePatterns[i], eta);
        if((fabs(tempPvalues.first) > 1.0) || (fabs(tempPvalues.second) > 1.0)) continue;

        //std::cout << "P_R at " << candidatePatterns[i].second << " " << descentContexts.second << " = " << flows(candidatePatterns[i].second, descentContexts.second) << endl;
        //std::cout << "right pvalue = " << 1.0 - tempPvalues.second << " for " << connections(candidatePatterns[i].second + 1, descentContexts.second).size() << " out of " << connections(candidatePatterns[i].second, descentContexts.second).size() << endl;

        //std::cout << "P_L at " << candidatePatterns[i].first << " " << descentContexts.first << " = " << flows(candidatePatterns[i].first, descentContexts.first) << endl;
        //std::cout << "left pvalue = " << 1.0 - tempPvalues.first << " for " << connections(candidatePatterns[i].first - 1, descentContexts.first).size() << " out of " << connections(candidatePatterns[i].first, descentContexts.first).size() << endl;

        // pattern IS significant
        if(isPatternSignificant(tempPvalues, alpha))
        {
            //std::cout << "Pattern is significant at [" << 1-tempPvalues.first << " --- " << 1-tempPvalues.second << "]" << endl;
            //std::cout << "Right descent context is [" << descentContexts.second << " -> " << candidatePatterns[i].second << "]" << endl;
            //std::cout << "Left descent context is [" << descentContexts.first << " -> " << candidatePatterns[i].first << "]" << endl;

            patterns.push_back(candidatePatterns[i]);
            pvalues.push_back(tempPvalues);

            //found a MORE significant pattern
            if((patterns.size() == 1) || (pvalues.back() < pvalues.front()))
            {
                swap(patterns.front(), patterns.back());
                swap(pvalues.front(), pvalues.back());
                //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!found a better SP" << endl;
            }
        }
        //std::cout << "END----------------------------------------------------------------------------------------------" << endl;
    }

    return patterns.size() > 0;
}

vector<Connection> RDSGraph::getRewirableConnections(const ConnectionMatrix &connections, const Range &bestSP, double alpha) const
{
    vector<Connection> validConnections = connections(bestSP.second, bestSP.first);

    return validConnections;
}

void RDSGraph::rewire(const std::vector<Connection> &connections, unsigned int ec)
{
    assert(nodes[ec].type == LexiconTypes::EC);

    for(unsigned int i = 0; i < connections.size(); i++)
        paths[connections[i].first][connections[i].second] = ec;

    updateAllConnections();
}

void RDSGraph::rewire(const vector<Connection> &connections, EquivalenceClass *ec)
{
    nodes.push_back(RDSNode(ec, LexiconTypes::EC));
    rewire(connections, nodes.size() - 1);
}

void RDSGraph::rewire(const vector<Connection> &connections, SignificantPattern *pattern)
{
    nodes.push_back(RDSNode(pattern, LexiconTypes::SP));

    // split connections into groups
    unsigned int patternLength = pattern->size();
    vector<vector<Connection> > separateConnections;
    for(unsigned int i = 0; i < connections.size(); i++)
    {
        bool found = false;
        for(unsigned int j = 0; j < separateConnections.size(); j++)
            if(separateConnections[j].front().first == connections[i].first)
            {
                bool inserted = false;
                for(unsigned int k = 0; k < separateConnections[j].size(); k++)
                    if(separateConnections[j][k].second > connections[i].second)
                    {
                        inserted = true;
                        separateConnections[j].insert(separateConnections[j].begin() + k, connections[i]);
                        break;
                    }

                if(!inserted) separateConnections[j].push_back(connections[i]);

                found = true;
                break;
            }

        if(!found)
        {
            separateConnections.push_back(vector<Connection>());
            separateConnections.back().push_back(connections[i]);
        }
    }

    // rewire each group
    for(unsigned int i = 0; i < separateConnections.size(); i++)
        paths[separateConnections[i].front().first].rewire(separateConnections[i], nodes.size() - 1, patternLength);

    updateAllConnections();
}

void RDSGraph::updateAllConnections()
{
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        nodes[i].setConnections(vector<Connection>());
        nodes[i].parents.clear();
    }

    corpusSize = 0;
    for(unsigned int i = 0; i < paths.size(); i++)
    {
        corpusSize += paths[i].size();
         for(unsigned int j = 0; j < paths[i].size(); j++)
             nodes[paths[i][j]].addConnection(Connection(i, j));
    }

    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        if(nodes[i].type == LexiconTypes::SP)
        {
            SignificantPattern *sp = static_cast<SignificantPattern *>(nodes[i].lexicon);
            for(unsigned int j = 0; j < sp->size(); j++)
                nodes[sp->at(j)].parents.push_back(Connection(i, sp->find(sp->at(j))));
        }
        else if(nodes[i].type == LexiconTypes::EC)
        {
            EquivalenceClass *ec = static_cast<EquivalenceClass *>(nodes[i].lexicon);
            for(unsigned int j = 0; j < ec->size(); j++)
                nodes[ec->at(j)].parents.push_back(Connection(i, 0));
        }
    }
}

bool RDSGraph::findRightDescentRow(unsigned int &descentRow, const NRMatrix<double> &descents, const Range &range, unsigned int column, double eta) const
{
    assert(range.first >= column);

    for(unsigned int i = range.first; i <= range.second; i++)
        if(descents(i, column) < eta)
        {
            std::cout << "searched for right descent in (" << range.first << ":" << i << ", " << column << ")" << endl;
            descentRow = i;
            return true;
        }
    std::cout << "searched for right descent in (" << range.first << ":" << range.second << ", " << column << ")" << endl;
    return false;
}

bool RDSGraph::findLeftDescentRow(unsigned int &descentRow, const NRMatrix<double> &descents, const Range &range, unsigned int column, double eta) const
{
    assert(range.second <= column);

    for(int i = range.second; i >= static_cast<int>(range.first); i--)
        if(descents(i, column) < eta)
        {
            std::cout << "searching for left descent in (" << i << ":" << range.second << ", " << column << ")" << endl;
            descentRow = i;
            return true;
        }
    std::cout << "searching for left descent in (" << range.first << ":" << range.second << ", " << column << ")" << endl;
    return false;
}

double RDSGraph::computeRightSignificance(const ConnectionMatrix &connections, const NRMatrix<double> &flows, const pair<unsigned int, unsigned int> &descentPoint, double eta) const
{
    unsigned int row = descentPoint.first;
    unsigned int col = descentPoint.second;
    assert(row > col);

    double significance = 0.0;
    unsigned int patternOccurences = connections(row - 1, col).size();
    unsigned int descentOccurences = connections(row, col).size();
    for(unsigned int i = 0; i <= descentOccurences; i++)
        significance += binom(i, patternOccurences, eta * flows(row - 1, col));

    return min(max(significance, 0.0), 1.0);
}

double RDSGraph::computeLeftSignificance(const ConnectionMatrix &connections, const NRMatrix<double> &flows, const pair<unsigned int, unsigned int> &descentPoint, double eta) const
{
    unsigned int row = descentPoint.first;
    unsigned int col = descentPoint.second;
    assert(row < col);

    double significance = 0.0;
    unsigned int patternOccurences = connections(row + 1, col).size();
    unsigned int descentOccurences = connections(row, col).size();
    for(unsigned int i = 0; i <= descentOccurences; i++)
        significance += binom(i, patternOccurences, eta * flows(row + 1, col));

    return min(max(significance, 0.0), 1.0);;
}

double RDSGraph::findBestRightDescentColumn(unsigned int &bestColumn, NRMatrix<double> &pvalueCache, const ConnectionMatrix &connections, const NRMatrix<double> &flows, const NRMatrix<double> &descents, const Range &pattern, double eta) const
{
    double pvalue = 2.0;
    pair<unsigned int, unsigned int> descentPoint(pattern.second + 1, bestColumn);
    for(unsigned int i = 0; i <= pattern.first; i++)
    {
        descentPoint.second = i;
        if(!(descents(descentPoint.first, descentPoint.second) < eta)) continue;
        if(pvalueCache(pattern.second + 1, i) > 1.0)
            pvalueCache(pattern.second + 1, i) = computeRightSignificance(connections, flows, descentPoint, eta);

        if(pvalueCache(pattern.second + 1, i) < pvalue)
        {
            bestColumn = i;
            pvalue = pvalueCache(pattern.second + 1, i);
        }
    }

    return pvalue;
}

double RDSGraph::findBestLeftDescentColumn(unsigned int &bestColumn, NRMatrix<double> &pvalueCache, const ConnectionMatrix &connections, const NRMatrix<double> &flows, const NRMatrix<double> &descents, const Range &pattern, double eta) const
{
    double pvalue = 2.0;
    pair<unsigned int, unsigned int> descentPoint(pattern.first - 1, bestColumn);
    for(unsigned int i = pattern.second; i < connections.numberOfColumns(); i++)
    {
        descentPoint.second = i;
        if(!(descents(descentPoint.first, descentPoint.second) < eta)) continue;
        if(pvalueCache(pattern.first - 1, i) > 1.0)
            pvalueCache(pattern.first - 1, i) = computeLeftSignificance(connections, flows, descentPoint, eta);

        if(pvalueCache(pattern.first - 1, i) < pvalue)
        {
            bestColumn = i;
            pvalue = pvalueCache(pattern.first - 1, i);
        }
    }

    return pvalue;
}

void RDSGraph::filterConnections(vector<Connection> &connections, const SearchPath &searchPath, const Range &range, unsigned int startOffset) const
{
    for(unsigned int i = range.first; i <= range.second; i++)
    {
        vector<Connection> tempConnections;
        for(unsigned int j = 0; j < connections.size(); j++)
        {
            unsigned int currentPath = connections[j].first;
            unsigned int currentStart = connections[j].second;
            unsigned int currentPosition = currentStart + i - startOffset;
            unsigned int testNode = searchPath[i];

            if(currentPosition < paths[currentPath].size())
            {
                unsigned int currentNode = paths[currentPath][currentPosition];

                if((currentNode == testNode) ||
                  ((nodes[testNode].type == LexiconTypes::EC) && static_cast<EquivalenceClass *>(nodes[testNode].lexicon)->has(currentNode)))
                    tempConnections.push_back(connections[j]);
            }
        }
        connections = tempConnections;
    }
}

vector<Connection> RDSGraph::getAllNodeConnections(unsigned int nodeIndex) const
{
    vector<Connection> connections = nodes[nodeIndex].getConnections();
    if(nodes[nodeIndex].type == LexiconTypes::EC)
    {
        vector<Connection> tempConnections = getEquivalenceConnections(*static_cast<EquivalenceClass *>(nodes[nodeIndex].lexicon));
        connections.insert(connections.end(), tempConnections.begin(), tempConnections.end());
    }

    return connections;
}

vector<Connection> RDSGraph::getEquivalenceConnections(const EquivalenceClass &ec) const
{
    vector<Connection> connections;
    for(unsigned int i = 0; i < ec.size(); i++)
    {
        vector<Connection> tempConnections = nodes[ec[i]].getConnections();
        connections.insert(connections.end(), tempConnections.begin(), tempConnections.end());
    }

    return connections;
}

SearchPath RDSGraph::encode(const vector<string> &sequence) const
{
    SearchPath path;
    for(unsigned int i = 0; i < sequence.size(); i++)
    {
        bool foundNode = false;
        for(unsigned int j = 0; j < nodes.size(); j++)
            if(nodes[j].type == LexiconTypes::Symbol)
            {
                BasicSymbol *symbol = static_cast<BasicSymbol *>(nodes[j].lexicon);
                if(symbol->getSymbol() == sequence[i])
                {
                    path.push_back(j);
                    foundNode = true;
                    break;
                }
            }
        assert(foundNode);
    }

    return path;
}
/*
unsigned int RDSGraph::predict(const SearchPath &searchPath) const
{
    unsigned int startPos = searchPath.size()-3;

    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    std::cout << printNode(searchPath[startPos]) << endl << endl;
    vector<ParseTree<unsigned int> > initialParseTrees;
    ParseTree<unsigned int> workingParseTree(searchPath[startPos], createParseNode(searchPath[startPos]));
    bottomUpSearchLeft(initialParseTrees, workingParseTree, searchPath[startPos]);
    for(unsigned int i = 0; i < initialParseTrees.size(); i++)
    {
        std::cout << "Initial Parse Tree " << i << endl;
        printParseTree(initialParseTrees[i].root, initialParseTrees[i], 0);
        std::cout << endl << endl << endl;
    }
    std::cout << "***************************************************************" << endl;

    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    std::cout << printNode(searchPath[startPos]) << endl << endl;
    vector<unsigned int> leftCorners;
    vector<ParseTree<unsigned int> > parseTrees;
    bottomUpSearchLeft(parseTrees, leftCorners, initialParseTrees);
    for(unsigned int i = 0; i < parseTrees.size(); i++)
    {
        ParseNode<unsigned int> corner = parseTrees[i].nodes[leftCorners[i]];
        std::cout << "Parse Tree " << i << endl;
        std::cout << "Tree Left Corner = [" << corner.value << " ----> " << corner.parent.first << " at " << corner.parent.second << "]" << endl;
        printParseTree(parseTrees[i].root, parseTrees[i], 0);
        std::cout << endl << endl << endl;
    }
    std::cout << "***************************************************************" << endl;


    for(int i = startPos-1; i >= 0; i--)
    {
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        std::cout << printNode(searchPath[i]) << endl;
        vector<ParseTree<unsigned int> > newParseTrees;
        workingParseTree = ParseTree<unsigned int>(searchPath[i], createParseNode(searchPath[i]));
        bottomUpSearchRight(newParseTrees, workingParseTree, searchPath[i], leftCorners, parseTrees);
        for(unsigned int j = 0; j < newParseTrees.size(); j++)
        {
            std::cout << "Parse Tree " << j << endl;
            printParseTree(newParseTrees[j].root, newParseTrees[j], 0);
            std::cout << endl << endl << endl;
        }

        parseTrees.clear();
        leftCorners.clear();
        bottomUpSearchLeft(parseTrees, leftCorners, newParseTrees);
        
        std::cout << "----------------------------------------------------------------" << endl;
        std::cout << printNode(searchPath[i]) << endl;
        for(unsigned int j = 0; j < parseTrees.size(); j++)
        {
            ParseNode<unsigned int> corner = parseTrees[j].nodes[leftCorners[j]];
            std::cout << "Parse Tree " << j << endl;
            std::cout << "Tree Left Corner = [" << corner.value << " ----> " << corner.parent.first << " at " << corner.parent.second << "]" << endl;
            printParseTree(parseTrees[j].root, parseTrees[j], 0);
            std::cout << endl << endl << endl;
        }
        std::cout << "***************************************************************" << endl;
    }

    return 0;
}
*/
void RDSGraph::bottomUpSearchLeft(vector<ParseTree<unsigned int> > &parseTrees, ParseTree<unsigned int> &workingParseTree, unsigned int unit) const
{
    // check pattern connections
    vector<Connection> parents = nodes[unit].parents;
    for(unsigned int i = 0; i < parents.size(); i++)
    {
        unsigned int parentPattern = parents[i].first;
        unsigned int parentPosition = parents[i].second;

        workingParseTree.pushOn(createParseNode(parentPattern), parents[i]);

        if(nodes[parentPattern].type == LexiconTypes::SP)
        {
            if(parentPosition > 0)
                parseTrees.push_back(workingParseTree);
            else
                bottomUpSearchLeft(parseTrees, workingParseTree, parentPattern);
        }
        else if(nodes[parentPattern].type == LexiconTypes::EC)
            bottomUpSearchLeft(parseTrees, workingParseTree, parentPattern);
        else
            assert(false);

        workingParseTree.popOff(parents[i]);
    }

    // check path connections
    for(unsigned int i = 0; i < nodes[unit].connections.size(); i++)
    {
        Connection con = nodes[unit].connections[i];
        con.first += nodes.size();
        workingParseTree.pushOn(createParseNode(con.first), con);
        parseTrees.push_back(workingParseTree);
        workingParseTree.popOff(con);
    }
}

void RDSGraph::bottomUpSearchRight(vector<ParseTree<unsigned int> > &parseTrees, ParseTree<unsigned int> &workingParseTree, unsigned int unit, const vector<unsigned int> &leftCorners, const vector<ParseTree<unsigned int> > &oldParseTrees) const
{
    assert(oldParseTrees.size() == leftCorners.size());
    assert(workingParseTree.nodes.size() > 0);

    vector<Connection> parentConnections = nodes[unit].parents;
    vector<Connection> pathConnections = nodes[unit].connections;

    // add path connections
    vector<Connection> testConnections = parentConnections;
    for(unsigned int i = 0; i < nodes[unit].connections.size(); i++)
        testConnections.push_back(Connection(pathConnections[i].first+nodes.size(), pathConnections[i].second));

    for(unsigned int i = 0; i < testConnections.size(); i++)
    {
        unsigned int testPattern = testConnections[i].first;
        unsigned int testPosition = testConnections[i].second;
        
        unsigned int oldWorkingRoot = workingParseTree.root;
        workingParseTree.pushOn(createParseNode(testPattern), testConnections[i]);

        bool isRootPattern = (testPattern >= nodes.size());
        if(isRootPattern || (nodes[testPattern].type == LexiconTypes::SP))  // check root patterns and other significant patterns
        {
            for(unsigned int j = 0; j < leftCorners.size(); j++)
            {
                const ParseNode<unsigned int> &leftCornerParseNode = oldParseTrees[j].getParseNode(leftCorners[j]);
                if((static_cast<int>(testPattern) == leftCornerParseNode.parent.first) && (testPosition == (leftCornerParseNode.parent.second-1)))
                {
                    parseTrees.push_back(oldParseTrees[j]);
                    parseTrees.back().attach(workingParseTree, oldWorkingRoot, testConnections[i]);
                }
            }

            if(!isRootPattern)   // no match is found and is not root pattern, so must be significant pattern
            {
                SignificantPattern *sp = static_cast<SignificantPattern *>(nodes[testPattern].lexicon);
                if(testPosition == (sp->size()-1))  // test position is right corner of test pattern then recursive call
                    bottomUpSearchRight(parseTrees, workingParseTree, testPattern, leftCorners, oldParseTrees);
            }
        }
        else if(nodes[testPattern].type == LexiconTypes::EC)
            bottomUpSearchRight(parseTrees, workingParseTree, testPattern, leftCorners, oldParseTrees);
        else
            assert(false);
            
        workingParseTree.popOff(testConnections[i]);
    }
}

void RDSGraph::bottomUpSearchLeft(vector<ParseTree<unsigned int> > &parseTrees, vector<unsigned int> &leftCorners, const vector<ParseTree<unsigned int> > &oldParseTrees) const
{
    for(unsigned int i = 0; i < oldParseTrees.size(); i++)
    {
        bool foundMatch = false;
        int matchedNode = oldParseTrees[i].bottomUpSearchLeft(oldParseTrees[i].leaves.front());

        if(matchedNode >= 0)
        {
            foundMatch = true;
            leftCorners.push_back(matchedNode);
            parseTrees.push_back(oldParseTrees[i]);
        }

        if(!foundMatch)
        {
            vector<unsigned int> newLeftCorners;
            vector<ParseTree<unsigned int> > newParseTrees;
            ParseTree<unsigned int> workingParseTree = oldParseTrees[i];
            bottomUpSearchLeft(newParseTrees, workingParseTree, workingParseTree.getRootNode().value);

            for(unsigned int j = 0; j < newParseTrees.size(); j++)
            {
                int newMatchedNode = newParseTrees[j].bottomUpSearchLeft(newParseTrees[j].leaves.front());
                if(newMatchedNode >= 0)
                {
                    leftCorners.push_back(newMatchedNode);
                    parseTrees.push_back(newParseTrees[j]);
                }
            }
        }
    }
}

ParseNode<unsigned int> RDSGraph::createParseNode(unsigned int theNode) const
{
    if(theNode < nodes.size())
    {
        if(nodes[theNode].type == LexiconTypes::SP)
        {
            SignificantPattern *sp = static_cast<SignificantPattern *>(nodes[theNode].lexicon);
            return ParseNode<unsigned int>(theNode, sp->size());
        }
        else if(nodes[theNode].type == LexiconTypes::EC)
            return ParseNode<unsigned int>(theNode, 1);
        else
            return ParseNode<unsigned int>(theNode, 0);
    }
    else
        return ParseNode<unsigned int>(theNode, paths[theNode-nodes.size()].size());
}

unsigned int RDSGraph::predict(const SearchPath &searchPath) const
{
    for(unsigned int i = 0; i < searchPath.size(); i++)
    {
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        std::cout << printNode(searchPath[i]) << endl << endl;
        vector<ParseTree<unsigned int> > initialParseTrees;
        ParseTree<unsigned int> workingParseTree(searchPath[i], createParseNode(searchPath[i]));
        initialiseParseTrees(initialParseTrees, workingParseTree, searchPath[i]);
        for(unsigned int j = 0; j < initialParseTrees.size(); j++)
        {
            std::cout << "Initial Parse Tree " << j << endl;
            printParseTree(initialParseTrees[j].root, initialParseTrees[j], 0);
            std::cout << endl << endl << endl;
        }
        std::cout << "***************************************************************" << endl;
    }

/*
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    std::cout << printNode(searchPath[startPos]) << endl << endl;
    vector<unsigned int> leftCorners;
    vector<ParseTree<unsigned int> > parseTrees;
    bottomUpSearchLeft(parseTrees, leftCorners, initialParseTrees);
    for(unsigned int i = 0; i < parseTrees.size(); i++)
    {
        ParseNode<unsigned int> corner = parseTrees[i].nodes[leftCorners[i]];
        std::cout << "Parse Tree " << i << endl;
        std::cout << "Tree Left Corner = [" << corner.value << " ----> " << corner.parent.first << " at " << corner.parent.second << "]" << endl;
        printParseTree(parseTrees[i].root, parseTrees[i], 0);
        std::cout << endl << endl << endl;
    }
    std::cout << "***************************************************************" << endl;


    for(int i = startPos-1; i >= 0; i--)
    {
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        std::cout << printNode(searchPath[i]) << endl;
        vector<ParseTree<unsigned int> > newParseTrees;
        workingParseTree = ParseTree<unsigned int>(searchPath[i], createParseNode(searchPath[i]));
        bottomUpSearchRight(newParseTrees, workingParseTree, searchPath[i], leftCorners, parseTrees);
        for(unsigned int j = 0; j < newParseTrees.size(); j++)
        {
            std::cout << "Parse Tree " << j << endl;
            printParseTree(newParseTrees[j].root, newParseTrees[j], 0);
            std::cout << endl << endl << endl;
        }

        parseTrees.clear();
        leftCorners.clear();
        bottomUpSearchLeft(parseTrees, leftCorners, newParseTrees);
        
        std::cout << "----------------------------------------------------------------" << endl;
        std::cout << printNode(searchPath[i]) << endl;
        for(unsigned int j = 0; j < parseTrees.size(); j++)
        {
            ParseNode<unsigned int> corner = parseTrees[j].nodes[leftCorners[j]];
            std::cout << "Parse Tree " << j << endl;
            std::cout << "Tree Left Corner = [" << corner.value << " ----> " << corner.parent.first << " at " << corner.parent.second << "]" << endl;
            printParseTree(parseTrees[j].root, parseTrees[j], 0);
            std::cout << endl << endl << endl;
        }
        std::cout << "***************************************************************" << endl;
    }
*/
    return 0;
}

void RDSGraph::initialiseParseTrees(std::vector<ParseTree<unsigned int> > &parseTrees, ParseTree<unsigned int> &workingParseTree, unsigned int unit) const
{     
    bool success = false;
    vector<Connection> parents = nodes[unit].parents;
    for(unsigned int i = 0; i < parents.size(); i++)
    {
        unsigned int parentPattern = parents[i].first;
        unsigned int parentPosition = parents[i].second;
        
        workingParseTree.pushOn(createParseNode(parentPattern), parents[i]);

        if(nodes[parentPattern].type == LexiconTypes::SP)
        {
            if(parentPosition == 0)    // since all SP are at least of length 2 
            {
                parseTrees.push_back(workingParseTree);
                success = true;
            }
        }
        else if(nodes[parentPattern].type == LexiconTypes::EC)   // if parent is EC, recursive call
        {
            initialiseParseTrees(parseTrees, workingParseTree, parentPattern);
            success = true;
        }
        else    // not possible to reach here
            assert(false);
            
        workingParseTree.popOff(parents[i]);
    }
    
    if(!success)
        parseTrees.push_back(workingParseTree);
}

string RDSGraph::printSignificantPattern(const SignificantPattern &sp) const
{
    ostringstream sout;

    for(unsigned int i = 0; i < sp.size(); i++)
    {
        unsigned tempIndex = sp[i];
        if(nodes[tempIndex].type == LexiconTypes::EC)
            sout << "E" << tempIndex;
        else if(nodes[tempIndex].type == LexiconTypes::SP)
            sout << "P" << tempIndex;
        else
            sout << *(nodes[tempIndex].lexicon);

        if(i < (sp.size() - 1)) sout << " - ";
    }

    return sout.str();
}

string RDSGraph::printEquivalenceClass(const EquivalenceClass &ec) const
{
    ostringstream sout;

    for(unsigned int i = 0; i < ec.size(); i++)
    {
        unsigned tempIndex = ec[i];
        if(nodes[tempIndex].type == LexiconTypes::EC)
            sout << "E" << tempIndex;
        else if(nodes[tempIndex].type == LexiconTypes::SP)
            sout << "P" << tempIndex;
        else
            sout << *(nodes[tempIndex].lexicon);

        if(i < (ec.size() - 1)) sout << " | ";
    }

    return sout.str();
}

string RDSGraph::printNode(unsigned int node) const
{
    ostringstream sout;

    if(nodes[node].type == LexiconTypes::EC)
    {
        EquivalenceClass *ec = static_cast<EquivalenceClass *>(nodes[node].lexicon);
        sout << "E[" << printEquivalenceClass(*ec) << "]";
    }
    else if(nodes[node].type == LexiconTypes::SP)
    {
        SignificantPattern *sp = static_cast<SignificantPattern *>(nodes[node].lexicon);
        sout << "P[" << printSignificantPattern(*sp) << "]";
    }
    else
        sout << *(nodes[node].lexicon);

    return sout.str();
}

string RDSGraph::printPath(const SearchPath &path) const
{
    ostringstream sout;

    sout << "[";
    for(unsigned int i = 0; i < path.size(); i++)
    {
        unsigned tempIndex = path[i];
        if(nodes[tempIndex].type == LexiconTypes::EC)
            sout << "E" << tempIndex;
        else if(nodes[tempIndex].type == LexiconTypes::SP)
            sout << "P" << tempIndex;
        else
            sout << *(nodes[tempIndex].lexicon);

        if(i < (path.size() - 1)) sout << " - ";
    }
    sout << "]";

    return sout.str();
}

void RDSGraph::printParseTree(unsigned int parseNode, const ParseTree<unsigned int> &parseTree, unsigned int tabLevel) const
{
    ParseNode<unsigned int> currentParseNode = (*parseTree.nodes.find(parseNode)).second;

    if(currentParseNode.children.size() == 0)
        std::cout << "*";

    for(unsigned int i = 0; i < tabLevel; i++)
        std::cout << "\t";

    if(currentParseNode.value < nodes.size())
        std::cout << "Pattern = [" << currentParseNode.value << " ----> " << currentParseNode.parent.first << " at " << currentParseNode.parent.second << "]" << endl;
    else
        std::cout << "Path = [" << currentParseNode.value << "(" << currentParseNode.value - nodes.size() << ") ----> " << currentParseNode.parent.first << " at " << currentParseNode.parent.second << "]" << endl;

    for(unsigned int i = 0; i < currentParseNode.children.size(); i++)
        if(currentParseNode.children[i] >= 0)
            printParseTree(currentParseNode.children[i], parseTree, tabLevel+1);
        else
        {
            for(unsigned int j = 0; j < tabLevel+1; j++)
                std::cout << "\t";
            std::cout << "{}" << endl;
        }
}

void printInfo(const ConnectionMatrix &connections, const NRMatrix<double> &flows, const NRMatrix<double> &descents)
{
    for(unsigned int i = 0; i < connections.numberOfRows(); i++)
    {
        for(unsigned int j = 0; j < connections.numberOfColumns(); j++)
            std::cout << connections(i, j).size() << "\t";
        std::cout << endl;
    }
    std::cout << endl << endl << endl;
    for(unsigned int i = 0; i < flows.numberOfRows(); i++)
    {
        for(unsigned int j = 0; j < flows.numberOfColumns(); j++)
            std::cout << flows(i, j) << "\t";
        std::cout << endl;
    }
    std::cout << endl << endl << endl;
    for(unsigned int i = 0; i < descents.numberOfRows(); i++)
    {
        for(unsigned int j = 0; j < descents.numberOfColumns(); j++)
            std::cout << descents(i, j) << "\t";
        std::cout << endl;
    }
    std::cout << endl << endl << endl;
}
