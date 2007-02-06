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
    // look possible significant pattern found with help of equivalence class
    ConnectionMatrix connections;
    NRMatrix<double> flows, descents;
    computeConnectionMatrix(connections, searchPath);
    computeDescentsMatrix(flows, descents, connections);

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
    // Significant Pattern Info
    bool patternFound = false;
    BootstrapInfo bestBootstrapInfo;
    SignificantPatternInfo bestPatternInfo;
    EquivalenceClassInfo bestECInfo;

    SearchPathInfo searchPathInfo;
    searchPathInfo.searchPath = searchPath;
    searchPathInfo.alreadyTested = false;
    for(unsigned int i = 0; (i+params.contextSize-1) < searchPath.size(); i++)
    {
        BootstrapInfo bootstrapInfo;
        bootstrapInfo.context = Range(i, i+params.contextSize-1);
        SearchPath boostedPath_part = bootstrap(bootstrapInfo, searchPath(i, i+params.contextSize-1), params.overlapThreshold);
        SearchPath boostedPath = searchPath.substitute(i, i+params.contextSize-1, boostedPath_part);

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
    bool patternFound = false;
    vector<SearchPath> previousPaths;
    for(unsigned int i = context.first + 1; i < context.second; i++)
    {
        // look for possible equivalence class
        EquivalenceClass ec;
        ConnectionMatrix connections;
        SearchPath generalPath = computeGeneralisedSubpaths(ec, connections, searchPath, context.first, i, context.second);

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
    // look possible significant pattern found with help of equivalence class
    vector<Range> patterns;
    vector<SignificancePair> pvalues;
    NRMatrix<double> flows, descents;
    computeDescentsMatrix(flows, descents, connections);
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

void RDSGraph::computeConnectionMatrix(ConnectionMatrix &connections, const SearchPath &searchPath) const
{
    // calculate subpath distributions, symmetrical matrix
    unsigned dim = searchPath.size();
    connections = ConnectionMatrix(dim, dim);
    for(unsigned int i = 0; i < dim; i++)
    {
        connections(i, i) = getAllNodeConnections(searchPath[i]);

        // compute the column from the diagonal
        for(unsigned int j = i + 1; j < dim; j++)
        {
            connections(j, i) = filterConnections(connections(j - 1, i), j-i, searchPath(j, j));
            connections(i, j) = connections(j, i);
        }
    }
}

SearchPath RDSGraph::computeGeneralisedSubpaths(EquivalenceClass &ec, ConnectionMatrix &connections, const SearchPath &searchPath, unsigned int prefixStart, unsigned int slotIndex, unsigned int postfixEnd)
{
    assert(prefixStart < slotIndex);
    assert(slotIndex < postfixEnd);

    // get the candidate connections
    vector<Connection> equivalenceConnections = getAllNodeConnections(searchPath[prefixStart]);
    equivalenceConnections = filterConnections(equivalenceConnections, 0,                       searchPath(prefixStart, slotIndex - 1));
    equivalenceConnections = filterConnections(equivalenceConnections, slotIndex+1-prefixStart, searchPath(slotIndex + 1, postfixEnd));

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
            computeConnectionMatrix(connections, generalPath);
            nodes.pop_back();
            updateAllConnections();
        }
        else
            computeConnectionMatrix(connections, generalPath);
    }
    else
        computeConnectionMatrix(connections, generalPath);

    return generalPath;
}

SearchPath RDSGraph::bootstrap(BootstrapInfo &bootstrapInfo, const SearchPath &searchPath, double overlapThreshold) const
{
    // find all possible connections
    vector<Connection> equivalenceConnections = getAllNodeConnections(searchPath[0]);
    equivalenceConnections = filterConnections(equivalenceConnections, searchPath.size()-1, searchPath(searchPath.size()-1, searchPath.size()-1));

    // find potential ECs
    bootstrapInfo.encounteredECs.clear();
    for(unsigned int i = 1; i < searchPath.size()-1; i++)
    {
        bootstrapInfo.encounteredECs.push_back(EquivalenceClass());
        for(unsigned int j = 0; j < equivalenceConnections.size(); j++)
        {
            unsigned int currentPath = equivalenceConnections[j].first;
            unsigned int currentStart = equivalenceConnections[j].second;

            bootstrapInfo.encounteredECs.back().add(paths[currentPath][currentStart+i]);
        }
    }

    // init bootstrap data
    bootstrapInfo.overlapECs = searchPath(1, searchPath.size()-2);
    bootstrapInfo.overlapRatios = vector<double>(searchPath.size()-2, 0.0);

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
        bootstrapPath[i + 1] = bootstrapInfo.overlapECs[i];
    }

    return bootstrapPath;
}

void RDSGraph::computeDescentsMatrix(NRMatrix<double> &flows, NRMatrix<double> &descents, const ConnectionMatrix &connections) const
{
    // calculate P_R and P_L
    unsigned dim = connections.numRows();
    flows = NRMatrix<double>(dim, dim, -1.0);
    for(unsigned int i = 0; i < dim; i++)
        for(unsigned int j = 0; j < dim; j++)
            if(i > j)
                flows(i, j) = static_cast<double>(connections(i, j).size()) / connections(i-1, j).size();
            else if(i < j)
                flows(i, j) = static_cast<double>(connections(i, j).size()) / connections(i+1, j).size();
            else
                flows(i, j) = static_cast<double>(connections(i, j).size()) / corpusSize;

    // calculate D_R and D_L
    descents = NRMatrix<double>(dim, dim, -1.0);
    for(unsigned int i = 0; i < dim; i++)
        for(unsigned int j = 0; j < dim; j++)
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
    unsigned int pathLength = descents.numRows();
    vector<unsigned int> candidateEndRows;
    vector<unsigned int> candidateStartRows;
    for(unsigned int i = 0; i < descents.numRows(); i++)
    {
        for(int j = i - 1; j >= 0; j--)
            if(descents(i, j) < eta)
            {
                candidateEndRows.push_back(i - 1);
                break;
            }

        for(unsigned int j = i + 1; j < descents.numCols(); j++)
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
    for(unsigned int i = pattern.second; i < connections.numCols(); i++)
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

vector<Connection> RDSGraph::filterConnections(const vector<Connection> &init_cons, unsigned int start_offset, const SearchPath &search_path) const
{
    vector<Connection> filtered_cons;
    for(unsigned int i = 0; i < init_cons.size(); i++)
    {
        unsigned int cur_path = init_cons[i].first;
        unsigned int cur_pos = init_cons[i].second;

        // discard current connection because the path is not long enough to match the search path (segment)
        if((cur_pos+start_offset+search_path.size()) > paths[cur_path].size())
            continue;

        unsigned int count = search_path.size();
        for(unsigned int j = 0; j < search_path.size(); j++)
        {
            unsigned int actual_pos = j+cur_pos+start_offset;
            if(nodes[search_path[j]].type == LexiconTypes::EC)
            {   // if node on search path is EC and it contains the node and temp path
                if(!(static_cast<EquivalenceClass *>(nodes[search_path[j]].lexicon)->has(paths[cur_path][actual_pos])))
                    break;
            }
            else// else just test if they are the same node (BasicSymbol)
                if(search_path[j] != paths[cur_path][actual_pos])
                    break;

            count--;
        }

        // 0 if search_path completely matches temp_path
        if(count == 0)
            filtered_cons.push_back(init_cons[i]);
    }

    return filtered_cons;
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

void printInfo(const ConnectionMatrix &connections, const NRMatrix<double> &flows, const NRMatrix<double> &descents)
{
    for(unsigned int i = 0; i < connections.numRows(); i++)
    {
        for(unsigned int j = 0; j < connections.numCols(); j++)
            std::cout << connections(i, j).size() << "\t";
        std::cout << endl;
    }
    std::cout << endl << endl << endl;
    for(unsigned int i = 0; i < flows.numRows(); i++)
    {
        for(unsigned int j = 0; j < flows.numCols(); j++)
            std::cout << flows(i, j) << "\t";
        std::cout << endl;
    }
    std::cout << endl << endl << endl;
    for(unsigned int i = 0; i < descents.numRows(); i++)
    {
        for(unsigned int j = 0; j < descents.numCols(); j++)
            std::cout << descents(i, j) << "\t";
        std::cout << endl;
    }
    std::cout << endl << endl << endl;
}
