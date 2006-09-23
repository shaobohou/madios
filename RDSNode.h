#ifndef RDSNODE_H
#define RDSNODE_H

#include "LexiconUnit.h"
#include "Definitions.h"
#include "ADIOSUtils.h"

class RDSNode
{
    public:
        LexiconUnit *lexicon;
        LexiconTypes::LexiconEnum type;
        std::vector<Connection> connections;
        std::vector<Connection> parents;

        RDSNode();
        explicit RDSNode(LexiconUnit *lexicon, LexiconTypes::LexiconEnum type);
        RDSNode(const RDSNode &other);
        ~RDSNode();

        RDSNode& operator=(const RDSNode &other);

        void addConnection(const Connection &con);
        const std::vector<Connection>& getConnections() const;
        void setConnections(const std::vector<Connection> &connections);

        bool addParent(const Connection &newParent);

    private:
        void deepCopy(const RDSNode &other);
};

#endif
