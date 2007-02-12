#ifndef PARSE_TREE_H
#define PARSE_TREE_H

#include <vector>
#include "Definitions.h"

template <class T>
class ParseNode
{
    friend class ParseTree;

    public:
        ParseNode()
        {
            parent = Connection(0, 0);
        };

        ParseNode(const T &value, const std::vector<unsigned int> &children)
        {
            parent = Connection(0, 0);
            this->value = value;
            this->children = children;
        };

    private:
        T value;
        Connection parent;
        std::vector<unsigned int> children;
};

template <class T>
class ParseTree
{
    public:
        ParseTree()
        {
            nodes.push_back(ParseNode());   // nodes[0] is always the root
        };

    private:
        std::vector<unsigned int> leaves;
        std::vector<ParseNode<T> > nodes;
};

#endif
