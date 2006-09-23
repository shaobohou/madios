#ifndef PARSE_TREE_H
#define PARSE_TREE_H

#include <map>

#include "Definitions.h"

typedef std::pair<int, unsigned int> ParseConnection;

template <class T>
class ParseNode
{
    public:
        T value;
        ParseConnection parent;
        std::vector<int> children;

        ParseNode()
        {
            parent.first = -1;
        };

        ParseNode(const T &value, unsigned int numberOfChildren)
        {
            parent.first = -1;
            this->value = value;
            for(unsigned int i = 0; i < numberOfChildren; i++)
                children.push_back(-1);
        };
};

template <class T>
class ParseTree
{
    public:
        unsigned int root;
        std::vector<unsigned int> leaves;
        std::map<unsigned int, ParseNode<T> > nodes;

        ParseTree()
        {
            root = 0;
        };

        ParseTree(unsigned int rootKey, const ParseNode<T> &rootNode)
        {
            root = rootKey;
            nodes.insert(std::pair<unsigned int, ParseNode<T> >(rootKey, rootNode));
        };

        void pushOn(const ParseNode<T> &newNode, const Connection &newConnection)
        {
            unsigned int newKey = newConnection.first;

            nodes[newKey] = newNode;
            nodes[root].parent = newConnection;
            nodes[newKey].children[newConnection.second] = root;
            root = newKey;

            updateLeaves(root);
        };

        void popOff(const Connection &newConnection)
        {
            assert(newConnection.first == root);

            unsigned int oldRoot = root;
            root = nodes[root].children[newConnection.second];
            nodes[root].parent.first = -1;

            ParseNode<T> tempNode = nodes[oldRoot];
            for(unsigned int i = 0; i < tempNode.children.size(); i++)
                if((tempNode.children[i] >= 0) && (tempNode.children[i] != static_cast<int>(root)))
                    eraseNode(tempNode.children[i]);
            nodes.erase(oldRoot);

            updateLeaves(root);
        };
        
        void attach(const ParseTree<unsigned int> &branch, unsigned int branchPoint, const Connection &branchConnection)
        {
             ParseNode<unsigned int> branchNode = branch.getParseNode(branchPoint);

             nodes[branchConnection.first].children[branchConnection.second] = branchPoint;
             nodes[branchPoint] = branch.getParseNode(branchPoint);
             nodes[branchPoint].parent = branchConnection;

             for(unsigned int i = 0; i < branchNode.children.size(); i++)
                 if(branchNode.children[i] >= 0)
                     attach(branch, branchNode.children[i], Connection(branchPoint, i));
                     
             updateLeaves(root);
        }

        int bottomUpSearchLeft(unsigned int testNode) const
        {
            const ParseNode<unsigned int> &tempNode = getParseNode(testNode);
            if(tempNode.parent.first >= 0)
            {
                if(tempNode.parent.second > 0)
                    return testNode;
                else
                    return bottomUpSearchLeft(tempNode.parent.first);
            }
            else
                return -1;
        };

        const ParseNode<T>& getRootNode() const
        {
            return getParseNode(root);
        };

        const ParseNode<T>& getParseNode(unsigned int key) const
        {
            return (*nodes.find(key)).second;
        };

    private:
        void eraseNode(unsigned int key)
        {
            ParseNode<T> tempNode = nodes[key];
            for(unsigned int i = 0; i < tempNode.children.size(); i++)
                if(tempNode.children[i] >= 0)
                    eraseNode(tempNode.children[i]);

            nodes.erase(key);
        };

        void updateLeaves(unsigned int testNode)
        {
            if(testNode == root)
                leaves.clear();

            ParseNode<T> tempNode = getParseNode(testNode);
            if(tempNode.children.size() == 0)
                leaves.push_back(testNode);
            else
                for(unsigned int i = 0; i < tempNode.children.size(); i++)
                    if(tempNode.children[i] >= 0)
                        updateLeaves(tempNode.children[i]);
        };
};

#endif
