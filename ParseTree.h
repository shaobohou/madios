#ifndef PARSE_TREE_H
#define PARSE_TREE_H

#include <vector>
#include "Definitions.h"

template <class T>
class ParseTree;

template <class T>
class ParseNode
{
    friend class ParseTree<T>;

    public:
        ParseNode()
        {
            parent = Connection(0, 0);
        }

        ParseNode(const T &value, const Connection &parent)
        {
            this->value = value;
            this->parent = parent;
        }

        const T& getValue() const
        {
            return value;
        }

        const std::vector<unsigned int>& getChildren() const
        {
            return children;
        }

        std::vector<unsigned int> rewireChildren(unsigned int start, unsigned int finish, unsigned new_node)
        {
            std::vector<unsigned int> subsumed_part(children.begin()+start, children.begin()+finish+1);
            children.erase( children.begin()+start, children.begin()+finish+1);
            children.insert(children.begin()+start, new_node);
            return subsumed_part;

        }

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
            nodes.push_back(ParseNode<T>());   // nodes[0] is always the root
        }

        ParseTree(const std::vector<T> &values)
        {
            nodes.push_back(ParseNode<T>());   // nodes[0] is always the root
            for(unsigned int i = 0; i < values.size(); i++)
            {
                nodes.front().children.push_back(nodes.size());
                nodes.push_back(ParseNode<T>(values[i], Connection(0, i)));
            }
        }

        const std::vector<ParseNode<T> >& getNodes() const
        {
            return nodes;
        }

        void rewire(unsigned int start, unsigned int finish, const T &new_node)
        {
            nodes.push_back(ParseNode<T>(new_node, Connection(0, 0)));
            nodes.back().children = nodes.front().rewireChildren(start, finish, nodes.size()-1);
            for(unsigned int i = 0; i < nodes.back().children.size(); i++)
                nodes[nodes.back().children[i]].parent = Connection(nodes.size()-1, i);
        }

        void print(unsigned int node, unsigned int tab_level) const
        {
            for(unsigned int i = 0; i < tab_level; i++)
                std::cout << "\t";
            std::cout << node << " ---> " << nodes[node].value << std::endl;
            for(unsigned int i = 0; i < nodes[node].children.size(); i++)
                print(nodes[node].children[i], tab_level+1);
        }

    private:
        std::vector<unsigned int> leaves;
        std::vector<ParseNode<T> > nodes;
};

#endif
