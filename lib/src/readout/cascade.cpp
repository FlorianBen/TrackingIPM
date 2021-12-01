#include <vector>

struct Node {
    Node* parent;
    std::vector<Node*> childs;

    Node();

    void createChild();
};