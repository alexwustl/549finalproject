#include "shared.h"

#define BRANCHING NUM
#define NODE_SIZE (BRANCHING - 1)

struct data_wrapper {
    int data;
    char value[VALUE_BYTES];
};

#define DATA_TYPE struct data_wrapper

struct tree_node {
    struct tree_node * parent;
    struct tree_node * children[BRANCHING];
    DATA_TYPE values[NODE_SIZE];
    int num_vals;
};