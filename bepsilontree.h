

#ifndef BEPSILONTREE_H
#define BEPSILONTREE_H

// #include <math.h>
#include <stdbool.h>
#include "shared.h"

#define BRANCHING (NUM-1)
#define EPSILON 0.5
// #define BEPSILON 3 //((int) floor(pow(BRANCHING,EPSILON)))
#define BEXCESS (BRANCHING-BEPSILON)

#define NODE_SIZE (BEPSILON-1)

#define timestamp int


struct data_wrapper {
    int data;
    char value[VALUE_BYTES];
    bool valid;
};

#define DATA_TYPE struct data_wrapper


enum message_type {DELETE = 0, ADD = 1};

struct message {
    DATA_TYPE key;
    timestamp time;
    enum message_type type;    
    struct message * parent;
    struct message * left;
    struct message * right;
    int color;
};

struct tree_node {
    struct node * parent;
    struct message buffer[BEXCESS];
    struct message * root;
    struct node * children[NODE_SIZE+1];
    DATA_TYPE values[NODE_SIZE];
    int num_buffer;
    int num_vals;
};

struct leaf_node {
    struct node * parent;
    DATA_TYPE values[BRANCHING];
    int num_vals;
};

struct node {
    bool leaf;
    union {
        struct tree_node t;
        struct leaf_node l;
    } n;
};

void check_avl_valid(struct message * root, int height);
bool flush_node(struct node * nodex, struct node ** root);
bool insert_message(struct message * m, struct node * nodex, struct node ** root);
bool delete_value_leaf(DATA_TYPE val, struct node * leafx, struct node ** root);
int check_tree_valid(struct node * tree, int min, struct node * parent);


#define MESSAGE_LESS_THAN(a,b) (a->key < b->key || (a->key == b->key && a->time < b->time))
#define MESSAGE_GREATER_THAN(a,b) (a->key > b->key || (a->key == b->key && a->time > b->time))
#define MESSAGE_EQUAL(a,b) (a->key == b->key && a->time == b->time)

#endif
