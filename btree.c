#include "btree.h"

//https://en.wikipedia.org/wiki/B-tree#Informal_description

#define IS_LEAF(node) (node->children[0] == NULL)

struct tree_node * create_tree() {
    return calloc(1,sizeof(struct tree_node));
}

void print_tree(struct tree_node * tree, int depth) {
    if(tree == NULL) {
        return;
    }
    printf("Depth: %d, Num Vals: %d Ptr %p Parent %p [",depth,tree->num_vals,tree,tree->parent);
    for(int i=0;i<tree->num_vals;++i) {
        printf("%d, ",tree->values[i].data);
    }
    printf("]\n");
    for(int i=0;i<tree->num_vals+1;++i) {
        print_tree(tree->children[i],depth+1);
    }
    
    
}

int check_tree_valid(struct tree_node * tree, int min, struct tree_node * parent) {
    
    if(tree == NULL) {
        return -10000;
    }
    if(IS_LEAF(tree)) {
        if(tree->parent != parent) {
            printf("Error 3 %p %p\n",parent,tree->parent);
            exit(1);
        }
        int start= min;
        for(int i=0;i<tree->num_vals;++i) {
            if(tree->values[i].data <= start) {
                printf("Errorx %d\n",start);
                exit(1);
            }
            start = tree->values[i].data;
        }
        return start;
    }
    if(tree->parent != parent) {
        printf("Error 3 %p %p %d\n",parent,tree->parent,tree->values[0].data);
        exit(1);
    }
    min = check_tree_valid(tree->children[0],min,tree);
    for(int i=0;i<tree->num_vals;++i) {
        if(tree->values[i].data < min) {
            printf("Error 2 %d\n",min);
            exit(1);
        }
        min = check_tree_valid(tree->children[i+1],tree->values[i].data,tree);
    }
    return min;
}

int find_less_loc(struct tree_node * tree, DATA_TYPE val) {
    int start = 0;
    int end = tree->num_vals;
    while(start < end) {
        int mid = (start+end) / 2;
        if(DATA_LESS_THAN(tree->values[mid],val)) {
            start = mid+1;
        } else {
            end = mid;
        }
    }
    return start;
}

/* Find node that contains this value, or a possible node */
struct tree_node * find_node(struct tree_node * tree, DATA_TYPE val, bool possible_node) {
    if(tree == NULL) {
        return NULL;
    }
    int start = 0;
    int end = tree->num_vals;
    while(start < end) {
        int mid = (start+end) / 2;
        if(DATA_EQUAL(tree->values[mid],val)) {
            return tree;
        } else if(DATA_LESS_THAN(tree->values[mid],val)) {
            start = mid+1;
        } else {
            end = mid;
        }
    }
    if (possible_node && tree->children[end] == NULL) {
        return tree;
    }
    return find_node(tree->children[end], val, possible_node);
}

DATA_TYPE * find_nodemk2(struct tree_node * tree, DATA_TYPE val) {
    if(tree == NULL) {
        return NULL;
    }
    int start = 0;
    int end = tree->num_vals;
    while(start < end) {
        int mid = (start+end) / 2;
        if(DATA_EQUAL(tree->values[mid],val)) {
            return &tree->values[mid];
        } else if(DATA_LESS_THAN(tree->values[mid],val)) {
            start = mid+1;
        } else {
            end = mid;
        }
    }
    return find_nodemk2(tree->children[end], val);
}

DATA_TYPE * find_pos(struct tree_node * tree, DATA_TYPE val) {
    int start = 0;
    int end = tree->num_vals;
    while(start < end) {
        int mid = (start+end) / 2;
        if(DATA_EQUAL(tree->values[mid],val)) {
            return &tree->values[mid];
        } else if(DATA_LESS_THAN(tree->values[mid],val)) {
            start = mid+1;
        } else {
            end = mid;
        }
    }
    return &tree->values[start];
}

void insert_val(struct tree_node * tree, DATA_TYPE val, struct tree_node * brought_child, struct tree_node ** root) {
    DATA_TYPE i;
    int index = find_less_loc(tree,val);
    DATA_TYPE * pos = find_pos(tree,val);
    if((pos-tree->values) < NODE_SIZE && DATA_EQUAL(*pos,val)) {
        //Update val
        memcpy(&(*pos->value),&val.value,sizeof(val.value));
        return;
    }
    if(tree->num_vals < NODE_SIZE) {
        memmove(&tree->values[index+1],&tree->values[index],sizeof(DATA_TYPE)*(NODE_SIZE-index-1));
        memmove(&tree->children[index+2],&tree->children[index+1],sizeof(struct tree_node *)*(BRANCHING-index-2));
        tree->values[index] = val;
        tree->children[index+1] = brought_child;
        if(brought_child != NULL) {
            brought_child->parent = tree;
        }
        tree->num_vals++;
    } else {
        struct tree_node * new_child = create_tree();
        new_child->parent = tree->parent;
        if(index == NODE_SIZE/2) {
            memcpy(&new_child->values,&tree->values[index],sizeof(DATA_TYPE)*(NODE_SIZE/2));
            memcpy(&new_child->children[1],&tree->children[index+1],sizeof(struct tree_node *)*(NODE_SIZE/2));
            new_child->children[0] = brought_child;
        } else if(index > NODE_SIZE/2) {
            memcpy(&new_child->values,&tree->values[NODE_SIZE/2+1],sizeof(DATA_TYPE)*(NODE_SIZE/2-1));
            memcpy(&new_child->children,&tree->children[NODE_SIZE/2+1],sizeof(struct tree_node *)*(NODE_SIZE/2));
            new_child->num_vals = NODE_SIZE/2-1;
            insert_val(new_child,val,brought_child, NULL);
            val = tree->values[NODE_SIZE/2];
        } else {
            memcpy(&new_child->values,&tree->values[NODE_SIZE/2],sizeof(DATA_TYPE)*(NODE_SIZE/2));
            memcpy(&new_child->children,&tree->children[NODE_SIZE/2],sizeof(struct tree_node *)*(NODE_SIZE/2+1));
            tree->num_vals = NODE_SIZE/2-1;
            i = tree->values[NODE_SIZE/2-1];
            insert_val(tree,val,brought_child, NULL);
            val = i;
        }
        
        tree->num_vals = NODE_SIZE/2;
        new_child->num_vals = NODE_SIZE/2;
        for(int i=0;i<=NODE_SIZE/2;++i) {
            if(new_child->children[i] != NULL) {
                new_child->children[i]->parent = new_child;
            }
        }
        
        if(tree == *root) {
            struct tree_node * new_root = create_tree();
            new_root->num_vals = 1;
            new_root->children[0] = tree;
            new_root->children[1] = new_child;
            tree->parent = new_root;
            new_child->parent = new_root;
            new_root->values[0] = val;
            *root = new_root;
            return;
        }
        insert_val(tree->parent,val,new_child,root);
    }
}

void insert_value(int valx, struct tree_node ** root) {
    DATA_TYPE val;
    val.data = valx;
    struct tree_node * possible = find_node(*root, val, true);
    insert_val(possible, val, NULL, root);
}

struct tree_node * get_far_right(struct tree_node * tree) {
    while(!IS_LEAF(tree)) {
        tree = tree->children[tree->num_vals];
    }
    return tree;
}

struct tree_node * get_largest_left(struct tree_node * tree,DATA_TYPE * val) {
    int pos = val - tree->values;
    if(!IS_LEAF(tree)) {
        tree = get_far_right(tree->children[pos]);
    }
    return tree;
}

struct tree_node * get_right_sibling(struct tree_node * tree) {
    if(tree->parent == NULL) {
        return NULL;
    }
    for(int i=0;i<tree->parent->num_vals;++i) {
        if(tree->parent->children[i] == tree) {
            return tree->parent->children[i+1];
        }
    }
    return NULL;
}

struct tree_node * get_left_sibling(struct tree_node * tree) {
    if(tree->parent == NULL) {
        return NULL;
    }
    for(int i=1;i<=tree->parent->num_vals;++i) {
        if(tree->parent->children[i] == tree) {
                return tree->parent->children[i-1];
        }
    }
    return NULL;
}


void rebalance(struct tree_node * tree, struct tree_node ** root) {
    if(tree->num_vals >= NODE_SIZE/2 || tree == *root) {
        return;
    }
    struct tree_node * left_sibling = get_left_sibling(tree);
    if(left_sibling && left_sibling->num_vals > NODE_SIZE/2) {
        DATA_TYPE * pos = find_pos(tree->parent,left_sibling->values[0]);
        memmove(&tree->values[1],tree->values,(NODE_SIZE-1)*(sizeof(DATA_TYPE)));
        memmove(&tree->children[1],tree->children,(NODE_SIZE)*(sizeof(struct tree_node *)));
        tree->values[0] = *pos;
        tree->children[0] = left_sibling->children[left_sibling->num_vals];

        if(tree->children[0]) {
            tree->children[0]->parent = tree;
        }
        tree->num_vals++;

        *pos = left_sibling->values[left_sibling->num_vals-1];
        left_sibling->num_vals--;
        return;
    }
    struct tree_node * right_sibling = get_right_sibling(tree);
    if(right_sibling && right_sibling->num_vals > NODE_SIZE/2) {
        DATA_TYPE * pos = find_pos(tree->parent,right_sibling->values[0]) - 1;
        tree->values[NODE_SIZE/2-1] = *pos;
        tree->children[NODE_SIZE/2] = right_sibling->children[0];
        if(right_sibling->children[0]) {
            right_sibling->children[0]->parent = tree;
        }
        *pos = right_sibling->values[0];
        tree->num_vals++;
        memmove(right_sibling->values,&right_sibling->values[1],(NODE_SIZE-1)*(sizeof(DATA_TYPE)));
        memmove(right_sibling->children,&right_sibling->children[1],(NODE_SIZE)*(sizeof(struct tree_node *)));
        right_sibling->num_vals--;
        return;
    }

    DATA_TYPE * pos;
    if(left_sibling == NULL) {
        left_sibling = tree;
        tree = right_sibling;
        pos = find_pos(tree->parent,right_sibling->values[0]) - 1;
    } else {
        pos = find_pos(tree->parent,left_sibling->values[0]);
    }
    //Merge left_sibling and tree
    left_sibling->values[left_sibling->num_vals] = *pos;
    
    memcpy(&left_sibling->values[left_sibling->num_vals+1],tree->values,sizeof(DATA_TYPE)*tree->num_vals);
    memcpy(&left_sibling->children[left_sibling->num_vals+1],tree->children,sizeof(struct tree_node *)*(tree->num_vals+1));
    int num = (pos-left_sibling->parent->values);
    memmove(pos,pos+1,sizeof(DATA_TYPE)*(NODE_SIZE-num-1));
    memmove(&left_sibling->parent->children[num+1],&left_sibling->parent->children[num+2],sizeof(struct tree_node *)*(NODE_SIZE-num-1));
    
    left_sibling->parent->num_vals--;
    left_sibling->num_vals += (1 + tree->num_vals);
    free(tree);

    for(int i=0;i<=left_sibling->num_vals;++i) {
        if(left_sibling->children[i]) {
            left_sibling->children[i]->parent = left_sibling;
        }
    }

    if(left_sibling->parent == *root && (*root)->num_vals == 0) {
        free(*root);
        *root = left_sibling;
        left_sibling->parent = NULL;
        return;
    }

    rebalance(left_sibling->parent,root);
}

void delete_value(int valx, struct tree_node ** root) {
    DATA_TYPE val;
    val.data = valx;
    struct tree_node * tree = find_node(*root,val,false);
    if(tree == NULL) {
        return;
    }
    DATA_TYPE * pos = find_pos(tree,val);
    struct tree_node * starter_node = get_largest_left(tree,pos);
    
    if(starter_node != tree) {
        *pos = starter_node->values[--starter_node->num_vals];        
    } else {
        memmove(pos,pos+1,sizeof(DATA_TYPE)*(NODE_SIZE-(pos-starter_node->values)-1));
        --starter_node->num_vals;
    }
    rebalance(starter_node,root);
}

void destroy_tree(struct tree_node * tree) {
    if(tree == NULL) {
        return;
    }
    for(int i=0;i<tree->num_vals+1;++i) {
        destroy_tree(tree->children[i]);
    }
    free(tree);
}



int main(int argc, char * argv[]) {
    srand(12345);
    struct tree_node * tree = create_tree();
    for(int i=0;i<NUM_ITERS;++i) {
        insert_value(i,&tree);
    }
    for(int i=1;i<NUM_ITERS;++i) {
        DATA_TYPE temp;
        int xxx=(rand() % (NUM_ITERS-1)) +1;
        temp.data = xxx;
        DATA_TYPE * x = find_nodemk2(tree,temp);
        if(!x || x->data !=xxx) {
            printf("AAAAxh\n");
            printf("%d\n",temp.data);
            exit(1);
        } 
    }
    destroy_tree(tree);
}

