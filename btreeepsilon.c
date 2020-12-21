#include "bepsilontree.h"
#include "red_black3.h"
//https://en.wikipedia.org/wiki/B-tree#Informal_description

#define IS_LEAF(node) (node->leaf || node->n.t.children[0] == NULL)
#define UPDATE_VAL_INSERT(val,new) ((*(val)).valid = true)
#define UPDATE_VAL_MARK_DEL(val) ((*(val)).valid = false)


int time = 0;

struct node * create_tree() {
    return calloc(1,sizeof(struct node));
}



void print_tree(struct node * treen, int depth) {
    if(treen == NULL) {
        return;
    }
    if(treen->leaf) {
        struct leaf_node * leaf = &treen->n.l;
        printf("Leaf Depth: %d, Num Vals: %d Ptr %p Parent %p [",depth,leaf->num_vals,treen,leaf->parent);
        for(int i=0;i<leaf->num_vals;++i) {
            printf("%d",leaf->values[i].data);
            if(!leaf->values[i].valid) {
                printf("i");
            }
            printf(", ");
        }
        printf("]\n");
    } else {
        struct tree_node *  tree = &treen->n.t;
        printf("Inner Depth: %d, Num Vals: %d Num Messages: %d Ptr %p Parent %p [",depth,tree->num_vals,tree->num_buffer,treen,tree->parent);
        for(int i=0;i<tree->num_vals;++i) {
            printf("%d",tree->values[i].data);
            if(!tree->values[i].valid) {
                printf("i");
            }
            printf(", ");
        }
        printf("] Messages %p[",tree->root);
        
        for(int i=0;i<tree->num_buffer;++i) {
            printf("%p %d %p %d %d, ",&tree->buffer[i],tree->buffer[i].key.data, tree->buffer[i].parent, tree->buffer[i].time,tree->buffer[i].key.valid);
        }
        printf("]\n");
        for(int i=0;i<tree->num_vals+1;++i) {
            print_tree(tree->children[i],depth+1);
        }
    }
}

int find_less_loc(DATA_TYPE * values, DATA_TYPE val, int num_vals) {
    int start = 0;
    int end = num_vals;
    while(start < end) {
        int mid = (start+end) / 2;
        if(DATA_LESS_THAN(values[mid],val)) {
            start = mid+1;
        } else {
            end = mid;
        }
    }
    return start;
}

struct message * find_message(struct tree_node * tree, DATA_TYPE val) {
    struct message * n = tree->root;
    while(n!=RBTREE_NULL) {
        if(DATA_EQUAL(n->key,val)) {
            return n;
        } else if(DATA_LESS_THAN(n->key,val)) {
            n = n->right;
        } else {
            n = n->left;
        }
    }
    return NULL;
}

int find_pos_index(struct tree_node * tree, DATA_TYPE val) {
    int start = 0;
    int end = tree->num_vals;
    while(start < end) {
        int mid = (start+end) / 2;
        if(DATA_EQUAL(tree->values[mid],val)) {
            return mid;
        } else if(DATA_LESS_THAN(tree->values[mid],val)) {
            start = mid+1;
        } else {
            end = mid;
        }
    }
    return start;
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

DATA_TYPE * find_parent_pos(struct node * tree) {
    struct tree_node * p = &tree->n.t.parent->n.t;
    for(int i=0;i<=p->num_vals;++i) {
        if(p->children[i] == tree) {
            return &p->values[i];
        }
    }
    return NULL;
}

DATA_TYPE * find_pos_leaf(struct leaf_node * tree, DATA_TYPE val) {
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

bool optional_flush(struct node * treex,struct node ** root) {
    if(treex->leaf) {
        printf("Error!\n");
        exit(1);
    }
    struct tree_node * tree = &treex->n.t;
    if(tree->num_buffer == BEXCESS) {
        flush_node(treex,root);
        return true;
    }
    return false;
}


static void swap_messages(rbnode_t * a, rbnode_t * b, rbtree_t * tree) {
    if(tree->root == a) {
        tree->root = b;
    } else if(tree->root == b) {
        tree->root = a;
    }
    rbnode_t * aleft = a->left;
    rbnode_t * bleft = b->left;
    rbnode_t * aright = a->right;
    rbnode_t * bright = b->right;
    bool skip =  false;
	if(a->parent != RBTREE_NULL) {
		if(a->parent->right == a) {
            if(a->parent->left == b) {
                a->parent->left = a;
                skip = true;
            }
			a->parent->right = b;
		} else {
            if(a->parent->right == b) {
                a->parent->right = a;
                skip = true;
            }
			a->parent->left = b;
		}
	}
	if(!skip && b->parent != RBTREE_NULL) {
		if(b->parent->right == b) {
			b->parent->right = a;
		} else {
			b->parent->left = a;
		}
	}
	if(aleft != RBTREE_NULL && aright != a) {
		aleft->parent = b;
	}
	if(aright != RBTREE_NULL && aright != a) {
		aright->parent = b;
	}
	if(bleft != RBTREE_NULL && bright != b) {
		bleft->parent = a;
	}
	if(bright != RBTREE_NULL && bright != b) {
		bright->parent = a;
	}
	struct message copy = *a;
	*a = *b;
	*b = copy;
}

void remove_message(struct message * m, struct tree_node * tree) {
    int key = m->key.data;
    if(&tree->buffer[tree->num_buffer-1] != m) {
        swap_messages(m,&tree->buffer[tree->num_buffer-1],tree);
    }
    rbtree_delete(tree,key);
    tree->num_buffer--;
    if(tree->num_buffer == 0) {
        rbtree_init(tree);
    }
}

bool insert_val(struct node * treex, DATA_TYPE val, struct node * brought_child, struct node ** root, struct node * caller) {
    
    if(treex == NULL) {
        struct node * new_rootx = create_tree();
        struct tree_node * new_root = &new_rootx->n.t;
        new_root->values[0] = val;
        new_root->num_vals = 1;
        new_root->children[0] = caller;
        new_root->children[1] = brought_child;
        caller->n.l.parent = new_rootx;
        brought_child->n.l.parent = new_rootx;
        *root = new_rootx;
        rbtree_init(new_root);
        return true;
    }
    if(treex->leaf) {
        if(!val.valid) {
            printf("Uhoh %d %d\n",val.data,val.valid);
            exit(1);
        }    
        struct leaf_node * tree = &treex->n.l;
        DATA_TYPE i;
        int index = find_less_loc(tree->values,val,tree->num_vals);
        if(tree->num_vals < BRANCHING) {
            memmove(&tree->values[index+1],&tree->values[index],sizeof(DATA_TYPE)*(BRANCHING-index-1));
            tree->values[index] = val;
            tree->num_vals++;
            return false;
        } else {
            
            struct node * new_childx = create_tree();
            new_childx->leaf = true;
            struct leaf_node * new_child = &new_childx->n.l;
            new_childx->leaf = true;
            new_child->parent = tree->parent;
            if(index == BRANCHING/2) {
                memcpy(&new_child->values,&tree->values[index],sizeof(DATA_TYPE)*(BRANCHING/2));
            } else if(index > BRANCHING/2) {
                memcpy(&new_child->values,&tree->values[BRANCHING/2+1],sizeof(DATA_TYPE)*(BRANCHING/2-1));
                new_child->num_vals = BRANCHING/2-1;
                insert_val(new_childx,val,brought_child, root,NULL);
                val = tree->values[BRANCHING/2];
            } else {
                memcpy(&new_child->values,&tree->values[BRANCHING/2],sizeof(DATA_TYPE)*(BRANCHING/2));
                tree->num_vals = BRANCHING/2-1;
                i = tree->values[BRANCHING/2-1];
                insert_val(treex,val,brought_child, root,NULL);
                val = i;
            }
            tree->num_vals = BRANCHING/2;
            new_child->num_vals = BRANCHING/2;
            insert_val(tree->parent,val,new_childx,root,treex);
            return true;
        }
    } else {
        
        struct tree_node * tree = &treex->n.t;
        //Find if message already exists
        struct message * pos_m = find_message(tree,val);
        if(pos_m && DATA_EQUAL(pos_m->key,val)) {
            if(pos_m->type == ADD) {
                UPDATE_VAL_INSERT(&val,pos_m->key);
            } else if(pos_m->type == DELETE) {
                UPDATE_VAL_MARK_DEL(&val);
            }
            remove_message(pos_m,tree);
            return false;
        }
        DATA_TYPE i;
        int index = find_less_loc(tree->values,val,tree->num_vals);
        if(tree->num_vals < NODE_SIZE) {
            memmove(&tree->values[index+1],&tree->values[index],sizeof(DATA_TYPE)*(NODE_SIZE-index-1));
            memmove(&tree->children[index+2],&tree->children[index+1],sizeof(struct tree_node *)*(NODE_SIZE-index-1));
            tree->values[index] = val;
            tree->children[index+1] = brought_child;
            if(brought_child != NULL) {
                if(brought_child->leaf) {
                    brought_child->n.l.parent = treex;
                } else {
                    brought_child->n.t.parent = treex;
                }
            }
            tree->num_vals++;
            return optional_flush(treex,root);
        } else {
            struct message buffer[BEXCESS];
            memcpy(&buffer,tree->buffer,sizeof(buffer));
            struct node * new_childx = create_tree();
            struct tree_node * new_child = &new_childx->n.t;
            rbtree_init(new_child);
            new_child->parent = tree->parent;
            if(index == NODE_SIZE/2) {
                memcpy(&new_child->values,&tree->values[index],sizeof(DATA_TYPE)*(NODE_SIZE/2));
                memcpy(&new_child->children[1],&tree->children[index+1],sizeof(struct tree_node *)*(NODE_SIZE/2));
                new_child->children[0] = brought_child;
            } else if(index > NODE_SIZE/2) {
                memcpy(&new_child->values,&tree->values[NODE_SIZE/2+1],sizeof(DATA_TYPE)*(NODE_SIZE/2-1));
                memcpy(&new_child->children,&tree->children[NODE_SIZE/2+1],sizeof(struct tree_node *)*(NODE_SIZE/2));
                new_child->num_vals = NODE_SIZE/2-1;
                insert_val(new_childx,val,brought_child, root,NULL);
                val = tree->values[NODE_SIZE/2];
            } else {
                memcpy(&new_child->values,&tree->values[NODE_SIZE/2],sizeof(DATA_TYPE)*(NODE_SIZE/2));
                memcpy(&new_child->children,&tree->children[NODE_SIZE/2],sizeof(struct tree_node *)*(NODE_SIZE/2+1));
                tree->num_vals = NODE_SIZE/2-1;
                i = tree->values[NODE_SIZE/2-1];
                insert_val(treex,val,brought_child, root,NULL);
                val = i;
            }
            
            tree->num_vals = NODE_SIZE/2;
            new_child->num_vals = NODE_SIZE/2;
            for(int i=0;i<=NODE_SIZE/2;++i) {
                if(new_child->children[i] != NULL) {
                    if(new_child->children[i]->leaf) {
                        new_child->children[i]->n.l.parent = new_childx;
                    } else {
                        new_child->children[i]->n.t.parent = new_childx;
                    }
                }
            }
            int current_buffer = tree->num_buffer;
            tree->num_buffer = 0;
            rbtree_init(tree);
            for(int i=0;i<current_buffer;++i) {
                if(DATA_LESS_THAN(buffer[i].key,val)) {
                    insert_message(&buffer[i],treex,NULL);
                } else {
                    insert_message(&buffer[i],new_childx,NULL);
                }
            }
            if(treex == *root) {
                struct node * new_rootx = create_tree();
                struct tree_node * new_root = &new_rootx->n.t;
                new_root->num_vals = 1;
                new_root->children[0] = treex;
                new_root->children[1] = new_childx;
                tree->parent = new_rootx;
                new_child->parent = new_rootx;
                new_root->values[0] = val;
                *root = new_rootx;
                rbtree_init(new_root);
                return true;
            }
            insert_val(tree->parent,val,new_childx,root,treex);
            return true;
        }
    }

}

bool insert_message(struct message * m, struct node * nodex, struct node ** root) {
    if(nodex->leaf) {

        struct leaf_node * leaf = &nodex->n.l;
        if(m->type == DELETE) {
            return delete_value_leaf(m->key,nodex,root);
        }
        DATA_TYPE * v = find_pos_leaf(leaf,m->key);
        
        if((v-leaf->values) < BRANCHING && DATA_EQUAL(*v,m->key)) {
            UPDATE_VAL_INSERT(v,m->key);
            return false;
        }
        if(m->type == ADD) {
            return insert_val(nodex,m->key,NULL,root,NULL);
        }
        return false;
    }
    struct tree_node * node = &nodex->n.t;
    DATA_TYPE * v = find_pos(node,m->key);
    if((v - node->values) < NODE_SIZE && DATA_EQUAL(*v,m->key)) {
        if(m->type == ADD){
            UPDATE_VAL_INSERT(v,m->key);
        } else {
            UPDATE_VAL_MARK_DEL(v);
        }
        return false;
    }
    //Find if message already exists
    struct message * pos_m = find_message(node,m->key);
    if(pos_m && DATA_EQUAL(pos_m->key,m->key)) {
        pos_m->time = m->time;
        pos_m->type = m->type;
        pos_m->key = m->key;
        return false;
    }
    struct message * mes = &(node->buffer[node->num_buffer++]);
    mes->time = m->time;
    mes->type = m->type; 
    mes->key = m->key;
    rbtree_insert(node,mes);
    return optional_flush(nodex, root);
}

bool flush_node(struct node * nodex, struct node ** root) {
    struct tree_node * node = &nodex->n.t;
    bool done = false;
    while(!done && node->num_buffer > 0) {
        struct message * m = &node->buffer[--node->num_buffer];
        rbtree_delete(node,m->key.data);
        int correct_child = find_pos_index(node,m->key);
        done = insert_message(m,node->children[correct_child],root);
    }
}


void insert_value(int valx, struct node ** root) {
    DATA_TYPE val;
    val.data = valx;
    val.valid = true; //As inserting
    if((*root)->leaf) {
        int index = find_less_loc((*root)->n.l.values,val,(*root)->n.l.num_vals);
        if(DATA_EQUAL((*root)->n.l.values[index],val)) {
            UPDATE_VAL_INSERT(&(*root)->n.l.values[index],val);
        } else if((*root)->n.l.num_vals == BRANCHING) {
            insert_val(*root,val,NULL,root,NULL);
        } else {
            memmove(&(*root)->n.l.values[index+1],&(*root)->n.l.values[index],sizeof(DATA_TYPE)*(BRANCHING-index-1));
            (*root)->n.l.values[index] = val;
            (*root)->n.l.num_vals++;
        }
        return;
    }

    //Non leaf
    
    DATA_TYPE * v = find_pos(&(*root)->n.t,val);

    if((v - (&(*root)->n.t.values[0])) < NODE_SIZE && DATA_EQUAL(*v,val)) {
        UPDATE_VAL_INSERT(v,val);
        return;
    }
    //Find if message already exists
    struct message * pos_m = find_message(&(*root)->n.t,val);
    if(pos_m && DATA_EQUAL(pos_m->key,val)) {
        pos_m->time = time++;
        pos_m->type = ADD;
        pos_m->key.valid = true;
        return;
    }
    struct message * mes = &((*root)->n.t.buffer[(*root)->n.t.num_buffer++]);
    mes->time = time++;
    mes->type = ADD; 
    mes->key = val;
    rbtree_insert(&(*root)->n.t,mes);
    if(!mes->key.valid) {
        printf("Error\n");
        exit(1);
    }
    if((*root)->n.t.num_buffer == BEXCESS) {
        flush_node(*root, root);
    }
}

struct node * get_far_right(struct node * tree) {
    while(!IS_LEAF(tree)) {
        tree = tree->n.t.children[tree->n.t.num_vals];
    }
    return tree;
}

struct node * get_largest_left(struct node * tree) {
    if(!IS_LEAF(tree)) {
        tree = get_far_right(tree->n.t.children[0]);
    }
    return tree;
}

struct node * get_right_sibling(struct node * tree) {
    if(tree->leaf) {
        for(int i=0;i<tree->n.l.parent->n.t.num_vals;++i) {
            if(tree->n.l.parent->n.t.children[i] == tree) {
                return tree->n.l.parent->n.t.children[i+1];
            }
        }
        return NULL;
    }
    if(tree->n.t.parent == NULL) {
        return NULL;
    }
    for(int i=0;i<tree->n.t.parent->n.t.num_vals;++i) {
        if(tree->n.t.parent->n.t.children[i] == tree) {
            return tree->n.t.parent->n.t.children[i+1];
        }
    }
    return NULL;
}

struct node * get_left_sibling(struct node * tree) {
    if(tree->leaf) {
        for(int i=1;i<=tree->n.l.parent->n.t.num_vals;++i) {
            if(tree->n.l.parent->n.t.children[i] == tree) {
                return tree->n.l.parent->n.t.children[i-1];
            }
        }
        return NULL;
    }
    if(tree->n.t.parent == NULL) {
        return NULL;
    }
    for(int i=1;i<=tree->n.t.parent->n.t.num_vals;++i) {
        if(tree->n.t.parent->n.t.children[i] == tree) {
            return tree->n.t.parent->n.t.children[i-1];
        }
    }
    return NULL;
}

void check_messages(struct tree_node * tree, DATA_TYPE * pos) {
    struct message * pos_m = find_message(tree,*pos);
    if(pos_m && DATA_EQUAL(pos_m->key,*pos)) {
        if(pos_m->type == ADD) {
            UPDATE_VAL_INSERT(pos,pos_m->key);
        } else if(pos_m->type == DELETE) {
            UPDATE_VAL_MARK_DEL(pos);
        }
        remove_message(pos_m,tree);
        return;
    }
}

bool rebalance(struct node * treex, struct node ** root) {
    
    if(treex->leaf) {
        struct leaf_node * tree = &treex->n.l;
        if(tree->num_vals >= BRANCHING/2 || treex==(*root)) {
            return false;
        }
        struct node * left_siblingx = get_left_sibling(treex);
        struct leaf_node * left_sibling = left_siblingx == NULL ? NULL : &left_siblingx->n.l;
        if(left_sibling && left_sibling->num_vals > BRANCHING/2) {
            DATA_TYPE * pos = find_pos(&(tree->parent->n.t),left_sibling->values[0]);
            memmove(&tree->values[1],tree->values,(BRANCHING-1)*(sizeof(DATA_TYPE)));
            tree->values[0] = *pos;
            tree->num_vals++;
            *pos = left_sibling->values[left_sibling->num_vals-1];
            left_sibling->num_vals--;
            check_messages(&(tree->parent->n.t),pos);
            if(!tree->values[0].valid) {
                return delete_value_leaf(tree->values[0],treex,root);
            }
            return false;
        }
        struct node * right_siblingx = get_right_sibling(treex);
        struct leaf_node * right_sibling = right_siblingx == NULL ? NULL : &right_siblingx->n.l;
        if(right_sibling && right_sibling->num_vals > BRANCHING/2) {
            DATA_TYPE * pos = find_pos(&(tree->parent->n.t),right_sibling->values[0]) - 1;
            tree->values[tree->num_vals] = *pos;
            *pos = right_sibling->values[0];
            tree->num_vals++;
            memmove(right_sibling->values,&right_sibling->values[1],(BRANCHING-1)*(sizeof(DATA_TYPE)));
            right_sibling->num_vals--;
            check_messages(&(tree->parent->n.t),pos);
            if(!tree->values[tree->num_vals-1].valid) {
                return delete_value_leaf(tree->values[tree->num_vals-1],treex,root);
            }
            return false;
        }
        if(left_siblingx == NULL && right_siblingx == NULL) {
            return false;
        }
        
        DATA_TYPE * pos;
        if(left_siblingx == NULL) {
            left_siblingx = treex;
            left_sibling = tree;
            treex = right_siblingx;
            tree = right_sibling;
            pos = find_pos(&(tree->parent->n.t),right_sibling->values[0]) - 1;
        } else {
            pos = find_pos(&(tree->parent->n.t),left_sibling->values[0]);
        }
        if(left_sibling->parent == *root) {
            return false;
        }
        //Merge left_sibling and tree
        bool valid = pos->valid;

        if(valid) {
            left_sibling->values[left_sibling->num_vals] = *pos;
            memcpy(&left_sibling->values[left_sibling->num_vals+1],tree->values,sizeof(DATA_TYPE)*tree->num_vals);
            int num = (pos-left_sibling->parent->n.t.values);
            memmove(pos,pos+1,sizeof(DATA_TYPE)*(NODE_SIZE-num-1));
            memmove(&left_sibling->parent->n.t.children[num+1],&left_sibling->parent->n.t.children[num+2],sizeof(struct tree_node *)*(NODE_SIZE-num-1));
            left_sibling->parent->n.t.num_vals--;
            left_sibling->num_vals += (1 + tree->num_vals);
        } else {
            memcpy(&left_sibling->values[left_sibling->num_vals],tree->values,sizeof(DATA_TYPE)*tree->num_vals);
            int num = (pos-left_sibling->parent->n.t.values);
            memmove(pos,pos+1,sizeof(DATA_TYPE)*(NODE_SIZE-num-1)); 
            memmove(&left_sibling->parent->n.t.children[num+1],&left_sibling->parent->n.t.children[num+2],sizeof(struct tree_node *)*(NODE_SIZE-num-1));
            left_sibling->parent->n.t.num_vals--;
            left_sibling->num_vals += (tree->num_vals);
        }
        free(treex);
        rebalance(left_sibling->parent,root);
        return true;
    }
    //Non-leaf
    struct tree_node * tree = &treex->n.t;
    if(tree->num_vals >= NODE_SIZE/2 || treex==*root) {
        return false;
    }
    struct node * left_siblingx = get_left_sibling(treex);
    struct tree_node * left_sibling = left_siblingx == NULL ? NULL : &left_siblingx->n.t;
    if(left_sibling && left_sibling->num_vals > NODE_SIZE/2 && (left_sibling->num_buffer + tree->num_buffer < BEXCESS)) {
        DATA_TYPE * pos = find_pos(&(tree->parent->n.t),left_sibling->values[0]);
        memmove(&tree->values[1],tree->values,(NODE_SIZE-1)*(sizeof(DATA_TYPE)));
        memmove(&tree->children[1],tree->children,(NODE_SIZE)*(sizeof(struct node *)));
        tree->values[0] = *pos;
        tree->children[0] = left_sibling->children[left_sibling->num_vals];
        if(tree->children[0]) {
            if(IS_LEAF(tree->children[0])) {
                tree->children[0]->n.l.parent = treex;
            } else {
                tree->children[0]->n.t.parent = treex;
            }
        }
        tree->num_vals++;
        *pos = left_sibling->values[left_sibling->num_vals-1];
        left_sibling->num_vals--;
        for(int i=left_sibling->num_buffer-1;i>=0;--i) {
            if(DATA_GREATER_THAN(left_sibling->buffer[i].key,*pos)) {
                tree->buffer[tree->num_buffer++] = left_sibling->buffer[i];
                rbtree_insert(tree,&tree->buffer[tree->num_buffer-1]);
                remove_message(&left_sibling->buffer[i],left_sibling);
                
            }
        }
        check_messages(&(tree->parent->n.t),pos);
        return false;
    } else if(left_sibling && (left_sibling->num_vals > NODE_SIZE/2 || left_sibling->num_buffer + tree->num_buffer >= BEXCESS)) {
        left_siblingx = NULL;
        left_sibling = NULL;
    }
    struct node * right_siblingx = get_right_sibling(treex);
    struct tree_node * right_sibling = right_siblingx == NULL ? NULL : &right_siblingx->n.t;
    if(right_sibling && right_sibling->num_vals > NODE_SIZE/2 &&(right_sibling->num_buffer + tree->num_buffer <BEXCESS)) {
        DATA_TYPE * pos = find_pos(&(tree->parent->n.t),right_sibling->values[0]) - 1;
        tree->values[tree->num_vals] = *pos;
        tree->children[tree->num_vals+1] = right_sibling->children[0];
        if(right_sibling->children[0]) {
            if(IS_LEAF(right_sibling->children[0])) {
                right_sibling->children[0]->n.l.parent = treex;
            } else {
                right_sibling->children[0]->n.t.parent = treex;
            }
        }
        *pos = right_sibling->values[0];
        tree->num_vals++;
        memmove(right_sibling->values,&right_sibling->values[1],(NODE_SIZE-1)*(sizeof(DATA_TYPE)));
        memmove(right_sibling->children,&right_sibling->children[1],(NODE_SIZE)*(sizeof(struct tree_node *)));
        right_sibling->num_vals--;
        for(int i=right_sibling->num_buffer-1;i>=0;--i) {
            if(DATA_LESS_THAN(right_sibling->buffer[i].key,*pos)) {
                tree->buffer[tree->num_buffer++] = right_sibling->buffer[i];
                rbtree_insert(tree,&tree->buffer[tree->num_buffer-1]);
                remove_message(&right_sibling->buffer[i],right_sibling);
            }
        }
        check_messages(&(tree->parent->n.t),pos);
        return false;
    } else if(right_sibling && (right_sibling->num_vals > NODE_SIZE/2 || right_sibling->num_buffer + tree->num_buffer >=BEXCESS)) {
        right_siblingx = NULL;
        right_sibling = NULL;
    }
    if(left_siblingx == NULL && right_siblingx == NULL) {
        return false;
    }
    DATA_TYPE * pos;
    if(left_sibling == NULL) {
        left_sibling = tree;
        left_siblingx = treex;
        tree = right_sibling;
        treex = right_siblingx;
        pos = find_parent_pos(right_siblingx) - 1;
    } else {
        pos = find_parent_pos(left_siblingx);
    }
    if(left_sibling->parent == (*root)) {
        return false;
    }
    //Merge left_sibling and tree
    left_sibling->values[left_sibling->num_vals] = *pos;
    
    memcpy(&left_sibling->values[left_sibling->num_vals+1],tree->values,sizeof(DATA_TYPE)*tree->num_vals);
    memcpy(&left_sibling->children[left_sibling->num_vals+1],tree->children,sizeof(struct tree_node *)*(tree->num_vals+1));
    int num = (pos-left_sibling->parent->n.t.values);
    memmove(pos,pos+1,sizeof(DATA_TYPE)*(NODE_SIZE-num-1));
    memmove(&left_sibling->parent->n.t.children[num+1],&left_sibling->parent->n.t.children[num+2],sizeof(struct tree_node *)*(NODE_SIZE-num-1));
    
    left_sibling->parent->n.t.num_vals--;
    left_sibling->num_vals += (1 + tree->num_vals);
    free(treex);

    for(int i=0;i<=left_sibling->num_vals;++i) {
        if(left_sibling->children[i]) {
            if(IS_LEAF(left_sibling->children[i])) {
                left_sibling->children[i]->n.l.parent = left_siblingx;
            } else {
                left_sibling->children[i]->n.t.parent = left_siblingx;
            }
        }
    }
    rebalance(left_sibling->parent,root);
    return true;
}

bool delete_value_leaf(DATA_TYPE val, struct node * leafx, struct node ** root) {
    
    struct leaf_node * leaf = &leafx->n.l;
    int index = find_less_loc(leaf->values,val,leaf->num_vals);
    if(index >= BRANCHING || !DATA_EQUAL(leaf->values[index],val)) {
        //Value did not exist, so noop
        return false;
    }
    memmove(&leaf->values[index],&leaf->values[index+1],sizeof(DATA_TYPE)*(BRANCHING-index-1));
    leaf->num_vals--;
    return rebalance(leafx,root);
}

void delete_value(int valx, struct node ** root) {
    DATA_TYPE val;
    val.data = valx;
    val.valid = false; //As deleting
    if((*root)->leaf) {
        delete_value_leaf(val,*root,root);
        return;
    }
    //Non leaf
    
    DATA_TYPE * v = find_pos(&(*root)->n.t,val);
    if(DATA_EQUAL(*v,val)) {
        UPDATE_VAL_MARK_DEL(v);
        return;
    }
    //Find if message already exists
    struct message * pos_m = find_message(&(*root)->n.t,val);

    if(pos_m && DATA_EQUAL(pos_m->key,val)) {
        pos_m->time = time++;
        pos_m->type = DELETE;
        pos_m->key.valid = false;
        return;
    }
    
    struct message * mes = &((*root)->n.t.buffer[(*root)->n.t.num_buffer++]);
    mes->time = time++;
    mes->type = DELETE; 
    mes->key = val;
    rbtree_insert(&(*root)->n.t,mes);
    optional_flush(*root,root);
}

DATA_TYPE * find_node(struct node * treex, DATA_TYPE val) {
    if(treex == NULL) {
        return NULL;
    }
    if(!treex->leaf) {
        struct tree_node * tree = &treex->n.t;
        struct message * m = find_message(tree,val);
        if(m && DATA_EQUAL(m->key,val)){
            if(m->type == ADD) {
                return &m->key;
            }
            return NULL;
        }
        int start = 0;
        int end = tree->num_vals;
        while(start < end) {
            int mid = (start+end) / 2;
            if(DATA_EQUAL(tree->values[mid],val)) {
                if(tree->values[mid].valid){
                    return &tree->values[mid];
                }
                return NULL;
            } else if(DATA_LESS_THAN(tree->values[mid],val)) {
                start = mid+1;
            } else {
                end = mid;
            }
        }
        return find_node(tree->children[end], val);
    }
    struct leaf_node * leaf = &treex->n.l;
    int start = 0;
    int end = leaf->num_vals;
    while(start < end) {
        int mid = (start+end) / 2;
        if(DATA_EQUAL(leaf->values[mid],val)) {
            return &leaf->values[mid];
        } else if(DATA_LESS_THAN(leaf->values[mid],val)) {
            start = mid+1;
        } else {
            end = mid;
        }
    }
    return NULL;
}

void destroy_tree(struct node * tree) {
    if(tree == NULL) {
        return;
    }
    if(!tree->leaf) {
        for(int i=0;i<tree->n.t.num_vals+1;++i) {
            destroy_tree(tree->n.t.children[i]);
        }
    }
    free(tree);
}

int check_tree_valid(struct node * tree, int min, struct node * parent) {
    
    if(tree == NULL) {
        return -10000;
    }
    if(tree->leaf) {
        if(tree->n.l.parent != parent) {
            printf("Error 3 %p %p\n",parent,tree->n.l.parent);
            exit(1);
        }
        int start= min;
        if(tree->n.l.num_vals > BRANCHING || tree->n.l.num_vals < 0) {
            printf("Error 9\n");
            exit(1);
        }
        for(int i=0;i<tree->n.l.num_vals;++i) {
            if(tree->n.l.values[i].data <= start) {
                printf("Errorx %d %d %p %d\n",start,tree->n.l.values[i].data,tree,i);
                print_tree(tree->n.l.parent->n.t.parent,0);
                exit(1);
            }
            start = tree->n.l.values[i].data;
        }
        return start;
    }
    struct tree_node * t = &tree->n.t;
    if(t->parent != parent) {
        printf("Error 3 %p %p %d\n",parent,t->parent,t->values[0].data);
        exit(1);
    }
    int min2 = min;
    if(t->num_buffer < 0 || t->num_buffer > BEXCESS) {
        printf("Error 7\n");
        int * x = 0;
        *x = 5;
        exit(1);
    }
    if(t->num_vals < 0 || t->num_vals > NODE_SIZE) {
        printf("Error 10\n");
        int * x = 0;
        *x = 5;
        exit(1);
    }
    if(t->num_buffer == 0 && t->root != RBTREE_NULL) {
        printf("Error 46 %p\n",t);
        print_tree(tree,0);
        exit(1);
    }
    for(int i=0;i<t->num_buffer;++i) {
        if(t->buffer[i].key.data <= min) {
            printf("Error 4 %d %p %d\n",i,tree,min);
            int * x = 0;
            *x = 5;
            exit(1);
        }
        if(t->buffer[i].key.valid != t->buffer[i].type) {
            printf("Error5 %d %d %p\n",t->buffer[i].key.data,t->buffer[i].key.valid,&t->buffer[i].type);
            exit(1);
        }
        if(t->buffer[i].key.data > min2) {
            min2 = t->buffer[i].key.data;
        }
        if(t->buffer[i].parent == RBTREE_NULL && t->root != &t->buffer[i])
        {
            printf("Error 40\n");
            print_tree(tree,0);
            exit(1);
        }
        if(t->buffer[i].parent!=RBTREE_NULL) {
            if(t->buffer[i].parent->left != &t->buffer[i] && t->buffer[i].parent->right != &t->buffer[i]) {
                printf("Error 41 %p\n",tree);
                exit(1);
            }
        }
        if(t->buffer[i].left != RBTREE_NULL && t->buffer[i].left->parent != &t->buffer[i]) {
            printf("Error 42\n");
            exit(1);
        }
        if(t->buffer[i].left != RBTREE_NULL && t->buffer[i].left->key.data >= t->buffer[i].key.data) {
            printf("Error 45 %p\n",tree);
            print_tree(tree,0);
            exit(1);
        }
        if(t->buffer[i].right != RBTREE_NULL && t->buffer[i].right->parent != &t->buffer[i]) {
            printf("Error 43\n");
            exit(1);
        }
        if(t->buffer[i].right != RBTREE_NULL && ((t->buffer[i].right->key.data) <= (t->buffer[i].key.data))) {
            printf("Error 45\n");
            exit(1);
        }
    }
    min = check_tree_valid(t->children[0],min,tree);
    for(int i=0;i<t->num_vals;++i) {
        if(t->values[i].data <= min) {
            printf("Error 2 %d %d\n",min,t->values[i].data);
            exit(1);
        }
        min = check_tree_valid(t->children[i+1],t->values[i].data,tree);
    }
    if(min < min2) {
        min = min2;
    }
    return min;
}

int main(int argc, char * argv[]) {
    struct node * tree = create_tree();
    tree->leaf = true;
    srand(12345);
    for(int i=1;i<NUM_ITERS;++i) {
        insert_value(i,&tree);
    }
    for(int i=1;i<NUM_ITERS;++i) {
        DATA_TYPE temp;
        int xxx = (rand() % (NUM_ITERS-1)) +1;
        temp.data = xxx;
        DATA_TYPE * x = find_node(tree,temp);
        if(!x || x->data !=xxx) {
            printf("AAAAh\n");
            exit(1);
        } 
    }
    destroy_tree(tree);
}