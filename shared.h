#ifndef SHARED_H
#define SHARED_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

// #define NUM 9 //This must be odd

// #define VALUE_BYTES 64
// #define NUM_ITERS 500000

#define DATA_EQUAL(a,b) ((a).data == (b).data)
#define DATA_LESS_THAN(a,b) ((a).data < (b).data)
#define DATA_GREATER_THAN(a,b) ((a).data > (b).data)

#endif