//#include <malloc.h>  // deprecated on Mac OSX ??
#include <stdlib.h>
#include <stdio.h>

#include "memory.h"
#include "tree.h"


/* extern  variables */
extern  StrucTree2  tree;

static void initNearMemory( void );

void *nearAlloc( unsigned int n )
{
   void *p;

   p = malloc(n);
   if (p == NULL) perror("nearAlloc");

   return( p );
} /*nearAlloc*/


char  *    ppMemory;
PStrucFace2  *ptreeface;

static void initNearMemory( void )
{
   tree.root = (PStrucNode2d)nearAlloc( S_StrucNode2d );
   tree.maxFace = 500000;
   tree.face = (StrucFace2**)nearAlloc( tree.maxFace*sizeof(PStrucFace2) );
   ptreeface = tree.face;
   tree.maxVicinityFace = 200;
   tree.vicinityFace = (StrucFace2**)nearAlloc( tree.maxVicinityFace*sizeof(StrucFace2) );

   return;
}/*initNearMemory*/


void initMemory( void )
{
   initNearMemory();

   return;
}/*initMemory*/

void freeMemory(void) {
    free(tree.vicinityFace),  tree.vicinityFace = NULL;
    while (tree.nFace > 0)  remFace(tree.face[0]);
    free(tree.root),  tree.root = NULL;
}
