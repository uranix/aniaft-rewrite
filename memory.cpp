//#include <malloc.h>  // deprecated on Mac OSX ??
#include <stdlib.h>
#include <stdio.h>

#include "memory.h"
#include "tree.h"


/* extern  variables */
extern  StrucMesh2  mesh;
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
    int i;
    free(tree.vicinityFace),  tree.vicinityFace = NULL;
    free(mesh.nRTria),  mesh.nRTria = NULL;
    free(mesh.nRPoint),  mesh.nRPoint = NULL;
    free(mesh.bCut),  mesh.bCut = NULL;
    free(mesh.region),  mesh.region = NULL;
    free(mesh.nVert),  mesh.nVert = NULL;
    for (i=0; i<mesh.iRLine; i++)  free(mesh.boundVert[i]);
    free(mesh.boundVert),  mesh.boundVert = NULL;
    while (tree.nFace > 0)  remFace(tree.face[0]);
    free(tree.root),  tree.root = NULL;
}
