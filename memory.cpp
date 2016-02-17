//#include <malloc.h>  // deprecated on Mac OSX ??
#include <stdlib.h>
#include <stdio.h>

#include "memory.h"
#include "tree.h"


/* extern  variables */
extern  StrucMesh2  mesh;
extern  StrucTree2  tree2;

static void initNearMemory( void );

void *nearAlloc( unsigned int n )
{
   void *p;

   p = malloc(n);
   if (p == NULL) perror("nearAlloc");

   return( p );
} /*nearAlloc*/


char  *    ppMemory;
PStrucFace2  *ptree2face;

static void initNearMemory( void )
{
   tree2.root = (PStrucNode2d)nearAlloc( S_StrucNode2d );
   tree2.maxFace = 500000;
   tree2.face = (StrucFace2**)nearAlloc( tree2.maxFace*sizeof(PStrucFace2) );
   ptree2face = tree2.face;
   tree2.maxVicinityFace = 200;
   tree2.vicinityFace = (StrucFace2**)nearAlloc( tree2.maxVicinityFace*sizeof(StrucFace2) );

   return;
}/*initNearMemory*/


void initMemory( void )
{
   initNearMemory();

   return;
}/*initMemory*/

void freeMemory(void) {
    int i;
    free(tree2.vicinityFace),  tree2.vicinityFace = NULL;
    free(mesh.nRTria),  mesh.nRTria = NULL;
    free(mesh.nRPoint),  mesh.nRPoint = NULL;
    free(mesh.bCut),  mesh.bCut = NULL;
    free(mesh.region),  mesh.region = NULL;
    free(mesh.nVert),  mesh.nVert = NULL;
    for (i=0; i<mesh.iRLine; i++)  free(mesh.boundVert[i]);
    free(mesh.boundVert),  mesh.boundVert = NULL;
    while (tree2.nFace > 0)  remFace(tree2.face[0]);
    free(tree2.root),  tree2.root = NULL;
}
