/*
 * Main developer: Nico Van Cleemput
 *
 * Copyright (C) 2015 Nico Van Cleemput.
 * Licensed under the GNU AFFERO GPL, read the file LICENSE for details.
 */


#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include "bitset.h"
#include "boolean.h"


#define MAXN 34            /* the maximum number of vertices */
#define MAXE (6*MAXN-12)    /* the maximum number of oriented edges */
#define MAXF (2*MAXN-4)      /* the maximum number of faces */
#define MAXVAL (MAXN-1)  /* the maximum degree of a vertex */
#define MAXCODELENGTH (MAXN+MAXE+3)

#define INFI (MAXN + 1)

typedef struct e /* The data type used for edges */ {
    int start; /* vertex where the edge starts */
    int end; /* vertex where the edge ends */
    int rightface; /* face on the right side of the edge
                          note: only valid if make_dual() called */
    struct e *prev; /* previous edge in clockwise direction */
    struct e *next; /* next edge in clockwise direction */
    struct e *inverse; /* the edge that is inverse to this one */
    int mark, index; /* two ints for temporary use;
                          Only access mark via the MARK macros. */
    
    bitset incident_faces;
} EDGE;

EDGE *firstedge[MAXN]; /* pointer to arbitrary edge out of vertex i. */
int degree[MAXN];

EDGE *facestart[MAXF]; /* pointer to arbitrary edge of face i. */
int faceSize[MAXF]; /* pointer to arbitrary edge of face i. */

bitset neighbours[MAXN];

bitset verticesInFace[MAXF];

EDGE edges[MAXE];

static int markvalue = 30000;
#define RESETMARKS {int mki; if ((markvalue += 2) > 30000) \
       { markvalue = 2; for (mki=0;mki<MAXE;++mki) edges[mki].mark=0;}}
#define MARK(e) (e)->mark = markvalue
#define MARKLO(e) (e)->mark = markvalue
#define MARKHI(e) (e)->mark = markvalue+1
#define UNMARK(e) (e)->mark = markvalue-1
#define ISMARKED(e) ((e)->mark >= markvalue)
#define ISMARKEDLO(e) ((e)->mark == markvalue)
#define ISMARKEDHI(e) ((e)->mark > markvalue)

int nv;
int ne;
int nf;

//data for the cycle we're building
bitset currentCycle;
int firstVertexCycle;
EDGE *firstEdgeCycle;

/* Returns a bitset containing the indices of all faces contained
 * between from and to in a clockwise directions.
 */
bitset facesBetween(EDGE *from, EDGE *to){
    EDGE *e, *elast;
    bitset faces = EMPTY_SET;
    
    e = elast = from;
    
    while(e != to){
        ADD(faces, e->rightface);
        e = e->next;
    }
    
    return faces;
}

boolean finishCycle(EDGE *newEdge, bitset saturatedFaces, bitset facesRight,
        bitset facesLeft, bitset emptyFaces){
    int i;
    
    if(IS_NOT_EMPTY(INTERSECTION(facesRight, facesLeft))){
        //a face cannot be on both sides of the cycles
        return FALSE;
    }
    
    ADD_ALL(saturatedFaces, newEdge->incident_faces);
    ADD_ALL(facesRight, facesBetween(firstEdgeCycle, newEdge->inverse));
    ADD_ALL(facesLeft, facesBetween(newEdge->inverse, firstEdgeCycle));
    
    for(i = 0; i < nf; i++){
        if(!CONTAINS(saturatedFaces, i)){
            //we found an empty face
            ADD(emptyFaces, i);
        }
    }
    if(IS_NOT_EMPTY(INTERSECTION(emptyFaces, facesLeft)) && 
            IS_NOT_EMPTY(INTERSECTION(emptyFaces, facesRight))){
        //there is an empty face on both sides of the cycle
        return FALSE;
    }
    
    return TRUE;
}

boolean continueCycle(EDGE *newEdge, int remainingVertices,
        bitset saturatedFaces, bitset facesRight, bitset facesLeft,
        bitset emptyFaces){
    //the end point of newEdge has not been added to the cycle yet
    //the faces between newEdge and the last edge have been added to left and right
    
    EDGE *e, *elast;
    int i;
    
    if(IS_NOT_EMPTY(INTERSECTION(facesRight, facesLeft))){
        //a face cannot be on both sides of the cycles
        return FALSE;
    }
    
    bitset lastVertexSingleton = SINGLETON(newEdge->start);
    for(i = 0; i < nf; i++){
        if(IS_NOT_EMPTY(INTERSECTION(verticesInFace[i], lastVertexSingleton)) &&
                !CONTAINS(saturatedFaces, i) &&
                CONTAINS_ALL(currentCycle, verticesInFace[i])){
            //we found a new empty face
            ADD(emptyFaces, i);
        }
    }
    if(IS_NOT_EMPTY(INTERSECTION(emptyFaces, facesLeft)) && 
            IS_NOT_EMPTY(INTERSECTION(emptyFaces, facesRight))){
        //there is an empty face on both sides of the cycle
        return FALSE;
    }
    
    ADD(currentCycle, newEdge->end);
    ADD_ALL(saturatedFaces, newEdge->incident_faces);
    
    if(remainingVertices == 1){
        //this was the last vertex: we try to close the cycle
        if(CONTAINS(neighbours[newEdge->end], firstVertexCycle)){
            //the cycle can be closed
            
            //find closing edge
            e = elast = firstedge[newEdge->end];
            while (e->end != firstVertexCycle){
                e = e->next;
            }
            
            finishCycle(e, saturatedFaces,
                    UNION(facesRight, facesBetween(e, newEdge->inverse)),
                    UNION(facesLeft, facesBetween(newEdge->inverse, e)),
                    emptyFaces);
        }
    } else {
        //just continue the cycle
        e = elast = firstedge[newEdge->end];
        do {
            if(!CONTAINS(currentCycle, e->end)){
                if(continueCycle(e, remainingVertices - 1, saturatedFaces,
                        UNION(facesRight, facesBetween(e, newEdge->inverse)),
                        UNION(facesLeft, facesBetween(newEdge->inverse, e)),
                        emptyFaces)){
                    return TRUE;
                }
            }
            e = e->next;
        } while (e!=elast);
    }
    
    REMOVE(currentCycle, newEdge->end);
    
    return FALSE;
}

boolean hasPathHamiltonianCycle(){
    int i;
    int minDegree;
    int minDegreeVertex;
    EDGE *e, *elast;
    EDGE *e2, *elast2;
    
    //look for the vertex with the smallest degree
    minDegree = 6; //there is always a vertex with degree at most 5 in a plane graph
    minDegreeVertex = 0; //just some vertex
    for(i = 0; i < nv; i++){
        if(degree[i] < minDegree){
            minDegree = degree[i];
            minDegreeVertex = i;
        }
    }
    
    //start looking for a cycle
    firstVertexCycle = minDegreeVertex;
    currentCycle = SINGLETON(minDegreeVertex);
    e = elast = firstedge[minDegreeVertex];
    do {
        ADD(currentCycle, e->end);
        firstEdgeCycle = e;
        bitset saturatedFaces = e->incident_faces;
        e2 = elast2 = firstedge[e->end];
        do {
            if(!CONTAINS(currentCycle, e2->end)){
                bitset facesRight = facesBetween(e2, e->inverse);
                bitset facesLeft = facesBetween(e->inverse, e2);
                if(continueCycle(e2, nv - 2, saturatedFaces,
                        facesRight, facesLeft, EMPTY_SET)){
                    return TRUE;
                }
            }
            e2 = e2->next;
        } while (e2 != elast2);
        REMOVE(currentCycle, e->end);
        e = e->next;
    } while (e != elast);
    
    return FALSE;
}

//=============== Output of binary graph code ===========================

void writeCode(FILE *f, unsigned short code[], int length){
    int i;
    for(i = 0; i < length; i++){
        putc(code[i], f);
    }
}

//=============== Reading and decoding planarcode ===========================

EDGE *findEdge(int from, int to) {
    EDGE *e, *elast;

    e = elast = firstedge[from];
    do {
        if (e->end == to) {
            return e;
        }
        e = e->next;
    } while (e != elast);
    fprintf(stderr, "error while looking for edge from %d to %d.\n", from, to);
    exit(0);
}

/* Store in the rightface field of each edge the number of the face on
   the right hand side of that edge.  Faces are numbered 0,1,....  Also
   store in facestart[i] an example of an edge in the clockwise orientation
   of the face boundary, and the size of the face in facesize[i], for each i.
   Returns the number of faces. */
void makeDual() {
    register int i, sz;
    register EDGE *e, *ex, *ef, *efx;

    RESETMARKS;

    nf = 0;
    for (i = 0; i < nv; ++i) {

        e = ex = firstedge[i];
        do {
            if (!ISMARKEDLO(e)) {
                facestart[nf] = ef = efx = e;
                sz = 0;
                do {
                    ef->rightface = nf;
                    MARKLO(ef);
                    ef = ef->inverse->prev;
                    ++sz;
                } while (ef != efx);
                faceSize[nf] = sz;
                ++nf;
            }
            e = e->next;
        } while (e != ex);
    }
}

void decodePlanarCode(unsigned short* code) {
    /* complexity of method to determine inverse isn't that good, but will have to satisfy for now
     */
    int i, j, codePosition;
    int edgeCounter = 0;
    EDGE *inverse;

    nv = code[0];
    codePosition = 1;

    for (i = 0; i < nv; i++) {
        degree[i] = 0;
        firstedge[i] = edges + edgeCounter;
        edges[edgeCounter].start = i;
        edges[edgeCounter].end = code[codePosition] - 1;
        edges[edgeCounter].next = edges + edgeCounter + 1;
        if (code[codePosition] - 1 < i) {
            inverse = findEdge(code[codePosition] - 1, i);
            edges[edgeCounter].inverse = inverse;
            inverse->inverse = edges + edgeCounter;
        } else {
            edges[edgeCounter].inverse = NULL;
        }
        neighbours[i] = SINGLETON(code[codePosition] - 1);
        edgeCounter++;
        codePosition++;
        for (j = 1; code[codePosition]; j++, codePosition++) {
            if (j == MAXVAL) {
                fprintf(stderr, "MAXVAL too small: %d\n", MAXVAL);
                exit(0);
            }
            edges[edgeCounter].start = i;
            edges[edgeCounter].end = code[codePosition] - 1;
            edges[edgeCounter].prev = edges + edgeCounter - 1;
            edges[edgeCounter].next = edges + edgeCounter + 1;
            if (code[codePosition] - 1 < i) {
                inverse = findEdge(code[codePosition] - 1, i);
                edges[edgeCounter].inverse = inverse;
                inverse->inverse = edges + edgeCounter;
            } else {
                edges[edgeCounter].inverse = NULL;
            }
            ADD(neighbours[i], code[codePosition] - 1);
            edgeCounter++;
        }
        firstedge[i]->prev = edges + edgeCounter - 1;
        edges[edgeCounter - 1].next = firstedge[i];
        degree[i] = j;

        codePosition++; /* read the closing 0 */
    }

    ne = edgeCounter;

    makeDual();

    // nv - ne/2 + nf = 2
    
    //store some additional data
    for(i = 0; i < MAXF; i++){
        verticesInFace[i] = EMPTY_SET;
    }
    
    for(i = 0; i < ne; i++){
        edges[i].incident_faces = UNION(SINGLETON(edges[i].rightface),
                SINGLETON(edges[i].inverse->rightface));
        ADD(verticesInFace[edges[i].rightface], edges[i].end);
    }
}

/**
 * 
 * @param code
 * @param laenge
 * @param file
 * @return returns 1 if a code was read and 0 otherwise. Exits in case of error.
 */
int readPlanarCode(unsigned short code[], int *length, FILE *file) {
    static int first = 1;
    unsigned char c;
    char testheader[20];
    int bufferSize, zeroCounter;
    
    int readCount;


    if (first) {
        first = 0;

        if (fread(&testheader, sizeof (unsigned char), 13, file) != 13) {
            fprintf(stderr, "can't read header ((1)file too small)-- exiting\n");
            exit(1);
        }
        testheader[13] = 0;
        if (strcmp(testheader, ">>planar_code") == 0) {

        } else {
            fprintf(stderr, "No planarcode header detected -- exiting!\n");
            exit(1);
        }
        //read reminder of header (either empty or le/be specification)
        if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
            return FALSE;
        }
        while (c!='<'){
            if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
                return FALSE;
            }
        }
        //read one more character
        if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
            return FALSE;
        }
    }

    /* possibly removing interior headers -- only done for planarcode */
    if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
        //nothing left in file
        return (0);
    }

    if (c == '>') {
        // could be a header, or maybe just a 62 (which is also possible for unsigned char
        code[0] = c;
        bufferSize = 1;
        zeroCounter = 0;
        code[1] = (unsigned short) getc(file);
        if (code[1] == 0) zeroCounter++;
        code[2] = (unsigned short) getc(file);
        if (code[2] == 0) zeroCounter++;
        bufferSize = 3;
        // 3 characters were read and stored in buffer
        if ((code[1] == '>') && (code[2] == 'p')) /*we are sure that we're dealing with a header*/ {
            while ((c = getc(file)) != '<');
            /* read 2 more characters: */
            c = getc(file);
            if (c != '<') {
                fprintf(stderr, "Problems with header -- single '<'\n");
                exit(1);
            }
            if (!fread(&c, sizeof (unsigned char), 1, file)) {
                //nothing left in file
                return (0);
            }
            bufferSize = 1;
            zeroCounter = 0;
        }
    } else {
        //no header present
        bufferSize = 1;
        zeroCounter = 0;
    }

    if (c != 0) /* unsigned chars would be sufficient */ {
        code[0] = c;
        if (code[0] > MAXN) {
            fprintf(stderr, "Constant N too small %d > %d \n", code[0], MAXN);
            exit(1);
        }
        while (zeroCounter < code[0]) {
            code[bufferSize] = (unsigned short) getc(file);
            if (code[bufferSize] == 0) zeroCounter++;
            bufferSize++;
        }
    } else {
        readCount = fread(code, sizeof (unsigned short), 1, file);
        if(!readCount){
            fprintf(stderr, "Unexpected EOF.\n");
            exit(1);
        }
        if (code[0] > MAXN) {
            fprintf(stderr, "Constant N too small %d > %d \n", code[0], MAXN);
            exit(1);
        }
        bufferSize = 1;
        zeroCounter = 0;
        while (zeroCounter < code[0]) {
            readCount = fread(code + bufferSize, sizeof (unsigned short), 1, file);
            if(!readCount){
                fprintf(stderr, "Unexpected EOF.\n");
                exit(1);
            }
            if (code[bufferSize] == 0) zeroCounter++;
            bufferSize++;
        }
    }

    *length = bufferSize;
    return (1);


}

 //====================== USAGE =======================

void help(char *name) {
    fprintf(stderr, "The program %s checks whether a plane triangulation contains a hamiltonian\ncycle with all faces that are missed on the same side.\n\n", name);
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "    -f, --filter\n");
    fprintf(stderr, "       Filter graphs that have the property.\n");
    fprintf(stderr, "    -i, --invert\n");
    fprintf(stderr, "       Invert the filter.\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    /*=========== commandline parsing ===========*/
    boolean invert = FALSE;
    boolean filter = FALSE;
    int c;
    char *name = argv[0];
    static struct option long_options[] = {
        {"invert", no_argument, NULL, 'i'},
        {"filter", no_argument, NULL, 'f'},
        {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "hif", long_options, &option_index)) != -1) {
        switch (c) {
            case 'i':
                invert = TRUE;
                break;
            case 'f':
                filter = TRUE;
                break;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case '?':
                usage(name);
                return EXIT_FAILURE;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }
    
    /*=========== process graphs ===========*/
    unsigned long long numberOfGraphs = 0ULL;
    unsigned long long numberOfGraphsWithPathHamiltonianCycle = 0ULL;
    unsigned long long numberOfGraphsWithoutPathHamiltonianCycle = 0ULL;
    unsigned short code[MAXCODELENGTH];
    int length;
    if(filter){
        fprintf(stdout, ">>planar_code<<");
    }
    while (readPlanarCode(code, &length, stdin)) {
        decodePlanarCode(code);
        if(hasPathHamiltonianCycle()){
            numberOfGraphsWithPathHamiltonianCycle++;
            if(filter && !invert){
                writeCode(stdout, code, length);
            }
        } else {
            numberOfGraphsWithoutPathHamiltonianCycle++;
            if(filter && invert){
                writeCode(stdout, code, length);
            }
        }
        numberOfGraphs++;
    }
    
    fprintf(stderr, "Read %llu graph%s.\n", numberOfGraphs, 
                numberOfGraphs==1 ? "" : "s");
    
    fprintf(stderr, "   %llu graph%s %s a path-hamiltonian cycle.\n",
                numberOfGraphsWithPathHamiltonianCycle, 
                numberOfGraphsWithPathHamiltonianCycle==1 ? "" : "s", 
                numberOfGraphsWithPathHamiltonianCycle==1 ? "has" : "have");
    fprintf(stderr, "   %llu graph%s %s not have a path-hamiltonian cycle.\n",
                numberOfGraphsWithoutPathHamiltonianCycle, 
                numberOfGraphsWithoutPathHamiltonianCycle==1 ? "" : "s", 
                numberOfGraphsWithoutPathHamiltonianCycle==1 ? "does" : "do");
    

    return EXIT_SUCCESS;
}
