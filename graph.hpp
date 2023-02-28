#include <limits.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h>
#include <fstream>
#include <set>
#include <sstream>

const int UNSET_ELEMENT = -1;

typedef struct Set {
    int elements;
    int* canonicalElements;
    int* rank;
} Set;

typedef struct WeightedGraph {
    int edges;
    int nodes;
    int* edgeList;
} WeightedGraph;

// initialize and allocate memory for the members of the graph
void newWeightedGraph(WeightedGraph* graph, const int edges) {
    graph->edges = edges;
    graph->edgeList = (int*)calloc(edges * 3, sizeof(int));
}

// read a previously generated maze file and store it in the graph
void readGraphFile(WeightedGraph* graph, const char filePath[]) {
    // open the file
    std::ifstream infile;
    infile.open(filePath);
    if (!infile) {
        std::cout << "Cannot open the file" << std::endl;
        exit(0);
    }
    std::string line;
    std::stringstream ss;
    int nodes = 0, edges = 0;
    edges = std::count(std::istreambuf_iterator<char>(infile),
                       std::istreambuf_iterator<char>(), '\n');
    infile.clear();
    infile.seekg(0);
    newWeightedGraph(graph, edges);

    edges = 0;
    std::set<int> nodeList;
    while (getline(infile, line)) {
        if (line[0] < '0' || line[0] > '9') {
            continue;
        }
        ss.clear();
        ss << line;

        int source;
        int destination;
        int weight = rand() % 100 + 1;
        if (ss >> source && ss >> destination) {
            graph->edgeList[3 * edges] = source;
            graph->edgeList[3 * edges + 1] = destination;
            graph->edgeList[3 * edges + 2] = weight;
            //printf("Added edge(%d, %d, %d)\n", source, destination, weight);
            nodeList.insert(source);
            nodeList.insert(destination);
            edges++;
        }
    }
    graph->edges = edges;
    graph->nodes = nodeList.size();
    printf("Made graph %d %d\n", graph->edges, graph->nodes);
    infile.close();
}

// print all edges of the graph in "from to weight" format
void printWeightedGraph(const WeightedGraph* graph) {
    printf("------------------------------------------------\n");
    for (int i = 0; i < graph->edges; i++) {
        printf("Edge %d : ", i);
        for (int j = 0; j < 3; j++) {
            printf("%d\t", graph->edgeList[i * 3 + j]);
        }
        printf("\n");
    }
    printf("------------------------------------------------\n");
}

void newSet(Set* set, const int elements) {
    set->elements = elements;
    set->canonicalElements =
        (int*)malloc(elements * sizeof(int));  // maintain parent
    memset(set->canonicalElements, UNSET_ELEMENT, elements * sizeof(int));
    set->rank = (int*)calloc(elements, sizeof(int));  //  maintain rank
}

// return the canonical element of a vertex with path compression
int findSet(const Set* set, const int vertex) {
    if (set->canonicalElements[vertex] == UNSET_ELEMENT) {
        return vertex;
    } else {
        set->canonicalElements[vertex] =
            findSet(set, set->canonicalElements[vertex]);
        return set->canonicalElements[vertex];
    }
}

// merge the set of parent1 and parent2 with union by rank
void unionSet(Set* set, const int parent1, const int parent2) {
    int root1 = findSet(set, parent1);
    int root2 = findSet(set, parent2);

    if (root1 == root2) {
        return;
    }
    // Attach smaller rank tree under root of high
    else if (set->rank[root1] < set->rank[root2]) {
        set->canonicalElements[root1] = root2;
    } else if (set->rank[root1] > set->rank[root2]) {
        set->canonicalElements[root2] = root1;
    }
    // If ranks are same, then make one as root and
    // increment its rank by one
    else {
        set->canonicalElements[root1] = root2;
        set->rank[root2] = set->rank[root1] + 1;
    }
}

// copy an edge
void copyEdge(int* to, int* from) { memcpy(to, from, 3 * sizeof(int)); }

// scatter the edge list of a graph
void scatterEdgeList(int* edgeList, int* edgeListPart, const int elements,
                     int* elementsPart) {
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Scatter(edgeList, *elementsPart * 3, MPI_INT, edgeListPart,
                *elementsPart * 3, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == size - 1 && elements % *elementsPart != 0) {
        // number of elements and processes isn't divisible without remainder
        *elementsPart = elements % *elementsPart;
    }

    if (elements / 2 + 1 < size && elements != size) {
        if (rank == 0) {
            fprintf(stderr, "Unsupported size/process combination, exiting!\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}

// cleanup set data
void deleteSet(Set* set) {
    free(set->canonicalElements);
    free(set->rank);
}

// cleanup graph data
void deleteWeightedGraph(WeightedGraph* graph) { free(graph->edgeList); }

// find a MST of the graph using Boruvka's algorithm
void mstBoruvka(const WeightedGraph* graph, WeightedGraph* mst) {
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    bool parallel = size != 1;

    // send number of edges and nodes
    int edges;
    int nodes;
    if (rank == 0) {
        edges = graph->edges;
        nodes = graph->nodes;
        MPI_Bcast(&edges, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Bcast(&edges, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    // scatter the edges to search in them
    int edgesPart = (edges + size - 1) / size;
    int* edgeListPart = (int*)malloc(edgesPart * 3 * sizeof(int));
    if (parallel) {
        scatterEdgeList(graph->edgeList, edgeListPart, edges, &edgesPart);
    } else {
        edgeListPart = graph->edgeList;
    }

    // create needed data structures
    Set set = {.elements = 0, .canonicalElements = NULL, .rank = NULL};
    newSet(&set, nodes);

    int edgesMST = 0;
    int* closestEdge = (int*)malloc(nodes * 3 * sizeof(int));
    int* closestEdgeRecieved;
    if (parallel) {
        closestEdgeRecieved = (int*)malloc(nodes * 3 * sizeof(int));
    }

    for (int i = 1; i < nodes && edgesMST < nodes - 1; i *= 2) {
        // reset all closestEdge
        for (int j = 0; j < nodes; j++) {
            closestEdge[j * 3 + 2] = INT_MAX;
        }

        // find closestEdge
        for (int j = 0; j < edgesPart; j++) {
            int* currentEdge = &edgeListPart[j * 3];
            int canonicalElements[2] = {findSet(&set, currentEdge[0]),
                                        findSet(&set, currentEdge[1])};

            // eventually update closestEdge
            if (canonicalElements[0] != canonicalElements[1]) {
                for (int k = 0; k < 2; k++) {
                    bool closestEdgeNotSet =
                        closestEdge[canonicalElements[k] * 3 + 2] == INT_MAX;
                    bool weightSmaller =
                        currentEdge[2] <
                        closestEdge[canonicalElements[k] * 3 + 2];
                    if (closestEdgeNotSet || weightSmaller) {
                        copyEdge(&closestEdge[canonicalElements[k] * 3],
                                 currentEdge);
                    }
                }
            }
        }

        if (parallel) {
            int from;
            int to;
            for (int step = 1; step < size; step *= 2) {
                if (rank % (2 * step) == 0) {
                    from = rank + step;
                    if (from < size) {
                        MPI_Recv(closestEdgeRecieved, nodes * 3, MPI_INT, from,
                                 0, MPI_COMM_WORLD, &status);

                        // combine all closestEdge parts
                        for (int i = 0; i < nodes; i++) {
                            int currentVertex = i * 3;
                            if (closestEdgeRecieved[currentVertex + 2] <
                                closestEdge[currentVertex + 2]) {
                                copyEdge(&closestEdge[currentVertex],
                                         &closestEdgeRecieved[currentVertex]);
                            }
                        }
                    }
                } else if (rank % step == 0) {
                    to = rank - step;
                    MPI_Send(closestEdge, nodes * 3, MPI_INT, to, 0,
                             MPI_COMM_WORLD);
                }
            }
            // publish all closestEdge parts
            MPI_Bcast(closestEdge, nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
        }

        // add new edges to MST
        for (int j = 0; j < nodes; j++) {
            if (closestEdge[j * 3 + 2] != INT_MAX) {
                int from = closestEdge[j * 3];
                int to = closestEdge[j * 3 + 1];

                // prevent adding the same edge twice
                if (findSet(&set, from) != findSet(&set, to)) {
                    if (rank == 0) {
                        copyEdge(&mst->edgeList[edgesMST * 3],
                                 &closestEdge[j * 3]);
                    }
                    edgesMST++;
                    unionSet(&set, from, to);
                }
            }
        }
    }

    // clean up
    deleteSet(&set);
    free(closestEdge);
    if (parallel) {
        free(closestEdgeRecieved);
        free(edgeListPart);
    }
}
