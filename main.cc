#include "graph.hpp"

void Compute_MST(const char filePath[]) {
	// MPI variables and initialization
	int rank;
	int size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	// graph Variables
	WeightedGraph graph = { .edges = 0, .nodes = 0,   .edgeList = NULL };
	WeightedGraph mst = { .edges = 0, .nodes = 0,	.edgeList = NULL };

	if (rank == 0) {

		// read the maze file and store it in the graph
        double start = MPI_Wtime();
		readGraphFile(&graph, filePath);
        double end = MPI_Wtime();

        printf("Parse time : %f s\n", end - start);
		// print the edges of the read graph
		//printf("Original Graph:\n");
		//printWeightedGraph(&graph);

		newWeightedGraph(&mst, graph.nodes - 1);
        mst.nodes = graph.nodes;
        //printf("Complete\n");
	}
	
    double start = MPI_Wtime();
	// use Boruvka's algorithm
	mstBoruvka(&graph, &mst);

    double end = MPI_Wtime();

	if (rank == 0) {

		// print the edges of the MST
		printf("Minimum Spanning Tree (Boruvka):\n");
		//printWeightedGraph(&mst);

		unsigned long weightMST = 0;
		for (int i = 0; i < mst.edges; i++) {
			weightMST += mst.edgeList[i * 3 + 2];
		}

		printf("MST weight: %lu\n", weightMST);
		printf("Time elapsed: %f s\n", end - start);
		// cleanup
		deleteWeightedGraph(&graph);
		deleteWeightedGraph(&mst);
	}

	MPI_Finalize();
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        // printf("Execute ./a.out input_graph_file numberOfProcesses\n");
        exit(0);
    }
    Compute_MST(argv[1]);
    return 0;
}
