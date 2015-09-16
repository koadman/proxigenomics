/*
 * MLR-MCL (Multi-Level Regularized Markov Clustering) - Version 1.2
 * 
 * Copyright 2010, The Ohio State University. 
 * Portions Copyright 1997, Regents of the University of Minnesota.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or
 * without modification, are
 * permitted provided that the following conditions are met:
 * 
 * 1. Redistributions, with or without modifications, of source 
 * code must retain the above
 * copyright notice, this list of conditions and the following
 * disclaimer.
 * 
 * 2. Redistributions, with or without modifications, in binary form 
 * must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer
 * in the documentation and/or other materials provided with the
 * distribution.
 * 
 * 3. The names of the Ohio
 * State University, the University of Minnesota and
 * their contributors may not be used to endorse or promote products
 * derived from this software without specific prior permission. 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 * TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 * 
 */

/*
 * Contributors: Venu Satuluri, Srinivasan Parthasarathy
 * 
 * Reference: "Scalable Graph Clustering using Stochastic Flows",
 * KDD 2009. http://doi.acm.org/10.1145/1557019.1557101
 */


/*
 *
 */

#include <metis.h>

void print_help(const char * program_name)
{
	fprintf(stderr, "Usage: %s -e <clusteringFile>",
	program_name); 
	fprintf(stderr, " [-m <mergeHeuristic>]");
	fprintf(stderr, " -n <numMerges> <GraphFile>");
	fprintf(stderr, " [-o <outputFile>]\n");
	fprintf(stderr, "MergeHeuristic: 1 - cut/maxvol (default),");
	fprintf(stderr, "2 - ncut optimal (can lead to high imbalance)");
	fprintf(stderr, "\n");

}

int main(int argc, char *argv[])
{
	GraphType graph;
	char filename[256], clusterFile[256], outputFile[256];
	int numMerges = 0;
	int clusterFileGiven = 0;
	int defaultOutputFile = 1;
	int mergeHeuristic = 1;

	if ( argc < 6 )
	{
		print_help(argv[0]);
		return -1;
	}

	for( argv++; *argv != NULL; argv++)
	{
		if ( (*argv)[0] == '-' )
		{
			switch((*argv)[1])
			{
				case 'm':
				case 'M':
				  mergeHeuristic = atoi(*(++argv));
				  break;
				case 'n':
				case 'N':
				  numMerges = atoi(*(++argv));
				  break;
				case 'e':
				case 'E':
				  strcpy(clusterFile, *(++argv));
				  clusterFileGiven = 1;
				  break;
				case 'o':
				case 'O':
				  defaultOutputFile=0;
				  strcpy(outputFile, *(++argv));
				  break;
				default:
				   	printf("Invalid switch %s\n", *argv);
			  		print_help(argv[0]);
	  				return -1;
			}
		}
		else
			strcpy(filename, *argv);
	}

	if ( numMerges == 0 || clusterFileGiven == 0 )
	{
		print_help(argv[0]);
		return -1;
	}

	if ( defaultOutputFile > 0 )
	{
		sprintf(outputFile, "%s.%dmerges", clusterFile,
		numMerges);
	}
				 
	timer iotmr, maintmr, totaltmr;
	cleartimer(iotmr); cleartimer(maintmr); cleartimer(totaltmr);

	starttimer(iotmr); starttimer(totaltmr);

	int wgtflag, txtFormat = 0;
	ReadGraph(&graph, filename, &wgtflag, 1, txtFormat);
	if (graph.nvtxs <= 0) {
	  printf("Empty graph. Nothing to do.\n");
	  exit(0);
	}

	int *part;
	int nparts;

	part=idxmalloc(graph.nvtxs,"main: part");
	nparts=readClustering(clusterFile, part, graph.nvtxs);
	nparts=mapPartition(part, graph.nvtxs);
	stoptimer(iotmr);

	starttimer(maintmr);

	ListGraph* cg = createClusterGraph(part, nparts, &graph);
	printf("Cluster graph: %d nodes, %d edges\n", cg->nvtxs,
			cg->nedges);
	idxtype *part_mapper = idxmalloc(nparts, "main:part_mapper");
	for ( int i=0; i<nparts; i++ )
		part_mapper[i] = i;

	for ( int i=0; i<numMerges; i++ )
		mergeBestClusters(cg, part_mapper, nparts, mergeHeuristic);
	
	for ( int i=0; i<graph.nvtxs; i++ )
		part[i] = part_mapper[part[i]];

	int new_nparts = mapPartition(part, graph.nvtxs);
	if ( new_nparts != nparts - numMerges )
	{
		printf("Yikes! new_nparts (%d) != ", new_nparts);
		printf("nparts (%d) - numMerges(%d)\n", nparts, numMerges);
	}

	stoptimer(maintmr);
	starttimer(iotmr);

	my_WritePartition(outputFile, part, graph.nvtxs, 0);
	printf("Output written to file: %s\n", outputFile);

	stoptimer(iotmr);
	stoptimer(totaltmr);

	printf("IO Time: %.2f\n", gettimer(iotmr));
	printf("Processing time: %.2f\n", gettimer(maintmr));
	printf("Total time: %.2f\n", gettimer(totaltmr));
}
