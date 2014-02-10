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


#include <metis.h>

void print_help(char program_name[])
{
	fprintf(stderr,"Usage: %s ",
			program_name);
	fprintf(stderr, " [-e <clusteringFile>] [-n <nSingletons>] <GraphFile>\n");
}

int main(int argc, char* argv[])
{
	int nparts,writeNcutVector=0;
	int	defaultOutputFile=1,wgtflag=0,i, txtFormat=0;
	int isDirected=0, hasPagerank=0, minClusterSize=20,
	maxClusterSize=50000;
	long numMemberships = 0;
	char clusterFile[256], filename[256], outputFile[256],
			pagerankFile[256],membershipsFile[256];
	int nSingletons = 0;
	float ncut=0;
	idxtype minEdgeWeight = 0;
	GraphType graph;
	idxtype* part;

	if ( argc < 4 )
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
				case 'e':
				case 'E':
					strcpy(clusterFile,*(++argv));
					break;
				case 'n':
				case 'N':
					nSingletons = atoi(*(++argv));
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

	if ( defaultOutputFile > 0 )
	{
		sprintf(outputFile,"ncutVector.%s",clusterFile);
		if ( numMemberships > 0 )
			sprintf(outputFile, "ncutVector");
	}
		
	ReadGraph(&graph, filename, &wgtflag, 1, txtFormat);
	if (graph.nvtxs <= 0) {
	  printf("Empty graph. Nothing to do.\n");
	  exit(0);
	}

	if ( minEdgeWeight > 0 && wgtflag&1 > 0 )
	{
		printf("Going to retain only edges with weight %d or above, in the graph\n"
					, minEdgeWeight ); 

		exactPruneGraph(&graph, 0, minEdgeWeight );
		printf("Pruned graph contains %d edges.\n", graph.nedges);
	}

	graph.isDirected = isDirected;
	if ( isDirected > 0 )
	{
		if ( hasPagerank > 0 )
		{
			graph.pagerank  = (wgttype*) GKmalloc(
			sizeof(wgttype)*graph.nvtxs, "ncut:pagerank");
			FILE *fpin;
			if ( (fpin = fopen(pagerankFile, "r")) == NULL )
			{
				printf("Could not open %s for reading\n",
				pagerankFile);
				abort();
			}

			int i;
			for ( i=0; i<graph.nvtxs; i++ )
			{
				fscanf(fpin, "%f", &graph.pagerank[i]);
				graph.pagerank[i] = pow(10.0, graph.pagerank[i]);
			}
		}
	}

	if ( numMemberships > 0 )
	{
		idxtype *members, *groups_xadj, *groupIds, ngroups=0;
		members = idxmalloc(numMemberships,	"readMemberships:members");
		readMemberships( membershipsFile, numMemberships,
			members, &groups_xadj, &groupIds, &ngroups);

		float *ncutVector = (float*) GKmalloc(
		sizeof(float)*ngroups, "ncut.c:ncutVector");
		float avgNcut
			= Overlap_ComputeNcutVector(&graph, members, groups_xadj,
				ngroups, ncutVector);
		
		FILE *fp;
		if ( (fp=fopen(outputFile,"w")) == NULL )
		{
			fprintf(stderr,"Could not open file:%s\n",outputFile);
			return -1;
		}

		int i, numClusters=0;
		wgttype wgtdAvgCut = 0, totalWgts = 0;
		avgNcut = 0;
		for ( i=0; i<ngroups; i++ )
		{
			int grpSize = (groups_xadj[i+1]-groups_xadj[i]);
			if ( grpSize >= minClusterSize && grpSize <=
				maxClusterSize )
			{
				fprintf( fp, "%d %.4f\n", groupIds[i],
				ncutVector[i]);
				avgNcut += ncutVector[i];
				numClusters++;
				wgtdAvgCut +=  grpSize * ncutVector[i];
				totalWgts += grpSize;
			}
		}
		wgtdAvgCut /= totalWgts;
		avgNcut /= numClusters++;
		fprintf(fp, "Avg. Conductance: %.4f\n", avgNcut);
		fprintf(fp, "Avg. Weighted Conductance: %.4f\n", wgtdAvgCut);
		fclose(fp);

		GKfree( (void**)&ncutVector, (void**)&groups_xadj,
			(void**)&groupIds, (void**)&members, LTERM);
		exit(0);	

	}

	part=idxmalloc(graph.nvtxs,"main: part");
	nparts=readClustering(clusterFile, part, graph.nvtxs);
	nparts=mapPartition(part, graph.nvtxs);
//	fprintf(stderr, "Num partitions:%d\n", nparts);
	if ( writeNcutVector > 0 )
	{
		FILE *fp;
		float* ncutVector = 
				(float*)malloc(sizeof(float)*graph.nvtxs);
		ncut=ComputeNCutVector(&graph,part,nparts,ncutVector);
		if ( (fp=fopen(outputFile,"w")) == NULL )
		{
			fprintf(stderr,"Could not open file:%s\n",outputFile);
			return -1;
		}
		for(i=0;i<nparts;i++)
			fprintf(fp,"%.3f\n",ncutVector[i]);
		fclose(fp);
		free(ncutVector);
	}
	else
	{
		int noOfSingletons = 0;
		idxtype* newIds = lookForSingletons(&graph,
		&noOfSingletons);
		ncut=ComputeNCut(&graph,part,nparts);
		if ( nSingletons > noOfSingletons )
			nSingletons -= noOfSingletons; // don't want double
			//counting
		ncut -= nSingletons;
		nparts -= nSingletons;

//		fprintf(stderr, "noOfSingletons:%d nSingletons:%d\n",
//		noOfSingletons, nSingletons);
		//ncut=ComputeConductance(&graph,part,nparts);
		if ( noOfSingletons > 0 )
		{
		
/*			GraphType* noSingletonGraph;
			getSubgraph(&graph, nodeMap, graph.nvtxs-noOfSingletons, 
						wgtflag, &noSingletonGraph);
			free(graph.gdata);
			mapIndices(part, nodeMap, graph.nvtxs,
			npart-noOfSingletons);
			float newNcut = ComputeNCut( noSingletonGraph, 
*/
			idxtype* newPart = idxmalloc( graph.nvtxs-noOfSingletons,
								"ncut:newPart");

			for ( i=0; i<graph.nvtxs; i++ )
			{
				if ( newIds[i] > -1 )
					newPart[newIds[i]] = part[i];
			}

			int nnewparts = mapPartition(newPart, graph.nvtxs-noOfSingletons);

		/*	fprintf(stderr, "nnewparts:%d, graph.nvtxs:%d,",
			nnewparts, graph.nvtxs);
			fprintf(stderr, "noOfSingletons:%d\n",
			noOfSingletons);
			fprintf(stderr, "nnewparts:%d\n", nnewparts);
*/
			idxtype *clusterSizes = histogram(newPart,
						graph.nvtxs-noOfSingletons, nnewparts);

			int maxSize = clusterSizes[idxamax(nnewparts, clusterSizes)];
			float avgClusterSize =
			//(graph.nvtxs-noOfSingletons-nSingletons)*1.0/(nnewparts);
			(graph.nvtxs-noOfSingletons-nSingletons)*1.0/(nparts-noOfSingletons);
			float balance =	(maxSize*1.0) / avgClusterSize;
			float stdDevn = stdDeviation(clusterSizes, nnewparts);
		//	float avgNcut = ncut * 1.0/nnewparts;
			float avgNcut = ncut * 1.0/(nparts-noOfSingletons);
			float normStdDevn = stdDevn/avgClusterSize;
	
		// Warning: This computation only works if the singletons
		// have been placed in their own clusters. This works for
		// MLR-MCL, in other words, because it is guaranteed to
		// place singletons in their own clusters.
			printf("For graph without singletons\n");
			printf("Clusters: %d N-Cut: %.3f", 
		//		nnewparts, ncut);
				nparts-noOfSingletons, ncut);
			printf(" AvgN-Cut: %.3f Balance: %.2f ",avgNcut,
						balance); 
			printf("Std_Deviation: %.2f ", stdDevn);
			printf("Coefficient_of_Variation: %.2f\n", normStdDevn);

			printf("For input graph\n");
	//		free( newPart );
	//		free( clusterSizes );
		}

	//	free(newIds);
	}

	idxtype* clusterSizes = histogram(part, graph.nvtxs, nparts);
	int maxSize = clusterSizes[idxamax(nparts, clusterSizes)];
	float avgClusterSize = (graph.nvtxs)*1.0/(nparts);
	float balance = (maxSize*1.0)/avgClusterSize;
	float stdDevn = stdDeviation(clusterSizes, nparts);
	float normStdDevn = stdDevn/avgClusterSize;
	float avgNcut = ncut * 1.0/nparts;
	
	printf("Clusters: %d N-Cut: %.3f AvgN-Cut: %.3f", nparts,
						ncut, avgNcut );
	printf(" Balance: %.2f Std.Deviation: %.2f ",
				 balance, stdDevn);
	printf("Coefficient_of_Variation: %.2f\n", normStdDevn);
}
