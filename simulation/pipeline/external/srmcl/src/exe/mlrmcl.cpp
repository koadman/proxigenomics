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


void print_help(char program_name[])
{
	FILE *fp = stderr;
	fprintf(fp, "Usage: %s [options] graph_file\n", program_name);
	fprintf(fp, "Options:\n");
	fprintf(fp, "\t-c <integer>\t\tsize of coarsest graph (default 1000)\n");
	fprintf(fp, "\t-b <float>\t\tbalance parameter (default 0.5)\n");
	fprintf(fp, "\t-i <float>\t\tinflation parameter (default 2.0)\n");
	fprintf(fp, "\t-o <string>\t\toutput file\n");
	return;
}

int main(int argc, char *argv[])
{
	int i, npart;
	idxtype *part;
	float ncut=0;
	GraphType graph;
	char filename[256],outputFile[256];
	int wgtflag = 0, addSelfLoop=1, outputFileGiven=0, txtFormat=0 ;
	int randomInit = 0;
	idxtype minEdgeWeight = 0;
	Options opt;
	timer TOTALTmr, METISTmr, IOTmr;

	initOptions(&opt);

	if (argc < 2) {
		print_help(argv[0]);
		exit(0);
	}
	
	for (argv++; *argv != NULL; argv++){
	    if ((*argv)[0] == '-')
		{
			int temp;
	      	switch ((*argv)[1])
			{
			case 'b':
			case 'B':
			  opt.penalty_power=atof(*(++argv));
			  break;
			case 'i':
			case 'I':
			  opt.gamma=atof(*(++argv));
			  break;
			case 'o':
			case 'O':
			  strcpy(outputFile,*(++argv));
			  outputFileGiven=1;
			  break;
			case 'D'://quality threshold. This is a post-processing step proposed in SR-MCL. If you dont want post-processing (this is what original MLR-MCL, R-MCL, MCL do, please set "-d 0"  
			case 'd':
			  opt.quality_threshold = atof(*(++argv));
			  break;
			case 'w':
			case 'W':
		          opt.weighted_density = true; 
			  break;

			case 'c':
			case 'C':
			  opt.coarsenTo= atoi(*(++argv));
			  break;
			default:
			  printf("Invalid option %s\n", *argv);
			  print_help(argv[0]);
			  exit(0);
			}
		}
	    else
		{
	      strcpy(filename, *argv);
	    }
	}  

	if ( randomInit > 0 )
		InitRandom(time(NULL));
	else
		InitRandom(-1);

	cleartimer(TOTALTmr);
	cleartimer(METISTmr);
	cleartimer(IOTmr);

	starttimer(TOTALTmr);
	starttimer(IOTmr);

	ReadGraph(&graph, filename, &wgtflag, addSelfLoop, txtFormat);

	if ( opt.matchType == MATCH_UNSPECIFIED )
	{
//		opt.matchType = (graph.nvtxs>50000) ? MATCH_POWERLAW_FC :
//							MATCH_SHEMN;
		opt.matchType = MATCH_SHEMN;
	}
	
	stoptimer(IOTmr);

	if (graph.nvtxs <= 0) {
	  printf("Empty graph. Nothing to do.\n");
	  exit(0);
	}

	int noOfSingletons = 0; 
	GraphType *noSingletonGraph ;
	idxtype* nodeMap = lookForSingletons(&graph, &noOfSingletons);
	if ( noOfSingletons > 0 )
	{
		getSubgraph(&graph, nodeMap, graph.nvtxs-noOfSingletons, 
						wgtflag, &noSingletonGraph);
		GKfree((void**)&(graph.xadj), (void**)&(graph.adjncy), LTERM);
		if ( wgtflag&1 > 0 )
			GKfree( (void**)&(graph.adjwgt), LTERM);
//		free(graph.gdata);
		printf("Found %d singleton nodes in the", noOfSingletons);
		printf(" input graph. Removing them.\n");
	}

	if ( !outputFileGiven )
	{
		strcpy(outputFile, filename);

		sprintf(outputFile,"%s.c%d.i%1.1f.b%1.1f",outputFile,opt.coarsenTo,opt.gamma,opt.penalty_power);
	}
	
	printf("Input graph information ---------------------------------------------------\n");
	printf("  Name: %s, #Vertices: %d, #Edges: %d\n", filename, graph.nvtxs, graph.nedges/2);
	printf("Output shall be placed in the file %s\n",
	outputFile);
	fflush(stdout);

	part = idxmalloc(graph.nvtxs, "main: part");

	printf("------------------------------------------------\n");
	printf("Clustering....\n");
	fflush(stdout);
	starttimer(METISTmr);         //YK: main algorithm starts here!

	if ( noOfSingletons > 0 )
	{
		
		mlmcl(&(noSingletonGraph->nvtxs), noSingletonGraph->xadj, noSingletonGraph->adjncy,
			noSingletonGraph->vwgt,noSingletonGraph->adjwgt, &wgtflag, part, opt ); 
	}
	else	
	{
		mlmcl(&graph.nvtxs, graph.xadj, graph.adjncy,graph.vwgt,
			graph.adjwgt, &wgtflag, part, opt ); 
	}

	stoptimer(METISTmr); 

	printf("------------------------------------------------\n");
	if ( noOfSingletons > 0 )
	{
		npart=mapPartition(part,noSingletonGraph->nvtxs);
		ncut=ComputeNCut(noSingletonGraph, part,npart);
//		printf("In graph that does not include singletons,");
//		printf("No. of Clusters: %d, N-Cut value: %.2f\n",npart,ncut);


		idxtype *clusterSizes = histogram(part,
					graph.nvtxs-noOfSingletons, npart);

		int maxSize = clusterSizes[idxamax(npart, clusterSizes)];
		float avgClusterSize =
						(graph.nvtxs-noOfSingletons)*1.0/(npart);
		float balance =	(maxSize*1.0) /
				((graph.nvtxs-noOfSingletons)*1.0/npart);
		float stdDevn = stdDeviation(clusterSizes, npart);
		float avgNcut = ncut * 1.0/npart;
		float normStdDevn = stdDevn/avgClusterSize;

	// Warning: This computation only works if the singletons
	// have been placed in their own clusters. This works for
	// MLR-MCL, in other words, because it is guaranteed to
	// place singletons in their own clusters.
		printf("Output statistics for graph without singletons\n");
		printf("Clusters: %d N-Cut: %.3f", 
					npart, ncut);
		printf(" AvgN-Cut: %.3f Balance in cluster sizes: %.2f ",avgNcut,
					balance); 
		printf("Std_Deviation in cluster sizes: %.2f ", stdDevn);
		printf("Coefficient_of_Variation: %.2f\n", normStdDevn);

		free( clusterSizes );

		npart += noOfSingletons;
	//	ncut += noOfSingletons;
		printf("Output statistics for original graph\n");

		mapIndices(part, nodeMap, graph.nvtxs, npart-noOfSingletons);
	}
	else
	{
		npart=mapPartition(part,graph.nvtxs);
		ncut=ComputeNCut(&graph, part,npart);
	}

	idxtype* clusterSizes = histogram(part, graph.nvtxs, npart);
	int maxSize = clusterSizes[idxamax(npart, clusterSizes)];
	float avgClusterSize = (graph.nvtxs)*1.0/(npart);
	float balance = (maxSize*1.0)/(graph.nvtxs*1.0/npart);
	float stdDevn = stdDeviation(clusterSizes, npart);
	float avgNcut = ncut * 1.0/npart;
	float normStdDevn = stdDevn/avgClusterSize;
	
	printf("Clusters: %d N-Cut: %.3f AvgN-Cut: %.3f", npart,
						ncut, avgNcut );
	printf(" Balance in cluster sizes: %.2f Std.Deviation in cluster sizes: %.2f ",
				 balance, stdDevn);
	printf("Coefficient_of_Variation: %.2f\n", normStdDevn);

	starttimer(IOTmr);
	my_WritePartition(outputFile, part, graph.nvtxs, opt.gamma); 
	if ( noOfSingletons > 0 )
	{
		free(nodeMap);
		nodeMap = NULL;
	}

	printf("\nOutput is written to file: %s\n", outputFile);
	stoptimer(IOTmr);
	stoptimer(TOTALTmr);
	
	printf("\nTiming Information --------------------------------------------------\n");
	printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
	printf("  Partitioning: \t\t %7.3f   (MLR-MCL time)\n", gettimer(METISTmr));
	printf("  Total:        \t\t %7.3f\n", gettimer(TOTALTmr));
	printf("**********************************************************************\n");


	GKfree((void**)&graph.xadj, (void**)&graph.adjncy, (void**)&graph.vwgt, 
				(void**)&graph.adjwgt, (void**)&part, LTERM);
}  


