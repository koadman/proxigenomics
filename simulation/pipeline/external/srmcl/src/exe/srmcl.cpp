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
#include <iostream>
#include <fstream>



void print_help(char program_name[])
{
	FILE *fp = stderr;
	fprintf(fp, "Usage: %s [options] graph_file\n", program_name);
	fprintf(fp, "Options:\n");
	//fprintf(fp, "\t-c <integer>\t\tsize of coarsest graph (default 1000)\n");
	fprintf(fp, "\t-b <float>\t\tbalance parameter in R-MCL (default 0.5) \n");
	fprintf(fp, "\t-i <float>\t\tinflation parameter in MCL (default 2.0) \n");
	fprintf(fp, "\t-w <string>\t\tquality threshold in SR-MCL (default 0 and the default function is density * sqrt(size) ) \n");
	fprintf(fp, "\t-p <string>\t\tredundacy threshold in SR-MCL (default 0.6) \n");
	fprintf(fp, "\t-q <string>\t\tpenalty ratio in SR-MCL (default 1.25) \n");
	fprintf(fp, "\t-t <string>\t\tnumber of times SR-MCL executing R-MCL (default: 30)\n");
	fprintf(fp, "\t-o <string>\t\toutput file\n");

	return;
}

int main(int argc, char *argv[])
{
	int i, npart;
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

			
			//the new parameters introduced in SR-MCL:
			case 'W':
			case 'w':
			  opt.quality_threshold = atof(*(++argv));
			  break;
			case 'p':
			case 'P':
			  opt.redundancy_threshold = atof(*(++argv));
			  break;
			//beta in SR-MCL
			case 'q':
			case 'Q':
			  opt.attractor_penalty_power = atof(*(++argv));
			  break;
			case 't':
			case 'T':
			  opt.time_rmcl = atoi(*(++argv));	
			  break;
			//case 'w':
			//case 'W':
		        //  opt.weighted_density = true; 
			//  break;



		
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
    
    // make a default output file name if none specified (mzd)
    if ( outputFileGiven == 0 ) {
        sprintf(outputFile, "%s.out", filename);
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

	ReadGraph(&graph, filename, &wgtflag, addSelfLoop, txtFormat);  //add self loop

	if ( opt.matchType == MATCH_UNSPECIFIED )
	{
//		opt.matchType = (graph.nvtxs>50000) ? MATCH_POWERLAW_FC :
//							MATCH_SHEMN;
		opt.matchType = MATCH_SHEMN;
		//printf("MATCH_SHEMN\n");
	}
	
	stoptimer(IOTmr);

	if (graph.nvtxs <= 0) {
	  printf("Empty graph. Nothing to do.\n");
	  exit(0);
	}

	int noOfSingletons = 0; 
	GraphType *noSingletonGraph ;
	idxtype* nodeMap = lookForSingletons(&graph, &noOfSingletons);

    // Added mzd - if the graph is trivial (completely unconnected)
    // just write out the trivial solution, which for metisCL format
    // is a blank line per node.
    if ( noOfSingletons == graph.nvtxs ) {
        printf("Graph is composed entirely of singletons, solution is trivial\n");
        std::ofstream outfs;
        outfs.open(outputFile);
        for (int i=0; i<graph.nvtxs; i++) {
            outfs << "" << std::endl;
        }
        outfs.close();
        exit(0);
    }

	if ( noOfSingletons > 0 )
	{
		getSubgraph(&graph, nodeMap, graph.nvtxs-noOfSingletons, 
						wgtflag, &noSingletonGraph);
		GKfree((void**)&(graph.xadj), (void**)&(graph.adjncy), LTERM);
		if ( (wgtflag&1) > 0 )
			GKfree( (void**)&(graph.adjwgt), LTERM);
//		free(graph.gdata);
		printf("Found %d singleton nodes in the", noOfSingletons);
		printf(" input graph. Removing them.\n");
	}


	printf("Input graph information ---------------------------------------------------\n");
	printf("  Name: %s, #Vertices: %d, #Edges: %d\n", filename, graph.nvtxs, graph.nedges/2);
	printf("Output shall be placed in the file %s\n", outputFile);
	fflush(stdout);




	std::vector<std::set<idxtype> > part(graph.nvtxs);

	for(int i=0; i<graph.nvtxs; i++)
		part[i] = std::set<idxtype>();
	//YK: end
	
	printf("------------------------------------------------\n");
	printf("Clustering....\n");
	fflush(stdout);
	starttimer(METISTmr);         //YK: main algorithm starts here!
	if ( noOfSingletons > 0 )
	{
		srmcl(&(noSingletonGraph->nvtxs), noSingletonGraph->xadj, noSingletonGraph->adjncy,
			noSingletonGraph->vwgt,noSingletonGraph->adjwgt, &wgtflag, part, opt);
        //part will be the solution. It contains cluster ids
	}
	else	
	{
		srmcl(&graph.nvtxs, graph.xadj, graph.adjncy,graph.vwgt,
			graph.adjwgt, &wgtflag, part, opt); 
	}

	stoptimer(METISTmr); 

	printf("------------------------------------------------\n");

	starttimer(IOTmr);

	my_WritePartition(outputFile, part, graph.nvtxs, opt.gamma);  //YK



	printf("\nOutput is written to file: %s\n", outputFile);
	stoptimer(IOTmr);
	stoptimer(TOTALTmr);
	
	printf("\nTiming Information --------------------------------------------------\n");
	printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
	printf("  Partitioning: \t\t %7.3f   (soft R-MCL time)\n", gettimer(METISTmr));
	printf("  Total:        \t\t %7.3f\n", gettimer(TOTALTmr));
	printf("**********************************************************************\n");


	//GKfree((void**)&graph.xadj, (void**)&graph.adjncy, (void**)&graph.vwgt, (void**)&graph.adjwgt, (void**)&part, LTERM);
}  


