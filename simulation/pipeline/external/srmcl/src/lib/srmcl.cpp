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
//#define TEST_OUTPUT
//#define CLUSTER_COEFFICIENT

class Cluster{
	public:
		double density;
		int cID;
		int size;
		Cluster(int, int, double);
};
Cluster::Cluster(int cID, int size, double density){
	this->cID = cID;
	this->density = density;
	this->size = size;
}

int compareCluster (const void * a, const void * b)
{
  double diff = (*(Cluster*)b).density*sqrt((double)(*(Cluster*)b).size)  -  (*(Cluster*)a).density * sqrt((double) (*(Cluster*)a).size);
//  double diff = (*(Cluster*)b).density  -  (*(Cluster*)a).density ;
  if( diff > 0)
	return true;
  else if( diff < 0)
	return false;
  else
	if((*(Cluster*)b).size  -  (*(Cluster*)a).size > 0)
		return true;
	else
		return false;
}

std::vector< std::set<idxtype> > srmclWithGraph(GraphType *graph,std::vector<std::set<idxtype> >& indices,Options opt )
{
	GraphType *cgraph,*t,*cgraphFinest, *cgraphCoarest;
	int nvtxs=graph->nvtxs;
	CtrlType ctrl;
	int levels=0;
	Matrix *M=NULL;

	int ct = ctrl.CoarsenTo = opt.coarsenTo, runMcl=0;
	ctrl.CType=opt.matchType;

	ctrl.optype=OP_KMETIS;
	ctrl.dbglvl=1;

	if ( ct < 5 )
		ctrl.maxvwgt = (int) ceil( (1.5*nvtxs)/ct );
	else if ( ct <= 100 )
		ctrl.maxvwgt = (int) ceil( (2.0*nvtxs)/ct );
	else if ( ct > 100 )
	{
		// we can allow some imbalance here.
		ctrl.maxvwgt = (int) ceil( 10.0 * nvtxs/ct );
	}

	my_AllocateWorkSpace(&ctrl,graph);

	cgraphCoarest = cgraph = Coarsen2Way(&ctrl,graph);
	

	FreeWorkSpace(&ctrl, graph);


	for(levels=0,t=cgraph ; t->finer!=NULL ; t=t->finer,levels++);

	int num_level = levels;
	printf("Coarsened %d levels\n", levels);
	printf("Number of vertices in coarsest graph:%d\n", cgraph->nvtxs);
	fflush(stdout);


	int times_rmcl= 1;	

	while ( cgraph != NULL )
	{
		Matrix *M0, *Mnew;
		int iters;

		M0=setupCanonicalMatrix(cgraph->nvtxs, cgraph->nedges, cgraph->xadj, cgraph->adjncy, cgraph->adjwgt, opt.ncutify);  //do weight normalization here
#ifdef TEST_OUTPUT
		//printf("level: %d, M0:\n",levels);
		//printMatrix(M0);
#endif


		if ( cgraph->coarser == NULL )
		{
		// if this is the coarsest graph, flow matrix is same as
		// the transition matrix. 
			M=M0;
	//		dumpMatrix(M);
		}

		if (cgraph->finer == NULL){
			iters=opt.num_last_iter;
			cgraphFinest = cgraph;
		}
		else{
			iters=opt.iter_per_level;
		}

		//printMatrix(M);
		Mnew=dprmcl(M,M0,cgraph, opt, iters, levels);        //YK: main process!
		M=Mnew;

		if ( cgraph->finer != NULL )  //YK: map nodes to the finer graph
		{
			int nnzRefinedGraph = 2*M->nnz;
			if ( ctrl.CType == MATCH_POWERLAW_FC )
			{
				int ii;
				nnzRefinedGraph = 0;
				for ( ii=0; ii<cgraph->finer->nvtxs; ii++ )
				{
					int tx=cgraph->finer->cmap[ii];
					nnzRefinedGraph +=
					M->xadj[tx+1]-M->xadj[tx];
				}
			}
			else if ( ctrl.CType == MATCH_HASH  && levels <= 3 )
			{
				int ii;
				nnzRefinedGraph=0;
				for ( ii=0; ii<cgraph->nvtxs; ii++ )
				{
					nnzRefinedGraph += cgraph->vwgt[ii] *
							(M->xadj[ii+1] - M->xadj[ii]);
				}
			}
//			dumpMatrix(M);
			Mnew=propagateFlow(M,cgraph,cgraph->finer,nnzRefinedGraph);
			//printf("times:%d, levels:%d Mnew\n",times_rmcl, levels);
			//printMatrix(Mnew);
			if ( M != NULL )
			{
				freeMatrix(M);
			}
			M=Mnew;
		}
		cgraph=cgraph->finer;  //change to the finer graph


		levels--;
		printf("Done level %d (%d times MLR-MCL)\n", levels+1,times_rmcl);
		fflush(stdout);
	}
	fflush(stdout);
	idxtype* firstIndices =  (idxtype*) malloc(sizeof(idxtype) * M->nvtxs);
	idxcopy(nvtxs, M->attractors, firstIndices);


	bool* isAttractor = (bool*) malloc(sizeof(bool) * M->nvtxs);  
	getAttractorsForAll(M, isAttractor);
	int** countAttractor = (int**) malloc(sizeof(int*) * (num_level+1));  //count the number of times the node being an attractor node
	for(int i=0; i <= num_level ; i++){
		countAttractor[i] = (int*) malloc (sizeof(int) * nvtxs);
		for(int j=0; j<M->nvtxs; j++){
			countAttractor[i][j] = 0;
		}
	}
#ifdef TEST_OUTPUT
	{
		for(int i=0; i<M->nvtxs; i++){
			int level = 1, nodeIdx = i;
			for(cgraph = cgraphFinest; cgraph->coarser!=NULL; cgraph = cgraph->coarser, level++){
				idxtype oldNodeIdx = nodeIdx;
				nodeIdx = cgraph->cmap[oldNodeIdx];
				printf("  cmap: level:%d  vid: %d to %d \n", level, oldNodeIdx+1, nodeIdx+1);
			}
		}
	}
#endif
	for(int i=0; i<M->nvtxs; i++){
		
		if(isAttractor[i]){
			countAttractor[0][i]++;
#ifdef TEST_OUTPUT
			printf("attractor vid:%d(%d times), %d-th times MLRMCL\n",i+1,countAttractor[0][i],times_rmcl);fflush(stdout);
#endif
			int level = 1, nodeIdx = i;
			for(cgraph = cgraphFinest; cgraph->coarser!=NULL; cgraph = cgraph->coarser, level++){
				idxtype oldNodeIdx = nodeIdx;
				nodeIdx = cgraph->cmap[nodeIdx];

				countAttractor[level][nodeIdx]++;
#ifdef TEST_OUTPUT
				//printf("countAttractor[%d][%d] = %d;\n", level, nodeIdx,countAttractor[level][nodeIdx]);
#endif
			}
			
		}
	}
	
	//construct hash table
	idxtype* hashtable = idxmalloc(nvtxs,"Hashtable");
	for(int i=0;i<nvtxs;i++)
		hashtable[i]=-1;  	// clear hashtable.
	int num_clusters = 0;  //include singleton
	for(int i=0;i<nvtxs;i++){
		if ( hashtable[firstIndices[i]] == -1){
			hashtable[firstIndices[i]] = num_clusters;
			num_clusters++;
		}
	}
	/*idxtype* invHashTable = idxmalloc(num_clusters,"inverse hashtable");  //map cluster id back to the attractor id
	for(int i=0;i<nvtxs;i++){
		if ( hashtable[firstIndices[i]] != -1){
			invHashTable[hashtable[firstIndices[i]]] = firstIndices[i];
		}
	}*/
	//construct clusters
	std::vector< std::set<idxtype> > clusters(0);   //clusters[i] contains nodes' numbers (in original graph) in cluster i
	for(int i=0; i<num_clusters ; i++){
		std::set<idxtype> newCluster;
		clusters.push_back(newCluster);
	}
	for(int i=0; i<nvtxs ; i++){
		if( firstIndices[i] != -1){
			idxtype cID = hashtable[firstIndices[i]];
			indices[i].insert(cID); 
			if(cID == -1)
				continue;
			clusters[ cID  ].insert( i );
		}
	}
	idxtype* lastIndices = firstIndices;
	bool* lastAttractor = isAttractor;
	idxtype* oldHashtable = hashtable;

	printf("number of clusters in 1st time MLR-MCL: %d\n", num_clusters);fflush(stdout);

	idxtype* thisIndices =  (idxtype*) malloc(sizeof(idxtype) * M->nvtxs);
	idxtype* newHashtable = idxmalloc(nvtxs,"newHashtable");
	bool* thisAttractor = (bool*) malloc(sizeof(bool) * M->nvtxs);  
	while(times_rmcl < opt.time_rmcl){  //do more than one times R-MCL but penalize attractor nodes (by assigning them higher inflation rate)

		times_rmcl++;
		cgraph = cgraphCoarest;
		for(levels=0,t=cgraph ; t->finer!=NULL ; t=t->finer,levels++);

		while ( cgraph != NULL )
		{
			Matrix *M0, *Mnew;
			int iters;

			M0=setupCanonicalMatrix(cgraph->nvtxs, cgraph->nedges, cgraph->xadj, cgraph->adjncy, cgraph->adjwgt, opt.ncutify);  //do weight normalization here

			if ( cgraph->coarser == NULL )
			{
				M=M0;
			}
			

			if (cgraph->finer == NULL)
				iters=opt.num_last_iter;
			else
				iters=opt.iter_per_level;

			//Mnew=dprmcl(M,M0,cgraph, opt, iters,levels);  //original
			Mnew=dprmcl_penalizeAttractors(M,M0,cgraph, opt, iters, levels, countAttractor);        //YK: main process!
			
			M=Mnew;

			if ( cgraph->finer != NULL )  //YK: map nodes to the finer graph
			{
				int nnzRefinedGraph = 2*M->nnz;
				if ( ctrl.CType == MATCH_POWERLAW_FC )
				{
					int ii;
					nnzRefinedGraph = 0;
					for ( ii=0; ii<cgraph->finer->nvtxs; ii++ )
					{
						int tx=cgraph->finer->cmap[ii];
						nnzRefinedGraph +=
						M->xadj[tx+1] - M->xadj[tx];
					}
				}
				else if ( ctrl.CType == MATCH_HASH  && levels <= 3 )
				{
					int ii;
					nnzRefinedGraph=0;
					for ( ii=0; ii<cgraph->nvtxs; ii++ )
					{
						nnzRefinedGraph += cgraph->vwgt[ii] *
								(M->xadj[ii+1] - M->xadj[ii]);
					}
				}
				Mnew=propagateFlow(M,cgraph,cgraph->finer,nnzRefinedGraph);
				//printf("times:%d, levels:%d Mnew\n",times_rmcl, levels);
				//printMatrix(Mnew);
				if ( M != NULL )
				{
					freeMatrix(M);
				}
				M=Mnew;
			}
			cgraph=cgraph->finer;  //change to the finer graph

			// These two didn't get freed earlier, when we freed
			// gdata.
			

			levels--;
			printf("Done level %d (%d-th times MLR-MCL)\n", levels+1,times_rmcl);fflush(stdout);

		}
		//update clusters
		
		idxcopy(nvtxs, M->attractors, thisIndices);

		
		getAttractorsForAll(M, thisAttractor);
		for(int i=0; i<M->nvtxs; i++){
			if(thisAttractor[i]){
				countAttractor[0][i]++;
	#ifdef TEST_OUTPUT
					printf("attractor vid:%d(%d times), %d-th times MLRMCL\n",i+1,countAttractor[0][i],times_rmcl);
	#endif
				int level = 1, nodeIdx = i;
				for(cgraph = cgraphFinest; cgraph->coarser!=NULL; cgraph = cgraph->coarser, level++){
					nodeIdx = cgraph->cmap[nodeIdx];
					countAttractor[level][nodeIdx]++;
				}
				//printf("attractor:%d, times:%d\n",i+1,times_rmcl);
			}
		}
	
		//construct hash table

		for(int i=0;i<nvtxs;i++)
			newHashtable[i]=-1;  	// clear newHashtable.

		int new_num_clusters = 0 ;
		for(int i=0;i<nvtxs;i++){
			if ( newHashtable[thisIndices[i]] == -1){
				newHashtable[thisIndices[i]] = num_clusters+new_num_clusters;  //cluster ID
				if(num_clusters+new_num_clusters == 22704 || num_clusters+new_num_clusters == 29925)
					printf("cluster #%d's attractor node %d\n",num_clusters+new_num_clusters,thisIndices[i]);
				new_num_clusters++;
			}
		}
		//construct clusters
		for(int i=0; i<new_num_clusters;  i++){
			std::set<idxtype>* newCluster = new std::set<idxtype>();
			clusters.push_back(*newCluster);
		}
		for(int i=0; i<nvtxs ; i++){
			if( thisIndices[i] != -1){
				idxtype cID = newHashtable[thisIndices[i]];
				indices[i].insert(cID); 
				if(cID == -1)
					continue;
				clusters[ cID ].insert( i );
			}
		}
		printf("number of clusters in %dth times MLR-MCL: %d\ttotal # clusters:%zu\n", times_rmcl, new_num_clusters,clusters.size());
		num_clusters += new_num_clusters;

		//detect last attractor is split to multiple clusters
		/*
		double* countToCluster = (double*) malloc(sizeof(double) * num_clusters) ; 
		for(idxtype vID=0; vID<nvtxs ; vID++){
			
			//printf("test:%d\n",vID+1);fflush(stdout);
			if(lastAttractor[vID] && !thisAttractor[vID]){
				for(int i=0; i<num_clusters; i++)
					countToCluster[i] = 0;
				int lastCluster = oldHashtable[lastIndices[vID]];
				//printf("lastCluster:cid%d (attractor:%d size:%zu), num clusters:%d\n",lastCluster,vID+1,clusters[lastCluster].size(),num_clusters);fflush(stdout);
				
				for(std::set<idxtype>::iterator nodeIterator = clusters[lastCluster].begin(); nodeIterator != clusters[lastCluster].end(); nodeIterator++){
					
					int currentCluster = newHashtable[thisIndices[*nodeIterator]];
					countToCluster[currentCluster] ++;
				}
				for(idxtype cID = 0; cID < num_clusters; cID++){
					if(countToCluster[cID] / clusters[lastCluster].size() >= opt.ratio_nodes_another_cluster && clusters[cID].find(vID)==clusters[cID].end()  ){
						//clusters[ cID ].insert(vID);
						//indices[ vID ].insert(cID);
#ifdef TEST_OUTPUT
						printf("^add extra cluster: vid:%d to cid:%d\n",vID+1, cID);fflush(stdout);
#endif
					}
				}
				
				
			}
		}	
		free(countToCluster);*/
		

		lastIndices = thisIndices;
		lastAttractor = thisAttractor;
		oldHashtable = newHashtable;
		
	}
	printf("total number of clusters after %d times MLR-MCL:%zu\n", opt.time_rmcl, clusters.size());fflush(stdout);

	free(hashtable);
	for(int i=0; i <= num_level ; i++)
		free(countAttractor[i]);
	free(countAttractor);
	free(isAttractor);
	free(firstIndices);
	free(thisIndices);
	free(newHashtable);
	free(thisAttractor);
	freeMatrix(M);
	return clusters;
}


void srmcl(int* nvtxs, idxtype* xadj, idxtype* adjncy, idxtype
*vwgt, idxtype* adjwgt, int* wgtflag, std::vector<std::set<idxtype> >& indices, Options opt)
{
	



	GraphType *graph = (GraphType*)malloc(sizeof(GraphType));
	my_SetUpGraph(graph, *nvtxs, xadj, adjncy, vwgt, adjwgt, *wgtflag, 1);


	std::vector< std::set<idxtype> > clusters = srmclWithGraph(graph, indices, opt);   //main procedure!


	//post processing: prune out clusters
	printf("start post-processing - prune clusters.\n");
	int num_clusters = clusters.size();

	double max_weight = 10000;
//printf("double max_weight = 10000;\n");
	//prune out clusters according to their clustering coefficient
#ifdef CLUSTER_COEFFICIENT
	double* cluster_coefficients = (double*) malloc(sizeof(double) * num_clusters);  
	//weighted cc is according to [B. Zhang and S. Horvath, Stat. App. Genet. Mol. Biol. 4, 17 2005.]
	for(int vIdx = 0; vIdx < *nvtxs; vIdx++){   //v1
		for(std::set<idxtype>::iterator cID = indices[vIdx].begin(); cID != indices[vIdx].end(); cID ++){ //for each cluster, calculate its clustering coefficient
			double numerator = 0;
			double denominator1 = 0;
			double denominator2 = 0;
			for(int adjIdx = xadj[vIdx]; adjIdx < xadj[vIdx+1]; adjIdx++){
				if(indices[adjncy[adjIdx]].find(*cID) != indices[adjncy[adjIdx]].end() && adjncy[adjIdx] != vIdx){  //if another node v1 is also in this cluster 
					double wij = (double)adjwgt[adjIdx] / max_weight;
					denominator1 += wij;
					denominator2 += pow(wij,2.0);
					for(int adjIdx2 = adjIdx+1; adjIdx2 < xadj[vIdx+1]; adjIdx2++){
						if(indices[adjncy[adjIdx2]].find(*cID) != indices[adjncy[adjIdx2]].end()  && adjncy[adjIdx2] != vIdx){ //if another node v2 is also in this cluster 
							double wik = (double)adjwgt[adjIdx2] / max_weight;
							for(int i=xadj[adjncy[adjIdx]]; i<xadj[adjncy[adjIdx]+1]; i++){ //find whether v1 and v2 are connected
								if(adjncy[i] == adjncy[adjIdx2]){  //if v1 and v2 are connected
									double wjk =  (double)adjwgt[i] / max_weight;
									numerator += wij*wik*wjk;
									//if(*cID == 3)
										
									break;
								}
							}
						}
					}
				}
			}
			denominator1 = pow(denominator1, 2.0);
			double denominator = denominator1 - denominator2;			
			numerator *= 2;
			double cluster_coefficient = 0; 
			if(numerator != 0)
				cluster_coefficient = numerator / denominator;
			//printf("cID:%d node:%d value:%.3f (%.3f/%.3f)\n",*cID, vIdx+1,cluster_coefficient,numerator,denominator );
			cluster_coefficients[*cID] += cluster_coefficient;
		}
	}
#endif

	//calculate weighted density


	double* num_internal_edges = (double*) malloc( sizeof(double) * num_clusters);
	for(int i=0; i<num_clusters; i++)
		num_internal_edges[i] = 0;
	for(int vID=0; vID<*nvtxs ; vID++){
		for(std::set<idxtype>::iterator cIterator = indices[vID].begin(); cIterator != indices[vID].end(); cIterator++){
			idxtype cID = *(cIterator);

			for( int j = xadj[vID]; j < xadj[vID+1]; j++ ){
				if( adjncy[j]!=vID && clusters[cID].find( adjncy[j] ) != clusters[cID].end()  ){  //contains it
					if(opt.weighted_density)
						num_internal_edges[ cID ] += (adjwgt[j] / max_weight);   //weighted version
					else
						num_internal_edges[ cID ] ++; 		
				}
			}
		}

	}

//printf("double* densities = (double*) malloc(sizeof(double) * (*nvtxs));  \n"); fflush(stdout);
	double* densities = (double*) malloc(sizeof(double) * num_clusters);  
	for(int cID = 0; cID < num_clusters ; cID++){
		int size = clusters[cID].size();
		if(size <= 1)
			densities[cID] = 0;
		else
			densities[cID] = num_internal_edges[ cID ] / size / (size-1);
	}
	free(num_internal_edges);


	int num_pruned_clusters_density = 0;	
	for(int cID = 0; cID < num_clusters ; cID++){
#ifdef CLUSTER_COEFFICIENT
		if(clusters[cID].size() > 1){
			cluster_coefficients[cID] /= clusters[cID].size();
			densities[cID] = cluster_coefficients[cID] ;
		}
#endif

#ifdef TEST_OUTPUT	
		printf("cluster %d: \tdensity:%.3f\tsize:%zu\n",cID, densities[cID], clusters[cID].size());
#endif
		if(clusters[cID].size()<=2 || densities[cID] * sqrt((double)clusters[cID].size()) < opt.quality_threshold){  //remove the cluster
		//if(clusters[cID].size()<=2 || densities[cID]  < opt.quality_threshold){  //remove the cluster
			for(std::set<idxtype>::iterator nodeIterator = clusters[cID].begin(); nodeIterator != clusters[cID].end(); nodeIterator++){
				//printf(" nodeIterator:%d\n",*nodeIterator);
				if(indices[*nodeIterator].find(cID) != indices[*nodeIterator].end()){
					//printf("  if nodeIterator:%d\n",*(indices[*nodeIterator].find(cID)));
					indices[*nodeIterator].erase(cID);
				}
				//printf(" nodeIterator:%d done\n",*nodeIterator);
			}
			num_pruned_clusters_density++;

			clusters[cID].clear();
			densities[cID] = -1;
#ifdef CLUSTER_COEFFICIENT
			cluster_coefficients [cID] = -1;
#endif
		}
	}
	printf("number of clusters pruned out since their density are smaller than %.3f:\t%d\n", opt.quality_threshold, num_pruned_clusters_density);

		
	//prune out clusters according to their redundancy (sort clusters by density * sqrt(size))
	int num_pruned_clusters_overlap = 0;
	int num_clusters_after_pruning_density = (num_clusters-num_pruned_clusters_density);
	printf("sort clusters\n");
//printf("malloc test: %d, sizeof(Cluster):%d (int:%d)\n",num_clusters_after_pruning_density,sizeof(Cluster),sizeof(int)); fflush(stdout);
	Cluster* cluster_array = (Cluster*) malloc(sizeof(Cluster)*(num_clusters_after_pruning_density ));
//printf("malloc test2: %d\n",num_clusters_after_pruning_density); fflush(stdout);
	int tempIdx = 0;
	for(int i=0; i < num_clusters; i++){
		if(densities[i] >= 0){
			Cluster c(i,clusters[i].size(),	densities[i] );
			cluster_array[tempIdx] = c;
			tempIdx++;
		}
	}

	double avg_cluster_size = 0;
	qsort(cluster_array, num_clusters_after_pruning_density  , sizeof(Cluster), compareCluster);
	for(int i=0; i<num_clusters_after_pruning_density   ; i++){
		int cID1 = cluster_array[i].cID;
		if(clusters[cID1].size() == 0)
			continue;
		//printf("examined cluster cID1: %d\n",cID1);fflush(stdout);	
		for(int j=i+1; j<num_clusters_after_pruning_density ; j++){
			
			//calculate overlap size
			int cID2 = cluster_array[j].cID;
			//printf(" examined cluster cID2: %d\n",cID2);fflush(stdout);	
			if(clusters[cID2].size() == 0)
				continue;	
			double overlap = 0;

			for ( std::set<idxtype>::iterator iterator = clusters[cID2].begin(); iterator != clusters[cID2].end(); iterator++ ){
				
				if ( clusters[cID1].find(*iterator) != clusters[cID1].end() ){
					overlap++;
				}
					
			}
			//printf("overlap: %1.0f\n",overlap);fflush(stdout);	
			//calculate neighbor affinity
			float overlapNA = pow(overlap,2) / clusters[cID1].size() / clusters[cID2].size();
			if(overlapNA >= opt.redundancy_threshold){  //remove the cluster with cID2
				
				for(std::set<idxtype>::iterator nodeIterator = clusters[cID2].begin(); nodeIterator != clusters[cID2].end(); nodeIterator++){

					indices[*nodeIterator].erase(cID2);
					
				}
				num_pruned_clusters_overlap++;
#ifdef TEST_OUTPUT
				printf(" overlapNA:%.3f, (keep cID:%d, remove cID:%d) test size:%1f, %d, %d\n",overlapNA,cID1,overlap, cID2,clusters[cID1].size(),clusters[cID2].size());fflush(stdout);
#endif
				clusters[cID2].clear();
				densities[cID2] = 0;
			}
			
		}
		avg_cluster_size += clusters[cID1].size();
	}
	printf("number of overlapped clusters pruned out since their NA are larger than %.3f:\t%d\n", opt.redundancy_threshold, num_pruned_clusters_overlap);fflush(stdout);

	num_clusters -=  (num_pruned_clusters_density + num_pruned_clusters_overlap);
	printf("*********afer prunning, total # clusters:\t%d**********\n",num_clusters);fflush(stdout);

	avg_cluster_size /= num_clusters;
	printf("*********average cluster size:\t%.3f**********\n",avg_cluster_size);fflush(stdout);

	int coverage = 0;
	for(int vID=0; vID<*nvtxs ; vID++){
		if(indices[vID].size() > 0)
			coverage++;
	}
	printf("*********coverage:\t%d**********\n",coverage);fflush(stdout);

	free(cluster_array);
	free(graph);
}



