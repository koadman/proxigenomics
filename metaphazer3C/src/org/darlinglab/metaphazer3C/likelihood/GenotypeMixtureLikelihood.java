package org.darlinglab.metaphazer3C.likelihood;

import java.util.Stack;

import org.darlinglab.metaphazer3C.data.SnvGraph;
import org.darlinglab.metaphazer3C.model.GenotypeMixture;

/**
 * Implements a likelihood calculator for a mixture of genotypes in HiC/3C data
 * Uses the multinomial density to calculate likelihood of observed data given model genotypes
 * @author koadman
 *
 */
public class GenotypeMixtureLikelihood {
	public GenotypeMixtureLikelihood(GenotypeMixture gm, SnvGraph snps){
		this.gm = gm;
		this.snvs = snps;
		for(int i=0; i < snps.edges.length; i++){
			dirtyList.push(new Integer(i));			
		}
		
		// pre-calculate the multinomial density coefficients for each pair of sites
		// these don't change unless the data changes, and are slow to calculate, so just do it once up-front
		log_md_coefficient = new double[gm.genotypes[0].states.length][gm.genotypes[0].states.length];
		for(SnvGraph.Edge e : snps.edges){
			double lognumer = 0;
			double logdenom = 0;
			for(int j=0; j<e.counts.length; j++){					
				for(int k=0; k<e.counts[j].length; k++){
					lognumer += e.counts[j][k];
					logdenom += logGamma(e.counts[j][k]+1);
				}
			}
			lognumer = logGamma((int)lognumer+1);
			log_md_coefficient[e.node_1][e.node_2] = lognumer - logdenom;
		}
	}

	/** calculate the natural log of the Gamma function.
	 *  consider this approximate due to accumulation of floating point error */
	private double logGamma(int n){
		double lf = 0;
		for (int i = 2; i < n; i++) {
		  lf += Math.log(i);
		}
		return lf;
	}
	
	public void updateLikelihood(){
		if(dirty){
			// complete recalculation
			ll = 0;
			for(int e=0; e<snvs.edges.length; e++){
				edge_ll[e] = multinomialDensity(snvs.edges[e]);
				ll += edge_ll[e];
			}
			dirtyList.clear();
			dirty = false;
		}else{
			while(!dirtyList.isEmpty()){
				Integer d = dirtyList.pop();
				ll -= edge_ll[d];
				edge_ll[d] = multinomialDensity(snvs.edges[d]);
				ll += edge_ll[d];
			}
		}
	}
	
	/** calculate the multinomial density for an edge
	 * 
	 * @param e  The edge
	 * @return   the log likelihood
	 */
	public double multinomialDensity(SnvGraph.Edge e){
		double density = 0;

		// count up the frequencies of each link type in the genotypes
		float[][] link_freqs = new float[SnvGraph.ALPHABET_SIZE][SnvGraph.ALPHABET_SIZE];
		double normalizer = 0;
		for(int g = 0; g < gm.genotypes.length; g++){
			link_freqs[gm.genotypes[g].states[e.node_1]][gm.genotypes[g].states[e.node_2]] += gm.abundances[g];
			// add in sequencing error 
			link_freqs[gm.genotypes[g].states[e.node_1]][gm.genotypes[g].states[e.node_2]] += gm.snv_error_rate;
			// add in chimera error. calculate as f_i * f_j * error_rate
			link_freqs[gm.genotypes[g].states[e.node_1]][gm.genotypes[g].states[e.node_2]] += gm.abundances[e.node_1] * gm.abundances[e.node_2] * gm.ligation_error_rate;
			// normalize back to a prob. distribution
			normalizer += link_freqs[gm.genotypes[g].states[e.node_1]][gm.genotypes[g].states[e.node_2]];
		}

		
		// calculate log like
		for(int j=0; j < e.counts.length; j++){
			for(int k=0; k < e.counts[j].length; k++){
				link_freqs[j][k] /= normalizer;
				density += log_md_coefficient[e.node_1][e.node_2] + Math.log(Math.pow(link_freqs[j][k], e.counts[j][k]));
			}
		}
		
		return density;
	}
	
	public double getLikelihood(){
		if(dirty || !dirtyList.isEmpty()){
			updateLikelihood();
		}
		return ll;
	}
	
	double[][] log_md_coefficient; /**< stores the ratio of Gammas coefficient in the multinomial density */
	double[] edge_ll; /**< the current log likelihood for each edge in the SNP graph */ 
	boolean dirty = false;
	Stack<Integer> dirtyList = new Stack<Integer>(); /**< things that have changed since */
	double ll = 0; /**< the log likelihood */
	GenotypeMixture gm;
	SnvGraph snvs;

}

