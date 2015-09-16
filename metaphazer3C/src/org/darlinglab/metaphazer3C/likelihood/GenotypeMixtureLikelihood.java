package org.darlinglab.metaphazer3C.likelihood;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Stack;

import org.darlinglab.metaphazer3C.data.SnvGraph;
import org.darlinglab.metaphazer3C.model.GenotypeMixture;
import org.darlinglab.metaphazer3C.model.ModelListener;

/**
 * Implements a likelihood calculator for a mixture of genotypes in HiC/3C data
 * Uses the multinomial density to calculate likelihood of observed data given model genotypes
 * @author koadman
 *
 */
public class GenotypeMixtureLikelihood implements ModelListener {
	public GenotypeMixtureLikelihood(GenotypeMixture gm, SnvGraph snvs){
		this.gm = gm;
		this.snvs = snvs;
		this.dirty = true;
		
		// pre-calculate the multinomial density coefficients for each pair of sites
		// these don't change unless the data changes, and are slow to calculate, so just do it once up-front
		this.log_md_coefficient = new double[snvs.nodes.length][snvs.nodes.length];
		for(SnvGraph.Edge e : snvs.edges){
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
	
	public void modelChanged(Object o){
		if(o instanceof SnvGraph.Node){
			for(SnvGraph.Edge e : ((SnvGraph.Node)o).edges){
				dirtyList.add(e);
			}
		}
	}
	
	public void updateLikelihood(){
		if(dirty){
			// do a complete recalculation
			ll = 0;
			for(SnvGraph.Edge e : snvs.edges){
				edge_ll.put(e, multinomialDensity(e));
				ll += edge_ll.get(e);
			}
			dirtyList.clear();
			dirty = false;
		}else{
			// just recalculate items on the dirty list
			for(SnvGraph.Edge e : dirtyList){
				ll -= edge_ll.get(e);
				edge_ll.put(e, multinomialDensity(e));
				ll += edge_ll.get(e);
			}
			dirtyList.clear();
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
			link_freqs[gm.genotypes[g].getState(e.node_1)][gm.genotypes[g].getState(e.node_2)] += gm.abundances[g];
			// add in sequencing error 
			link_freqs[gm.genotypes[g].getState(e.node_1)][gm.genotypes[g].getState(e.node_2)] += gm.snv_error_rate;
			// add in chimera error. calculate as f_i * f_j * error_rate
			link_freqs[gm.genotypes[g].getState(e.node_1)][gm.genotypes[g].getState(e.node_2)] += gm.abundances[e.node_1] * gm.abundances[e.node_2] * gm.ligation_error_rate;
			// normalize back to a prob. distribution
			normalizer += link_freqs[gm.genotypes[g].getState(e.node_1)][gm.genotypes[g].getState(e.node_2)];
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
		updateLikelihood();
		return ll;
	}
	
	double[][] log_md_coefficient; /**< stores the ratio of Gammas coefficient in the multinomial density */
	HashMap<SnvGraph.Edge,Double> edge_ll; /**< the current log likelihood for each edge in the SNP graph */ 
	boolean dirty = false;
	HashSet<SnvGraph.Edge> dirtyList = new HashSet<SnvGraph.Edge>(); /**< things that have changed since */
	double ll = 0; /**< the log likelihood */
	GenotypeMixture gm;
	SnvGraph snvs;

}

