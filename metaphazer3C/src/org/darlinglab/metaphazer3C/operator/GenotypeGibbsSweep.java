package org.darlinglab.metaphazer3C.operator;

import org.darlinglab.metaphazer3C.Random;
import org.darlinglab.metaphazer3C.data.SnvGraph;
import org.darlinglab.metaphazer3C.likelihood.GenotypeMixtureLikelihood;
import org.darlinglab.metaphazer3C.model.Genotype;

/**
 * An operator to change the genotype in the current model
 * @author koadman
 *
 */
public class GenotypeGibbsSweep {
	GenotypeMixtureLikelihood gml;
	SnvGraph snvs;
	
	public GenotypeGibbsSweep(GenotypeMixtureLikelihood gml, SnvGraph snvs){
		this.gml = gml;
		this.snvs = snvs;
	}

	public void move(Genotype g){
		for(int i=0; i < snvs.nodes.length; i++){
			// Gibbs sample a new state for site i
			// calculate likelihood for each nucleotide
			double[] state_ll = new double[SnvGraph.ALPHABET_SIZE];
			double max_ll = -1e-300;
			for(int a=0; a<SnvGraph.ALPHABET_SIZE; a++){
				g.setState(i, a);
				state_ll[a] = gml.getLikelihood();
				max_ll = state_ll[a] > max_ll ? state_ll[a] : max_ll;
			}
			// now convert to a sampling distribution
			double ll_sum = 0;
			for(int a=0; a<SnvGraph.ALPHABET_SIZE; a++){
				state_ll[a] += max_ll;
				state_ll[a] = Math.exp(state_ll[a]);
				ll_sum += state_ll[a];
			}
			// normalize
			for(int a=0; a<SnvGraph.ALPHABET_SIZE; a++){
				state_ll[a] /= ll_sum;
			}
			// sample
			int new_state = Random.sampleCategorical(state_ll);
			g.setState(i, new_state);
			// done!
		}
	}
	
}
