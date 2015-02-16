package org.darlinglab.metaphazer3C.operator;

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
		for(int i=0; i < g.states.length; i++){
			// Gibbs sample a new state for site i
			// calculate likelihood for each nucleotide
			double[] state_ll = new double[SnvGraph.ALPHABET_SIZE];
			for(int a=0; a<SnvGraph.ALPHABET_SIZE; a++){
				g.states[i] = a;
				state_ll[a] = gml.getLikelihood();
			}
			// now sample one
		}
	}
	
}
