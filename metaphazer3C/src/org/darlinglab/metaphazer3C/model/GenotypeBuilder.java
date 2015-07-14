package org.darlinglab.metaphazer3C.model;

import org.darlinglab.metaphazer3C.Random;
import org.darlinglab.metaphazer3C.data.SnvGraph;

public class GenotypeBuilder {
	public static Genotype randomGenotype(SnvGraph snps){
		Genotype g = new Genotype();
		g.states = new int[snps.nodes.length];
		
		for(int i=0; i<g.states.length; i++){
			g.states[i] = Random.sampleCategorical(snps.nodes[i].freqs);
		}
		return g;
	}
}
