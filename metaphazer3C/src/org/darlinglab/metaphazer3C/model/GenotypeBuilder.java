package org.darlinglab.metaphazer3C.model;

import org.darlinglab.metaphazer3C.Random;
import org.darlinglab.metaphazer3C.data.SnvGraph;

public class GenotypeBuilder {
	public static Genotype randomGenotype(SnvGraph snvs){
		Genotype g = new Genotype();
		g.states = new int[snvs.nodes.length];
		
		final float[] equalFrequencies = new float[SnvGraph.ALPHABET_SIZE];
		for(int i=0; i<equalFrequencies.length; i++){
			equalFrequencies[i] = 1.0f / ((float)SnvGraph.ALPHABET_SIZE);
		}
		for(int i=0; i<g.states.length; i++){
			g.states[i] = Random.sampleCategorical(equalFrequencies);
		}
		return g;
	}
}
