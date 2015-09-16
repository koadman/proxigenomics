package org.darlinglab.metaphazer3C;

import java.io.File;
import java.io.IOException;

import org.darlinglab.metaphazer3C.data.SnvGraph;
import org.darlinglab.metaphazer3C.data.SnvGraphBuilder;
import org.darlinglab.metaphazer3C.likelihood.GenotypeMixtureLikelihood;
import org.darlinglab.metaphazer3C.model.Genotype;
import org.darlinglab.metaphazer3C.model.GenotypeBuilder;
import org.darlinglab.metaphazer3C.model.GenotypeMixture;
import org.darlinglab.metaphazer3C.operator.GenotypeGibbsSweep;

public class MetaPhazer3C {

	public static void main(String[] args) {
		// parse the SNV linkage data
		File graphFile = new File(args[1]);
		SnvGraph snvs = null;
		try {
			snvs = SnvGraphBuilder.build(graphFile);
		} catch (IOException e) {
			System.err.println("Error parsing input file " + graphFile);
			e.printStackTrace();
			return;
		}
		
		// create an initial model 
		final int GENOTYPE_COUNT = 2;
		GenotypeMixture gm = new GenotypeMixture();
		gm.genotypes = new Genotype[GENOTYPE_COUNT];
		gm.abundances = new float[GENOTYPE_COUNT];
		for(int g=0; g < gm.genotypes.length; g++){
			// initialize the genotype with a random state
			gm.genotypes[g] = GenotypeBuilder.randomGenotype(snvs);
			// initialize to uniform abundances
			gm.abundances[g] = 1.0f/((float)GENOTYPE_COUNT);
		}
		GenotypeMixtureLikelihood gml = new GenotypeMixtureLikelihood(gm, snvs);
		GenotypeGibbsSweep ggs = new GenotypeGibbsSweep(gml, snvs);

		// loop over genotype updates until the end of time
		// TODO: add operators to change genotype abundances
		while(true){
			for(Genotype g : gm.genotypes){
				ggs.move(g);
			}
			System.out.println("Current LL: " + gml.getLikelihood());
		}
	}

}
