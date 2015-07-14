package org.darlinglab.metaphazer3C;

import java.io.File;

import org.darlinglab.metaphazer3C.data.SnvGraph;
import org.darlinglab.metaphazer3C.data.SnvGraphBuilder;
import org.darlinglab.metaphazer3C.likelihood.GenotypeMixtureLikelihood;
import org.darlinglab.metaphazer3C.model.Genotype;
import org.darlinglab.metaphazer3C.model.GenotypeMixture;
import org.darlinglab.metaphazer3C.operator.GenotypeGibbsSweep;

public class MetaPhazer3C {

	public static void main(String[] args) {
		File graph_file = new File(args[1]);
		SnvGraph snvs = SnvGraphBuilder.build(graph_file);
		// initialize the model with a random state
		
		// loop over updates until the end of time
		GenotypeMixture gm = new GenotypeMixture();
		GenotypeMixtureLikelihood gml = new GenotypeMixtureLikelihood(gm, snvs);
		GenotypeGibbsSweep ggs = new GenotypeGibbsSweep(gml, snvs);

		// TODO: add operators to change genotype abundances
		while(true){
			for(Genotype g : gm.genotypes){
				ggs.move(g);
			}
		}
	}

}
