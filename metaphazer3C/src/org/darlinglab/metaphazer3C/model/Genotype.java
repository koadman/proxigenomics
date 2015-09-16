package org.darlinglab.metaphazer3C.model;

import org.darlinglab.metaphazer3C.data.SnvGraph;

/**
 *  a genotype specifies one of the possible states (e.g. nucleotides) at every site
 * @author koadman
 */
public class Genotype {
	public int getState(int i){
		return states[i];
	}
	public void setState(int i, int state){
		states[i] = state;
		for(ModelListener m : ml){
			m.modelChanged(snvs.nodes[i]);
		}
	}
	int[] states;
	ModelListener[] ml;
	SnvGraph snvs;
}
