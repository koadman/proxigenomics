package org.darlinglab.metaphazer3C.model;

/**
 * A class to represent a mixture of genotypes, each with their own relative abundance
 * @author koadman
 */
public class GenotypeMixture {
	public Genotype[] genotypes;
	public float[] abundances;
	public double snv_error_rate = 0.0001; /**< the rate at which sites are misread by the sequencer */
	public double ligation_error_rate = 0.001; /**< the rate at which chimeric ligation products form */
}
