package org.darlinglab.metaphazer3C;

public class Random {
	public static java.util.Random randy = new java.util.Random();

	public static int sampleCategorical(int[] counts){
		float[] freqs = new float[counts.length];
		for(int i=0; i<freqs.length; i++){
			freqs[i] = counts[i];
		}
		return sampleCategorical(freqs);
	}
	public static int sampleCategorical(float[] freqs){
		double[] d = new double[freqs.length];
		for(int i=0; i<freqs.length; i++){
			d[i] = freqs[i];
		}
		return sampleCategorical(d);
	}
	public static int sampleCategorical(double[] freqs){
		double r = randy.nextDouble();
		double sum = 0;
		for( double f : freqs){ sum += f; }
		r /= sum;
		sum = 0;
		for(int i=0; i<freqs.length; i++){
			sum += freqs[i];
			if(sum > r) return i;
		}
		return freqs.length;
	}
}
