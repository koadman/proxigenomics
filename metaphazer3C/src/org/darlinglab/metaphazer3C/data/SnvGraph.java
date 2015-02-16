package org.darlinglab.metaphazer3C.data;

/**
 * an undirected graph where nodes are variant sites and edges link nodes that were observed together
 * in a sequence read or read pair
 * @author koadman
 *
 */
public class SnvGraph {

		static public final int ALPHABET_SIZE = 4; // assume nucleotides for now (not mC, etc.)

		public class Node{
			public int site;
			public float[] freqs;
			Edge[] edges; // all edges (undirected)
			public Node(){
				freqs = new float[ALPHABET_SIZE];
			}
		}

		public class Edge{
			public int node_1;
			public int node_2;
			public int[][] counts;
			public Edge(){
				counts = new int[ALPHABET_SIZE][ALPHABET_SIZE];
			}
		}

		public Node[] nodes;
		public Edge[] edges;
}
