package org.darlinglab.metaphazer3C.data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.Vector;
import java.lang.String;

/**
 * builds an SNV graph from a file in the MCL "abc" format
 * weights are expected to be link counts, with nodes as Site.Allele
 * @author koadman
 *
 */
public class SnvGraphBuilder {
	public static SnvGraph build(java.io.File f) throws IOException {
		SnvGraph snvs = new SnvGraph();
		BufferedReader reader = new BufferedReader(new FileReader(f));
		String line;
		HashMap<Integer,Integer> siteMap = new HashMap<>();
		int siteId = 0;
		HashMap<Integer,SnvGraph.Node> nodeMap = new HashMap<>();
		HashMap<SnvGraph.Node, HashMap<SnvGraph.Node,SnvGraph.Edge>> edgeMap = new HashMap<>();
		while(true){
			line = reader.readLine();
			if(line == null)
				break;
			StringTokenizer tok = new StringTokenizer(line);
			String a = tok.nextToken();
			String b = tok.nextToken();
			String c = tok.nextToken();
			int aSite = parseSite(a);
			int bSite = parseSite(b);
			int aAllele = parseAllele(a);
			int bAllele = parseAllele(b);
			int count = Integer.parseInt(c);
			if(!nodeMap.containsKey(aSite)){
				nodeMap.put(aSite, snvs.new Node());
				siteMap.put(aSite, siteId++);
			}
			if(!nodeMap.containsKey(bSite)){
				nodeMap.put(bSite, snvs.new Node());
				siteMap.put(bSite, siteId++);
			}
			// now find or create the edge and add the genotype link weight to it
			SnvGraph.Node aNode = nodeMap.get(aSite);
			SnvGraph.Node bNode = nodeMap.get(bSite);
			SnvGraph.Edge e;
			if(!edgeMap.containsKey(aNode)){
				edgeMap.put(aNode, new HashMap<SnvGraph.Node, SnvGraph.Edge>());
			}
			if(!edgeMap.containsKey(bNode)){
				edgeMap.put(bNode, new HashMap<SnvGraph.Node, SnvGraph.Edge>());
			}
			if(!edgeMap.get(aNode).containsKey(bNode)){
				e = snvs.new Edge();
				e.node_1 = siteMap.get(aSite);
				e.node_2 = siteMap.get(bSite);
				edgeMap.get(aNode).put(bNode,e);
				edgeMap.get(bNode).put(aNode,e);
			}
			e = edgeMap.get(aNode).get(bNode);
			if(e.node_1 == aSite && e.node_2 == bSite)
				e.counts[aAllele][bAllele] = count;
			else if(e.node_1 == bSite && e.node_2 == aSite)
				e.counts[bAllele][aAllele] = count;
			else
				System.err.println("Mismatched node ID and reference\n");
		}
		reader.close();
		return null;
	}
	
	private static int parseSite(String a){
		StringTokenizer aTok = new StringTokenizer(a, ":");
		String aSiteStr = aTok.nextToken();
		return Integer.parseInt(aSiteStr);
	}
	
	private static int parseAllele(String a){
		StringTokenizer aTok = new StringTokenizer(a, ":");
		String str = aTok.nextToken();
		str = aTok.nextToken();
		if(str.equalsIgnoreCase("A")){
			return 0;
		}else if(str.equalsIgnoreCase("C")){
			return 1;
		}else if(str.equalsIgnoreCase("G")){
			return 2;
		}else if(str.equalsIgnoreCase("T")){
			return 3;
		}
		throw new RuntimeException("Error parsing allele from MCL ABC format file");
	}
}
