/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tdij;

import tetrad.data.*;
import tetrad.graph.*;
import tetrad.util.*;
import tetrad.stat.*;
import tetrad.search.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

//import java.io.ObjectInputStream;
//import java.text.NumberFormat;
import java.util.*;

/**
 *
 * @author XIM33
 */
public class TdiAppOnFGS {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String readInFile = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\SgaDegTumor_test2.csv";
        String writeOutFile = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\strParentsAndChildren_test2.txt";

        Graph TDIgraph = readData(readInFile);

        //test node parents and children and count
        try (BufferedWriter bf = new BufferedWriter(new FileWriter(writeOutFile))) {
            //write out parents and childrens to a file
            for (Node eachNode : TDIgraph.getNodes()) {
                if (eachNode.getGeneType() == "SGA") {
                    List<Node> nodeChildren = TDIgraph.getChildren(eachNode);
                    for (Node child : nodeChildren) {
                        bf.write("Node " + eachNode.getName() + "'child is " + child + "\r\n");
                    }
                } else {
                    List<Node> nodeParents = TDIgraph.getParents(eachNode);
                    for (Node parent : nodeParents) {
                        bf.write("Node " + eachNode.getName() + "'parent is " + parent + "\r\n");
                    }
                }
            }
            //test count
            for (Edge edge : TDIgraph.getEdges()) {
                int scores[] = edge.getEdgeScores();
                bf.write("Edge " + edge.toString() + " score are " + scores[3] + "," + scores[2] + "," + scores[1] + "," + scores[0] + "\r\n");
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        Node node = TDIgraph.getNode("D2");
        if (node.getGeneType()== "SGA"){
            List<Node> children = TDIgraph.getChildren(node);
            for (Node child : children) {
                System.out.println("Node " + node.getName() + "'child is " + child + "\n");
            }
        }else{
            List<Node> parents = TDIgraph.getParents(node);
            for (Node parent : parents) {
                System.out.println("Node " + node.getName() + "'parent is " + parent + "\n");
            }
        }

    }

    public static Graph readData(String fileIn) {

        int totalTumors = 8;
        Map<String, Integer> mapSGAcount = new HashMap<String, Integer>();
        Map<String, Integer> mapDEGcount = new HashMap<String, Integer>();
        Map<String, Integer> mapSGAtoDEGcount = new HashMap<String, Integer>();
        Set<String> pairSGATumor = new HashSet<String>();
        Set<String> pairDEGTumor = new HashSet<String>();

        //read in file and save strs, edges and counts
        try (BufferedReader br = new BufferedReader(new FileReader(fileIn))) {
            String sCurrentLine;
            int lineCounter = 0;
            String[] items;

            while ((sCurrentLine = br.readLine()) != null) {
                lineCounter++;
                if (lineCounter == 1) {
                    continue;
                } else {
                    items = sCurrentLine.split(",");
                    String strSGA = items[0];
                    String strDEG = items[1];
                    String strSGADEG = items[0] + "," + items[1];

                    //SGA and DEG both not new -> tumor must new, so edge+1, new SGA pair, new DEG pair, SGA count+1, DEG count+1
                    if (mapSGAcount.containsKey(strSGA) && mapDEGcount.containsKey(strDEG)) {
                        mapSGAtoDEGcount.put(strSGADEG, mapSGAtoDEGcount.get(strSGADEG) + 1);
                        pairSGATumor.add(items[0] + items[2]);
                        pairDEGTumor.add(items[1] + items[2]);
                        mapSGAcount.put(strSGA, mapSGAcount.get(strSGA) + 1);
                        mapDEGcount.put(strDEG, mapDEGcount.get(strDEG) + 1);
                    } //SGA and DEG both new -> add edge, add SGApair, add DEGpair, add SGA str, add DEG str
                    else if (!mapSGAcount.containsKey(strSGA) && !mapDEGcount.containsKey(strDEG)) {
                        mapSGAtoDEGcount.put(strSGADEG, 1);
                        pairSGATumor.add(items[0] + items[2]);
                        pairDEGTumor.add(items[1] + items[2]);
                        mapSGAcount.put(strSGA, 1);
                        mapDEGcount.put(strDEG, 1);
                    } //new SGA str and old DEG str  -> new SGAstrtumor pair, new edge
                    // if DEGpair not exit -> add pair, DEG str count + 1
                    else if (!mapSGAcount.containsKey(strSGA)) {
                        //SGA related
                        mapSGAcount.put(strSGA, 1);
                        pairSGATumor.add(items[0] + items[2]);
                        mapSGAtoDEGcount.put(strSGADEG, 1);
                        //DEG related
                        if (!pairDEGTumor.contains(items[1] + items[2])) {
                            pairDEGTumor.add(items[1] + items[2]);
                            mapDEGcount.put(strDEG, mapDEGcount.get(strDEG) + 1);
                        }

                    } //new DEG str and old SGA str -> new DEGstrtumor pair, new edge
                    // if SGApair not exit -> add pair, SGA str count + 1
                    else if (!mapDEGcount.containsKey(strDEG)) {
                        //DEG related
                        mapDEGcount.put(strDEG, 1);
                        pairDEGTumor.add(items[1] + items[2]);
                        mapSGAtoDEGcount.put(strSGADEG, 1);
                        //SGA related
                        if (!pairSGATumor.contains(items[0] + items[2])) {
                            pairSGATumor.add(items[0] + items[2]);
                            mapSGAcount.put(strSGA, mapSGAcount.get(strSGA) + 1);
                        }
                    }

                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        //build graph
        EdgeListGraphSingleConnections TDIgraph = new EdgeListGraphSingleConnections();

        //add node
        for (String itemSGA : mapSGAcount.keySet()) {
            Node nodeSGA = new DiscreteVariable(itemSGA);
            nodeSGA.setGeneType("SGA");
            TDIgraph.addNode(nodeSGA);

        }
        for (String itemDEG : mapDEGcount.keySet()) {
            Node nodeDEG = new DiscreteVariable(itemDEG);
            nodeDEG.setGeneType("DEG");
            TDIgraph.addNode(nodeDEG);
        }
        //add edge
        for (String itemEdge : mapSGAtoDEGcount.keySet()) {
            String items[] = itemEdge.split(",");
            Node nodeSGA = TDIgraph.getNode(items[0]);
            Node nodeDEG = TDIgraph.getNode(items[1]);
            Edge edgeSGAtoDEG = new Edge(nodeSGA, nodeDEG, Endpoint.TAIL, Endpoint.ARROW);
            int c11 = mapSGAtoDEGcount.get(itemEdge);
            int c10 = mapSGAcount.get(items[0]) - c11;
            int c01 = mapDEGcount.get(items[1]) - c11;
            int c00 = totalTumors - c11 - c10 - c01;
            edgeSGAtoDEG.setEdgeScores(c00, c01, c10, c11);
            TDIgraph.addEdge(edgeSGAtoDEG);
        }

        return TDIgraph;
    }
}
