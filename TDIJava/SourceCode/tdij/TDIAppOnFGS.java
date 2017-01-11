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
        if (node.getGeneType() == "SGA") {
            List<Node> children = TDIgraph.getChildren(node);
            for (Node child : children) {
                System.out.println("Node " + node.getName() + "'child is " + child + "\n");
            }
        } else {
            List<Node> parents = TDIgraph.getParents(node);
            for (Node parent : parents) {
                System.out.println("Node " + node.getName() + "'parent is " + parent + "\n");
            }
        }

    }

    public static Graph readData(String fileIn) {

        int totalTumors = 4;
        Map<String, Set<String> > mapSGATumors = new HashMap<String, Set<String> >();
        Map<String, Set<String> > mapDEGTumors = new HashMap<String, Set<String> >();
        Map<String, Set<String> > mapSGAtoDEGTumors = new HashMap<String, Set<String> >();
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
                    String strTumor = items[2];
                    String strSGADEG = strSGA + "," + strDEG;
                    

                    if(mapSGAtoDEGTumors.containsKey(strSGADEG)){
                        mapSGAtoDEGTumors.get(strSGADEG).add(strTumor);
                    }else{
                        Set<String> setTumors = new HashSet<String>();
                        setTumors.add(strTumor);
                        mapSGAtoDEGTumors.put(strSGADEG, setTumors);
                    }

                    if(mapSGATumors.containsKey(strSGA)){
                        mapSGATumors.get(strSGA).add(strTumor);
                    }else{
                        Set<String> setTumors = new HashSet<String>();
                        setTumors.add(strTumor);
                        mapSGATumors.put(strSGA, setTumors);
                    }
                    
                    if(mapDEGTumors.containsKey(strDEG)){
                        mapDEGTumors.get(strDEG).add(strTumor);
                    }else{
                        Set<String> setTumors = new HashSet<String>();
                        setTumors.add(strTumor);
                        mapDEGTumors.put(strDEG, setTumors);
                    }
                }
            }
                    //SGA and DEG both not new 
        
        } catch (IOException e) {
            e.printStackTrace();
        }

        //build graph
        EdgeListGraphSingleConnections TDIgraph = new EdgeListGraphSingleConnections();

        //add nodes
        for (String itemSGA : mapSGATumors.keySet()) {
            Node nodeSGA = new DiscreteVariable(itemSGA);
            nodeSGA.setGeneType("SGA");
            TDIgraph.addNode(nodeSGA);

        }
        for (String itemDEG : mapDEGTumors.keySet()) {
            Node nodeDEG = new DiscreteVariable(itemDEG);
            nodeDEG.setGeneType("DEG");
            TDIgraph.addNode(nodeDEG);
        }
        //add edges
        for (String itemEdge : mapSGAtoDEGTumors.keySet()) {
            String items[] = itemEdge.split(",");
            String SGAname = items[0];
            String DEGname = items[1];
            Node nodeSGA = TDIgraph.getNode(SGAname);
            Node nodeDEG = TDIgraph.getNode(DEGname);
            Edge edgeSGAtoDEG = new Edge(nodeSGA, nodeDEG, Endpoint.TAIL, Endpoint.ARROW);
            int c11 = mapSGAtoDEGTumors.get(itemEdge).size();
            int c10 = mapSGATumors.get(SGAname).size() - c11;
            int c01 = mapDEGTumors.get(DEGname).size() - c11;
            //Need to search if SGA DEG 10 and 01 have the same tumor
            
            Set<String> intersectionSGAandDEG = new HashSet<String>(mapSGATumors.get(SGAname));
            intersectionSGAandDEG.retainAll(mapDEGTumors.get(DEGname));
            
            Set<String> setSGADEG = mapSGAtoDEGTumors.get(itemEdge);
            intersectionSGAandDEG.removeAll(setSGADEG);
            
            int c00 = totalTumors - c11 - c10 - c01 + intersectionSGAandDEG.size();
            edgeSGAtoDEG.setEdgeScores(c00, c01, c10, c11);
            TDIgraph.addEdge(edgeSGAtoDEG);
        }

        return TDIgraph;
    }
}
