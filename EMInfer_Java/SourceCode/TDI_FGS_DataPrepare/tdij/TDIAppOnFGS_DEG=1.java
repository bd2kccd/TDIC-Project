/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tdij;

import tetrad.data.*;
import tetrad.graph.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.exp;
import static java.lang.Math.log;

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

        String SgaDegTumorFile = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\SgaDegTumor_PanCancerAtlas.csv";
//        String graphOutputFile = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\strParentsAndChildren_test2.txt";
        String DEGInputFile = "C:\\Users\\XIM33\\Documents\\NetBeansProjects\\DataSource\\PanCancer.Atlas.tumorDeginput.txt";
        String SGAStateFile = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\SGAStatesFile_PanCancerAtlas.csv";
//        String SgaDegTumorFile = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\SgaDegTumor_test2.csv";
//        String graphOutputFile = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\strParentsAndChildren_test2.txt";
//        String DEGInputFile = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\DEGInput.txt";
//        String SGAStateFile = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\SGAStatesFile.csv";

        TdiAppOnFGS app = new TdiAppOnFGS();
        EdgeListGraphSingleConnections TDIGraph = app.readData(SgaDegTumorFile);
        app.inferSGAState(DEGInputFile, TDIGraph, SGAStateFile);

//        //test node parents and children and count
//        try (BufferedWriter bf = new BufferedWriter(new FileWriter(graphOutputFile))) {
//            //write out parents and childrens to a file
//            for (Node eachNode : TDIGraph.getNodes()) {
//                if (eachNode.getGeneType() == "SGA") {
//                    List<Node> nodeChildren = TDIGraph.getChildren(eachNode);
//                    for (Node child : nodeChildren) {
//                        bf.write("Node " + eachNode.getName() + "'child is " + child + "\r\n");
//                    }
//                } else {
//                    List<Node> nodeParents = TDIGraph.getParents(eachNode);
//                    for (Node parent : nodeParents) {
//                        bf.write("Node " + eachNode.getName() + "'parent is " + parent + "\r\n");
//                    }
//                }
//            }
//            //test count
//            for (Edge edge : TDIGraph.getEdges()) {
//                double scores[] = edge.getEdgeScores();
//                bf.write("Edge " + edge.toString() + " score are " + scores[3] + "," + scores[2] + "," + scores[1] + "," + scores[0] + "\r\n");
//            }
//
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//
//        Node node = TDIGraph.getNode("D2");
//        if (node.getGeneType() == "SGA") {
//            List<Node> children = TDIGraph.getChildren(node);
//            for (Node child : children) {
//                System.out.println("Node " + node.getName() + "'child is " + child + "\n");
//            }
//        } else {
//            List<Node> parents = TDIGraph.getParents(node);
//            for (Node parent : parents) {
//                System.out.println("Node " + node.getName() + "'parent is " + parent + "\n");
//            }
//        }

    }

    public List<List<String>> readDegFile(String DegFile) {
        List<String> items = new ArrayList<String>();
        List<List<String>> listOfDegs = new ArrayList<List<String>>();
        List<String> tumorNames = new ArrayList<String>();
        int lineCount = 0;
        String sCurrentLine;
        try (BufferedReader br = new BufferedReader(new FileReader(DegFile))) {
            while ((sCurrentLine = br.readLine()) != null) {
                lineCount++;
                if (lineCount % 2 == 1) {
                    tumorNames.add(sCurrentLine);
                } else {
                    listOfDegs.add(items);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return listOfDegs;
    }

    public EdgeListGraphSingleConnections readData(String fileIn) {

        int totalTumors;
        Set<String> setTotalTumors = new HashSet<String>();
        Map<String, Set<String>> mapSGATumors = new HashMap<String, Set<String>>();
        Map<String, Set<String>> mapDEGTumors = new HashMap<String, Set<String>>();
        Map<String, Set<String>> mapSGAtoDEGTumors = new HashMap<String, Set<String>>();
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

                    setTotalTumors.add(strTumor);
                    if (mapSGAtoDEGTumors.containsKey(strSGADEG)) {
                        mapSGAtoDEGTumors.get(strSGADEG).add(strTumor);
                    } else {
                        Set<String> setTumors = new HashSet<String>();
                        setTumors.add(strTumor);
                        mapSGAtoDEGTumors.put(strSGADEG, setTumors);
                    }

                    if (mapSGATumors.containsKey(strSGA)) {
                        mapSGATumors.get(strSGA).add(strTumor);
                    } else {
                        Set<String> setTumors = new HashSet<String>();
                        setTumors.add(strTumor);
                        mapSGATumors.put(strSGA, setTumors);
                    }

                    if (mapDEGTumors.containsKey(strDEG)) {
                        mapDEGTumors.get(strDEG).add(strTumor);
                    } else {
                        Set<String> setTumors = new HashSet<String>();
                        setTumors.add(strTumor);
                        mapDEGTumors.put(strDEG, setTumors);
                    }
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        
       
        //build graph
        totalTumors = setTotalTumors.size();
        
        EdgeListGraphSingleConnections TDIGraph = new EdgeListGraphSingleConnections();
        
        //add nodes
        for (String itemSGA : mapSGATumors.keySet()) {
            Node nodeSGA = new DiscreteVariable(itemSGA);
            nodeSGA.setGeneType("SGA");
            TDIGraph.addNode(nodeSGA);

        }
        for (String itemDEG : mapDEGTumors.keySet()) {
            Node nodeDEG = new DiscreteVariable(itemDEG);
            nodeDEG.setGeneType("DEG");
            TDIGraph.addNode(nodeDEG);
        }
        //add edges
        for (String itemEdge : mapSGAtoDEGTumors.keySet()) {
            String items[] = itemEdge.split(",");
//            System.out.println(itemEdge);
            String SGAname = items[0];
            String DEGname = items[1];
            Node nodeSGA = TDIGraph.getNode(SGAname, "SGA");
            Node nodeDEG = TDIGraph.getNode(DEGname, "DEG");
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
//            System.out.println(edgeSGAtoDEG.toString() + "\n");
            TDIGraph.addEdge(edgeSGAtoDEG);

        }

        return TDIGraph;
    }

    public void inferSGAState(String DegFile, EdgeListGraphSingleConnections TDIGraph, String SGAStateFile) {

        List<String[]> listOfDEGs = new ArrayList<String[]>();
        List<String> tumorNames = new ArrayList<String>();
        List<String> SGANames = new ArrayList<String>();
        List<List<Double>> tumorSGAStates = new ArrayList<List<Double>>();

        int lineCount = 0;
        String sCurrentLine;
        String[] items;
        try (BufferedReader br = new BufferedReader(new FileReader(DegFile))) {
            while ((sCurrentLine = br.readLine()) != null) {
                lineCount++;
                items = sCurrentLine.split(",");
                if (lineCount % 2 == 1) {
                    tumorNames.add(items[0]);
                } else {
                    listOfDEGs.add(items);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        List<Node> nodes = TDIGraph.getNodes();

        for (Node node : nodes) {

            if (node.getGeneType() == "SGA") {
                //get SGA ame
                String nameSGA = node.getName();
                SGANames.add(nameSGA);
                
                double pS1_DEGs;
                double pDEG_S0;
                double pDEG_S1;

                List<Double> DEGStates = new ArrayList<Double>();

                //Loop DEG list
                for (String[] Degs : listOfDEGs) {

                    pDEG_S0 = 0.0;
                    pDEG_S1 = 0.0;
 
                    //for each DEG, combine with SGA, get edge scores, then calculate pDEG_S0,pDEG_S1
                    for (String nameDEG : Degs) {
                        Node node1 = TDIGraph.getNode(nameSGA, "SGA");
                        Node node2 = TDIGraph.getNode(nameDEG, "DEG");
                        if (node2 == null) {
                            System.out.println("DEG " + nameDEG + " is not in our dataset\n");
                            continue;
                        }
                        Edge edge = TDIGraph.getEdge(node1, node2);
                        if (edge == null) {
                            System.out.println("Edge " + nameSGA + "-->" + nameDEG + " is not in our dataset\n");
                            continue;
                        }
                        double[] scores = edge.getEdgeScores();
                        //scores[3] S1D1; scores[2] S1D0; scores[1] S0D1; scores[0] S0D0
                        pDEG_S0 += log((scores[2]) / (scores[2] + scores[0]));
                        pDEG_S1 += log((scores[3]) / (scores[3] + scores[1]));
                    }
                    
                    pS1_DEGs = 1 / (1 + exp(pDEG_S0 - pDEG_S1));
                    
                    DEGStates.add(pS1_DEGs);
                }

                tumorSGAStates.add(DEGStates);

            } else {
                break;//SGA nodes was added in sequence. Fist SGA, then DEG
            }
        }

        //Output results to a file
        try (BufferedWriter bf = new BufferedWriter(new FileWriter(SGAStateFile))) {
            //write in SGA names
            for (String SGAName : SGANames) {
                bf.write(", " + SGAName);
            }
            bf.write("\n");

            //Write in SGA probabily of each tumor
            for (int i = 0; i < tumorNames.size(); i++) {
                bf.write(tumorNames.get(i));
                for (int j = 0; j < SGANames.size(); j++) {
                    bf.write("," + tumorSGAStates.get(j).get(i));
                }
                bf.write("\n");
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
