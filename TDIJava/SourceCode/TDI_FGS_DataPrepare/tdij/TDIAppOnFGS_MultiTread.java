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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author XIM33
 */
public class TdiAppOnFGS {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.out.println("Working Directory = " + System.getProperty("user.dir"));

        String inputSgaDegTumor = "../DataSource/SgaDegTumor_PanCancerAtlas.csv";
//        String graphOutputFile = "..\\DataSource\\strParentsAndChildren_test2.txt";
        String inputTumorDEG = "../DataSource/PanCancer.Atlas.tumorDeginput.txt";
        String outputSGAState = "../DataSource/SGAStatesFile_PanCancerAtlas.csv";
//        String inputSgaDegTumor = "..\\DataSource\\SgaDegTumor_test2.csv";
//        String OutputGraph = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\strParentsAndChildren_test2.txt";
//        String inputTumorDEG = "..\\DataSource\\DEGInput.txt";
//        String outputSGAState = "..\\DataSource\\SGAStatesFile.csv";

        TdiAppOnFGS app = new TdiAppOnFGS();
        Object[] tripletsInfo = app.readTriplets(inputSgaDegTumor);
        Object[] tumorDegs = app.readDegFile(inputTumorDEG);
        Object[] SGAState = app.inferSGAState(tripletsInfo, tumorDegs);
        app.output(SGAState, outputSGAState);

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
//                    }u709
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

    /**
     * Read in triplet and return mapEdgeScores and setSGA;
     *
     * @param fileIn triple file name
     * @return object[0]: Map<String, double[]> mapEdgeScores object[1]: Set
     * setSGA
     */
    public Object[] readTriplets(String fileIn) {

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

        totalTumors = setTotalTumors.size();

        //build map mapEdgeScores
        Map<String, double[]> mapEdgeScores = new HashMap<String, double[]>();
        for (String itemEdge : mapSGAtoDEGTumors.keySet()) {
            String items[] = itemEdge.split(",");
            String SGAname = items[0];
            String DEGname = items[1];
            int c11 = mapSGAtoDEGTumors.get(itemEdge).size();
            int c10 = mapSGATumors.get(SGAname).size() - c11;
            int c01 = mapDEGTumors.get(DEGname).size() - c11;

            //find the common tumors of C10 and C01
            Set<String> intersectionSGAandDEG = new HashSet<String>(mapSGATumors.get(SGAname));
            intersectionSGAandDEG.retainAll(mapDEGTumors.get(DEGname));
            //This intersection contains C11 tumors, so need to remove
            Set<String> setSGADEG = mapSGAtoDEGTumors.get(itemEdge);
            intersectionSGAandDEG.removeAll(setSGADEG);

            int c00 = totalTumors - c11 - c10 - c01 + intersectionSGAandDEG.size();
            double[] scores = {c00, c01, c10, c11};

            mapEdgeScores.put(itemEdge, scores);
        }

        //convert set to listSGA to gurantee the order
        List listSGA = new ArrayList<String>();
        for(String nameSGA:mapSGATumors.keySet()){
            listSGA.add(nameSGA);
        }
        
        //define return objects
        Object[] returnValues = new Object[2];
        returnValues[0] = mapEdgeScores;
        returnValues[1] = listSGA;

        return returnValues;
    }

    /**
     * Read in Tumor DEG list and return listOfTumors and listOfDegs
     *
     * @param DegFile: file name
     * @return Object[0] List<String> ListOfTumors
     * @return Object[1] List<String[]> listOfDegs
     */
    public Object[] readDegFile(String DegFile) {

        List<String[]> listOfDegs = new ArrayList<String[]>();
        List<String> listOfTumors = new ArrayList<String>();
        int lineCount = 0;
        String[] items;
        String sCurrentLine;
        try (BufferedReader br = new BufferedReader(new FileReader(DegFile))) {
            while ((sCurrentLine = br.readLine()) != null) {
                lineCount++;
                if (lineCount % 2 == 1) {
                    listOfTumors.add(sCurrentLine);
                } else {
                    items = sCurrentLine.split(",");
                    listOfDegs.add(items);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new Object[]{listOfTumors, listOfDegs};
    }

    /**
     *
     * @param triplets
     * @param tumorDegs
     * @return
     */
    public Object[] inferSGAState(Object[] triplets, Object[] tumorDegs) {
        Map<String, double[]> mapEdgeScores = (Map<String, double[]>) triplets[0];
        List<String> listSGA = (List<String>) triplets[1];
        List<String> listOfTumors = (List<String>) tumorDegs[0];
        List<String[]> listOfDEGs = (List<String[]>) tumorDegs[1];
        
//        List<String> SGANames = new ArrayList<String>();
        int sizeSGA = listSGA.size();
        List<List<Double>>  tumorSGAStates = new ArrayList<List<Double>>(Collections.nCopies(sizeSGA,null));

        ExecutorService executor = Executors.newFixedThreadPool(10);
//        for (String SGA : setSGA) {
          for (int i=0; i<listSGA.size(); i++){

            final int idx = i;
              
            executor.submit(new Runnable() {
                @Override
                public void run() {
                    double pS1_DEGs;

                    double pDEG_S0;
                    double pDEG_S1;
                    System.out.println("Processing " + idx + "/" + sizeSGA + "th SGA\n");
                    String nameSGA = listSGA.get(idx);
                    List<Double> SGAStates = new ArrayList<Double>();

                    //Loop DEG list
                    for (String[] Degs : listOfDEGs) {

                        pDEG_S0 = 0.0;
                        pDEG_S1 = 0.0;

                        //for each DEG, combine with SGA, get edge scores, then calculate pDEG_S0,pDEG_S1
                        for (String nameDEG : Degs) {
                            double[] scores = mapEdgeScores.get(nameSGA + "," + nameDEG);

                            if (scores == null) {
//                                System.out.println("Edge " + nameSGA + "-->" + nameDEG + " is not in the triplets\n");
                                continue;
                            }

                            //scores[3] S1D1; scores[2] S1D0; scores[1] S0D1; scores[0] S0D0
                            pDEG_S0 += log((scores[2] + 0.0000001) / (scores[2] + scores[0]));
                            pDEG_S1 += log((scores[3] + 0.0000001) / (scores[3] + scores[1]));
                        }

                        pS1_DEGs = 1 / (1 + exp(pDEG_S0 - pDEG_S1));

                        SGAStates.add(pS1_DEGs);

                    }
                    tumorSGAStates.set(idx,SGAStates);

                }
            }
            
            );
        }
        executor.shutdown();

        System.out.println(
                "All tasks submitted.");
        try {
            executor.awaitTermination(5, TimeUnit.DAYS);
        } catch (InterruptedException ignored) {
        }
        System.out.println("All tasks completed.");

        return new Object[]{listSGA, listOfTumors, tumorSGAStates};
    }

    /**
     * Output results to a file
     *
     * @param Object SGA_tumor_SGAStates[0]: Set<String> setSGA Object
     * SGA_tumor_SGAStates[1]: List<String> listOfTumors Object
     * SGA_tumor_SGAStates[1]: List<List<Double>> tumorSGAStates
     * @param String SGAStateFile :output file name
     */
    public void output(Object[] SGA_tumor_SGAStates, String SGAStateFile) {
        List<String> listSGA = (List<String>) SGA_tumor_SGAStates[0];
        List<String> listOfTumors = (List<String>) SGA_tumor_SGAStates[1];
        List<List<Double>> tumorSGAStates = (List<List<Double>>) SGA_tumor_SGAStates[2];
        try (BufferedWriter bf = new BufferedWriter(new FileWriter(SGAStateFile))) {
            //write in SGA names
            for (String SGAName : listSGA) {
                bf.write(", " + SGAName);
            }
            bf.write("\n");

            //Write in SGA probabily of each tumor
            for (int i = 0; i < listOfTumors.size(); i++) {
                bf.write(listOfTumors.get(i));
                for (int j = 0; j < listSGA.size(); j++) {
                    bf.write("," + tumorSGAStates.get(j).get(i));
                }
                bf.write("\n");
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
