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
        String inputTumorDEG = "../DataSource/PanCancer.Atlas.tumorDEGinput.txt";
        String inputTumorSGA = "../DataSource/PanCancer.Atlas.tumorSGAinput.txt";
        String outputSGAState = "../DataSource/SGAStatesFile_PanCancerAtlas.csv";
//        String inputSgaDegTumor = "..\\DataSource\\SgaDegTumor_test2.csv";
//        String OutputGraph = "C:\\Users\\XIM33\\Documents\\JavaApplications\\TDIAppOnFGS\\data\\strParentsAndChildren_test2.txt";
//        String inputTumorDEG = "..\\DataSource\\DEGInput.txt";
//        String outputSGAState = "..\\DataSource\\SGAStatesFile.csv";

        TdiAppOnFGS app = new TdiAppOnFGS();
        Object[] tripletsInfo = app.readTriplets(inputSgaDegTumor);
        List listMapDEGs = app.readDegFile(inputTumorDEG);
        List listMapSGAs = app.readDegFile(inputTumorSGA);
        Object[] SGAState = app.inferSGAState(tripletsInfo, listMapSGAs, listMapDEGs);
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

        Map<String, Set<String>> mapSgaDegs = new HashMap<String, Set<String>>();
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

                    if (mapSgaDegs.containsKey(strSGA)) {
                        mapSgaDegs.get(strSGA).add(strDEG);
                    } else {
                        Set<String> setDEGs = new HashSet<String>();
                        setDEGs.add(strDEG);
                        mapSgaDegs.put(strSGA, setDEGs);
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

        //define return objects
        Object[] returnValues = new Object[2];
        returnValues[0] = mapEdgeScores;
        returnValues[1] = mapSgaDegs;

        return returnValues;
    }

    /**
     * Read in Tumor DEG list and return listOfMapDegs
     *
     * @param DegFile: file name
     * @return List[Map<String:int>] listOfMapDegs
     */
    public List<HashMap<String, Integer>> readDegFile(String DegFile) {

        List<HashMap<String, Integer>> listMapDEGs = new ArrayList<HashMap<String, Integer>>();
//        List<String[]> listDEGs= new ArrayList<String[]>();
        int lineCount = 0;
        String[] items;
        String sCurrentLine;
        try (BufferedReader br = new BufferedReader(new FileReader(DegFile))) {
            while ((sCurrentLine = br.readLine()) != null) {
                lineCount++;
                if (lineCount % 2 == 1) {
//                    listOfTumors.add(sCurrentLine);
                    continue;
                } else {
                    items = sCurrentLine.split(",");
                    HashMap<String, Integer> mapDEG = new HashMap<String, Integer>();
                    for (String item : items) {
                        mapDEG.put(item, 1);
                    }
                    listMapDEGs.add(mapDEG);
//                    listDEGs.add(items);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return listMapDEGs;
    }

    /**
     * Read in Tumor SGA list and return listOfMapDegs
     *
     * @param SgaFile: file name
     * @return List[Map<String:int>] listOfMapDegs
     */
    public List<HashMap<String, Integer>> readSgaFile(String SgaFile) {

        List<HashMap<String, Integer>> listMapSGAs = new ArrayList<HashMap<String, Integer>>();
//        List<String[]> listSGAs = new ArrayList<String[]>();

        int lineCount = 0;
        String[] items;
        String sCurrentLine;
        try (BufferedReader br = new BufferedReader(new FileReader(SgaFile))) {
            while ((sCurrentLine = br.readLine()) != null) {
                lineCount++;
                if (lineCount % 2 == 1) {
//                    listOfTumors.add(sCurrentLine);
                    continue;
                } else {
                    items = sCurrentLine.split(",");
                    HashMap<String, Integer> mapSGA = new HashMap<String, Integer>();
                    for (String item : items) {
                        mapSGA.put(item, 1);
                    }
                    listMapSGAs.add(mapSGA);
//                      listSGAs.add(items);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return listMapSGAs;
    }

    /**
     *
     * @param triplets
     * @param tumorDegs
     * @return
     */
    public Object[] inferSGAState(Object[] triplets, List<HashMap<String, Integer>> listTumorSGAs, List<HashMap<String, Integer>> listTumorDEGs) {
//    public Object[] inferSGAState(Object[] triplets, List<String[]> listTumorSGAs, List<String[]> listTumorDEGs) {
        if (listTumorSGAs.size() != listTumorDEGs.size()) {
            System.out.println("Tumor numbers from Driver SGA Table are not the same with the numbers from Taget DEG table");
            System.exit(1);
        }

        Map<String, double[]> mapEdgeScores = (Map<String, double[]>) triplets[0];
        Map<String, Set<String>> mapSGAtoDEGs = (Map<String, Set<String>>) triplets[1];

        //get listSGA to gurantee the order when do multithread
        ArrayList<String> listDriverSGAs = new ArrayList<String>();
        for (String nameSGA : mapSGAtoDEGs.keySet()) {
            listDriverSGAs.add(nameSGA);
        }

        int sizeDriver = listDriverSGAs.size();
        int sizeTumors = listTumorSGAs.size();
        System.out.println("sizeTumors = " + sizeTumors);
        
        List<List<Double>> tumorSGAStates = new ArrayList<List<Double>>(Collections.nCopies(sizeTumors, null));
        
        for (int i = 0; i < sizeTumors; i++) {
            System.out.println("Processing tumor " + i);
            List<Double> SGAStates = new ArrayList<Double>();

            for (int j = 0; j < sizeDriver; j++) {

                String driverSGA = listDriverSGAs.get(j);
                //if driver SGA is in the tumor SGA list, that means that dirver state is 1
                if (listTumorSGAs.get(i).containsKey(driverSGA)) 
                {
                    SGAStates.add(1.0);
                } 
                else//calculatte the SGA state
                {
                    //get the DEGs of drive SGA
                    Set<String> driverDEGs = mapSGAtoDEGs.get(driverSGA);

                    double logDEG_S0 = 0.0;
                    double logDEG_S1 = 0.0;
                    //check target DEG table of that tumor to see if each driver DEG is 0/1, then calculate the infer SGA state
                    for (String driverDEG : driverDEGs) 
                    {
                        //get the edge score
                        double[] scores = mapEdgeScores.get(driverSGA + "," + driverDEG);
                        if (scores == null){ 
                            System.out.println("Edge " + driverSGA + "-->" + driverDEG + " is not in the triplets\n");
                            System.exit(1);
                        }
                        //check DEGgState
                        int DEGState;
                        if (listTumorDEGs.get(i).containsKey(driverDEG)) {
                            DEGState = 1;
                        } else {
                            DEGState = 0;
                        }

                        //calculste infer SGA state                        
                        if (DEGState == 1) {
                            //scores[3] S1D1; scores[2] S1D0; scores[1] S0D1; scores[0] S0D0
                            logDEG_S0 += log((scores[1] + 0.0000001) / (scores[1] + scores[0]));
                            logDEG_S1 += log((scores[3] + 0.0000001) / (scores[3] + scores[2]));
                        } else {
                            logDEG_S0 += log((scores[0] + 0.0000001) / (scores[1] + scores[0]));
                            logDEG_S1 += log((scores[2] + 0.0000001) / (scores[3] + scores[2]));
                        }
                    }
                    double S1_DEGs = 1 / (1 + exp(logDEG_S0 - logDEG_S1));
                    SGAStates.add(S1_DEGs);
                }
            }
            tumorSGAStates.set(i, SGAStates);
        }
        return new Object[]{listDriverSGAs, tumorSGAStates};
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
//        List<String> listOfTumors = (List<String>) SGA_tumor_SGAStates[1];
        List<List<Double>> tumorSGAStates = (List<List<Double>>) SGA_tumor_SGAStates[1];
        
        try (BufferedWriter bf = new BufferedWriter(new FileWriter(SGAStateFile))) {
            //write in SGA names
            String strSGA = "";
            for (String SGAName : listSGA) {
//                bf.write(", " + SGAName);
                strSGA += SGAName + ",";
            }

            strSGA = strSGA.substring(0, strSGA.length() - 1);
            bf.write(strSGA + "\n");

            //Write in SGA probabily of each tumor
            for (int i = 0; i < tumorSGAStates.size(); i++) {
//                bf.write(listOfTumors.get(i));
                String strStates = "";
                for (int j = 0; j < listSGA.size(); j++) {
                    strStates += tumorSGAStates.get(i).get(j) + ",";
                }
                //remove the last ","
                strStates = strStates.substring(0, strStates.length() - 1);
                bf.write(strStates + "\n");
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
