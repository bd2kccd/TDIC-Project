/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eminfer;

/**
 *
 * @author XIM33
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
//import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
//import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.exp;
import static java.lang.Math.log;
//import static java.lang.Math.exp;
//import static java.lang.Math.log;
import java.util.*;

public class InferData {
    private Map<String, Set<String>> mapSgaDegs = new HashMap<String, Set<String>>(); ;
    private ArrayList<String> solidEdgeList = new ArrayList<String>(); ;
    private ArrayList<String> driverSGAs = new ArrayList<String>(); ;
    private ArrayList<String> targetDEGs = new ArrayList<String>(); ;

    private    Map<String, Integer> mapSGA = new HashMap<String, Integer>();
    private    Map<String, Integer> mapDEG = new HashMap<String, Integer>();

    private Map<String, Double[]> mapEdgeParam;
    private Map<String, Double[]> mapSGAParam;
    
    
    
    
    private ArrayList<ArrayList<Integer>> targetDEGTable ;
    private ArrayList<String> tumorNames ;//tumor names was read in when get targetDEGTable
    private ArrayList<ArrayList<Integer>> driverSGATable ;
    private ArrayList<ArrayList<Integer>> driverSGATable_ori;

    private ArrayList<ArrayList<Double>> driverActivationTable;
    private ArrayList<ArrayList<Integer>> inferDriverTable;
    

    private ArrayList<String> edgeSGAs ;
    private ArrayList<String> edgeDEGs ;
    private ArrayList<String> edgeList ;
    private ArrayList<String> GtMatrixSGAs ;
    private ArrayList<String> GeMatrixDEGs ;

    
   public InferData(String fileEdgeList, String fileGtMatrix, String fileGeMatrix){
           ReadData(fileEdgeList,fileGtMatrix,fileGeMatrix);
    }


    public void ReadData(String fileEdgeList, String fileGtMatrix, String fileGeMatrix) {
        GtMatrixSGAs = new ArrayList<String>();
        GeMatrixDEGs = new ArrayList<String>();
        System.out.println("Reading GtMatrix to get SGA names...");
        GtMatrixSGAs = getMatrixColumns(fileGtMatrix);
        System.out.println("Reading GeMatrix to get DEG names...");
        GeMatrixDEGs = getMatrixColumns(fileGeMatrix);
        readInEdgeList(fileEdgeList); //=>save to tripletSGAs, tripletDEGs, tripletEdgeList
        updateSGAsDEGsEdges();
        readInGtMatrix(fileGtMatrix);
        readInGeMatrix(fileGeMatrix);
        

// for test purpose        
//        for (String item : solidEdgeList){
//            System.out.println(item + "\n");
//        }
//        System.out.println("\n");
//        for (String item : driverSGAs){
//            System.out.println(item + "\n");
//        }
//        System.out.println("\n");
//        for (String item : targetDEGs){
//            System.out.println(item + "\n");
//        }
//        System.out.println("\n");
//        for (String item:mapSgaDegs.keySet()){
//            System.out.println(item + ":");
//            for (String value: mapSgaDegs.get(item)){
//                System.out.println(value+",");
//            }
//            System.out.println("\n");
//        }
//        
//        for (int i = 0; i < driverSGATable.size(); i++) {
//            for (int j = 0; j < driverSGATable.get(i).size(); j++) {
//                System.out.print(driverSGATable.get(i).get(j));
//                System.out.print(',');
//            }
//            System.out.println('\n');
//
//        }
//        
//        for(int i=0; i<targetDEGTable.size(); i++){
//            String line="";
//            for(int j =0; j<targetDEGTable.get(i).size(); j++)
//                line += targetDEGTable.get(i).get(j);
//            System.out.println(line);
//            
//        }

    }

    
    public ArrayList<String> getMatrixColumns(String fileGtMatrix) {
        
        ArrayList<String> columnNames = new ArrayList<String>();
        try (BufferedReader br = new BufferedReader(new FileReader(fileGtMatrix))) {
            String firstLine = br.readLine();
            String[] items = firstLine.split(",");

            for (int i = 1; i < items.length; i++) {
                columnNames.add(items[i]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return columnNames;
    }

    public void readInEdgeList(String fileEdgeList) {
        System.out.println("Reading EdgeList data...");
        edgeSGAs = new ArrayList<String>();
        edgeDEGs = new ArrayList<String>();
        edgeList = new ArrayList<String>();
        try (BufferedReader br = new BufferedReader(new FileReader(fileEdgeList))) {
            String sCurrentLine;
            int lineCounter = 0;
            String[] items;

            while ((sCurrentLine = br.readLine()) != null) {
                lineCounter++;
//                System.out.println("linecounter = " + lineCounter);
                if (lineCounter == 1) {
                    continue;//This is the title
                } else {
                    items = sCurrentLine.split(",");
//                    String strSGA = items[1].trim();
//                    String strDEG = items[2].trim();
                    String strSGA = items[0].trim();
                    String strDEG = items[1].trim();

                    if (!edgeList.contains(strSGA + "," + strDEG)) {
                        edgeList.add(strSGA + "," + strDEG);
                    }

                    if (!edgeSGAs.contains(strSGA)) {
                        edgeSGAs.add(strSGA);
                    }
                    if (!edgeDEGs.contains(strDEG)) {
                        edgeDEGs.add(strDEG);
                    }
                }
            }
            System.out.println("Finish reading edgelist.");
        } catch (IOException e) {
            e.printStackTrace();
        }
        
    }

    public void updateSGAsDEGsEdges(){
    //update solidEdgeList
    //get intersection of SGAs DEGs of Matrix and triplets and save to driverSGAs and targetDEGs 
    //build mapSgaDegs
        System.out.println("Combining edgelist and matrix info...");
    
        //get new edgelist
        String[] items;
        for (String strEdge : edgeList) {
            items = strEdge.split(",");
            String SGA = items[0];
            String DEG = items[1];
            if (GtMatrixSGAs.contains(SGA) && GeMatrixDEGs.contains(DEG)){
                solidEdgeList.add(strEdge);
            }
        }
        
        //get edgeSGAs, edgeDEGs and build mapSGaDegs
        ArrayList<String> edgeSGAs = new ArrayList<String>();
        ArrayList<String> edgeDEGs = new ArrayList<String>();
        for (String strEdge : solidEdgeList) {
            items = strEdge.split(",");
            String strSGA = items[0];
            String strDEG = items[1];
            
            if (!edgeSGAs.contains(strSGA)) {
                edgeSGAs.add(strSGA);
            }
            if (!edgeDEGs.contains(strDEG)) {
                edgeDEGs.add(strDEG);
            }
            
            if (mapSgaDegs.containsKey(strSGA)) {
                mapSgaDegs.get(strSGA).add(strDEG);
            } else {
                Set<String> setDEGs = new HashSet<String>();
                setDEGs.add(strDEG);
                mapSgaDegs.put(strSGA, setDEGs);
            }
        }
        //based on edgeSGAs, edgeDEGs create driverSGAs, targetDEGs
        //driverSGAs, targetDEGs have the same item with edgeSGAs, edgeDEGs, but have different order 
        for (String SGA : GtMatrixSGAs){
            if (edgeSGAs.contains(SGA)){
                driverSGAs.add(SGA);
            }
        }
        for (String DEG : GeMatrixDEGs){
            if (edgeDEGs.contains(DEG)){
                targetDEGs.add(DEG);
            }
        }
        //create driverSGAs and targetDEGs name to index maps
        for (int i = 0; i < driverSGAs.size(); i++) {
            mapSGA.put(driverSGAs.get(i), i);
        }
        for (int i = 0; i < targetDEGs.size(); i++) {
            mapDEG.put(targetDEGs.get(i), i);
        }


    }

    public void readInGtMatrix(String fileGtMatrix) {
        System.out.println("Reading GtMatrix...");
  
        driverSGATable = new ArrayList<ArrayList<Integer>>();
        driverSGATable_ori = new ArrayList<ArrayList<Integer>>();
        try (BufferedReader br = new BufferedReader(new FileReader(fileGtMatrix))) {
            String sCurrentLine;
            int lineCounter = 0;
            String[] items;

            ArrayList<Integer> SGAIndx = new ArrayList<Integer>();
            while ((sCurrentLine = br.readLine()) != null) {
                lineCounter++;
                items = sCurrentLine.split(",");
                if (lineCounter == 1) {//this is the first line
                    for (int i = 0; i < items.length; i++) {
                        if (driverSGAs.contains(items[i])) {
                            SGAIndx.add(i);
                        }
                    }
                } else {
                    //tumorNames.add(items[0]); tumor name should be obtained in reading GeMatrix, since crossinfer read testDEGMatrix and infer SGA state of those tumors
                    ArrayList<Integer> SGAValueRow = new ArrayList<Integer>();
                    for (int j=0; j<SGAIndx.size(); j++){
                        SGAValueRow.add(Integer.parseInt(items[SGAIndx.get(j)]));
                    }

                    driverSGATable_ori.add(SGAValueRow);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        //copy driverSGATable_ori to driverSGATable
        int tumorSize = driverSGATable_ori.size();
        int SGASize = driverSGATable_ori.get(0).size();
        for (int i = 0; i < tumorSize; i++) {
            ArrayList<Integer> SGAValueRow = new ArrayList<Integer>();
            for (int j = 0; j < SGASize; j++) {
                SGAValueRow.add(driverSGATable_ori.get(i).get(j));
            }  
            driverSGATable.add(SGAValueRow);
        }
//        //for test
//        for (int i = 0; i < tumorSize; i++) {
//            for (int j = 0; j < SGASize; j++) {
//                System.out.print(driverSGATable.get(i).get(j));
//                System.out.print(",");
//            }
//            System.out.println("\n");
//        }
    }

    public void readInGeMatrix(String fileGeMatrix) {
        System.out.println("Reading GeMatrix...");

        tumorNames = new ArrayList<String>();
        targetDEGTable = new ArrayList<ArrayList<Integer>>();
        try (BufferedReader br = new BufferedReader(new FileReader(fileGeMatrix))) {
            String sCurrentLine;
            int lineCounter = 0;
            String[] items;
            ArrayList<Integer> DEGIndx = new ArrayList<Integer>();
            while ((sCurrentLine = br.readLine()) != null) {
                lineCounter++;
                items = sCurrentLine.split(",");

                if (lineCounter == 1) {//this is the first line
                    for (int i = 0; i < items.length; i++) {
                        if (targetDEGs.contains(items[i])) {
                            DEGIndx.add(i);
                        }
                    }
                } else {
                    tumorNames.add(items[0]);
                    ArrayList<Integer> DEGValueRow = new ArrayList<Integer>();
                    for (int j=0; j<DEGIndx.size(); j++){
                        DEGValueRow.add(Integer.parseInt(items[DEGIndx.get(j)]));
                    }

                    targetDEGTable.add(DEGValueRow);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    public void updateDriverSGATable() {
        System.out.println("Update driverSGATable...");
        int tumorSize = driverSGATable.size();
        int SGASize = driverSGATable.get(0).size();
        for (int i = 0; i < tumorSize; i++) {
            for (int j = 0; j < SGASize; j++) {
                Integer v = inferDriverTable.get(i).get(j);
                driverSGATable.get(i).set(j,v);
            }
            
        }
//        for test purpose
//        for (int i = 0; i < tumorSize; i++) {
//            for (int j = 0; j < SGASize; j++) {
//                System.out.print(driverSGATable.get(i).get(j));
//                System.out.print(",");
//            }
//            System.out.println("\n");
//        }
        
    }
    
    public void estimateParams() {
        
        mapEdgeParam = new HashMap<String, Double[]>();
        mapSGAParam = new HashMap<String, Double[]>();
        mapSGAParam.clear();
        
        System.out.println("Estimate parameter...");
        //covert SGAs from list to map(SGAname:index)
        //covert DEGs from list to map(DEGname:index)
        
        Double[] edgeCount = new Double[]{1.0,1.0,1.0,1.0};//initialized with 1 to avoid divided by 0 in calculating probability
        Double[] SGACount = new Double[]{1.0,1.0};//initialized with 1 to avoid divided by 0 in calculating probability
        String[] items;
        for (String strEdge : solidEdgeList) {
            items = strEdge.split(",");
            String SGA = items[0];
            String DEG = items[1];
            //Some SGA, DEG in edge does not exist in driver SGA, or target DEG
            if( (!targetDEGs.contains(DEG))  || (!driverSGAs.contains(SGA) ) ){
                continue;
            }
            int SGAIndx = mapSGA.get(SGA);
            int DEGIndx = mapDEG.get(DEG);
            Arrays.fill(edgeCount, 1.0);//re-initialization
            Arrays.fill(SGACount, 1.0);
            //Need to check if current SGA has been counted before
            boolean SGACounted = false;
            if (mapSGAParam.containsKey(SGA)) {
                SGACounted = true;
            }
            for (int t = 0; t < driverSGATable.size(); t++) {
                int SGAValue = driverSGATable.get(t).get(SGAIndx);
                int DEGValue = targetDEGTable.get(t).get(DEGIndx);
//                edgeCount[SGAValue * 2 + DEGValue] += 1;//If table has missing value then this count is not appropriate
                if (SGAValue == 0 && DEGValue == 0)
                    edgeCount[0] += 1;
                else if (SGAValue == 0 && DEGValue == 1)
                    edgeCount[1] += 1;
                else if (SGAValue == 1 && DEGValue == 0)
                    edgeCount[2] += 1;
                else if (SGAValue == 1 && DEGValue == 1)
                    edgeCount[3] += 1;
                
                if (!SGACounted) {
//                    SGACount[SGAValue] += 1;//If table has missing value then this count is not appropriate
                    if (SGAValue == 0)
                        SGACount[0] += 1;
                    else if (SGAValue == 1)
                        SGACount[1] += 1;
                }
            }
            //normalize count to propability
            Double[] paramOfEdge = new Double[4];
            Double[] paramOfSGA = new Double[2];

            paramOfEdge[0] = edgeCount[0] / (edgeCount[0] + edgeCount[1]); //p(d=0|s=0)
            paramOfEdge[1] = edgeCount[1] / (edgeCount[0] + edgeCount[1]); //p(d=1|s=0)
            paramOfEdge[2] = edgeCount[2] / (edgeCount[2] + edgeCount[3]);//p(d=0|s=1)
            paramOfEdge[3] = edgeCount[3] / (edgeCount[2] + edgeCount[3]);//p(d=1|s=1)
            mapEdgeParam.put(strEdge, paramOfEdge);

            if (!SGACounted) {
                paramOfSGA[0] = SGACount[0] / (SGACount[0] + SGACount[1]);
                paramOfSGA[1] = SGACount[1] / (SGACount[0] + SGACount[1]);
                mapSGAParam.put(SGA, paramOfSGA);
             }
        }
    }
    
    public void inferDriverActivation(){
        System.out.println("Infer driver activation...");
        driverActivationTable = new ArrayList<ArrayList<Double>>();
        
        int tumorSize = targetDEGTable.size();
        int SGASize = driverSGAs.size();

        for (int i = 0; i < tumorSize; i++) {

            ArrayList<Double> driverActivationRow = new ArrayList<Double>();

            for (int j = 0; j < SGASize; j++) {

                if (driverSGATable_ori.get(i).get(j) == 1) {
                    driverActivationRow.add(1.0);
                } else {
                    String SGA = driverSGAs.get(j);
                    Set<String> DEGs = mapSgaDegs.get(SGA);
                    Double logPds0 = 0.0;
                    Double logPds1 = 0.0;
                    for (String DEG : DEGs) {
                        String edge = SGA + "," + DEG;
                        if (!mapEdgeParam.keySet().contains(edge)) {
                            continue;
                        }
                        Double[] edgeProb = new Double[4];

                        edgeProb = mapEdgeParam.get(edge);

                        Double Pd0s0 = edgeProb[0]; //p(d=0|s=0)
                        Double Pd1s0 = edgeProb[1]; //p(d=1|s=0)
                        Double Pd0s1 = edgeProb[2]; //p(d=0|s=1)
                        Double Pd1s1 = edgeProb[3]; //p(d=1|s=1)

                        int DEGIndx = mapDEG.get(DEG);
                        int DEGState = targetDEGTable.get(i).get(DEGIndx);
                        if (DEGState == 1) {
                            logPds0 += log(Pd1s0);
                            logPds1 += log(Pd1s1);
                        } 
                        else if(DEGState == 0) {
                            logPds0 += log(Pd0s0);
                            logPds1 += log(Pd0s1);
                        }
                    }
                    Double[] SGAProb = new Double[2];
                    SGAProb = mapSGAParam.get(SGA);
                    Double Ps0 = SGAProb[0];
                    Double Ps1 = SGAProb[1];
                    Double Ps1d = 1 / (1 + exp(logPds0 - logPds1) * Ps0 / Ps1);
                    driverActivationRow.add(Ps1d);

                }
            }

            driverActivationTable.add(driverActivationRow);
        }

    }
    
    public void thresholding(double T) {
        System.out.println("Thresholding...");
        inferDriverTable = new ArrayList<ArrayList<Integer>>();
        int tumorSize = driverActivationTable.size();
        int SGASize = driverActivationTable.get(0).size();
        for (int i = 0; i < tumorSize; i++) {
            ArrayList<Integer> inferDriverRow = new ArrayList<Integer>();
            for (int j = 0; j < SGASize; j++) {
                if (driverActivationTable.get(i).get(j) >= T) {
                    inferDriverRow.add(1);
                } else {
                    inferDriverRow.add(0);
                }
            }
            inferDriverTable.add(inferDriverRow);
        }
    }
    
    public double compareMatrix() {
        System.out.println("Comare Matrix....");
        int tumorSize = inferDriverTable.size();
        int SGASize = inferDriverTable.get(0).size();
        long totalNum = tumorSize * SGASize;
        long totalChange = 0;
//        add for output value one compare
        long numOldof1 = 0;
        long numNewof1 = 0;
        long numOld1New0 = 0;
        long numNew1Old0 = 0;
//      add ended
        for (int i = 0; i < tumorSize; i++) {
            for (int j = 0; j < SGASize; j++) {
                double vNew = inferDriverTable.get(i).get(j);
                double vOld = driverSGATable.get(i).get(j);
                if (vNew != vOld) { 
                    totalChange++;
                }
                //added begin
                if (vNew == 1){
                    numNewof1 += 1;
                    if(vOld == 0){
                        numNew1Old0 += 1;
                    }
                }
                if (vOld == 1){
                    numOldof1 += 1;
                    if(vNew == 0){
                        numOld1New0 += 1;
                    
                    }
                }
                //added end
              
                
            }
        }
        System.out.println("numOld1, numOld1New0, numNew1, numNew1Old0 :" + numOldof1 + "," + numOld1New0 + "," + numNewof1 + "," + numNew1Old0);
        return (double)totalChange/totalNum;
    }

    public void outputInferActivation( String FileInferDriver){
        try (BufferedWriter bf = new BufferedWriter(new FileWriter(FileInferDriver))) {
            //write in SGA names
            String strSGA = "";
            for (String SGAName : driverSGAs) {
                strSGA +=  "," + SGAName ;
            }
//            System.out.println("strSGA="+strSGA);
//            strSGA = strSGA.substring(0, strSGA.length() - 1);
            bf.write(strSGA + "\n");

            //Write in SGA probabily of each tumor
            for (int i = 0; i < driverActivationTable.size(); i++) {
                String strStates = "";
                for (int j = 0; j < driverSGAs.size(); j++) {
                    strStates +=  "," + driverActivationTable.get(i).get(j) ;
                }
                //remove the last ","
//                strStates = strStates.substring(0, strStates.length() - 1);
                bf.write(tumorNames.get(i) + strStates + "\n");
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    public void outputDriverSGATable(String fileDriverSGATable){
         try (BufferedWriter bf = new BufferedWriter(new FileWriter(fileDriverSGATable))) {
            //write in SGA names
            String strSGA = "";
            for (String SGAName : driverSGAs) {
                strSGA +=  SGAName + ",";
            }
            strSGA = strSGA.substring(0, strSGA.length() - 1);//remove last ,
            bf.write(strSGA + "\n");

            //Write in SGA probabily of each tumor
            for (int i = 0; i < driverSGATable.size(); i++) {
                String strStates = "";
                for (int j = 0; j < driverSGAs.size(); j++) {
                    strStates +=  driverSGATable.get(i).get(j) + ",";
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


          

    
