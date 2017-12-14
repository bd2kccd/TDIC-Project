/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package emforinferdriverstate;

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
//import static java.lang.Math.exp;
//import static java.lang.Math.log;
import java.util.*;

class DataReader {

    public Map<String, Set<String>> mapSgaDegs = new HashMap<String, Set<String>>();
    public ArrayList<String> edgeList = new ArrayList<String>();
    public ArrayList<String> driverSGAs = new ArrayList<String>();
    public ArrayList<String> targetDEGs = new ArrayList<String>();
    public ArrayList<String> tumorNames = new ArrayList<String>();
    public ArrayList<ArrayList<Integer>> targetDEGTable = new ArrayList<ArrayList<Integer>>();
    public ArrayList<ArrayList<Integer>> driverSGATable = new ArrayList<ArrayList<Integer>>();

    private ArrayList<String> tripletSGAs = new ArrayList<String>();
    private ArrayList<String> tripletDEGs = new ArrayList<String>();
    private ArrayList<String> tripletEdgeList = new ArrayList<String>();
    private ArrayList<String> GtMatrixSGAs = new ArrayList<String>();
    private ArrayList<String> GeMatrixDEGs = new ArrayList<String>();
   
    
    public DataReader(String fileTriplets, String fileGtMatrix, String fileGeMatrix) {
        GtMatrixSGAs = getMatrixColumns(fileGtMatrix);
        GeMatrixDEGs = getMatrixColumns(fileGeMatrix);
        readInTriples(fileTriplets);
        UpdateData();
        readInGtMatrix(fileGtMatrix);
        readInGeMatrix(fileGeMatrix);
        

// for test purpose        
//        for (String item : edgeList){
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
        System.out.println("Reading GtMatrix to get SGA names...");
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
    
    public void readInTriples(String fileTriplets) {
        System.out.println("Reading triplet data...");
        try (BufferedReader br = new BufferedReader(new FileReader(fileTriplets))) {
            String sCurrentLine;
            int lineCounter = 0;
            String[] items;

            while ((sCurrentLine = br.readLine()) != null) {
                lineCounter++;
                if (lineCounter == 1) {
                    continue;//This is the title
                } else {
                    items = sCurrentLine.split(",");
//                    String strSGA = items[1].trim();
//                    String strDEG = items[2].trim();
                    String strSGA = items[0].trim();
                    String strDEG = items[1].trim();

                    if (!tripletEdgeList.contains(strSGA + "," + strDEG)) {
                        tripletEdgeList.add(strSGA + "," + strDEG);
                    }

                    if (!tripletSGAs.contains(strSGA)) {
                        tripletSGAs.add(strSGA);
                    }
                    if (!tripletDEGs.contains(strDEG)) {
                        tripletDEGs.add(strDEG);
                    }
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void UpdateData(){
    //update edgeList
    //get intersection of SGAs DEGs of Matrix and triplets and save to driverSGAs and targetDEGs 
    //build mapSgaDegs
    
        //get new edgelist
        String[] items;
        for (String strEdge : tripletEdgeList) {
            items = strEdge.split(",");
            String SGA = items[0];
            String DEG = items[1];
            if (GtMatrixSGAs.contains(SGA) && GeMatrixDEGs.contains(DEG)){
                edgeList.add(strEdge);
            }
        }
        
        //get edgeSGAs, edgeDEGs and build mapSGaDegs
        ArrayList<String> edgeSGAs = new ArrayList<String>();
        ArrayList<String> edgeDEGs = new ArrayList<String>();
        for (String strEdge : edgeList) {
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
    }

    public void readInGtMatrix(String fileGtMatrix) {
        System.out.println("Reading GtMatrix...");
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
                    tumorNames.add(items[0]);
                    ArrayList<Integer> SGAValueRow = new ArrayList<Integer>();
                    for (int j=0; j<SGAIndx.size(); j++){
                        SGAValueRow.add(Integer.parseInt(items[SGAIndx.get(j)]));
                    }

                    driverSGATable.add(SGAValueRow);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void readInGeMatrix(String fileGeMatrix) {
        System.out.println("Reading GeMatrix...");
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
    

        
    public void updateDriverSGATable(ArrayList<ArrayList<Integer>> newSGATable) {
        System.out.println("Update driverSGATable...");
        int tumorSize = driverSGATable.size();
        int SGASize = driverSGATable.get(0).size();
        for (int i = 0; i < tumorSize; i++) {
            for (int j = 0; j < SGASize; j++) {
                Integer v = newSGATable.get(i).get(j);
                driverSGATable.get(i).set(j,v);
            }
            
        }
//        for test purpose
//         for (int i = 0; i < tumorSize; i++) {
//
//
//                System.out.print(driverSGATable.get(i).get(j));
//                System.out.print(" , ");
//            }
//            System.out.println("\n");
//        }

        
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