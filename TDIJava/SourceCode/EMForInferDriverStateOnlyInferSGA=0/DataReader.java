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
//import java.io.BufferedWriter;
import java.io.FileReader;
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

    private ArrayList<String> tripletSGAs = new ArrayList<String>();//tripletSGAs, driverSGAs have same SGAs, but different order
    private ArrayList<String> tripletDEGs = new ArrayList<String>();

    public DataReader(String fileTriplets, String fileGtMatrix, String fileGeMatrix) {
        readInTriples(fileTriplets);
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

                    if (!edgeList.contains(strSGA + "," + strDEG)) {
                        edgeList.add(strSGA + "," + strDEG);
                    }

                    if (!tripletSGAs.contains(strSGA)) {
                        tripletSGAs.add(strSGA);
                    }
                    if (!tripletDEGs.contains(strDEG)) {
                        tripletDEGs.add(strDEG);
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
                        if (tripletSGAs.contains(items[i])) {
                            driverSGAs.add(items[i]);
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
                        if (tripletDEGs.contains(items[i])) {
                            targetDEGs.add(items[i]);
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
//            for (int j = 0; j < SGASize; j++) {
//
//                System.out.print(driverSGATable.get(i).get(j));
//                System.out.print(" , ");
//            }
//            System.out.println("\n");
//        }

        
    }
}