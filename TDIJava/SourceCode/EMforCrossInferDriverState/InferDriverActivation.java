/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package emforcrossinferdriverstate;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author XIM33
 */
public class InferDriverActivation {

    public ArrayList<ArrayList<Double>> driverActivationTable;
    public ArrayList<ArrayList<Integer>> inferDriverTable;
    
    private ArrayList<String> driverSGAs = new ArrayList<String>();
    
    public InferDriverActivation(Map<String, Double[]> mapEdgeParam, Map<String, Double[]> mapSGAParam,
             ArrayList<ArrayList<Integer>> targetDEGTable,
            ArrayList<String> driverSGAs, ArrayList<String> targetDEGs, Map<String, Set<String>> mapSgaDegs) {
        this.driverSGAs = driverSGAs;
        System.out.println("Infer driver activation...");
        driverActivationTable = new ArrayList<ArrayList<Double>>();
        Map<String, Integer> mapDEGIndx = new HashMap<String, Integer>();
        for (int i = 0; i < targetDEGs.size(); i++) {
            mapDEGIndx.put(targetDEGs.get(i), i);
        }

        int tumorSize = targetDEGTable.size();
        int SGASize = driverSGAs.size();

        for (int i = 0; i < tumorSize; i++) {

            ArrayList<Double> driverActivationRow = new ArrayList<Double>();

            for (int j = 0; j < SGASize; j++) {

//                if (driverSGATable.get(i).get(j) == 1) {
//                    driverActivationRow.add(1.0);
//                } else {
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

                        Double Pd0s0 = edgeProb[0]; //s0d0
                        Double Pd1s0 = edgeProb[1]; //s0d1
                        Double Pd0s1 = edgeProb[2]; //s1d0
                        Double Pd1s1 = edgeProb[3]; //s1d1

                        int DEGIndx = mapDEGIndx.get(DEG);
                        int DEGState = targetDEGTable.get(i).get(DEGIndx);
                        if (DEGState == 1) {
                            logPds0 += log(Pd1s0);
                            logPds1 += log(Pd1s1);
                        } else {
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

//                }
            }

            driverActivationTable.add(driverActivationRow);
        }

//        //for test purpose
//        for (int i = 0; i < tumorSize; i++) {
//
//            for (int j = 0; j < SGASize; j++) {
////                System.out.print(driverSGAs.get(j) + "," );
//                System.out.print(driverActivationTable.get(i).get(j));
//                System.out.print(" , ");
//            }
//            System.out.println("\n");
//        }

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
//        //for test purpose
//
//        for (int i = 0; i < tumorSize; i++) {
//
//            for (int j = 0; j < SGASize; j++) {
//
//                System.out.print(inferDriverTable.get(i).get(j));
//                System.out.print(" , ");
//            }
//            System.out.println("\n");
//        }

    }
    
    public void updateInferDriverTable(ArrayList<ArrayList<Integer>> driverSGATable){
        System.out.println("Update InferDriverTable with original SGA = 1");
        int tumorSize = inferDriverTable.size();
        int SGASize = inferDriverTable.get(0).size();
          for (int i = 0; i < tumorSize; i++) {
            for (int j = 0; j < SGASize; j++) {
                if(driverSGATable.get(i).get(j) == 1){
                    inferDriverTable.get(i).set(j, 1);
                }
            }
          }
//        //for test purpose
//
//        for (int i = 0; i < tumorSize; i++) {
//            for (int j = 0; j < SGASize; j++) {
//
//                System.out.print(inferDriverTable.get(i).get(j));
//                System.out.print(" , ");
//            }
//            System.out.println("\n");
//        }
          
    }

    public double compareMatrix(ArrayList<ArrayList<Integer>> SGATable) {
        System.out.println("Comare Matrix....");
        int tumorSize = inferDriverTable.size();
        int SGASize = inferDriverTable.get(0).size();
        long totalNum = tumorSize * SGASize;
        long totalChange = 0;
        for (int i = 0; i < tumorSize; i++) {
            for (int j = 0; j < SGASize; j++) {
                double vNew = inferDriverTable.get(i).get(j);
                double vOld = SGATable.get(i).get(j);
                if (vNew != vOld) { 
                    totalChange++;
                }
            }
        }
        return (double)totalChange/totalNum;
    }
    
    public void outputInferActivation( String FileInferDriver, ArrayList<String> tumorNames){
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
    
}
