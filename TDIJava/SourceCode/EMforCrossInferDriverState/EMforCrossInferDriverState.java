/*
 This version, infer SGA state only on Gt=0
 */
package emforcrossinferdriverstate;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

class EMforCrossInferDriverState {

    public static void main(String[] args) {

//        String fileEdgeList = "../DataSource/tSgaDegTumor.csv";
//        String fileGtTrainingMatrix = "../DataSource/CrossMatrix/tGtTraining.10.csv";
//        String fileGtTestingMatrix = "../DataSource/CrossMatrix/tGtTesting.10.csv";
//        String fileGeTrainingMatrix = "../DataSource/CrossMatrix/tGeTraining.10.csv";
//        String fileGeTestingMatrix = "../DataSource/CrossMatrix/tGeTesting.10.csv";
//        String fileInferDriver = "../DataSource/CrossMatrix/tInferDriver.10.csv";
//        String fileDriverSGATable = "../DataSource/CrossMatrix/tDriverSGATable.10.csv";

        String fileEdgeList = "../DataSource/PANCANsigEdgeList.csv";
        String fileGtTrainingMatrix = "../DataSource/CrossMatrix/PANCAN.GtM.4468tumorsTraining.10.csv";
        String fileGtTestingMatrix = "../DataSource/CrossMatrix/PANCAN.GtM.4468tumorsTesting.10.csv";
        String fileGeTrainingMatrix = "../DataSource/CrossMatrix/PANCAN.GeM.4468tumorsTraining.10.csv";
        String fileGeTestingMatrix = "../DataSource/CrossMatrix/PANCAN.GeM.4468tumorsTesting.10.csv";
        String fileInferDriver = "../DataSource/CrossMatrix/PANCAN.InferDriver.10.csv";
//        String fileDriverSGATable = "../DataSource/CrossMatrix/PANCAN.DriverSGATable.10.csv";
        
        
        DataReader dataObj = new DataReader(fileEdgeList, fileGtTrainingMatrix, fileGeTrainingMatrix);
//        System.out.println("Original driverSGATable");
//        for (int i = 0; i < dataObj.driverSGATable.size(); i++) {
//            for (int j = 0; j < dataObj.driverSGATable.get(i).size(); j++) {
//                System.out.print(dataObj.driverSGATable.get(i).get(j));
//                System.out.print(',');
//            }
//            System.out.println('\n');
//
//        }

        int reRun = 0;
        double T = 0.5;
        do {
// System.out.println("newSGATable");
//for (int i = 0; i < dataObj.driverSGATable.size(); i++) {
//            for (int j = 0; j < dataObj.driverSGATable.get(i).size(); j++) {
//                System.out.print(dataObj.driverSGATable.get(i).get(j));
//                System.out.print(',');
//            }
//            System.out.println('\n');
//
//        }
            reRun += 1;

            EstimateParams paramObj = new EstimateParams(dataObj.edgeList, dataObj.driverSGAs, dataObj.targetDEGs, dataObj.driverSGATable, dataObj.targetDEGTable);
//        //test purpose
//        for (String edge : paramObj.mapEdgeParam.keySet()) {
//            System.out.print(edge + ":");
//            for (Double d : paramObj.mapEdgeParam.get(edge)) {
//                System.out.print(d);
//                System.out.print(',');
//            }
//            System.out.println("\n");
//        }
//
//        for (String SGA : dataObj.driverSGAs) {
//            System.out.print(SGA + ":");
//            for (Double d :paramObj.mapSGAParam.get(SGA)) {
//                System.out.print(d);
//                System.out.print(',');
//            }
//            System.out.println("\n");
//        }            

            InferDriverActivation actObj = new InferDriverActivation(paramObj.mapEdgeParam, paramObj.mapSGAParam,
                    dataObj.targetDEGTable, dataObj.driverSGAs, dataObj.targetDEGs, dataObj.mapSgaDegs);

//        //for test purpose
//        System.out.println("driverActivationTable");
//        for (int i = 0; i < actObj.driverActivationTable.size(); i++) {
//
//            for (int j = 0; j < actObj.driverActivationTable.get(i).size(); j++) {
//
//                System.out.print(actObj.driverActivationTable.get(i).get(j));
//                System.out.print(" , ");
//            }
//            System.out.println("\n");
//        }            


            actObj.thresholding(T);
            
            actObj.updateInferDriverTable( dataObj.driverSGATable);
            
            double change = actObj.compareMatrix(dataObj.driverSGATable);

            System.out.println("Current T is " + T);
            System.out.println("Change is " + change);
            System.out.println("This is the " + reRun + "th run");
            if (change < 0.001 || T > 1) {
                System.out.println("Total times of run is " + reRun + ". Final cut shreshold is " + T);
                
                dataObj.readInGtMatrix(fileGtTestingMatrix);
                dataObj.readInGeMatrix(fileGeTestingMatrix);
                InferDriverActivation actObjTesting = new InferDriverActivation(paramObj.mapEdgeParam, paramObj.mapSGAParam,
                    dataObj.targetDEGTable, dataObj.driverSGAs, dataObj.targetDEGs, dataObj.mapSgaDegs);
                
                actObjTesting.outputInferActivation(fileInferDriver);
//                dataObj.outputDriverSGATable(fileDriverSGATable);

                break;    
            
            } else {
                dataObj.updateDriverSGATable(actObj.inferDriverTable);
                T += 0.05;
                
            }

        } while (true);
    }
}
