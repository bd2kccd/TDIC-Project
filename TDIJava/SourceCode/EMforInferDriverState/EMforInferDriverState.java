/*
 This version, infer SGA state only on Gt=0
 */
package emforinferdriverstate;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

class EMforInferDriverState {

    public static void main(String[] args) {
        String fileTriplets = "../DataSource/PANCANsigEdgeList.csv";
        String fileGtMatrix = "../DataSource/Drivercallpertumor.csv";
        String fileGeMatrix = "../DataSource/PANCAN.GeM.4468tumors.csv";
        String fileInferState = "../DataSource/Drivercallpertumor.EmInferState.csv";
        String fileDriverSGATable = "../DataSource/Drivercallpertumor.driverSGATable.csv";

        DataReader dataObj = new DataReader(fileTriplets, fileGtMatrix, fileGeMatrix);

        int reRun = 0;
        double T = 0.5;
        do {
            reRun += 1;

            EstimateParams paramObj = new EstimateParams(dataObj.edgeList, dataObj.driverSGAs, dataObj.targetDEGs, dataObj.driverSGATable, dataObj.targetDEGTable);

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
                actObj.outputInferActivation(fileInferState);
                dataObj.outputDriverSGATable(fileDriverSGATable);

                break;    
            
            } else {
                dataObj.updateDriverSGATable(actObj.inferDriverTable);
                T += 0.05;
                
            }

        } while (true);
    }
}
