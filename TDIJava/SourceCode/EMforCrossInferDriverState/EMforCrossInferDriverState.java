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
//        String fileGtTrainingMatrix = "../DataSource/CrossMatrix/tGtTraining.50.csv";
//        String fileGtTestingMatrix = "../DataSource/CrossMatrix/tGtTesting.50.csv";
//        String fileGeTrainingMatrix = "../DataSource/CrossMatrix/tGeTraining.50.csv";
//        String fileGeTestingMatrix = "../DataSource/CrossMatrix/tGeTesting.50.csv";
//        String fileInferDriver = "../DataSource/CrossMatrix/tInferDriver.50.csv";
//        String fileDriverSGATable = "../DataSource/CrossMatrix/tDriverSGATable.50.csv";

        String fileEdgeList = "../DataSource/PANCANsigEdgeList.csv";
        String fileGtTrainingMatrix = "../DataSource/CrossMatrixDriver/PANCAN.Drivercallpertumor.4468tumors.training.10.csv";
//        String fileGtTestingMatrix = "../DataSource/CrossMatrixDriver/PANCAN.Drivercallpertumor.4468tumors.testing.10.csv";
        String fileGeTrainingMatrix = "../DataSource/CrossMatrixDriver/PANCAN.GeM.4468tumors.training.10.csv";
        String fileGeTestingMatrix = "../DataSource/CrossMatrixDriver/PANCAN.GeM.4468tumors.testing.10.csv";
        String fileInferDriver = "../DataSource/CrossMatrixDriver/PANCAN.Drivercallpertumor.4468tumors.InferDriver.10.withTumorID.csv";
//        String fileDriverSGATable = "../DataSource/CrossMatrixDriver/PANCAN.Drivercallpertumor.4468tumors.DriverSGATable.10.csv";
        
        
        DataReader dataObj = new DataReader(fileEdgeList, fileGtTrainingMatrix, fileGeTrainingMatrix);

        int reRun = 0;
        double T = 0.5;
        do {
            reRun += 1;

            EstimateParams paramObj = new EstimateParams(dataObj.edgeList, dataObj.driverSGAs, dataObj.targetDEGs, dataObj.driverSGATable, dataObj.targetDEGTable);
            InferDriverActivation actObj = new InferDriverActivation(paramObj.mapEdgeParam, paramObj.mapSGAParam,
                    dataObj.targetDEGTable, dataObj.driverSGAs, dataObj.targetDEGs, dataObj.mapSgaDegs);


            actObj.thresholding(T);
            
            actObj.updateInferDriverTable( dataObj.driverSGATable);
            
            double change = actObj.compareMatrix(dataObj.driverSGATable);

            System.out.println("Current T is " + T);
            System.out.println("Change is " + change);
            System.out.println("This is the " + reRun + "th run");
            if (change < 0.001 || T > 1) {
                System.out.println("Total times of run is " + reRun + ". Final cut shreshold is " + T);
                
//                dataObj.readInGtMatrix(fileGtTestingMatrix);
                dataObj.readInGeMatrix(fileGeTestingMatrix);
                InferDriverActivation actObjTesting = new InferDriverActivation(paramObj.mapEdgeParam, paramObj.mapSGAParam,
                    dataObj.targetDEGTable, dataObj.driverSGAs, dataObj.targetDEGs, dataObj.mapSgaDegs);
                
                actObjTesting.outputInferActivation(fileInferDriver, dataObj.tumorNames);
//                dataObj.outputDriverSGATable(fileDriverSGATable); //output contains tumorName, since tumor name is defined in dataObj, so no need to pass in

                break;    
            
            } else {
                dataObj.updateDriverSGATable(actObj.inferDriverTable);
                T += 0.05;
                
            }

        } while (true);
    }
}