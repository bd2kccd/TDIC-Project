/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package eminfer;

import java.util.Random;

/**
 *
 * @author XIM33
 */
public class CrossEMinferDriverState {
    public static void main(String[] args) {

//        String fileEdgeList = "../DataSource/tEdgeList.csv";
//        String fileGtTrainingMatrix = "../DataSource/CrossMatrix/tGtTraining.1.csv";
////        String fileGtTestingMatrix = "../DataSource/CrossMatrix/tGtTesting.1.csv";
//        String fileGeTrainingMatrix = "../DataSource/CrossMatrix/tGeTraining.1.csv";
//        String fileGeTestingMatrix = "../DataSource/CrossMatrix/tGeTesting.1.csv";
//        String fileInferDriver = "./result/tInferDriver.1.withTumorId.run3.csv";
////        String fileDriverSGATable = "./result/tDriverSGATableTesting.1.run2.csv";

//        String fileEdgeList = "../DataSource/PANCANsigEdgeList.csv";
//        String fileGtTrainingMatrix = "../DataSource/PANCAN.GtM.4468tumors.noCan.csv";
//        String fileGeTrainingMatrix = "../DataSource/PANCAN.GeM.4468tumors.csv";
////      String fileGtTestingMatrix = "../DataSource/CrossMatrixDriver/PANCAN.Drivercallpertumor.4468tumors.testing.10.csv";
//        String fileGeTestingMatrix = "../DataSource/pdx_brca_luad_after_normalization_abs_rm.csv";
//        String fileInferDriver = "../DataSource/pdx_brca_luad_after_normalization_abs_inferDriverState.csv";
////        String fileDriverSGATable = "../DataSource/CrossMatrixDriver/PANCAN.Drivercallpertumor.4468tumors.DriverSGATable.10.csv";
////        String fileGeTestingMatrix = "../DataSource/CrossMatrix/PANCAN.GeM.4468tumorsTesting.1.csv";
////        String fileInferDriver = "../DataSource/CrossMatrix/PANCAN.GeM.4468tumorsTesting.1.inferDriverState.csv";
        
//        String fileEdgeList = "C:/Users/XIM33/Documents/EMcrossInferSGAstate/sigEdgeList.csv";
//        String fileGtTrainingMatrix ="C:/Users/XIM33/Documents/EMcrossInferSGAstate/crossSGAnDEGmatrices/PanCancer13tts.SGAmatrix.5097.training.2.csv" ;
//        String fileGeTrainingMatrix = "C:/Users/XIM33/Documents/EMcrossInferSGAstate/crossSGAnDEGmatrices/PanCancer13tts.DEGmatrix.5097.training.2.csv";
//      //String fileGtTestingMatrix = "../DataSource/CrossMatrixDriver/PANCAN.Drivercallpertumor.4468tumors.testing.10.csv";
//        String fileGeTestingMatrix = "C:/Users/XIM33/Documents/EMcrossInferSGAstate/crossSGAnDEGmatrices/PanCancer13tts.DEGmatrix.5097.testing.2.csv";
//        String fileInferDriver = "C:/Users/XIM33/Documents/EMcrossInferSGAstate/crossInferState/PanCancer13tts.inferDriverState.2.csv";

//        String fileEdgeList = "./sigEdgeList.csv";
//        String fileGtTrainingMatrix ="./crossSGAnDEGmatrices/PanCancer13tts.SGAmatrix.5097.training.2.csv" ;
//        String fileGeTrainingMatrix = "./crossSGAnDEGmatrices/PanCancer13tts.DEGmatrix.5097.training.2.csv";
////      String fileGtTestingMatrix = "../DataSource/CrossMatrixDriver/PANCAN.Drivercallpertumor.4468tumors.testing.10.csv";
//        String fileGeTestingMatrix = "./crossSGAnDEGmatrices/PanCancer13tts.DEGmatrix.5097.testing.2.csv";
//        String fileInferDriver = "./crossInferState/PanCancer13tts.inferDriverState.2.csv";
        
//        System.out.println(args[0]);
//        System.out.println(args[1]);
//        System.out.println(args[2]);
//        System.out.println(args[3]);
//        System.out.println(args[4]);
        
        String fileEdgeList = args[0];
        String fileGtTrainingMatrix = args[1];
        String fileGeTrainingMatrix = args[2];
        String fileGeTestingMatrix = args[3];
        String fileInferDriver = args[4];
        
      
        InferData obj = new InferData(fileEdgeList, fileGtTrainingMatrix, fileGeTrainingMatrix);
        
        int reRun = 0;
        do {
            reRun += 1;
            if (reRun == 1)//first run
                obj.estimateParams_firsttime();
            else
                obj.estimateParams_afterfirsttime();
            
            obj.inferDriverActivation();
          
            double change = obj.compareMatrix();

            System.out.println("Change is " + change);
            System.out.println("This is the " + reRun + "th run");
            if (change < 0.005) {
//              if (reRun > 2){
                System.out.println("Total times of run is " + reRun );
                
                obj.readInGeMatrix(fileGeTestingMatrix);
                obj.inferDriverActivation_testingSet();
                obj.outputInferActivation(fileInferDriver);
               
//for test use
//                obj.outputInferActivation(fileInferDriver);
//                obj.outputDriverSGATable(fileDriverSGATable); //output contains tumorName, since tumor name is defined in dataObj, so no need to pass in

                break;    
            
            } else {
                obj.updateDriverSGATable();
                
            }

        } while (true);
    }
    
}
