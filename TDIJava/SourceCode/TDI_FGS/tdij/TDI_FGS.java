/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tdij;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import tetrad.data.*;
//
//import java.io.BufferedReader;
//import java.io.BufferedWriter;
//import java.io.FileReader;
//import java.io.FileWriter;
//import java.io.IOException;
//import static java.lang.Math.exp;
//import static java.lang.Math.log;
//
////import java.io.ObjectInputStream;
////import java.text.NumberFormat;
//import java.util.*;
//import java.util.concurrent.ExecutorService;
//import java.util.concurrent.Executors;
//import java.util.concurrent.TimeUnit;

//import tetrad.cli.data.ContinuousDataReader;
import tetrad.io.TabularContinuousDataReader;
import tetrad.io.DataReader;
import tetrad.cli.validation.DataValidation;
import tetrad.cli.validation.TabularContinuousData;
import tetrad.data.DataSet;
import tetrad.data.CovarianceMatrixOnTheFly;
import tetrad.graph.Graph;
import tetrad.search.Fgs;
import tetrad.search.SemBicScore;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Author : Jeremy Espino MD Created 2/12/16 9:44 AM
 */
public class TDI_FGS {

    public static void main(String[] args) throws Exception {

        // set path to Retention data
        Path dataFile = Paths.get("../DataSource/", "Retention.txt");

        Character delimiter = '\t';

        // perform data validation
        // note: assuming data has unique variable names and does not contain zero covariance pairs
        DataValidation dataValidation = new TabularContinuousData(dataFile, delimiter);
        if (!dataValidation.validate(System.err, true)) {
            System.exit(-128);
        }

        DataReader dataReader = new TabularContinuousDataReader(dataFile, delimiter);
        DataSet dataSet = dataReader.readInData();
        ICovarianceMatrix covariances = new CovarianceMatrixOnTheFly(dataSet);
        Fgs fgs = new Fgs(new SemBicScore(covariances));
        fgs.setOut(System.out);
//        fgs.setDepth(-1);
//        fgs.setIgnoreLinearDependent(false);
        fgs.setPenaltyDiscount(4.0);
        fgs.setNumPatternsToStore(0);  // always set to zero
        fgs.setFaithfulnessAssumed(true);
        fgs.setVerbose(true);

        Graph graph = fgs.search();
        System.out.println();
        System.out.println(graph.toString().trim());
        System.out.flush();

        try {
            
            FileOutputStream fout = new FileOutputStream("./output.txt", true);
            PrintStream out = new PrintStream(fout);
            out.println(graph.toString().trim());
        } catch (IOException ex) {
            System.out.println("There was a problem creating/writing to the temp file");
            ex.printStackTrace();
        }
        
    }

}
