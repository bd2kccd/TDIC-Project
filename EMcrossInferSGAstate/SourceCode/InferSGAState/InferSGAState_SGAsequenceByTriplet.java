/**
 * This version the output SGA sequence generated according to triplets, not the same with the original SGA matrix.
 * So will have a new version to the output SGA has the same squence with the GtMatrix.
 */
package infersgastate;

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
 * Data prepared using python code RemoveRandomTumorsOfEachCanTypeFromTriplets.py
 */
public class InferSGAState {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.out.println("Working Directory = " + System.getProperty("user.dir"));

        String inputSgaDegTumor = "../DataSource/PANCANsig.rm3tumorPerCan.SgaDegTumor2.csv";
        String inputTumorDEG = "../DataSource/PANCANsig.3tumorPerCan.TargetDEG2.csv";
        String outputSGAState = "../DataSource/PANCANsig.rm3tumorPerCan.SGAState2.csv";
//        String inputSgaDegTumor = "../DataSource/tSgaDegTumor.rm1tumorPerCan.csv";
//        String inputTumorDEG = "../DataSource/tTargetDegs.1tumorPerCan.csv";
//        String outputSGAState = "../DataSource/tSGAState.rm1tumorPerCan.csv";

        InferSGAState app = new InferSGAState();
        Object[] tripletsInfo = app.readTriplets(inputSgaDegTumor);
        Object[] DegInfo = app.readDegFile(inputTumorDEG);
        Object[] SGAStateInfo = app.inferSGAState(tripletsInfo, DegInfo);
        app.output(SGAStateInfo, outputSGAState);

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
        Object[] returnValues = new Object[3];
        returnValues[0] = mapEdgeScores;
        returnValues[1] = mapSgaDegs;
        returnValues[2] = totalTumors;

        return returnValues;
    }

    /**
     * Read in Tumor DEG list and return listOfMapDegs Tumor DEG list is
     * generated from python code generateTumorDegInputFromDegMatrix.py
     *
     * @param DegFile: file name
     * @return List[Map<String:int>] listOfMapDegs
     */
    public Object[] readDegFile(String DegFile) {
        //using fake map instead of list to speed to look up of the DEG  map(DEG:1)
        List<HashMap<String, Integer>> listMapDEGs = new ArrayList<HashMap<String, Integer>>();
        int lineCount = 0;
        String[] items;
        String[] DegNames = {};
        String sCurrentLine;
        List<String> tumorNames = new ArrayList<String>();
        try (BufferedReader br = new BufferedReader(new FileReader(DegFile))) {
            while ((sCurrentLine = br.readLine()) != null) {
                lineCount++;
                if (lineCount == 1) {
                    DegNames = sCurrentLine.split(",");
                    continue;
                }
                items = sCurrentLine.split(",");
                HashMap<String, Integer> mapDEG = new HashMap<String, Integer>();
                tumorNames.add(items[0]);
                for (int i = 1; i < items.length; i++) {
                    if (Integer.parseInt(items[i]) == 1) {
                        mapDEG.put(DegNames[i], 1);
                    }
                }
                listMapDEGs.add(mapDEG);
//                    listDEGs.add(items);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        return new Object[]{listMapDEGs, tumorNames};
    }

    /**
     *
     * @param triplets
     * @param tumorDegs
     * @return
     */
    public Object[] inferSGAState(Object[] triplets, Object[] DegInfo ) {
        Map<String, double[]> mapEdgeScores = (Map<String, double[]>) triplets[0];
        Map<String, Set<String>> mapSGAtoDEGs = (Map<String, Set<String>>) triplets[1];
        int numTripletTumors = (int) triplets[2];
        
        List<HashMap<String, Integer>> listTumorDEGs = (List<HashMap<String, Integer>>) DegInfo[0];
        List<String> tumorNames = (List<String>)DegInfo[1];
        //get listSGA to gurantee the order when do multithread
        ArrayList<String> listDriverSGAs = new ArrayList<String>();
        for (String nameSGA : mapSGAtoDEGs.keySet()) {
            listDriverSGAs.add(nameSGA);
        }

        int sizeDriver = listDriverSGAs.size();
        int sizeTumors = listTumorDEGs.size();
        System.out.println("sizeTumors = " + sizeTumors);

        List<List<Double>> tumorSGAStates = new ArrayList<List<Double>>(Collections.nCopies(sizeTumors, null));

        for (int i = 0; i < sizeTumors; i++) {
            System.out.println("Processing tumor " + i);
            List<Double> SGAStates = new ArrayList<Double>();

            for (int j = 0; j < sizeDriver; j++) {

                String driverSGA = listDriverSGAs.get(j);
                //get the DEGs of drive SGA
                Set<String> driverDEGs = mapSGAtoDEGs.get(driverSGA);

                double logDEG_S0 = 0.0;
                double logDEG_S1 = 0.0;
                int numSGA1 = 0;
                int numSGA0 = 0;
                //check target DEG table of that tumor to see if each driver DEG is 0/1, then calculate the infer SGA state
                for (String driverDEG : driverDEGs) {
                    //get the edge score
                    double[] scores = mapEdgeScores.get(driverSGA + "," + driverDEG);
                    if (scores == null) {
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
                    numSGA1 += scores[3] + scores[2];
                    numSGA0 += scores[1] + scores[0];
                }
                double pS0 = numSGA0 / (double) numTripletTumors;
                double pS1 = numSGA1 / (double) numTripletTumors;
//                double S1_DEGs = 1 / (1 + (exp(logDEG_S0 - logDEG_S1) * pS0 / pS1));
                double S1_DEGs = 1 / (1 + exp(logDEG_S0 - logDEG_S1));
                SGAStates.add(S1_DEGs);
            }
            tumorSGAStates.set(i, SGAStates);
        }
        return new Object[]{listDriverSGAs, tumorSGAStates, tumorNames};
    }

    /**
     * Output results to a file
     *
     * @param Object SGA_tumor_SGAStates[0]: Set<String> setSGA Object
     * SGA_tumor_SGAStates[1]: List<String> listOfTumors Object
     * SGA_tumor_SGAStates[1]: List<List<Double>> tumorSGAStates
     * @param String SGAStateFile :output file name
     */
    public void output(Object[] SGAStateInfo, String SGAStateFile) {

        List<String> listSGA = (List<String>) SGAStateInfo[0];
        List<List<Double>> tumorSGAStates = (List<List<Double>>) SGAStateInfo[1];
        List<String> tumorNames = (List<String>)SGAStateInfo[2];
        try (BufferedWriter bf = new BufferedWriter(new FileWriter(SGAStateFile))) {
            //write in SGA names
            String strSGA = "";
            for (String SGAName : listSGA) {
                strSGA += "," + SGAName;
            }
//            strSGA = strSGA.substring(0, strSGA.length() - 1);
            bf.write(strSGA + "\n");

            //Write in SGA probabily of each tumor
            for (int i = 0; i < tumorSGAStates.size(); i++) {
                String strStates = tumorNames.get(i);
                for (int j = 0; j < listSGA.size(); j++) {
                    strStates +=  "," +tumorSGAStates.get(i).get(j);
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
