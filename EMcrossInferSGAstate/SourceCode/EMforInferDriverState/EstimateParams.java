/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package emforinferdriverstate;

//import java.io.IOException;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import java.util.*;

/**
 *
 * @author XIM33
 */
public class EstimateParams {

    Map<String, Double[]> mapEdgeParam = new HashMap<String, Double[]>();
    Map<String, Double[]> mapSGAParam = new HashMap<String, Double[]>();

    public EstimateParams(List<String> edgeList, List<String> driverSGAs,
            List<String> targetDEGs, ArrayList<ArrayList<Integer>> SGATable, ArrayList<ArrayList<Integer>> DEGTable) {
        System.out.println("Estimate parameter...");
        //covert SGAs from list to map(SGAname:index)
        //covert DEGs from list to map(DEGname:index)
        Map<String, Integer> mapSGA = new HashMap<String, Integer>();
        Map<String, Integer> mapDEG = new HashMap<String, Integer>();

        for (int i = 0; i < driverSGAs.size(); i++) {
            mapSGA.put(driverSGAs.get(i), i);
        }
        for (int i = 0; i < targetDEGs.size(); i++) {
            mapDEG.put(targetDEGs.get(i), i);
        }

        Double[] edgeCount = new Double[]{1.0,1.0,1.0,1.0};//initialized with 1 to avoid divided by 0 in calculating probability
        Double[] SGACount = new Double[]{1.0,1.0};//initialized with 1 to avoid divided by 0 in calculating probability
        String[] items;
        for (String strEdge : edgeList) {
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
            for (int t = 0; t < SGATable.size(); t++) {
                int SGAValue = SGATable.get(t).get(SGAIndx);
                int DEGValue = DEGTable.get(t).get(DEGIndx);
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

            paramOfEdge[0] = edgeCount[0] / (edgeCount[0] + edgeCount[1]);
            paramOfEdge[1] = edgeCount[1] / (edgeCount[0] + edgeCount[1]);
            paramOfEdge[2] = edgeCount[2] / (edgeCount[2] + edgeCount[3]);
            paramOfEdge[3] = edgeCount[3] / (edgeCount[2] + edgeCount[3]);
            mapEdgeParam.put(strEdge, paramOfEdge);

            if (!SGACounted) {
                paramOfSGA[0] = SGACount[0] / (SGACount[0] + SGACount[1]);
                paramOfSGA[1] = SGACount[1] / (SGACount[0] + SGACount[1]);
                mapSGAParam.put(SGA, paramOfSGA);
             }
        }
//        //test purpose
//        for (String edge : mapEdgeParam.keySet()) {
//            System.out.print(edge + ":");
//            for (Double d : mapEdgeParam.get(edge)) {
//                System.out.print(d);
//                System.out.print(',');
//            }
//            System.out.println("\n");
//        }
//
//        for (String SGA : driverSGAs) {
//            System.out.print(SGA + ":");
//            for (Double d : mapSGAParam.get(SGA)) {
//                System.out.print(d);
//                System.out.print(',');
//            }
//            System.out.println("\n");
//        }
    }

}
