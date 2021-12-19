package com.company;

import java.io.*;
import java.util.HashMap;

public class Main {

    public static void main(String[] args) {
        String annotFile = "../plottingData/ebna_hisat.annot";
        String geneFile = "../plottingData/ebna_hisat.tsv";
        String plottingFile = "../plottingData/ebna_hisat.plotting";
        int counts = 53898152;
        // counts = (#totalReads - (2 * #readPairs)) + #readPairs
        // nookaew_cm: 5493225, hes_star: 23477686, ebna_hisat: 53898152

        double perMillionScalingFactor = counts / 1000000.0;

        HashMap<String, Integer> geneToLength = new HashMap<>();
        HashMap<String, Integer> geneToCounts = new HashMap<>();
        HashMap<String, Integer> geneToCountsPcrZero = new HashMap<>();

        try {
            BufferedReader geneReader = new BufferedReader(new FileReader(geneFile));
            String line;
            while ((line = geneReader.readLine()) != null) {
                if (!line.startsWith("genes")) {
                    line = line.replace("\n", "");
                    String[] info = line.split("\t");
                    geneToLength.put(info[0], Integer.parseInt(info[1]));
                    geneToCounts.put(info[0], 0);
                    geneToCountsPcrZero.put(info[0], 0);
                }
            }
            geneReader.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

        try {
            BufferedReader annotReader = new BufferedReader(new FileReader(annotFile));
            String line;
            while ((line = annotReader.readLine()) != null) {
                if (!line.startsWith("read_id")) {
                    line = line.replace("\n", "");
                    String[] info = line.split("\t");
                    String[] geneIds = info[2].split("\\|");
                    for (String gid : geneIds) {
                        geneToCounts.replace(gid, geneToCounts.get(gid) + 1);
                        if (info[3].equals("0")) {
                            geneToCountsPcrZero.replace(gid, geneToCountsPcrZero.get(gid) + 1);
                        }
                    }
                }
            }
            annotReader.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

        // write plotting file
        try{
            BufferedWriter bw = new BufferedWriter(new FileWriter(plottingFile));
            bw.write("gene_id\trpkm\trpkm_pcrZero\tlog2rpkm\tlog2rpkm_pcrZero\n");
            for(String gene:geneToLength.keySet()){
                double len = geneToLength.get(gene)/1000.0;
                double rpkm = (geneToCounts.get(gene)/perMillionScalingFactor)/len;
                double rpkm_pcrZero = (geneToCountsPcrZero.get(gene)/perMillionScalingFactor)/len;
                double log2rpkm = Math.log(rpkm)/Math.log(2);
                double log2rpkm_pcrZero = Math.log(rpkm_pcrZero)/Math.log(2);
                // round to 2 decimal points
                //rpkm = Math.round(rpkm*100)/100.0;
                //rpkm_pcrZero = Math.round(rpkm_pcrZero*100)/100.0;
                bw.write(gene+"\t"+rpkm+"\t"+rpkm_pcrZero+"\t"+log2rpkm+"\t"+log2rpkm_pcrZero+"\n");
            }
            bw.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
