package com.company;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class Main {
    private static ArrayList<String> genes;
    private static ArrayList<ArrayList<Double>> samples;
    private static String sampleNames;

    public static void main(String[] args) {
        Parser parser = new Parser("-gtf", "-counts", "-out");
        parser.parseArguments(args);

        String gtfPath = parser.getValue("-gtf");
        String countsFile = parser.getValue("-counts");
        String outPath = parser.getValue("-out");

        if (gtfPath == null || countsFile == null || outPath == null) {
            System.out.println("Usage:");
            System.out.println("-gtf <file.gtf>\tgtf-file for desired organism");
            System.out.println("-counts <counts.csv>\tcsv-file containing the counts (ENSEMBL-IDs!)");
            System.out.println("-out <output-file>\tpath to output-file (csv)");
            System.exit(0);
        }

        HashMap<String, Double> geneToLengthKb = parseGtf(gtfPath);
        computeTpm(countsFile, geneToLengthKb);
        writeOutput(outPath);
    }

    // read gtf and map gene-IDs to the length of the corresponding genes (IN KILOBASES!)
    private static HashMap<String, Double> parseGtf(String path) {
        HashMap<String, Double> geneToLength = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("#")) {
                    String[] parts = line.split("\t");
                    if (parts[2].equals("gene")) {
                        double length = Integer.parseInt(parts[4]) - Integer.parseInt(parts[3]);
                        // kilobases
                        length = length / 1000.0;
                        String geneId = parts[8].substring(parts[8].indexOf("gene_id \"") + 9);
                        geneId = geneId.substring(0, geneId.indexOf("\""));
                        geneToLength.put(geneId, length);
                    }
                }
            }

        } catch (FileNotFoundException e) {
            System.out.println("The gtf-file " + path + " was not found!");
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return geneToLength;
    }

    // read count file, copmute TPM values and save them
    private static void computeTpm(String path, HashMap<String, Double> geneToLengthKb) {
        samples = new ArrayList<>();
        genes = new ArrayList<>();

        // read count file and save in data structure (already divide by gene-length!)
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            // header (1st line)
            String line = br.readLine();
            line = line.replace("\n", "");
            line = line.replace("\"", "");
            String[] parts = line.split(",");
            StringBuilder sb = new StringBuilder();
            for (int i = 1; i < parts.length; i++) {
                samples.add(new ArrayList<>());
                sb.append(",").append(parts[i]);
            }
            sampleNames = sb.toString();

            while ((line = br.readLine()) != null) {
                line = line.replace("\n", "");
                line = line.replace("\"", "");
                parts = line.split(",");

                // count line
                Double geneLength = geneToLengthKb.getOrDefault(parts[0].split("\\.")[0], null);
                if (geneLength == null) {
                    System.out.println("Warning: The gene-ID " + parts[0].split("\\.")[0] + " from the count file is not in the gtf (line is omitted).");
                } else {
                    genes.add(parts[0]);
                    for (int i = 1; i < parts.length; i++) {
                        double count = Integer.parseInt(parts[i]);
                        samples.get(i - 1).add(count / geneLength);
                    }
                }

            }

        } catch (FileNotFoundException e) {
            System.out.println("The counts file " + path + " was not found!");
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // compute TPM values
        // for each sample: build sum of values, divide by 1,000,000 => "per-million-scaling-factor"
        // then: divide the numbers by the per-million-scaling-factor and round to 2 digits
        for (int i = 0; i < samples.size(); i++) {
            double pmsf = 0.0;
            for (double d : samples.get(i)) {
                pmsf += d;
            }
            pmsf = pmsf / 1000000.0;
            for (int j = 0; j < samples.get(i).size(); j++) {
                double val = samples.get(i).get(j) / pmsf;
                samples.get(i).set(j, roundToTwoDigits(val));
            }
        }
    }

    public static double roundToTwoDigits(double num) {
        return Math.round(num * 100) / 100.0;
    }

    // write output file
    private static void writeOutput(String path) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(path));
            bw.write("ENSEMBL" + sampleNames + "\n");
            for (int i = 0; i < genes.size(); i++) {
                StringBuilder line = new StringBuilder(genes.get(i));
                for (int j = 0; j < samples.size(); j++) {
                    ArrayList<Double> sample = samples.get(j);
                    line.append(",").append(sample.get(i));
                }
                line.append("\n");
                bw.write(line.toString());
            }
            bw.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
