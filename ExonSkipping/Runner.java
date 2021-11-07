package com.company;


import javax.swing.plaf.synth.Region;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class Runner {
    public static void main(String[] args) {
        // argument parser
        Parser parser = new Parser("-gtf", "-o", "-h");
        parser.setNonParameterOptions("-h");
        parser.parseArguments(args);

        if (parser.getValue("-gtf") == null || parser.getValue("-o") == null || parser.isSet("-h")) {
            System.out.println("Please specify the gtf and output file:");
            System.out.println("-gtf <path to input gtf-file>");
            System.out.println("-o <path to desired output file>");
            System.exit(0);
        }

        // TODO: test paths?
        readGtf(parser.getValue("-gtf"), parser.getValue("-o"));
    }

    // reads gtf and detects ES-SE gene-wise
    private static void readGtf(String inPath, String outPath) {
        try {
            // open reader
            BufferedReader br = new BufferedReader(new FileReader(inPath));
            BufferedWriter bw = new BufferedWriter(new FileWriter(outPath));
            bw.write("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_prots\tSV_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases\n");

            // read content of gtf-file line by line
            String line;
            String chr = "";
            String strand = "";
            String id = "";
            String name = "";
            int numTranscripts = 0;
            ArrayList<String> transcriptIds = new ArrayList<>();
            ArrayList<String> prots = new ArrayList<>();
            HashMap<String, String> information;
            HashMap<String, Protein> idToProtein = new HashMap<>();

            while ((line = br.readLine()) != null) {
                // ignore comments
                if (!line.startsWith("#")) {
                    String[] info = line.split("\t");
                    information = splitInfoColumn(info[8]);

                    // save Protein-ID (if it has one)
                    String protId = information.get("protein_id");
                    if (information.containsKey("protein_id") && !prots.contains(information.get("protein_id"))) {
                        prots.add(protId);
                        idToProtein.put(protId, new Protein(protId));
                    }

                    // new gene starts at current line
                    if (information.containsKey("gene_id") && !id.equals(information.get("gene_id"))) {
                        // work on previous gene & write results to output file
                        if (!chr.equals("")) {
                            String outputStart = id + "\t" + name + "\t" + chr + "\t" + strand + "\t" + prots.size() + "\t" + numTranscripts + "\t";
                            bw.write(handleGene(idToProtein, outputStart));
                        }
                        // clear/refill data structures
                        idToProtein.clear();
                        numTranscripts = 0;
                        prots.clear();
                        strand = info[6];
                        chr = info[0];
                        name = information.get("gene_name");
                        id = information.get("gene_id");
                        transcriptIds.clear();
                    }
                    // new transcript
                    if (information.containsKey("transcript_id") && !transcriptIds.contains(information.get("transcript_id"))) {
                        // count the transcripts
                        numTranscripts++;
                        transcriptIds.add(information.get("transcript_id"));

                    }

                    // add the coordinates of the CDS to the Protein
                    if (info[2].equals("CDS") && protId != null) {
                        idToProtein.get(protId).addCds(new int[]{Integer.parseInt(info[3]), Integer.parseInt(info[4])});
                    }

                }
            }
            // handle very last gene
            if (!chr.equals("")) {
                String outputStart = id + "\t" + name + "\t" + chr + "\t" + strand + "\t" + prots.size() + "\t" + numTranscripts + "\t";
                bw.write(handleGene(idToProtein, outputStart));
            }
            bw.close();
            br.close();

        } catch (FileNotFoundException e) {
            System.out.println("The gtf input file was not found!");
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private static HashMap<String, String> splitInfoColumn(String info) {
        // create Hash Map out of the last gtf-Column
        String[] pieces = info.split(";");
        HashMap<String, String> features = new HashMap<>();
        String feature;
        String value;
        for (String piece : pieces) {
            piece = piece.strip();
            feature = piece.split(" ")[0];
            value = piece.split(" ")[1];
            value = value.replace("\"", "");
            features.put(feature, value);
        }
        return features;
    }

    private static String handleGene(HashMap<String, Protein> idToProtein, String outputStart) {
        // get the information for one gene and write into output file
        // @outputStart contains the first 6 columns of the output file (identical for the whole gene)
        StringBuilder output = new StringBuilder();

        // create a Map that maps all the introns to the protein IDs that use them
        HashMap<String, String> intronToId = new HashMap<>();
        for (Protein p : idToProtein.values()) {
            ArrayList<String> introns = p.getIntrons().getRegionsAsString();
            for (String i : introns) {
                if (!intronToId.containsKey(i)) {
                    intronToId.put(i, p.getId());
                } else {
                    String old = intronToId.get(i);
                    intronToId.replace(i, old, old + "|" + p.getId());
                }
            }
        }

        // Set computations: iterate over all introns
        for (String intron : intronToId.keySet()) {
            // Set SV: all proteins that contain the candidate intron
            HashSet<String> sv = new HashSet<>(Arrays.asList(intronToId.get(intron).split("\\|")));

            // wt start and end sets: iterate over all proteins to find the ones that have an intron with the specific start/end
            int start = Integer.parseInt(intron.split(":")[0]);
            int end = Integer.parseInt(intron.split(":")[1]);
            HashSet<String> wt_start = new HashSet<>();
            HashSet<String> wt_end = new HashSet<>();

            for (Protein p : idToProtein.values()) {
                if (p.checkIntronStart(start)) wt_start.add(p.getId());
                if (p.checkIntronEnd(end)) wt_end.add(p.getId());
            }

            // intersect wt_start and wt_end
            HashSet<String> intersect = new HashSet<>(wt_start);
            intersect.retainAll(wt_end);

            // remove all values of sv from intersect
            intersect.removeAll(sv);

            // evaluate the output
            /* Output
             * start:end is the sv
             * get subregions to get the wt (introns within sv intron: more than 1!)
             * intersect contains all the wt_prots
             * sv contains all the sv_prots
             * compute with wt: min/max skipped exon/bases
             * */
            if (intersect.size() > 0) {
                // exon skipping event!
                output.append(outputStart);
                // sv
                output.append(start).append(":").append(end).append("\t");

                // wt
                StringBuilder wt_builder = new StringBuilder();
                int min_skipped_exons = end - start;
                int max_skipped_exons = 0;
                int min_skipped_bases = end - start;
                int max_skipped_bases = 0;

                for (Protein p : idToProtein.values()) {
                    RegionVector rv = p.getIntronsWithin(new int[]{start, end});
                    String s = null;
                    if (rv != null) s = rv.regionsAsPrintString();
                    // only add new introns that are not the big one and ignore cases with only one intron
                    if (s != null && !wt_builder.toString().contains(s) && !s.equals(start + ":" + end) && s.split(":").length > 2) {
                        wt_builder.append("|").append(s);

                        // compute min/max skipped exon/bases (and check if they change min/max over all proteins)
                        RegionVector introns = new RegionVector(new ArrayList<>());
                        String[] intronArray = s.split("\\|");
                        for (String i : intronArray) {
                            String[] startEnd = i.split(":");
                            int[] region = new int[]{Integer.parseInt(startEnd[0]), Integer.parseInt(startEnd[1])};
                            introns.addRegion(region);
                        }
                        int skipped_exons = introns.size() - 1;
                        min_skipped_exons = Math.min(skipped_exons, min_skipped_exons);
                        max_skipped_exons = Math.max(skipped_exons, max_skipped_exons);

                        int skipped_bases = introns.getInverseLength();
                        min_skipped_bases = Math.min(min_skipped_bases, skipped_bases);
                        max_skipped_bases = Math.max(max_skipped_bases, skipped_bases);
                    }
                }
                output.append(wt_builder.toString().substring(1)).append("\t");

                // wt_prots
                StringBuilder wt_prots = new StringBuilder();
                for (String s : intersect) {
                    wt_prots.append("|").append(s);
                }
                output.append(wt_prots.toString().substring(1)).append("\t");

                // sv_prots
                StringBuilder sv_prots = new StringBuilder();
                for (String s : sv) {
                    sv_prots.append("|").append(s);
                }
                output.append(sv_prots.toString().substring(1)).append("\t");

                output.append(min_skipped_exons).append("\t").append(max_skipped_exons).append("\t").append(min_skipped_bases).append("\t").append(max_skipped_bases).append("\n");
            }
        }
        return output.toString();
    }
}
