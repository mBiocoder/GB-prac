package com.company;

import augmentedTree.Interval;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Gene implements Interval {
    private int start;
    private int stop;
    private ArrayList<String> transcriptIds;
    private HashMap<String, RegionVector> transcripts;
    private String geneId;
    private String biotype;
    private String chr;
    private boolean strand;

    private HashMap<String, Integer> transcriptToCount;
    // saves skipped exons and maps them to the transcript that INCLUDES them
    private HashMap<String, String> skippedExons;
    private HashMap<String, int[]> skippedExonCounts;
    private HashMap<String, RegionVector> cds;
    private HashMap<String, int[][]> skippedExonToNeighbors;


    public Gene(int start, int stop, String id, String biotype, String chromosome, boolean strand) {
        this.start = start;
        this.stop = stop;
        this.geneId = id;
        this.biotype = biotype;
        this.chr = chromosome;
        this.transcriptIds = new ArrayList<>();
        this.transcripts = new HashMap<>();
        this.strand = strand;

        transcriptToCount = new HashMap<>();
        this.cds = new HashMap<>();
    }

    @Override
    public int getStart() {
        return this.start;
    }

    @Override
    public int getStop() {
        return this.stop;
    }

    public String getChromosome() {
        return this.chr;
    }

    public void countRead(RegionVector sr1, RegionVector sr2) {
        ArrayList<String> foundExons = new ArrayList<>();
        ArrayList<String> countedExclusive = new ArrayList<>();
        ArrayList<String> countedInclusive = new ArrayList<>();

        if (this.matchPlausibilityLevel(sr1, sr2).equals("3")) {
            for (String skipped : skippedExons.keySet()) {
                int a = Integer.parseInt(skipped.split("-")[0]);
                int b = Integer.parseInt(skipped.split("-")[1]) - 1;

                // only count a read as inclusive when it covers a skipped exon!
                if (sr1.containsRegionPartially(a, b) || sr2.containsRegionPartially(a, b)) {
                    foundExons.add(skipped);
                }
            }

            //########## pre-computations ##############################################################################
            String tr1 = transcriptIds.get(0);
            String tr2 = transcriptIds.get(1);

            HashMap<String, HashMap<String, RegionVector>> trToExonNeighborMap = new HashMap<>();
            HashMap<String, HashMap<String, RegionVector>> trToExonNeighborMapInclusive = new HashMap<>();
            HashMap<String, ArrayList<String>> trToContainedSkippedExons = new HashMap<>();
            // map the transcripts to all the skipped exons that they contain
            for (String exon : skippedExons.keySet()) {
                if (!trToContainedSkippedExons.containsKey(skippedExons.get(exon))) {
                    trToContainedSkippedExons.put(skippedExons.get(exon), new ArrayList<>());
                }
                trToContainedSkippedExons.get(skippedExons.get(exon)).add(exon);
            }
            // find the neighbors of the skipped exons on the transcripts that contain them and that skip them
            for (String tr : trToContainedSkippedExons.keySet()) {
                trToExonNeighborMap.put(tr, transcripts.get(tr).getNeighborsOfExons(trToContainedSkippedExons.get(tr)));
                String otherTr = tr1;
                if (tr.equals(tr1)) otherTr = tr2;
                // map: transcript that skips some exons to hash map that maps exon to neighboring regions on the other transcript
                trToExonNeighborMapInclusive.put(tr, transcripts.get(otherTr).getNeighborsOfSkippedExons(trToContainedSkippedExons.get(tr)));
            }

            // pre-computations
            boolean evidenceAgainstTr1 = false;
            boolean evidenceAgainstTr2 = false;
            boolean evidenceForTr1 = false;
            boolean evidenceForTr2 = false;
            HashMap<String, Boolean> exonSidesCovered = new HashMap<>();
            for (String exon : skippedExons.keySet()) {
                // get the neighboring regions of the skipped exon on the transcript that CONTAINS the exon
                RegionVector neighboring = trToExonNeighborMap.get(skippedExons.get(exon)).get(exon);
                // and the neigboring regions of the skipped exon on the transcript that SKIPS the exon
                RegionVector neighborsOfSkipped = trToExonNeighborMapInclusive.get(skippedExons.get(exon)).get(exon);
                // check if read covers the neighboring regions but exceeds them
                int readLeftStart = sr1.getLastRegion()[0];
                int readLeftEnd = sr1.getEnd();
                int readRightEnd = sr2.getFirstRegion()[1];
                int readRightStart = sr2.getStart();
                if (sr2.getStart() < sr1.getStart()) {
                    readLeftStart = sr2.getLastRegion()[0];
                    readLeftEnd = sr2.getEnd();
                    readRightEnd = sr1.getFirstRegion()[1];
                    readRightStart = sr1.getStart();
                }

                // is there evidence against a transcript that INCLUDES a skipped exon?
                boolean coversLeft = (neighboring.getStart() <= readLeftEnd) && (readLeftEnd <= neighboring.getFirstRegion()[1]); // read covers left neighbor
                boolean leftProblem = coversLeft && neighboring.getStart() > readLeftStart; // read exceeds left neighbor
                boolean coversRight = (neighboring.getEnd() >= readRightStart) && (readRightStart >= neighboring.getLastRegion()[0]); // read covers right neighbor
                boolean rightProblem = coversRight && neighboring.getEnd() < readRightEnd; // read exceeds right neighbor
                if ((coversLeft && coversRight) && (leftProblem || rightProblem)) {
                    if (skippedExons.get(exon).equals(tr1)) evidenceAgainstTr1 = true;
                    else evidenceAgainstTr2 = true;
                }
                exonSidesCovered.put(exon, (coversLeft && coversRight));

                // check if there is evidence for a transcript that INCLUDES a skipped exon
                coversLeft = (neighborsOfSkipped.getStart() <= readLeftEnd) && (readLeftEnd <= neighborsOfSkipped.getFirstRegion()[1]); // read covers left neighbor
                leftProblem = coversLeft && neighborsOfSkipped.getStart() > readLeftStart; // read exceeds left neighbor
                coversRight = (neighborsOfSkipped.getEnd() >= readRightStart) && (readRightStart >= neighborsOfSkipped.getLastRegion()[0]); // read covers right neighbor
                rightProblem = coversRight && neighborsOfSkipped.getEnd() < readRightEnd; // read exceeds right neighbor
                if ((coversLeft && coversRight) && (leftProblem || rightProblem)) {
                    if (skippedExons.get(exon).equals(tr1)) evidenceForTr2 = true;
                    else evidenceForTr1 = true;
                }
            }


            // ######################## count the read according to the precomputations ################################
            for (String exon : skippedExons.keySet()) {
                int a = Integer.parseInt(exon.split("-")[0]);
                int b = Integer.parseInt(exon.split("-")[1]) - 1;
                if (foundExons.contains(exon)) {
                    // read is counted as INCLUSIVE for that exon
                    if (!countedInclusive.contains(exon)) {
                        skippedExonCounts.get(exon)[0]++;
                        countedInclusive.add(exon);
                    }
                    // check: is there another exon in the direction of the other read that is not covered? => counted as well
                    boolean sr1Smaller = sr1.getStart() <= sr2.getStart();
                    ArrayList<String> between;
                    if (sr1.containsRegionPartially(a, b)) {
                        if (sr1Smaller) {
                            // look to the right
                            between = transcripts.get(skippedExons.get(exon)).completeRegionsInRange(sr1.getEnd(), sr2.getStart());
                        } else {
                            // look to the left
                            between = transcripts.get(skippedExons.get(exon)).completeRegionsInRange(sr2.getEnd(), sr1.getStart());
                        }
                    } else {
                        if (sr1Smaller) {
                            // look to the left
                            between = transcripts.get(skippedExons.get(exon)).completeRegionsInRange(sr1.getEnd(), sr2.getStart());
                        } else {
                            // look to the right
                            between = transcripts.get(skippedExons.get(exon)).completeRegionsInRange(sr2.getEnd(), sr1.getStart());
                        }
                    }
                    // if between contains a skipped exon => count it as inclusive
                    if (between != null) {
                        between.retainAll(skippedExons.keySet());
                        if (between.size() > 0) {
                            for (String e : between) {
                                if (!countedInclusive.contains(e)) {
                                    skippedExonCounts.get(e)[0]++;
                                    countedInclusive.add(e);
                                }
                            }
                        }
                    }

                } else if ((skippedExons.get(exon).equals(tr1) && evidenceForTr2) || (skippedExons.get(exon).equals(tr2) && evidenceForTr1)) {
                    // maybe INCLUSIVE?
                    // check: does the read deliver evidence by extending neighboring regions of transcript that does not contain the exon?
                    skippedExonCounts.get(exon)[0]++;
                    countedInclusive.add(exon);

                } else {
                    // Read might be EXCLUSIVE for the corresponding exon?
                    if (sr1.skipsRegion(a, b) || sr2.skipsRegion(a, b)) {
                        if (!countedExclusive.contains(exon)) {
                            countedExclusive.add(exon);
                            skippedExonCounts.get(exon)[1]++;
                        }
                        // if exon is skipped, also check other exons on same transcript
                        for (String e : skippedExons.keySet()) {
                            // if there is a different skipped exon on the same transcript
                            if (!e.equals(exon) && skippedExons.get(e).equals(skippedExons.get(exon)) && !countedExclusive.contains(e)) {
                                int c = Integer.parseInt(e.split("-")[0]);
                                int d = Integer.parseInt(e.split("-")[1]) - 1;
                                if (transcripts.get(skippedExons.get(e)).coversNeighborsOfRegion(c, d, sr1.merge(sr2))) {
                                    countedExclusive.add(e);
                                    skippedExonCounts.get(e)[1]++;
                                }
                            }
                        }
                    } else {
                        // check if there is other evidence for a skip: read uses intronic region of one transcript
                        if (!countedExclusive.contains(exon)) {
                            if (evidenceAgainstTr1 && skippedExons.get(exon).equals(tr1)) {
                                countedExclusive.add(exon);
                                skippedExonCounts.get(exon)[1]++;
                            } else if (evidenceAgainstTr2 && skippedExons.get(exon).equals(tr2)) {
                                countedExclusive.add(exon);
                                skippedExonCounts.get(exon)[1]++;
                            } else {
                                if (exonSidesCovered.get(exon)) {
                                    // check case where an exon is skipped, that is no exon skipping event
                                    HashSet<String> sr1Tr1Skipped = transcripts.get(tr1).whichSkipped(sr1);
                                    sr1Tr1Skipped.removeAll(skippedExons.keySet());
                                    HashSet<String> sr1Tr2Skipped = transcripts.get(tr2).whichSkipped(sr1);
                                    sr1Tr2Skipped.removeAll(skippedExons.keySet());
                                    HashSet<String> sr2Tr1Skipped = transcripts.get(tr1).whichSkipped(sr2);
                                    sr2Tr1Skipped.removeAll(skippedExons.keySet());
                                    HashSet<String> sr2Tr2Skipped = transcripts.get(tr2).whichSkipped(sr2);
                                    sr2Tr2Skipped.removeAll(skippedExons.keySet());
                                    if (sr1Tr1Skipped.size() > 0 || sr1Tr2Skipped.size() > 0 || sr2Tr1Skipped.size() > 0 || sr2Tr2Skipped.size() > 0) {
                                        if (!countedExclusive.contains(exon)) {
                                            countedExclusive.add(exon);
                                            skippedExonCounts.get(exon)[1]++;
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }
    }

    public ArrayList<String> outputLines() {
        // return output lines for all exons
        ArrayList<String> out = new ArrayList<>();

        // get the output line for every skipped exon
        for (String exon : skippedExons.keySet()) {
            int inclReads = skippedExonCounts.get(exon)[0];
            int exclReads = skippedExonCounts.get(exon)[1];
            int totalReads = inclReads + exclReads;

            double psi = inclReads / (double) totalReads;
            // only add exons with total count > 0
            if (totalReads > 0)
                out.add(this.geneId + "\t" + exon + "\t" + inclReads + "\t" + exclReads + "\t" + totalReads + "\t" + doubleToString(psi));
        }

        return out;
    }

    public void addExon(int start, int end, String trId) {
        this.transcripts.get(trId).addRegion(new int[]{start, end});
    }

    public void addCds(int start, int end, String trId) {
        this.cds.get(trId).addRegion(new int[]{start, end});
    }

    public void addTranscript(String trId) {
        this.transcriptIds.add(trId);
        this.transcripts.put(trId, new RegionVector(new ArrayList<>()));
        transcriptToCount.put(trId, 0);
        this.cds.put(trId, new RegionVector(new ArrayList<>()));
    }

    // compute the exon skipping events using the CDS!
    public void computeExonSkippings() {
        // the gtf only contains genes with two transcripts each!
        skippedExonToNeighbors = new HashMap<>();

        // little pre-check: method is only called when no trancripts are added anymore => check if there are really only 2
        if(this.transcripts.keySet().size() != 2){
            System.out.println("The provided gtf-file contains a gene with not exactly two transcripts. This case cannot be handled.");
            System.exit(0);
        }

        String tr1 = null;
        String tr2 = null;
        for (String tr : cds.keySet()) {
            if (tr1 == null) tr1 = tr;
            else if (tr2 == null) tr2 = tr;
        }

        HashMap<String, String> skippedExons = new HashMap<>();
        this.skippedExonCounts = new HashMap<>();

        for (int[] exon : cds.get(tr1).skippedInOther(cds.get(tr2))) {
            String e = exon[0] + "-" + (exon[1] + 1);
            skippedExons.put(e, tr1);
            skippedExonCounts.put(e, new int[]{0, 0, 0});
            skippedExonToNeighbors.put(e, transcripts.get(tr1).getNeighbors(exon[0], exon[1]));
        }
        for (int[] exon : cds.get(tr2).skippedInOther(cds.get(tr1))) {
            String e = exon[0] + "-" + (exon[1] + 1);
            skippedExons.put(e, tr2);
            skippedExonCounts.put(e, new int[]{0, 0, 0});
            skippedExonToNeighbors.put(e, transcripts.get(tr2).getNeighbors(exon[0], exon[1]));
        }
        this.skippedExons = skippedExons;
    }

    public String matchPlausibilityLevel(RegionVector sr1, RegionVector sr2) {
        // returns match plausibility level as String: either "1" or "2" or "3transcriptId1,transcriptId2,..."
        HashSet<String> Subregion = new HashSet<>();

        // check all transcripts: are the region vectors sub-RVs?
        for (String id : transcripts.keySet()) {
            boolean a = transcripts.get(id).hasAsSubRV(sr1);
            boolean b = transcripts.get(id).hasAsSubRV(sr2);
            if (a && b) Subregion.add(id);
        }

        // if Subregion contains at least one transcript => level 3 (transcriptomic)!
        if (Subregion.size() > 0) {
            return "3";
        }

        // check if its merged-transcriptomic
        // merge all transcript RVs to one RV
        RegionVector mergedTranscripts = transcripts.get(transcriptIds.get(0));
        for (int i = 1; i < transcriptIds.size(); i++) {
            mergedTranscripts = mergedTranscripts.merge(transcripts.get(transcriptIds.get(i)));
        }

        if (mergedTranscripts.Contains(sr1) && mergedTranscripts.Contains(sr2)) {
            // level 2 (merged-transcriptomic)
            return "2";
        }
        // level 1 (intronic)
        return "1";
    }

    // rounds to 3 decimal points and returns string with exactly 3 decimal points
    public String doubleToString(double d) {
        d = d * 1000;
        d = Math.round(d);
        d = d / 1000.0;
        String result = "" + d;
        for (int x = result.split("\\.")[1].length(); x < 3; x++) {
            result += "0";
        }
        return result;
    }
}
