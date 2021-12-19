package com.company;

import net.sf.samtools.*;

import java.io.*;
import java.util.*;

import augmentedTree.*;

public class Runner {
    public static void main(String[] args) {
        Long time;
        ArrayList<Long> times = new ArrayList<>();
        ArrayList<Long> writingTimes = new ArrayList<>();

        Parser parser = new Parser("-gtf", "-bam", "-o", "-frstrand");
        parser.parseArguments(args);


        String gtf = parser.getValue("-gtf");
        String bam = parser.getValue("-bam");
        String out = parser.getValue("-o");
        String fr = parser.getValue("-frstrand");

        // check if no parameter is set
        if (gtf == null || bam == null || out == null) {
            System.out.println("Usage:");
            System.out.println("-gtf <gtf_file>");
            System.out.println("-bam <bam_file>");
            System.out.println("-o <output_tsv>");
            System.out.println("[ -frstrand <true/false> ]");
            System.exit(0);
        }

        // check frstrand
        boolean specific = (fr != null);
        boolean frstrand = true;
        if (specific) {
            if (fr.equals("false")) frstrand = false;
            else if (!fr.equals("true")) {
                System.out.println("-frstrand must be <true> or <false>!");
                System.exit(0);
            }
        }

        // open file writer for output file
        BufferedWriter outputWriter = null;
        try {
            outputWriter = new BufferedWriter(new FileWriter(out));
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(0);
        }

        SAMFileReader samReader = new SAMFileReader(new File(bam));
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        Iterator<SAMRecord> it = samReader.iterator();
        SAMRecord sr;

        // ### read gtf file ###########################################################################################
        HashMap<String, IntervalTree<Gene>> chrToIntTree = new HashMap<>();
        HashMap<String, IntervalTree<Gene>> chrToIntTreePos = new HashMap<>();
        HashMap<String, IntervalTree<Gene>> chrToIntTreeNeg = new HashMap<>();
        HashSet<String> chromosomesSeen = new HashSet<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(gtf));
            String line;
            Gene g = null;
            String transcript = null;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("#")) {
                    String[] columns = line.split("\t");
                    String features = columns[8].strip();

                    // if new chromosome is seen: add to HashMaps (with interval trees)
                    if (!chromosomesSeen.contains(columns[0])) {
                        chromosomesSeen.add(columns[0]);
                        if (!specific) {
                            chrToIntTree.put(columns[0], new IntervalTree<>());
                        } else {
                            chrToIntTreePos.put(columns[0], new IntervalTree<>());
                            chrToIntTreeNeg.put(columns[0], new IntervalTree<>());
                        }
                    }

                    if (columns[2].equals("gene")) {
                        // new gene: add old one to Interval Tree
                        if (g != null) {
                            if (!specific) chrToIntTree.get(g.getChromosome()).add(g);
                            else {
                                if (g.getStrand()) chrToIntTreePos.get(g.getChromosome()).add(g);
                                else chrToIntTreeNeg.get(g.getChromosome()).add(g);
                            }
                        }
                        String gId = features.substring((features.indexOf("gene_id") + 9));
                        gId = gId.substring(0, gId.indexOf("\""));
                        boolean strand = columns[6].equals("+");
                        g = new Gene(Integer.parseInt(columns[3]), Integer.parseInt(columns[4]), gId, columns[1], columns[0], strand);
                    }

                    if (columns[2].equals("transcript")) {
                        // new transcript for current gene
                        transcript = features.substring(features.indexOf("transcript_id") + 15);
                        transcript = transcript.substring(0, transcript.indexOf("\""));
                        g.addTranscript(transcript);
                    }

                    if (columns[2].equals("exon")) {
                        // new exon for current transcript
                        g.addExon(Integer.parseInt(columns[3]), Integer.parseInt(columns[4]), transcript);
                    }
                }
            }
            // add last gene
            if (g != null) {
                if (!specific) chrToIntTree.get(g.getChromosome()).add(g);
                else {
                    if (g.getStrand()) chrToIntTreePos.get(g.getChromosome()).add(g);
                    else chrToIntTreeNeg.get(g.getChromosome()).add(g);
                }
            }

        } catch (FileNotFoundException e) {
            System.out.println("gtf-file \"" + gtf + "\" was not found!");
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // ### go through BAM-files and process reads ##################################################################
        HashMap<String, SAMRecord> lookup = new HashMap<>();
        // to keep track of current chromosome
        String chromosome = null;
        PriorityQueue<RegionVector> pcrIndexLookup = new PriorityQueue<>();

        while ((sr = it.next()) != null) {
            // clear lookup (for read mates) and pcrIndexLookup if we are on new chromosome
            if (!sr.getReferenceName().equals(chromosome)) {
                lookup.clear();
                pcrIndexLookup.clear();
                chromosome = sr.getReferenceName();
            }
            // remove regions that cannot be seen anymore from pcrIndexLookup
            int pos = sr.getAlignmentStart();
            RegionVector first = pcrIndexLookup.peek();
            while (first != null && first.getEnd() < pos) {
                pcrIndexLookup.poll();
                first = pcrIndexLookup.peek();
            }

            // check if read can be ignored
            IntervalTree<Gene> readSpecificTree;
            boolean readSenseToTranscript = (frstrand && sr.getFirstOfPairFlag()) || (!frstrand && !sr.getFirstOfPairFlag());
            boolean transcriptStrand = !sr.getReadNegativeStrandFlag();
            // get the strand where the transcript was on (depends on if the read is sense to transcript and the strand where the read matches)
            if (!readSenseToTranscript) transcriptStrand = sr.getReadNegativeStrandFlag();

            if (!specific) readSpecificTree = chrToIntTree.get(sr.getReferenceName());
            else {
                // get the right trees
                if (!transcriptStrand) readSpecificTree = chrToIntTreeNeg.get(sr.getReferenceName());
                else readSpecificTree = chrToIntTreePos.get(sr.getReferenceName());
            }

            if (ignorePreCheck(sr)) continue;

            String id = sr.getReadName();
            SAMRecord mate = lookup.get(id);
            String outLine = "";
            // if mate already is in lookup => HANDLE READ PAIR
            if (mate != null && !ignoreReadPair(sr, mate, readSpecificTree)) {
                time = System.nanoTime();
                // check split inconsistency
                boolean splitInconsistent = checkSplitInconsistent(sr, mate);

                if (splitInconsistent) {
                    outLine = sr.getReadName() + "\tsplit-inconsistent:true\n";
                } else {
                    // compute/look up pcrIndex:
                    RegionVector merged = samRecordToRegionVector(sr).merge(samRecordToRegionVector(mate));
                    // set transcript strand
                    if (!specific) merged.setTranscriptStrand(true);
                    else merged.setTranscriptStrand(transcriptStrand);

                    Iterator<RegionVector> iterator = pcrIndexLookup.iterator();
                    RegionVector element = null;
                    if (pcrIndexLookup.size() > 0) element = iterator.next();
                    Integer pcrIndex = null;
                    while (element != null && pcrIndex == null) {
                        if (element.equals(merged) && element.sameStrand(merged)) {
                            element.incNumber();
                            pcrIndex = element.getNumber();
                        }
                        try {
                            element = iterator.next();
                        } catch (NoSuchElementException e) {
                            element = null;
                        }
                    }
                    if (pcrIndex == null) {
                        pcrIndexLookup.add(merged);
                        pcrIndex = 1;
                    }

                    // some strandness handling...
                    SAMRecord firstRead;
                    SAMRecord secondRead;
                    if (sr.getFirstOfPairFlag()) {
                        firstRead = sr;
                        secondRead = mate;
                    } else {
                        firstRead = mate;
                        secondRead = sr;
                    }

                    IntervalTree<Gene> tree;
                    IntervalTree<Gene> other;
                    if (!specific) {
                        tree = chrToIntTree.get(sr.getReferenceName());
                        other = chrToIntTree.get(sr.getReferenceName());
                    } else {
                        // get the right trees (we already know on which strand the transcript was)
                        if (transcriptStrand) {
                            tree = chrToIntTreePos.get(sr.getReferenceName());
                            other = chrToIntTreeNeg.get(sr.getReferenceName());
                        } else {
                            tree = chrToIntTreeNeg.get(sr.getReferenceName());
                            other = chrToIntTreePos.get(sr.getReferenceName());
                        }
                    }

                    // actually handle read pair
                    outLine = handleReadPair(firstRead, secondRead, tree, other, pcrIndex);
                    times.add(System.nanoTime() - time);
                }
                try {
                    time = System.nanoTime();
                    outputWriter.write(outLine);
                    writingTimes.add(System.nanoTime() - time);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                // remove from lookup afterwards
                // TODO: maybe remove bc its inefficient?
                lookup.remove(id);
            } else {
                lookup.put(id, sr);
            }
        }
        try {
            outputWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // print time info
        long sum = 0;
        for (long t : times) sum += t;
        float mean = sum / (float) times.size();
        System.out.println("Mean time per read pair in ns: " + mean);
        System.out.println("That equals " + (1000000000 / mean) + " read pairs in 1 second");

        sum = 0;
        for (long t : writingTimes) sum += t;
        mean = sum / (float) writingTimes.size();
        System.out.println("Mean time per line that is wrote in ns: " + mean);
        System.out.println("That equals " + (1000000000 / mean) + " lines in 1 second");
    }

    private static String handleReadPair(SAMRecord sr1, SAMRecord sr2, IntervalTree tree, IntervalTree treeOtherDir, int pcrIndex) {
        StringBuilder outLine = new StringBuilder(sr1.getReadName() + "\t");

        // mismatches
        Integer mm1 = (Integer) sr1.getAttribute("NM");
        Integer mm2 = (Integer) sr2.getAttribute("NM");
        outLine.append("mm:").append(mm1 + mm2).append("\t");

        // clipping: calculate the differences between alignment start/end and unclipped start/end
        int clip = (sr1.getAlignmentStart() - sr1.getUnclippedStart()) + (sr1.getUnclippedEnd() - sr1.getAlignmentEnd()) +
                (sr2.getAlignmentStart() - sr2.getUnclippedStart()) + (sr2.getUnclippedEnd() - sr2.getAlignmentEnd());
        outLine.append("clipping:").append(clip).append("\t");
        // nsplit
        outLine.append("nsplit:").append(getNumSplit(sr1, sr2)).append("\t");

        // matching genes
        HashSet<Gene> matchingGenes = genesMatchingReadPair(sr1, sr2, tree);

        // check the match plausibility levels
        ArrayList<Gene> highestLevelGenes = new ArrayList<>();
        ArrayList<String> transcriptIds = new ArrayList<>();
        int highestLevel = 0;
        RegionVector rv1 = samRecordToRegionVector(sr1);
        RegionVector rv2 = samRecordToRegionVector(sr2);

        for (Gene g : matchingGenes) {
            String level = g.matchPlausibilityLevel(rv1, rv2);
            int lev = Integer.parseInt(level.substring(0, 1));
            if (lev > highestLevel) {
                highestLevel = lev;
                highestLevelGenes.clear();
            }
            if (lev == highestLevel) {
                highestLevelGenes.add(g);
                if (highestLevel == 3) {
                    // already remember the transcript IDs, needed later
                    transcriptIds.add(level.substring(1));
                }
            }
        }
        int gcount = highestLevelGenes.size();
        outLine.append("gcount:").append(gcount).append("\t");

        // in case level > 0: print rest information
        for (int i = 0; i < highestLevelGenes.size(); i++) {
            Gene g = highestLevelGenes.get(i);
            outLine.append(g.getId()).append(",").append(g.getBiotype()).append(":");
            if (highestLevel == 1) outLine.append("INTRON");
            else if (highestLevel == 2) outLine.append("MERGED");
            else if (highestLevel == 3) outLine.append(transcriptIds.get(i));
            if (i < highestLevelGenes.size() - 1) outLine.append("|");
        }
        if (highestLevelGenes.size() > 0) outLine.append("\t");

        // in case gcount == 0 (= highest level == 0): print gdist and antisense
        if (gcount == 0) {
            HashSet<Gene> leftNeighbors = new HashSet<>();
            HashSet<Gene> rightNeighbors = new HashSet<>();
            tree.getIntervalsLeftNeighbor(Math.min(sr1.getAlignmentStart(), sr2.getAlignmentStart()), Math.max(sr1.getAlignmentEnd(), sr2.getAlignmentEnd()), leftNeighbors);
            tree.getIntervalsRightNeighbor(Math.min(sr1.getAlignmentStart(), sr2.getAlignmentStart()), Math.max(sr1.getAlignmentEnd(), sr2.getAlignmentEnd()), rightNeighbors);

            // okay: 1 or more neighbors (genes can overlap each other) (overlap genes => gdist = 0) No neighbors also okay, but should only occur on one side!
            Integer shortestDist = null;
            int leftEnd = Math.min(sr1.getAlignmentStart(), sr2.getAlignmentStart());
            int rightEnd = Math.max(sr1.getAlignmentEnd(), sr2.getAlignmentEnd());
            for (Gene left : leftNeighbors) {
                int dist = leftEnd - left.getStop();
                if (shortestDist == null) shortestDist = dist;
                if (dist < shortestDist) {
                    shortestDist = dist;
                }
            }
            for (Gene right : rightNeighbors) {
                int dist = right.getStart() - rightEnd;
                if (shortestDist == null) shortestDist = dist;
                if (dist < shortestDist) {
                    shortestDist = dist;
                }
            }
            if (shortestDist != null) shortestDist--;
            if (shortestDist != null && shortestDist < 0) shortestDist = 0;
            outLine.append("gdist:").append(shortestDist).append("\t");

            // antisense
            HashSet<Gene> matchingGenesAntisense = genesMatchingReadPair(sr1, sr2, treeOtherDir);
            if (matchingGenesAntisense.size() > 0) outLine.append("antisense:true\t");
            else outLine.append("antisense:false\t");
        }

        // pcrIndex
        outLine.append("pcrindex: ").append(pcrIndex - 1);

        outLine.append("\n");
        return outLine.toString();
    }

    private static boolean ignorePreCheck(SAMRecord sr) {
        // Ignore reads if they are secondary or supplementary, if they are unmapped, if their mate is
        //unmapped, or if mate pairs are mapped to the same strand of the same chromosome.
        boolean ignore = sr.getNotPrimaryAlignmentFlag();
        ignore = ignore || sr.getReadUnmappedFlag();
        ignore = ignore || sr.getMateUnmappedFlag();
        return ignore || (sr.getReadNegativeStrandFlag() == sr.getMateNegativeStrandFlag() && sr.getReferenceName().equals(sr.getMateReferenceName()));
    }

    private static boolean ignoreReadPair(SAMRecord sr1, SAMRecord sr2, IntervalTree<Gene> genes) {
        // check if pair is intergenic and contains at least 1 whole gene
        boolean ignoreIntergenic = false;
        HashSet<Gene> containingSr1 = new HashSet<>();
        HashSet<Gene> containingSr2 = new HashSet<>();
        HashSet<Gene> genesSpannedByPair = new HashSet<>();
        genes.getIntervalsSpanning(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), containingSr1);
        genes.getIntervalsSpanning(sr2.getAlignmentStart(), sr2.getAlignmentEnd(), containingSr2);
        containingSr1.retainAll(containingSr2);

        int start = Math.min(sr1.getAlignmentStart(), sr2.getAlignmentStart());
        int end = Math.max(sr1.getAlignmentEnd(), sr2.getAlignmentEnd());
        genes.getIntervalsSpannedBy(start, end, genesSpannedByPair);

        if (genesSpannedByPair.size() > 0 && containingSr1.size() == 0) ignoreIntergenic = true;

        return ignoreIntergenic;
    }

    private static int getNumSplit(SAMRecord sr1, SAMRecord sr2) {
        // should only be called when read pair is not inconsitent!
        RegionVector a = samRecordToRegionVector(sr1).getReverse();
        RegionVector b = samRecordToRegionVector(sr2).getReverse();
        if (a.getNumRegions() == 0) return b.getNumRegions();
        if (b.getNumRegions() == 0) return a.getNumRegions();
        return a.numUniqueRegions(b);
    }

    private static boolean checkSplitInconsistent(SAMRecord sr1, SAMRecord sr2) {
        // create RegionVector Objects
        RegionVector rv1 = samRecordToRegionVector(sr1);
        RegionVector rv2 = samRecordToRegionVector(sr2);

        RegionVector overlap = rv2.overlap(rv1.getReverse());
        RegionVector overlap2 = rv1.overlap(rv2.getReverse());

        return (overlap.getNumRegions() > 0) || (overlap2.getNumRegions() > 0);
    }

    // return the region vector implied by a SAMRecord object
    private static RegionVector samRecordToRegionVector(SAMRecord sr) {
        RegionVector rv = new RegionVector(new ArrayList<>());
        for (net.sf.samtools.AlignmentBlock ab : sr.getAlignmentBlocks()) {
            rv.addRegion(new int[]{ab.getReferenceStart(), ab.getReferenceStart() + ab.getLength() - 1});
        }
        return rv;
    }

    // returns HashSet with all geneIDs of genes where read pair lies within
    private static HashSet<Gene> genesMatchingReadPair(SAMRecord sr1, SAMRecord sr2, IntervalTree tree) {
        HashSet<Gene> genesWithSr1Start = new HashSet<>();
        HashSet<Gene> genesWithSr1End = new HashSet<>();
        HashSet<Gene> genesWithSr2Start = new HashSet<>();
        HashSet<Gene> genesWithSr2End = new HashSet<>();
        tree.getIntervalsSpanning(sr1.getAlignmentStart(), genesWithSr1Start);
        tree.getIntervalsSpanning(sr1.getAlignmentEnd(), genesWithSr1End);
        tree.getIntervalsSpanning(sr2.getAlignmentStart(), genesWithSr2Start);
        tree.getIntervalsSpanning(sr2.getAlignmentEnd(), genesWithSr2End);

        genesWithSr1Start.retainAll(genesWithSr1End);
        genesWithSr2Start.retainAll(genesWithSr2End);

        genesWithSr1Start.retainAll(genesWithSr2Start);

        return genesWithSr1Start;
    }
}