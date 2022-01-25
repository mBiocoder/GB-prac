package com.company;

import net.sf.samtools.*;

import java.io.*;
import java.util.*;

import augmentedTree.*;

public class Runner {
    public static void main(String[] args) {
        Long totalTime = System.currentTimeMillis();
        Long time;
        ArrayList<Long> times = new ArrayList<>();

        Parser parser = new Parser("-gtf", "-bam", "-o");
        parser.parseArguments(args);

        String gtf = parser.getValue("-gtf");
        String bam = parser.getValue("-bam");
        String out = parser.getValue("-o");

        // check if no parameter is set
        if (gtf == null || bam == null || out == null) {
            System.out.println("Usage:");
            System.out.println("-gtf <gtf_file>");
            System.out.println("-bam <bam_file>");
            System.out.println("-o <output_tsv>");
            System.exit(0);
        }

        // ### read gtf file ###########################################################################################
        HashMap<String, IntervalTree<Gene>> chrToIntTree = new HashMap<>();
        HashSet<String> chromosomesSeen = new HashSet<>();

        time = System.currentTimeMillis();

        try {
            BufferedReader br = new BufferedReader(new FileReader(gtf));
            String line;
            Gene g = null;
            String transcript = null;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("#")) {
                    String[] columns = line.split("\t");
                    String features = columns[8].replace("\n", "");

                    // if new chromosome is seen: add to HashMaps (with interval trees)
                    if (!chromosomesSeen.contains(columns[0])) {
                        chromosomesSeen.add(columns[0]);
                        chrToIntTree.put(columns[0], new IntervalTree<>());

                    }

                    if (columns[2].equals("gene")) {
                        // new gene: add old one to Interval Tree
                        if (g != null) {
                            g.computeExonSkippings();
                            chrToIntTree.get(g.getChromosome()).add(g);
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

                    if (columns[2].equals("CDS")) {
                        // new exon for current transcript
                        g.addCds(Integer.parseInt(columns[3]), Integer.parseInt(columns[4]), transcript);
                    }
                }
            }
            // add last gene
            if (g != null) {
                g.computeExonSkippings();
                chrToIntTree.get(g.getChromosome()).add(g);
            }

        } catch (FileNotFoundException e) {
            System.out.println("gtf-file \"" + gtf + "\" was not found!");
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("GTF reading time in ms: " + (System.currentTimeMillis() - time));

        // ### go through BAM-files and count reads for every exon #####################################################
        // open file writer for output file
        BufferedWriter outputWriter = null;
        try {
            outputWriter = new BufferedWriter(new FileWriter(out));
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(0);
        }

        // open SAMFileReader
        SAMFileReader samReader = new SAMFileReader(new File(bam));
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        Iterator<SAMRecord> it = samReader.iterator();
        SAMRecord sr;

        HashMap<String, SAMRecord> lookup = new HashMap<>();
        // to keep track of current chromosome
        String chromosome = null;

        while ((sr = it.next()) != null) {
            // clear lookup (for read mates) and pcrIndexLookup if we are on new chromosome
            if (!sr.getReferenceName().equals(chromosome)) {
                lookup.clear();
                chromosome = sr.getReferenceName();
            }

            // check if read can be ignored
            if (ignorePreCheck(sr)) continue;

            IntervalTree<Gene> readSpecificTree = chrToIntTree.get(sr.getReferenceName());

            String id = sr.getReadName();
            SAMRecord mate = lookup.get(id);
            // if mate already is in lookup => Count Read Pair
            if (mate != null && !ignoreReadPair(sr, mate, readSpecificTree)) {
                time = System.nanoTime();
                // find genes where read pair fits, then count it
                HashSet<Gene> matchingGenes = genesMatchingReadPair(sr, mate, chrToIntTree.get(sr.getReferenceName()));
                RegionVector rv1 = samRecordToRegionVector(sr);
                RegionVector rv2 = samRecordToRegionVector(mate);
                for (Gene g : matchingGenes) {
                    g.countRead(rv1, rv2);
                }

                lookup.remove(id);
                times.add(System.nanoTime() - time);
            } else {
                lookup.put(id, sr);
            }
        }

        //############ go through genome data structure, compute psi values and write them to output file ##############
        time = System.currentTimeMillis();
        try {
            outputWriter.write("gene\texon\tnum_incl_reads\tnum_excl_reads\tnum_total_reads\tpsi\n");
            for (String chr : chrToIntTree.keySet()) {
                IntervalTree<Gene> tree = chrToIntTree.get(chr);
                for (Gene g : tree) {
                    for (String line : g.outputLines()) {
                        outputWriter.write(line + "\n");
                    }
                }
            }
            outputWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


        Long writingTime = System.currentTimeMillis() - time;
        System.out.println("Writig time: " + writingTime);
        // print time info
        System.out.println("total time in ms: " + (System.currentTimeMillis() - totalTime));

        long sum = 0;
        for (long t : times) sum += t;
        float mean = sum / (float) times.size();
        System.out.println("Mean time per read pair in ns: " + mean);
        System.out.println("That equals " + (1000000000 / mean) + " read pairs in 1 second");
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