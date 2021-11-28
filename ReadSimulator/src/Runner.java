import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class Runner {
    public static void main(String[] args) {
        int readLength = 0;
        int frLengthMean = 0;
        int SD = 0;
        File readcounts = null;
        double mutationrate = 0;
        File fasta = null;
        File fidx = null;
        File gtf = null;
        File od = null;
        if (args.length != 18) {
            String consoleParameters = "";
            for (int i = 0; i < args.length; i++)
                consoleParameters = consoleParameters.concat(" args[" + i + "]=" + args[i]);
            System.err.println("Required 18 parameters of following parameters: -length <int>, -frlength <int>, -SD <int>, -readcounts <table>, -mutationrate <double>, -fasta <file>, -fidx <file>, -gtf <file>, -od <file> ");
            System.exit(0);
        } else {
            for (int i = 0; i < args.length; i = i + 2) {
                switch (args[i]) {
                    case "-length":
                        readLength = Integer.parseInt(args[i + 1]);
                        continue;
                    case "-frlength":
                        frLengthMean = Integer.parseInt(args[i + 1]);
                        continue;
                    case "-SD":
                        SD = Integer.parseInt(args[i + 1]);
                        continue;
                    case "-readcounts":
                        readcounts = new File(args[i + 1]);
                        continue;
                    case "-mutationrate":
                        mutationrate = Double.parseDouble(args[i + 1]);
                        continue;
                    case "-fasta":
                        fasta = new File(args[i + 1]);
                        continue;
                    case "-fidx":
                        fidx = new File(args[i + 1]);
                        continue;
                    case "-gtf":
                        gtf = new File(args[i + 1]);
                        continue;
                    case "-od":
                        od = new File(args[i + 1]);
                }
            }
            System.out.println("|UserInput|\nreadLength:\t" + readLength + "\nfrLengthMean:\t" + frLengthMean + "\nSD:\t" + SD + "\nreadcounts:\t" + readcounts.toString() + "\nmutationrate:\t" + mutationrate + "\nfasta:\t" + fasta.toString() + "\nfidx:\t" + fidx.toString() + "\ngtf:\t" + gtf.toString() + "\nod:\t" + od.toString());
        }

        // Object[0] = String gene_id, Object[1] = String transcript_id, Object[2] = Integer count
        ArrayList<Object[]> allReadcounts = parseReadcounts(readcounts);

        // ---Parse GTF---
        GSE gse = new GSE(fasta, fidx, gtf);

        //Try for testsert
        //testCases(gse);

        // ---For every ReadcountLine, there's one outputLine containing readsFw, readsRw, readMappingInfo--
        ArrayList<OutputLines> outLines = new ArrayList<>();
        for (Object[] obj : allReadcounts) {
            System.out.println("Begin: getOutputLines for readcountLine GeneId:" + obj[0] + " transcriptId: " + obj[1] + " counts: " + obj[2]);
            String gene_id = (String) obj[0];
            String transcript_id = (String) obj[1];
            int counts = (Integer) obj[2];
            outLines.add(getOutputLines(gene_id, transcript_id, counts, gse, frLengthMean, SD, readLength, mutationrate));
        }

        // ---Write OutputLines---
        Boolean append = false;
        int counter = 0;
        System.out.println("Writing Output Lines...");
        for (OutputLines line : outLines) {
            // append to files
            counter = writeOutputLines(line, od, append, counter);
            append = true;
        }
        System.out.println("End of main method.");
    }

    private static OutputLines getOutputLines(String gene_id, String transcript_id, int counts, GSE gse, int frLengthMean,
                                              int SD, int readLength, double mutationrate) {
        ArrayList<String> readsFw = new ArrayList<>();
        ArrayList<String> readsRw = new ArrayList<>();
        ArrayList<String> readMappingInfo = new ArrayList<>();

        // get transcript sequence
        TreeSet<Exon> exonsInTranscript = gse.getExonsByTranscriptAndGene(transcript_id, gene_id);
        String transcriptSequence = "";
        Exon referenceExon = null;
        for (Exon e : exonsInTranscript) {
            String s = gse.getSequence(e.chr, e.start, e.end, e.strand);
            transcriptSequence = transcriptSequence.concat(s);
            referenceExon = new Exon(e);
        }
        for (int j = 0; j < counts; j++) {
            // Get fragment length from max(readlength, ND(mean, SD))
            int fragmentLength = GSE.getFragmentLength(frLengthMean, SD, readLength, transcriptSequence.length());
            // random fragmentPosition from 0 to length(t) - FL
            int fragmentStart = new Random().nextInt(transcriptSequence.length() - fragmentLength);
            // get fragment sequence
            String fragmentSeq = transcriptSequence.substring(fragmentStart, fragmentStart + fragmentLength);
            // get read sequences of length readlength (second read is reverse complement)
            String readFw = GSE.getReadSequence(fragmentSeq, readLength, false);
            String readRw = GSE.getReadSequence(fragmentSeq, readLength, true);
            // simulate mutations with required rate
            String[] mutatedReadFw = GSE.simulateMutations(readFw, mutationrate);
            String[] mutatedReadRw = GSE.simulateMutations(readRw, mutationrate);

            readsFw.add(mutatedReadFw[0]);
            readsRw.add(mutatedReadRw[0]);

            // Without mutation
//            String[] tempReadFw = {readFw, ""};
//            String[] tempReadRw = {readRw, ""};
//            readsFw.add(tempReadFw[0]);
//            readsRw.add(tempReadRw[0]);

            // readid chr gene transcript t_fw_regvec t_rw_regvec fw_regvec rw_regvec fw_mut rw_mut 0 19 ENSG00000104870 ENST00000221466 810-885 854-929 50017390-50017391|50017468-50017542 50017511-50017586 55
            int tFwRegVecStart = fragmentStart;
            int tFwRegVecEnd = fragmentStart + readLength;
            int tRwRegVecStart = fragmentStart + fragmentLength - readLength;
            int tRwRegVecEnd = fragmentStart + fragmentLength;

            int fwRegVecStart = gse.getGenomicPosition(transcript_id, gene_id, tFwRegVecStart, true, referenceExon.strand);
            int fwRegVecEnd = gse.getGenomicPosition(transcript_id, gene_id, tFwRegVecEnd, false, referenceExon.strand);
            int rwRegVecStart = gse.getGenomicPosition(transcript_id, gene_id, tRwRegVecStart, true, referenceExon.strand);
            int rwRegVecEnd = gse.getGenomicPosition(transcript_id, gene_id, tRwRegVecEnd, false, referenceExon.strand);

            String fwRegVec = gse.getGenomicRegionPosition(transcript_id, gene_id, fwRegVecStart, fwRegVecEnd);
            String rwRegVec = gse.getGenomicRegionPosition(transcript_id, gene_id, rwRegVecStart, rwRegVecEnd);

            String mappingInfo = referenceExon.chr + "\t" + referenceExon.gene_id + "\t" + referenceExon.transcript_id
                    + "\t" + tFwRegVecStart + "-" + tFwRegVecEnd + "\t" + tRwRegVecStart + "-" + tRwRegVecEnd + "\t"
                    + fwRegVec + "\t" + rwRegVec + "\t" + mutatedReadFw[1] + "\t" + mutatedReadRw[1];
            readMappingInfo.add(mappingInfo);
            // System.out.println("readFw: " + readFw + "(" + readFw.length() + ") mutated to " + mutatedReadFw[0] + "(" + mutatedReadFw[0].length() + ")");
            // System.out.println("readRw: " + readRw + "(" + readRw.length() + ") mutated to " + mutatedReadRw[0] + "(" + mutatedReadRw[0].length() + ")");
            // System.out.println("transcriptSeq (" + transcriptSequence.length()/* + ") seq: " + transcriptSequence */);
        }
        return new OutputLines(readsFw, readsRw, readMappingInfo);
    }

    public static int writeOutputLines(OutputLines line, File od, Boolean append, int counter) {
        // System.out.println("Begin: writeOutputLines. Append: " + append);
        // System.out.println("\tWriting fw.fastq...");
        try {
            int fwCounter = counter;
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(od + "/fw.fastq", append)));
            for (String string : line.readsFw) {
                out.println("@" + fwCounter);
                out.println(string);
                out.println("+" + fwCounter);
                // write string.length() times 'I'
                char[] chars = new char[string.length()];
                Arrays.fill(chars, 'I');
                out.println(new String(chars));
                fwCounter++;
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        // System.out.println("\tWriting rw.fastq...");
        try {
            int rwCounter = counter;
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(od + "/rw.fastq", append)));
            for (String string : line.readsRw) {
                out.println("@" + rwCounter);
                out.println(string);
                out.println("+" + rwCounter);
                // write string.length() times 'I'
                char[] chars = new char[string.length()];
                Arrays.fill(chars, 'I');
                out.println(new String(chars));
                rwCounter++;
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        // System.out.println("\tWriting read.mappinginfo...");
        int readCounter = counter;
        try {
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(od + "/read.mappinginfo", append)));
            if (!append)
                out.println("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut");
            for (String string : line.readMappingInfo) {
                out.println(readCounter + "\t" + string);
                readCounter++;
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return readCounter;
    }

    private static ArrayList<Object[]> parseReadcounts(File readcounts) {
        System.out.println("|parseReadcounts| Started...");
        ArrayList<Object[]> result = new ArrayList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(readcounts));
            String line = br.readLine(); // skip header
            while ((line = br.readLine()) != null) {
                // gather readcounts line information
                Object[] obj = new Object[3];
                obj[0] = line.split("\\t")[0]; // gene_id
                obj[1] = line.split("\\t")[1]; // transcript_id
                obj[2] = Integer.parseInt(line.split("\\t")[2]); // counts
                result.add(obj.clone());
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("|parseReadcounts| Finished. Found " + result.size() + " entries.");
        return result;
    }

    private static void testCases(GSE gse) {
        // ------------Test Cases----------------
        String geneId = "ENSG00000241978";
        String transcriptId = "ENST00000374525";
        /*
        TreeSet<Exon> exonsInTranscript = gse.getExonsByTranscriptAndGene(transcriptId, geneId);
        String transcriptSequence = "";
        Exon referenceExon = null;
        for (Exon e : exonsInTranscript) {
            String s = gse.getSequence(e.chr, e.start, e.end, e.strand);
            transcriptSequence = transcriptSequence.concat(s);
            referenceExon = new Exon(e);
        }
        String fragmentSeq = transcriptSequence.substring(2652, 2777);
        String readFw = GSE.getReadSequence(fragmentSeq, 75, false);
        String readRw = GSE.getReadSequence(fragmentSeq, 75, true);
        */

        int posStartTest1 = gse.getGenomicPosition(transcriptId, geneId, 2677, true, "+");
        int posEndTest1 = gse.getGenomicPosition(transcriptId, geneId, 2752, false, "+");
        String genomicRegionTest1 = gse.getGenomicRegionPosition(transcriptId, geneId, posStartTest1, posEndTest1);
        // Expected Output: 112.918.703-112.918.778 (t.end = e.end)
        if (!genomicRegionTest1.equals("112918703-112918778"))
            System.err.println("Test1 Found:\t" + genomicRegionTest1 + "\n\texpected:\t112918703-112918778 (t.end = e.end)");

        int posStartTest2 = gse.getGenomicPosition(transcriptId, geneId, 2573, true, "+");
        int posEndTest2 = gse.getGenomicPosition(transcriptId, geneId, 2648, false, "+");
        String genomicRegionTest2 = gse.getGenomicRegionPosition(transcriptId, geneId, posStartTest2, posEndTest2);
        // Expected Output: 112.918.599-112.918.674 (t.start=e.end)
        if (!genomicRegionTest2.equals("112918599-112918674"))
            System.err.println("Test2 Found:\t" + genomicRegionTest2 + "\n\texpected:\t112918599-112918674 (t.start=e.end)");

        geneId = "ENSG00000160767";
        transcriptId = "ENST00000361361";
        int posStartTest3 = gse.getGenomicPosition(transcriptId, geneId, 2682, true, "-");
        int posEndTest3 = gse.getGenomicPosition(transcriptId, geneId, 2757, false, "-");
        String genomicRegionTest3 = gse.getGenomicRegionPosition(transcriptId, geneId, posStartTest3, posEndTest3);
        // Expected Output: 155.217.333-155.217.408
        if (!genomicRegionTest3.equals("155217333-155217408"))
            System.err.println("Test3 Found:\t" + genomicRegionTest3 + "\n\texpected:\t155217333-155217408 (- strand)");

        int posStartTest4 = gse.getGenomicPosition(transcriptId, geneId, 739, true, "-");
        int posEndTest4 = gse.getGenomicPosition(transcriptId, geneId, 814, false, "-");
        String genomicRegionTest4 = gse.getGenomicRegionPosition(transcriptId, geneId, posStartTest4, posEndTest4);
        // Expected Output: 155.223.964-155.223.986|155.224.191-155.224.244
        if (!genomicRegionTest4.equals("155223964-155223986|155224191-155224244"))
            System.err.println("Test4 Found:\t" + genomicRegionTest4 + "\n\texpected:\t155223964-155223986|155224191-155224244");

        int posStartTest5 = gse.getGenomicPosition(transcriptId, geneId, 2020, true, "-");
        int posEndTest5 = gse.getGenomicPosition(transcriptId, geneId, 2095, false, "-");
        String genomicRegionTest5 = gse.getGenomicRegionPosition(transcriptId, geneId, posStartTest5, posEndTest5);
        // Expected Output: 155.218.089-155.218.095|155.218.196-155.218.265
        if (!genomicRegionTest5.equals("155218089-155218095|155218196-155218265"))
            System.err.println("Test5 Found:\t" + genomicRegionTest5 + "\n\texpected:\t155218089-155218095|155218196-155218265");

        geneId = "ENSG00000162946";
        transcriptId = "ENST00000439617";
        int posStartTest6 = gse.getGenomicPosition(transcriptId, geneId, 1680, true, "+");
        int posEndTest6 = gse.getGenomicPosition(transcriptId, geneId, 1755, false, "+");
        String genomicRegionTest6 = gse.getGenomicRegionPosition(transcriptId, geneId, posStartTest6, posEndTest6);
        // Expected Output: 231906810-231906817|231930988-231931043|231935854-231935867 (3 exons)
        if (!genomicRegionTest6.equals("231906810-231906817|231930988-231931043|231935854-231935867"))
            System.err.println("Test6 Found:\t" + genomicRegionTest6 + "\n\texpected:\t231906810-231906817|231930988-231931043|231935854-231935867 (3 exons)");
    }

}
