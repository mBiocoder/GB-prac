import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.*;

public class GSE {
    // Genome Sequence Extractor
    File fasta;
    File idx;
    HashMap<String, long[]> indices = new HashMap<>();
    File gtf;
    HashMap<String, Gene> genes = new HashMap<>();

    // Constructor
    public GSE(File fasta, File idx, File gtf) {
        super();
        this.fasta = fasta;
        this.idx = idx;
        this.indices = parseIdx(idx);
        this.gtf = gtf;
        this.genes = parseGtf(gtf);
    }


    public String getSequence(String chr, int seqStart, int seqEnd, String strand) {
        String result = "";
        //chrStart and chrLength are filled
        long chrLength = indices.get(chr)[0];
        long chrStart = indices.get(chr)[1];
        int lineSize = (int) indices.get(chr)[2];

        // Error checking
        if (seqEnd - seqStart > chrLength)
            System.err.println("Required seq length is longer than chr itself.");
        if (seqStart > seqEnd)
            System.err.println("seqStart > seqEnd.");
        if (seqStart < 0 || seqEnd < 1)
            System.err.println("seqStart or seqEnd < 0");

        try {
            StringBuilder sb = new StringBuilder();
            // ("r"=readOnly)
            RandomAccessFile raf = new RandomAccessFile(fasta, "r");
            // read chr from begin to end
            long rafStart = chrStart + seqStart + (seqStart / lineSize);
            // Jump to rafStart location
            raf.seek(rafStart);
            while (sb.length() < seqEnd - seqStart) {
                byte b;
                b = raf.readByte();
                if ((char) b == 'A' || (char) b == 'C' || (char) b == 'T' || (char) b == 'G' || (char) b == 'N') {
                    sb.append((char) b);
                }
            }
            result = sb.toString();
            //System.out.println(result);
            raf.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        if (strand.equals("-")) {
            // Inverse result string
            result = new StringBuilder(result).reverse().toString();
            result = invertSequence(result);
        }

        // Error checking
        if (result.isEmpty())
            System.err.println("No such sequence on chr:" + chr + "(" + seqStart + "," + seqEnd + ")");
        if (result.length() != (seqEnd - seqStart))
            System.err.println("Result length (" + result.length() + " does not correspond to seqStart-SeqEnd+1: " + seqStart + " - " + seqEnd + " = " + (seqEnd - seqStart));

        // System.out.println("|getSequence| Input chr" + chr + "(" + seqStart + "-" + seqEnd + ") strand: " + strand + ". Result string length: " + result.length() + "/" + (seqEnd - seqStart));
        return result;
    }

    public TreeSet<Exon> getExonsByTranscriptAndGene(String transcriptId, String geneId) {
        // Returns sorted TreeSet of Exons belonging to transcript-/geneId
        TreeSet<Exon> exons = this.genes.get(geneId).transcripts.get(transcriptId).exons;
        if (exons.size() == 0) {
            System.err.println("No exon, which belongs to transcriptID" + transcriptId + "was found!");
            return null;
        } else
            return exons;
    }

    public static int getFragmentLength(int mean, int standardDeviation, int readLength, int transcriptLength) {
        // returns a value from a normal distribution with a given mean and standardDeviation
        Random r = new Random();
        double number = r.nextGaussian() * (standardDeviation + 10) + mean;
        // re-draws number, if smaller then readLength or larger then transcriptLength
        while (number < readLength || number > transcriptLength)
            number = getFragmentLength(mean, standardDeviation, readLength, transcriptLength);
        return (int) number;
    }

    public static String getReadSequence(String fragment, int readlength, boolean reverseComplement) {
        // returns ReadSequence of a certain length for fragment
        String result;
        if (!reverseComplement)
            result = fragment.substring(0, readlength);
        else {
            // second read is reverse complement
            result = fragment.substring(fragment.length() - readlength, fragment.length());
            result = new StringBuilder(result).reverse().toString();
            result = invertSequence(result);
        }
        return result;
    }

    public static String[] simulateMutations(String sequence, double mutationratePercent) {
        // get mutated sequence and positions of mutations as string
        // result[0] = mutatedSeq, result[1] = cs positions of mutated chars
        String[] result = {"", ""};

        // go through all chars in sequence. select a random number [0,100]. If this number <= mutationrate: mutate the current char
        StringBuilder mutatedSeq = new StringBuilder(sequence);
        for (int i = 0; i < sequence.length(); i++) {
            int randomNumber = new Random().nextInt(100);
            if (randomNumber <= (mutationratePercent - 1)) {
                result[1] = result[1].concat(i + ",");
                if (mutatedSeq.charAt(i) == 'A') {
                    mutatedSeq.setCharAt(i, 'G');
                } else if (mutatedSeq.charAt(i) == 'T') {
                    mutatedSeq.setCharAt(i, 'C');
                } else if (mutatedSeq.charAt(i) == 'G') {
                    mutatedSeq.setCharAt(i, 'T');
                } else if (mutatedSeq.charAt(i) == 'C') {
                    mutatedSeq.setCharAt(i, 'A');
                } else {
                    System.err.println("|simulateMutations| Found inappropriate char: " + mutatedSeq.charAt(i) + " in "
                            + sequence + ".");
                }
            }
        }
        if (!result[1].isEmpty())
            result[1] = result[1].substring(0, result[1].length() - 1); // trimm comma
        result[0] = new String(mutatedSeq);
        return result;
    }

    public int getGenomicPosition(String transcriptId, String geneId, int positionInTranscript, boolean start, String strand) {
        // return the position on the gene of a position on a transcript. Go through all exons and reduce positionInTranscript every time we pass an exon
        // if boolean start = true --> positionInTranscript differs by 1
        TreeSet<Exon> exons = this.genes.get(geneId).transcripts.get(transcriptId).exons;

        if (strand.equals("+")) {
            for (Exon exon : exons) {
                int exonLength = exon.end - exon.start;
                if (exonLength < 1)
                    System.err.println("Length cannot be < 1");
                if (exonLength < positionInTranscript)
                    positionInTranscript -= exonLength;
                else if (exonLength == positionInTranscript && !start)
                    return exon.end + 1;
                else if (exonLength == positionInTranscript && start)
                    positionInTranscript = 0;
                else
                    return exon.start + positionInTranscript + 1;
            }
        } else if (strand.equals("-")) {
            for (Exon exon : exons) {
                int exonLength = exon.end - exon.start;
                if (exonLength < 1)
                    System.err.println("Length cannot be < 1");
                if (exonLength < positionInTranscript)
                    positionInTranscript -= exonLength;
                else if (exonLength == positionInTranscript && !start)
                    return exon.start + 1;
                else if (exonLength == positionInTranscript && start)
                    positionInTranscript = 0;
                else
                    return exon.end - positionInTranscript + 1;
            }
        } else
            System.err.println("|getGenomicPosition| Strand: " + strand + " not found.");

        System.err.println("|getGenomicPosition| positionInTranscript " + positionInTranscript + " not found. geneId: " + geneId + " transcriptId: " + transcriptId);
        return 0;
    }

    public String getGenomicRegionPosition(String transcriptId, String geneId, int genomicStart, int genomicEnd) {
        if (genomicStart > genomicEnd) { // eventually, swap start<->end
            int temp = genomicStart;
            genomicStart = genomicEnd;
            genomicEnd = temp;
        }
        String res = "";
        // TreeSet<Exon> exons = this.genes.get(geneId).transcripts.get(transcriptId).exons;
        Set<Exon> exons = new TreeSet<>((e1, e2) -> { // Sorts Exons ascending to their start-points
            if (e1.start > e2.start) {
                return 1;
            } else return -1;
        });
        exons.addAll(this.genes.get(geneId).transcripts.get(transcriptId).exons);

        for (Exon exon : exons) {
            if (exon.end >= genomicStart) {
                if (exon.end < genomicEnd - 1) {
                    res += (exon.start + 1) + "-" + (exon.end + 1) + "|";
                } else if (exon.end == genomicEnd - 1) {
                    res += (exon.start + 1) + "-" + (exon.end + 1);
                    break;
                } else {
                    res += (exon.start + 1) + "-" + genomicEnd;
                    break;
                }
            }
        }

        res = res.substring(res.indexOf("-"), res.length());
        res = genomicStart + res;
        return res;
    }

    private HashMap<String, long[]> parseIdx(File idx) {
        // Go through idx file & fetch start and length for a chromosome
        // Key = chr_id, long[0] = length, long[1] = start point, long[2] = line size.
        HashMap<String, long[]> result = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(idx));
            String line;
            while ((line = br.readLine()) != null) {
                String[] tabSeparated = line.split("\\t");
                String chr_id = tabSeparated[0];
                long length = Long.parseLong(tabSeparated[1]);
                long start = Long.parseLong(tabSeparated[2]);
                long lineSize = Long.parseLong(tabSeparated[3]);
                long[] temp = {length, start, lineSize};
                result.put(chr_id, temp.clone());
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return result;
    }

    private static HashMap<String, Gene> parseGtf(File gtfInput) {
        // Key is the always the respective ID
        HashMap<String, Gene> allGenes = new HashMap<>();
        int exonCounter = 0;
        int transcriptCounter = 0;
        int geneCounter = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(gtfInput));
            String line;
            int linecounter = 0; // for error printing only
            System.out.println("Begin: Parsing gtf file.");
            while ((line = br.readLine()) != null) {
                linecounter++;
                // ignore comments (beginning with "#")
                if (!line.startsWith("#")) {
                    // For every line in input
                    String[] tabSeparated = line.split("\\t");
                    String seqname = tabSeparated[0];
                    // String source = tabSeparated[1];
                    String feature = tabSeparated[2];
                    String start = tabSeparated[3];
                    String end = tabSeparated[4];
                    // String score = tabSeparated[5];
                    String strand = tabSeparated[6];
                    // String frame = tabSeparated[7];
                    String attribute = tabSeparated[8];
                    // -------For lines, which are exons:-------
                    if (feature.equalsIgnoreCase("exon")) {
                        exonCounter++;
                        // parameters needed to construct new Exon()
                        String exon_id = "";
                        String exon_number = "";
                        String transcript_id = "";
                        String transcript_name = "";
                        String gene_id = "";
                        String gene_name = "";
                        // gather parameters from String "attribute"
                        String[] attributeSeparated = attribute.split(";");
                        for (int i = 0; i < attributeSeparated.length; i++) {
                            if (attributeSeparated[i].contains("exon_id")) {
                                // get only value between quotation marks
                                exon_id = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            } else if (attributeSeparated[i].contains("transcript_id"))
                                transcript_id = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            else if (attributeSeparated[i].contains("transcript_name"))
                                transcript_name = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            else if (attributeSeparated[i].contains("exon_number"))
                                exon_number = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            else if (attributeSeparated[i].contains("gene_id"))
                                gene_id = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            else if (attributeSeparated[i].contains("gene_name"))
                                gene_name = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                        }
                        Exon e = new Exon(seqname, Integer.parseInt(start) - 1, Integer.parseInt(end), strand, exon_id, Integer.parseInt(exon_number), transcript_id, transcript_name, gene_id, gene_name);
                        // Check, if all exon values are actually filled. Note: e.gene_name can be empty, is not checked
                        if (e.start == 0 || e.end == 0 || e.chr.isEmpty() || e.strand.isEmpty() || e.exon_id.isEmpty() || e.exon_number == 0 || e.transcript_id.isEmpty() || e.gene_id.isEmpty())
                            System.err.println("Exon in line " + linecounter + " has an empty value!");
                        // Create data structure allGenes, contains transcripts, contains exons
                        if (!allGenes.containsKey(e.gene_id)) {
                            geneCounter++;
                            allGenes.put(e.gene_id, new Gene());
                        }
                        if (!allGenes.get(e.gene_id).transcripts.containsKey(e.transcript_id)) {
                            transcriptCounter++;
                            allGenes.get(e.gene_id).transcripts.put(e.transcript_id, new Transcript());
                        }
                        exonCounter++;
                        allGenes.get(e.gene_id).transcripts.get(e.transcript_id).exons.add(new Exon(e));
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Finished: Parsing file.");
        System.out.println("\texonCounter: " + exonCounter);
        System.out.println("\ttranscriptCounter:" + transcriptCounter);
        System.out.println("\tgeneCounter:" + geneCounter);
        return allGenes;
    }

    public static String invertSequence(String transcriptSeq) {
        // invert given Sequence: A<->T, G<->C
        char[] result = new char[transcriptSeq.length()];
        for (int i = 0; i < transcriptSeq.length(); i++) {
            if (transcriptSeq.charAt(i) == (char) 'A')
                result[i] = 'T';
            else if (transcriptSeq.charAt(i) == (char) 'T')
                result[i] = 'A';
            else if (transcriptSeq.charAt(i) == (char) 'G')
                result[i] = 'C';
            else if (transcriptSeq.charAt(i) == (char) 'C')
                result[i] = 'G';
            else
                System.err.println("Invalid character found in transcriptSeq (" + (transcriptSeq.charAt(i) + ")"));
        }
        return new String(result);
    }

}
