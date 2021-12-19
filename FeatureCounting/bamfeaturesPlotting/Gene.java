package com.company;

import augmentedTree.Interval;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Gene implements Interval {
    private int start;
    private int stop;
    private ArrayList<String> transcriptIds;
    private String geneId;
    private String biotype;
    private String chr;
    private HashMap<String, RegionVector> transcripts;
    private boolean strand;

    public Gene(int start, int stop, String id, String biotype, String chromosome, boolean strand){
        this.start = start;
        this.stop = stop;
        this.geneId = id;
        this.biotype = biotype;
        this.chr = chromosome;
        this.transcriptIds = new ArrayList<>();
        this.transcripts = new HashMap<>();
        this.strand = strand;
    }

    public boolean getStrand(){
        return this.strand;
    }

    @Override
    public int getStart() {
        return this.start;
    }

    @Override
    public int getStop() {
        return this.stop;
    }

    public String getId(){
        return this.geneId;
    }

    public void addExon(int start, int end, String trId){
        this.transcripts.get(trId).addRegion(new int[]{start, end});
    }

    public String getChromosome(){
        return this.chr;
    }

    public String getBiotype(){
        return this.biotype;
    }

    public void addTranscript(String trId){
        this.transcriptIds.add(trId);
        this.transcripts.put(trId, new RegionVector(new ArrayList<>()));
    }

    public String matchPlausibilityLevel(RegionVector sr1, RegionVector sr2){
        // returns match plausibility level as String: either "1" or "2" or "3transcriptId1,transcriptId2,..."
        HashSet<String> Subregion = new HashSet<>();

        // check all transcripts: are the region vectors sub-RVs?
        for(String id:transcripts.keySet()){
            boolean a = transcripts.get(id).hasAsSubRV(sr1);
            boolean b = transcripts.get(id).hasAsSubRV(sr2);
            if(a && b) Subregion.add(id);
        }

        // if Subregion contains at least one transcript => level 3 (transcriptomic)!
        if(Subregion.size() > 0){
            StringBuilder sb = new StringBuilder();
            sb.append("3");
            for(String t:Subregion){
                sb.append(t).append(",");
            }
            String s = sb.toString();
            s = s.substring(0, s.length()-1);
            return s;
        }

        // check if its merged-transcriptomic
        // merge all transcript RVs to one RV
        RegionVector mergedTranscripts = transcripts.get(transcriptIds.get(0));
        for(int i = 1; i < transcriptIds.size(); i++){
            mergedTranscripts = mergedTranscripts.merge(transcripts.get(transcriptIds.get(i)));
        }

        if(mergedTranscripts.Contains(sr1) && mergedTranscripts.Contains(sr2)) {
            // level 2 (merged-transcriptomic)
            return "2";
        }
        // level 1 (intronic)
        return "1";
    }

    public int getMergedTranscriptsLength(){
        RegionVector merged = null;
        for(RegionVector rv:transcripts.values()){
            if(merged == null) merged = rv.copy();
            else {
                merged = merged.merge(rv);
            }
        }
        return merged.getNumBases();
    }
}
