package com.company;

import java.util.ArrayList;

public class Protein {
    private String id;
    private RegionVector cds;
    private RegionVector introns;

    public Protein(String id){
        this.id = id;
        this.cds = new RegionVector(new ArrayList<int[]>());
    }

    public String getId(){
        return this.id;
    }

    public RegionVector getCds(){
        return this.cds;
    }

    public RegionVector getIntrons() {
        // first compute the introns from the set of exons, then return them
        this.introns = cds.getInverse();
        return this.introns;
    }

    public void addCds(int[] newCds){
        this.cds.addRegion(newCds);
    }

    // checks if protein has an intron with specific start
    public boolean checkIntronStart(int start){
        if(this.introns==null) this.introns = cds.getInverse();
        return this.introns.searchRegion(start, null);
    }

    // checks if protein has an intron with specific end
    public boolean checkIntronEnd(int end){
        if(this.introns==null) this.introns = cds.getInverse();
        return this.introns.searchRegion(null, end);
    }

    public RegionVector getIntronsWithin(int[] otherIntron){
        if(this.introns==null) this.introns = cds.getInverse();
        return this.introns.findRegionsIn(otherIntron);
    }
}
