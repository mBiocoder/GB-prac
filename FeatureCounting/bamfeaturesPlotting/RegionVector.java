package com.company;

import java.util.ArrayList;

public class RegionVector implements Comparable<RegionVector> {
    private ArrayList<int[]> regions;
    // to count how often a region vector was found in the genome
    private int number;
    private boolean transcriptStrand;

    public RegionVector(ArrayList<int[]> regions) {
        this.regions = regions;
        this.number = 1;
    }

    public int getNumber() {
        return this.number;
    }

    public void incNumber() {
        // increase the number of this RegionVector by one
        this.number++;
    }

    // returns the regions between its own regions.
    public RegionVector getReverse() {
        ArrayList<int[]> reverseRegions = new ArrayList<>();
        Integer prev = null;
        for (int[] region : this.regions) {
            // numbers are inclusive!!!
            if (prev != null && region[0] != prev) reverseRegions.add(new int[]{prev, region[0] - 1});
            prev = region[1] + 1;
        }
        return new RegionVector(reverseRegions);
    }

    public ArrayList<int[]> getRegions() {
        return this.regions;
    }

    // calculates the overlapping regions of 2 RegionVector Objects
    public RegionVector overlap(RegionVector other) {
        ArrayList<int[]> overlaps = new ArrayList<>();
        ArrayList<int[]> otherRegions = other.getRegions();
        if (this.regions.size() == 0 || otherRegions.size() == 0) return new RegionVector(overlaps);

        int j = 0;
        for (int[] thisReg : this.regions) {
            int[] otherReg = otherRegions.get(j);
            // do nothing if other region is completely bigger than this region
            if (otherReg[0] <= thisReg[1]) {
                while (otherReg[1] < thisReg[0]) {
                    j++;
                    if (j >= otherRegions.size()) return new RegionVector(overlaps);
                    otherReg = otherRegions.get(j);
                }
                // check for overlap:
                if (otherReg[0] <= thisReg[1]) {
                    overlaps.add(new int[]{Math.max(thisReg[0], otherReg[0]), Math.min(thisReg[1], otherReg[1])});
                }
            }
        }
        return new RegionVector(overlaps);
    }

    public int getNumRegions() {
        return this.regions.size();
    }

    public void addRegion(int[] region) {
        // add region either in front or back (so that region vector stays sorted)
        // don't need to check all possible positions, bc regions occur either forward or backward (in the gtf)
        // returns if region was added in back or not
        if (this.regions.size() == 0) {
            this.regions.add(region);
        } else {
            if (this.regions.get(this.regions.size() - 1)[1] < region[0]) {
                // add in the back
                if (region[0] - 1 == this.regions.get(this.getNumRegions() - 1)[1]) {
                    this.regions.get(this.getNumRegions() - 1)[1] = region[1];
                } else {
                    this.regions.add(region);
                }
            } else {
                // add in the front
                if (region[1] + 1 == this.regions.get(0)[0]) {
                    this.regions.get(0)[0] = region[0];
                } else {
                    this.regions.add(0, region);
                }
            }
        }
    }

    // returns true if other RV is sub RV of this one
    public boolean hasAsSubRV(RegionVector other) {
        ArrayList<int[]> otherRegions = other.getRegions();
        int i = 0;
        int j = 0;
        while (i < this.regions.size() && j < otherRegions.size()) {
            int[] thisReg = this.regions.get(i);
            int[] otherReg = otherRegions.get(j);
            if (thisReg[1] < otherReg[0]) {
                if (j != 0) return false;
                i++;
            } else if (otherReg[1] < thisReg[0]) return false;
            else {
                // check if start and end fit
                if ((otherReg[0] < thisReg[0]) || (otherReg[0] > thisReg[0] && j != 0)) return false;
                if ((otherReg[1] > thisReg[1]) || (otherReg[1] < thisReg[1] && j != otherRegions.size() - 1))
                    return false;
                i++;
                j++;
            }
        }
        // check if all regions of other were checked
        return j == otherRegions.size();
    }

    // returns true if this RV contains the other one
    public boolean Contains(RegionVector other) {
        ArrayList<int[]> otherRegions = other.getRegions();
        int i = 0;
        int j = 0;
        while (i < this.regions.size() && j < otherRegions.size()) {
            int[] thisReg = this.regions.get(i);
            int[] otherReg = otherRegions.get(j);
            if (thisReg[1] < otherReg[0]) {
                i++;
            } else if (otherReg[1] < thisReg[0]) return false;
            else {
                // check if other region lies within this region
                if (otherReg[0] >= thisReg[0] && otherReg[1] <= thisReg[1]) {
                    j++;
                } else return false;
            }
        }
        // check if all regions of other were checked
        return j == otherRegions.size();
    }

    @Override
    public int compareTo(RegionVector other) {
        // compare 2 region vectors (according to their end (and if thats the same starting) position if they are not equal)
        // 0 => both are the same, 1 => this has greater end, -1 => other has greater end
        int dist = this.getEnd() - other.getEnd();
        if (dist == 0) {
            // check if they are exactly the same, else return 1
            if (this.equals(other)) return 0;
            else return 1;
        }
        if (dist > 0) return 1;
        return -1;
    }

    public int getEnd() {
        if (this.regions.size() > 0) {
            return this.regions.get(regions.size() - 1)[1];
        } else return 0;
    }

    public boolean equals(Object other) {
        if (!(other instanceof RegionVector)) return false;
        ArrayList<int[]> otherRegions = ((RegionVector) other).getRegions();
        if (otherRegions.size() != this.regions.size()) return false;
        for (int i = 0; i < this.regions.size(); i++) {
            if (this.regions.get(i)[0] != otherRegions.get(i)[0] || this.regions.get(i)[1] != otherRegions.get(i)[1])
                return false;
        }
        return true;
    }

    public RegionVector merge(RegionVector other) {
        // return resulting RegionVector when this is merged with another one
        ArrayList<int[]> newRegions = new ArrayList<>();
        ArrayList<int[]> otherRegions = other.getRegions();
        int i = 0;
        int j = 0;
        while (i < this.regions.size() && j < otherRegions.size()) {
            int[] thisReg = this.regions.get(i);
            int[] otherReg = otherRegions.get(j);
            if (thisReg[1] < otherReg[0]) {
                newRegions.add(new int[]{thisReg[0], thisReg[1]});
                i++;
            } else if (otherReg[1] < thisReg[0]) {
                newRegions.add(new int[]{otherReg[0], otherReg[1]});
                j++;
            } else {
                // there is an overlap. If regions don't agree, use bigger alternative
                newRegions.add(new int[]{Math.min(otherReg[0], thisReg[0]), Math.max(otherReg[1], thisReg[1])});
                i++;
                j++;
            }
        }
        while (i < this.regions.size()) {
            newRegions.add(new int[]{this.regions.get(i)[0], this.regions.get(i)[1]});
            i++;
        }
        while (j < otherRegions.size()) {
            newRegions.add(new int[]{otherRegions.get(j)[0], otherRegions.get(j)[1]});
            j++;
        }

        // unite neighboring regions that overlap/have no space in between
        if (newRegions.size() > 0) {
            ArrayList<int[]> delete = new ArrayList<>();
            int[] prev = newRegions.get(0);
            for (int k = 1; k < newRegions.size(); k++) {
                int[] reg = newRegions.get(k);
                if (reg[0] <= prev[1] + 1) {
                    prev[1] = Math.max(prev[1], reg[1]);
                    delete.add(reg);
                } else prev = reg;

            }
            for (int[] r : delete) {
                newRegions.remove(r);
            }
        }
        RegionVector result = new RegionVector(newRegions);
        return result;
    }

    public int numUniqueRegions(RegionVector other) {
        int unique = 0;
        int i = 0;
        int j = 0;
        while (i < this.getNumRegions() && j < other.getNumRegions()) {
            int[] thisReg = this.regions.get(i);
            int[] otherReg = other.getRegions().get(j);
            if (thisReg[1] < otherReg[0]) {
                unique++;
                i++;
            } else if (otherReg[1] < thisReg[0]) {
                unique++;
                j++;
            } else {
                if (thisReg[0] == otherReg[0] && thisReg[1] == otherReg[1]) unique++;
                else unique += 2;
                i++;
                j++;
            }
        }
        while (i < this.getNumRegions()) {
            unique++;
            i++;
        }
        while (j < other.getNumRegions()) {
            unique++;
            j++;
        }
        return unique;
    }

    public void setTranscriptStrand(boolean strand) {
        this.transcriptStrand = strand;
    }

    private boolean getTranscriptStrand() {
        return this.transcriptStrand;
    }

    public boolean sameStrand(RegionVector other) {
        return this.transcriptStrand == other.getTranscriptStrand();
    }

    public void setNumber(int number) {
        this.number = number;
    }

    public RegionVector copy() {
        ArrayList<int[]> newRegions = new ArrayList<>();
        for (int[] region : this.regions) {
            newRegions.add(new int[]{region[0], region[1]});
        }
        RegionVector copy = new RegionVector(newRegions);
        copy.setTranscriptStrand(this.transcriptStrand);
        copy.setNumber(this.number);
        return copy;
    }

    public int getNumBases() {
        int num = 0;
        for (int[] region : this.regions) {
            num += (region[1] - region[0] + 1);
        }
        return num;
    }
}
