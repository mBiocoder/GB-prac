package com.company;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class RegionVector implements Comparable<RegionVector> {
    private ArrayList<int[]> regions;
    private boolean transcriptStrand;

    public RegionVector(ArrayList<int[]> regions) {
        this.regions = regions;
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

    public RegionVector copy() {
        ArrayList<int[]> regions = new ArrayList<>();
        for (int[] reg : this.regions) {
            regions.add(new int[]{reg[0], reg[1]});
        }
        return new RegionVector(regions);
    }

    public ArrayList<int[]> trim(int start, int end) {
        // returns copy of regions within a certain range
        ArrayList<int[]> regions = new ArrayList<>();
        for (int[] reg : this.regions) {
            regions.add(new int[]{reg[0], reg[1]});
        }

        ArrayList<int[]> remove = new ArrayList<>();
        for (int[] reg : regions) {
            if (reg[1] < start) remove.add(reg);
            else if (reg[0] > end) remove.add(reg);
            else {
                if (reg[0] < start) reg[0] = start;
                if (reg[1] > end) reg[1] = end;
            }
        }
        regions.removeAll(remove);
        return regions;
    }

    public int getStart() {
        return this.regions.get(0)[0];
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

    public ArrayList<int[]> skippedInOther(RegionVector other) {
        // return all the regions from this RV, that are skipped in another one
        ArrayList<int[]> skipped = new ArrayList<>();
        ArrayList<int[]> otherRegions = other.getRegions();

        int i = 0;
        int j = 0;
        ArrayList<int[]> multipleSkipped = new ArrayList<>();
        while (i < this.regions.size() && j < otherRegions.size()) {
            if (this.regions.get(i)[1] < otherRegions.get(j)[0]) {
                if (j > 0 && i > 0 && i < this.regions.size() - 1) {
                    if (this.regions.get(i - 1)[1] == otherRegions.get(j - 1)[1] && this.regions.get(i + 1)[0] == otherRegions.get(j)[0]) {
                        skipped.add(this.regions.get(i));
                    } else if (this.regions.get(i - 1)[1] == otherRegions.get(j - 1)[1]) {
                        multipleSkipped.add(this.regions.get(i));
                    } else if (multipleSkipped.size() > 0) {
                        if (this.regions.get(i + 1)[0] == otherRegions.get(j)[0]) {
                            skipped.add(this.regions.get(i));
                            skipped.addAll(multipleSkipped);
                            multipleSkipped.clear();
                        } else multipleSkipped.add(this.regions.get(i));
                    }
                }
                i++;
            } else if (this.regions.get(i)[0] > otherRegions.get(j)[1]) {
                j++;
            } else {
                i++;
            }
        }

        return skipped;
    }

    public boolean containsRegionPartially(int start, int end) {
        // return  true if RV contains a region that is either start-end or smaller (any part of start-end) but not bigger!
        for (int[] reg : this.regions) {
            if (end < reg[0]) return false;
            else if (start <= reg[0] && end >= reg[1]) return true;
        }
        return false;
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

    // returns true if no part of the specified region is in the RV but it has regions on both sides
    public boolean skipsRegion(int start, int end) {
        boolean cuts = false;
        boolean left = false;
        boolean right = false;

        for (int[] reg : this.regions) {
            if (reg[1] < start) left = true;
            else if (reg[0] > end) right = true;
            else cuts = true;
        }
        return !cuts && left && right;
    }

    // get region to the left and to the right of a specified region
    public int[][] getNeighbors(int start, int end) {
        int[][] neighbors = new int[2][];

        for (int i = 0; i < this.regions.size(); i++) {
            if (this.regions.get(i)[0] == start && this.regions.get(i)[1] == end) {
                if (i > 0) {
                    neighbors[0] = this.regions.get(i - 1);
                }
                if (i < this.regions.size() - 1) {
                    neighbors[1] = this.regions.get(i + 1);
                }
            }
        }
        return neighbors;
    }

    // get region to the left and to the right of a specified region that is not included
    public int[][] getNeighborsOfNotIncludedRegion(int start, int end) {
        int[][] neighbors = new int[2][];

        for (int i = 0; i < this.regions.size(); i++) {
            if (this.regions.get(i)[1] < start && i < this.regions.size() - 1 && this.regions.get(i + 1)[0] > end) {
                neighbors[0] = this.regions.get(i);
                neighbors[1] = this.regions.get(i + 1);
                return neighbors;
            }
        }
        return null;
    }

    public boolean coversNeighborsOfRegion(int start, int end, RegionVector other) {
        int[][] neighbors = this.getNeighbors(start, end);
        ArrayList<int[]> left = new ArrayList<>();
        left.add(neighbors[0]);
        ArrayList<int[]> right = new ArrayList<>();
        right.add(neighbors[1]);

        RegionVector leftNeighbor = new RegionVector(left);
        RegionVector rightNeighbor = new RegionVector(right);

        return (other.overlap(leftNeighbor).getNumRegions() > 0 && other.overlap(rightNeighbor).getNumRegions() > 0);
    }

    // gets the neighboring Regions of a set of exons (if 2 skipped exons are next to each other, the neighboring exon cannot be a neighbor!)
    public HashMap<String, RegionVector> getNeighborsOfExons(ArrayList<String> exons) {
        HashMap<String, RegionVector> exonToNeighbors = new HashMap<>();
        // compute neighbors
        for (String exon : exons) {
            int a = Integer.parseInt(exon.split("-")[0]);
            int b = Integer.parseInt(exon.split("-")[1]) - 1;
            int[][] neighbors = getNeighbors(a, b);
            ArrayList<int[]> reg = new ArrayList<>();
            reg.add(neighbors[0]);
            reg.add(neighbors[1]);
            exonToNeighbors.put(exon, new RegionVector(reg));
        }
        // check if one of the computed neighbors is a skipped exon
        for (String exon : exonToNeighbors.keySet()) {
            int[] left = exonToNeighbors.get(exon).getRegions().get(0);
            String leftStr = left[0] + "-" + (left[1] + 1);
            int[] right = exonToNeighbors.get(exon).getRegions().get(1);
            String rightStr = right[0] + "-" + (right[1] + 1);

            // replace left neighbor with left neighbor of skipped exon to our left
            while (exons.contains(leftStr)) {
                exonToNeighbors.get(exon).getRegions().remove(0);
                exonToNeighbors.get(exon).getRegions().add(0, exonToNeighbors.get(leftStr).getRegions().get(0));
                left = exonToNeighbors.get(exon).getRegions().get(0);
                leftStr = left[0] + "-" + (left[1] + 1);
            }
            // same to the right
            while (exons.contains(rightStr)) {
                exonToNeighbors.get(exon).getRegions().remove(1);
                exonToNeighbors.get(exon).getRegions().add(exonToNeighbors.get(rightStr).getRegions().get(1));
                right = exonToNeighbors.get(exon).getRegions().get(1);
                rightStr = right[0] + "-" + (right[1] + 1);
            }
        }
        return exonToNeighbors;
    }

    // gets the neighboring Regions of a set of exons that are not in the RV (if 2 skipped exons are next to each other, the neighboring exon cannot be a neighbor!)
    public HashMap<String, RegionVector> getNeighborsOfSkippedExons(ArrayList<String> exons) {
        HashMap<String, RegionVector> exonToNeighbors = new HashMap<>();
        // compute neighbors
        for (String exon : exons) {
            int a = Integer.parseInt(exon.split("-")[0]);
            int b = Integer.parseInt(exon.split("-")[1]) - 1;
            int[][] neighbors = getNeighborsOfNotIncludedRegion(a, b);
            ArrayList<int[]> reg = new ArrayList<>();
            if (neighbors != null) {
                reg.add(neighbors[0]);
                reg.add(neighbors[1]);
            }
            exonToNeighbors.put(exon, new RegionVector(reg));
        }
        // check if one of the computed neighbors is a skipped exon
        /*for (String exon : exonToNeighbors.keySet()) {
            int[] left = exonToNeighbors.get(exon).getRegions().get(0);
            String leftStr = left[0] + "-" + (left[1] + 1);
            int[] right = exonToNeighbors.get(exon).getRegions().get(1);
            String rightStr = right[0] + "-" + (right[1] + 1);

            // replace left neighbor with left neighbor of skipped exon to our left
            while (exons.contains(leftStr)) {
                exonToNeighbors.get(exon).getRegions().remove(0);
                exonToNeighbors.get(exon).getRegions().add(0, exonToNeighbors.get(leftStr).getRegions().get(0));
                left = exonToNeighbors.get(exon).getRegions().get(0);
                leftStr = left[0] + "-" + (left[1] + 1);
            }
            // same to the right
            while (exons.contains(rightStr)) {
                exonToNeighbors.get(exon).getRegions().remove(1);
                exonToNeighbors.get(exon).getRegions().add(exonToNeighbors.get(rightStr).getRegions().get(1));
                right = exonToNeighbors.get(exon).getRegions().get(1);
                rightStr = right[0] + "-" + (right[1] + 1);
            }
        }*/
        return exonToNeighbors;
    }

    public int[] getFirstRegion() {
        return this.regions.get(0);
    }

    public int[] getLastRegion() {
        return this.regions.get(this.regions.size() - 1);
    }

    public HashSet<String> whichSkipped(RegionVector other) {
        HashSet<String> skipped = new HashSet<>();
        ArrayList<int[]> otherRegions = other.getRegions();
        ArrayList<int[]> thisRegions = this.trim(other.getStart(), other.getEnd());
        int i = 0;
        int j = 0;
        while (i < thisRegions.size() && j < otherRegions.size()) {
            int[] thisReg = thisRegions.get(i);
            int[] otherReg = otherRegions.get(j);
            if (thisReg[0] != otherReg[0] || thisReg[1] != otherReg[1]) {
                if (j > 0) {
                    if (otherReg[0] > thisReg[1] && otherRegions.get(j - 1)[1] < thisReg[0]) {
                        String e = thisReg[0] + "-" + (thisReg[1] + 1);
                        skipped.add(e);
                        i++;
                    } else return skipped;
                } else return skipped;
            } else {
                i++;
                j++;
            }
        }
        return skipped;
    }

    // returns regions that lie between start and end but no region is cut (only complete regions)
    public ArrayList<String> completeRegionsInRange(int start, int end) {
        if (start >= end) return null;
        ArrayList<String> result = new ArrayList<>();
        boolean take = false;
        for (int[] reg : this.regions) {
            String e = reg[0] + "-" + (reg[1] + 1);
            if (take) {
                if (reg[1] < end) result.add(e);
                else break;
            } else {
                if (reg[1] >= start) take = true;
                if (take) {
                    if (reg[0] >= start) result.add(e);
                    else if (end <= reg[1]) break;
                }
            }
        }
        return result;
    }
}
