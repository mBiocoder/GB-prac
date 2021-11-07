package com.company;

import java.util.ArrayList;

public class RegionVector {
    // in gtf: sequence numbering starts at 1, both start and end position are INCLUDED
    // here: start and end position are INCLUDED
    private ArrayList<int[]> regions;

    public RegionVector(ArrayList regions) {
        this.regions = regions;
    }

    public ArrayList<int[]> getRegions() {
        return this.regions;
    }

    // returns the regions of a vector as Strings: "<start>:<end>"
    public ArrayList<String> getRegionsAsString() {
        ArrayList<String> reg = new ArrayList<>();
        for (int[] r : this.regions) reg.add(r[0] + ":" + r[1]);
        return reg;
    }

    public void addRegion(int[] newRegion) {
        // regions should be sorted so insert new region at the right place
        // special cases: regions are empty or new region should be added in last position
        if (this.regions.size() == 0 || newRegion[0] > regions.get(regions.size() - 1)[1]) this.regions.add(newRegion);
        else {
            int i = 0;
            int[] r = this.regions.get(0);
            while (newRegion[0] > r[0]) {
                i++;
                r = this.regions.get(i);
            }
            this.regions.add(i, newRegion);
        }
    }

    public int size() {
        return regions.size();
    }

    @Override
    public String toString() {
        StringBuilder str = new StringBuilder();
        for (int[] region : regions) {
            str.append(region[0]).append("\t").append(region[1]).append("\n");
        }
        return str.toString();
    }

    public boolean equals(RegionVector other) {
        if (this.size() != other.size()) return false;

        for (int i = 0; i < this.size(); i++) {
            int[] a = this.regions.get(i);
            int[] b = other.getRegions().get(i);
            if (a[0] != b[0] || a[1] != b[1]) return false;
        }
        return true;
    }

    // returns the inverse RegionVector
    public RegionVector getInverse() {
        ArrayList<int[]> invRegions = new ArrayList<>();

        // first inverse region starts at end point of first region
        int invStart = regions.get(0)[1] + 1;
        // iterate over all regions and add the spaces between to an array list
        for (int i = 1; i < regions.size(); i++) {
            int[] exon = regions.get(i);
            invRegions.add(new int[]{invStart, exon[0]});
            invStart = exon[1] + 1;
        }

        return new RegionVector(invRegions);
    }

    // returns lngth of parts between the regions
    public int getInverseLength() {
        int length = 0;
        int temp = regions.get(0)[1];
        for (int i = 1; i < regions.size(); i++) {
            length += (regions.get(i)[0] - temp);
            temp = regions.get(i)[1];
        }
        return length;
    }

    // check if a specific region is in the vector (binary search)
    // if one parameter is null it is ignored
    public boolean searchRegion(Integer start, Integer end) {
        if(this.regions.size() == 0) return false;
        return binarySearch(0, this.regions.size(), start, end);
    }

    private boolean binarySearch(int a, int b, Integer start, Integer end) {
        if (a >= b - 1) {
            boolean startMatches = (start == null || start == this.regions.get(a)[0]);
            boolean endMatches = (end == null || end == this.regions.get(a)[1]);
            return startMatches && endMatches;
        }
        int middle = (a + b) / 2;
        boolean startMatches = (start == null || start == this.regions.get(middle)[0]);
        boolean endMatches = (end == null || end == this.regions.get(middle)[1]);
        if (startMatches && endMatches) return true;
        else {
            if ((start != null && start < this.regions.get(middle)[0]) || (end != null && end < this.regions.get(middle)[1]))
                return binarySearch(a, middle, start, end);
            else return binarySearch(middle, b, start, end);
        }
    }

    // find all the regions in a specific bigger region (subregions have to cover exact start and end)
    public RegionVector findRegionsIn(int[] biggerRegion) {
        ArrayList<int[]> subregions = new ArrayList<>();
        int start = biggerRegion[0];
        int end = biggerRegion[1];
        boolean add = false;
        for (int[] reg : this.regions) {
            if (reg[0] == start) add = true;
            if (add && reg[1] <= end) subregions.add(reg);
            if (reg[1] >= end) add = false;
        }
        // return null if end doesn't fit
        if(subregions.size() == 0 || subregions.get(subregions.size()-1)[1] != end) return null;
        return new RegionVector(subregions);
    }

    // return a string that contains all the regions of the vector in the format start1:end1|start2:end2|...
    public String regionsAsPrintString() {
        if (this.regions.size() == 0) return "";
        StringBuilder s = new StringBuilder();
        for (int[] reg : this.regions) {
            s.append("|").append(reg[0]).append(":").append(reg[1]);
        }
        return s.toString().substring(1);
    }
}
