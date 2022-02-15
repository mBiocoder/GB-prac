package com.company;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

public class GoEntry {
    private String id;
    private String name;
    private HashSet<String> parents;
    private HashSet<String> genes;
    private boolean signifTruth;
    private int numSignificantGenes;
    private ArrayList<Double> measuredFCs;
    private HashMap<String, Path> nodeToPath;

    private double hgPval;
    private double hgPadj;
    private double fejPval;
    private double fejPadj;
    private double ksPval;
    private double ksStat;
    private double ksPadj;

    public GoEntry(String id, String name, HashSet<String> parents) {
        this.id = id;
        this.name = name;
        this.parents = parents;
        this.genes = new HashSet<>();
        this.signifTruth = false;
        this.numSignificantGenes = 0;
        measuredFCs = new ArrayList<>();
    }

    public String getId() {
        return this.id;
    }

    public String getName() {
        return this.name;
    }

    public HashSet<String> getParents() {
        return this.parents;
    }

    // call it when the GoEntry is in the ground truth
    public void setSignifTrue() {
        this.signifTruth = true;
    }

    public boolean getSignifTruth() {
        return this.signifTruth;
    }

    public void addGene(String geneId, HashMap<String, GoEntry> idToEntry, HashSet<String> geneGos) {
        // add gene to the ArrayList
        this.genes.add(geneId);
        // this HashSet is needed in runner (all the GoEntry IDs where the gene belongs to)
        geneGos.add(this.id);
        // propagate: add gene to all parents (and they add it to their parents and so on...)
        for (String parent : parents) {
            GoEntry parentEntry = idToEntry.get(parent);
            if (!parentEntry.containsGene(geneId)) parentEntry.addGene(geneId, idToEntry, geneGos);
        }
    }

    public boolean containsGene(String geneId) {
        return this.genes.contains(geneId);
    }

    // definition of size: # assigned genes (whether measured or not!)
    public int getSize() {
        return this.genes.size();
    }

    public void incrementNumSignifGenes() {
        this.numSignificantGenes += 1;
    }

    public int getNumSignificantGenes() {
        return this.numSignificantGenes;
    }

    public void addMeasured(double fc) {
        this.measuredFCs.add(fc);
    }

    // perform hypergeometric test on the GOEntry
    public void hgTest(int numGenes, int numDiffGenes) {
        int k = numSignificantGenes;
        int n = this.getSize();

        // parameters: population size (all genes), # successes in population (all DE genes), sample size
        HypergeometricDistribution hg = new HypergeometricDistribution(numGenes, numDiffGenes, n);
        this.hgPval = hg.upperCumulativeProbability(k);
    }

    public void fisherTest(int numGenes, int numDiffGenes) {
        int a = this.numSignificantGenes;
        int c = this.getSize() - a;
        // jacknife: use all values -1 (but make sure the value is not zero!)
        int populationSize = numGenes - 1;
        if (populationSize < 0) {
            populationSize++;
        }
        int succInPopulation = numDiffGenes - 1;
        if (succInPopulation < 0) {
            succInPopulation++;
        }
        int sampleSize = a + c - 1;
        if (sampleSize < 0) {
            sampleSize++;
        }

        // parameters: population size (all genes), # successes in population (all DE genes), sample size
        HypergeometricDistribution hg = new HypergeometricDistribution(populationSize, succInPopulation, sampleSize);
        this.fejPval = hg.upperCumulativeProbability(a - 1);
    }

    // perform Kolmogorov-Smirnov test
    public void ksTest(ArrayList<Double> allFc) {
        KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();

        // background: all fold-changes that are not in this set
        double[] in_set_distrib = this.measuredFCs.stream().mapToDouble(Double::doubleValue).toArray();
        ArrayList<Double> bg_set = new ArrayList<>(allFc);
        this.measuredFCs.forEach(bg_set::remove);
        double[] bg_distrib = bg_set.stream().mapToDouble(Double::doubleValue).toArray();

        if (in_set_distrib.length >= 2 && bg_distrib.length >= 2) {
            this.ksPval = ks.kolmogorovSmirnovTest(in_set_distrib, bg_distrib);
            this.ksStat = ks.kolmogorovSmirnovStatistic(in_set_distrib, bg_distrib);
        } else {
            this.ksPval = 1.0;
            this.ksStat = 0.0;
        }
    }

    public int numOverlappingGenes(GoEntry other) {
        HashSet<String> overlap = new HashSet<>(this.genes);
        overlap.retainAll(other.getGenes());
        return overlap.size();
    }

    private HashSet<String> getGenes() {
        return this.genes;
    }

    public void preComputeAncestors(HashMap<String, GoEntry> goToEntry) {
        nodeToPath = new HashMap<>();
        this.ancestorHelper(0, goToEntry, this.nodeToPath, null);
    }

    private void ancestorHelper(int dist, HashMap<String, GoEntry> goToEntry, HashMap<String, Path> result, String child) {
        if (child == null) {
            result.put(this.id, new Path());
        } else if (!result.containsKey(this.id)) {
            Path p = result.get(child).copy();
            p.addNodeFront(child);
            result.put(this.id, p);
        } else {
            // check if this path is shorter than the one in the HashMap
            Path p = result.get(child).copy();
            p.addNodeFront(child);
            if (p.getLength() < result.get(this.id).getLength()) result.put(this.id, p);
        }


        // call function for every parent
        for (String parent : this.parents) {
            goToEntry.get(parent).ancestorHelper(dist + 1, goToEntry, result, this.id);
        }
    }

    public HashMap<String, Path> getAncestorMap() {
        return this.nodeToPath;
    }

    public void sethgPadj(double padj) {
        this.hgPadj = padj;
    }

    public void setfejPadj(double padj) {
        this.fejPadj = padj;
    }

    public void setksPadj(double padj) {
        this.ksPadj = padj;
    }

    public double getHgPval() {
        return hgPval;
    }

    public double getHgPadj() {
        return hgPadj;
    }

    public double getFejPval() {
        return fejPval;
    }

    public double getFejPadj() {
        return fejPadj;
    }

    public double getKsPval() {
        return ksPval;
    }

    public double getKsStat() {
        return ksStat;
    }

    public double getKsPadj() {
        return ksPadj;
    }
}
