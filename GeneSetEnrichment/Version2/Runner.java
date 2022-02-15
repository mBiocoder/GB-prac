package com.company;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;

import java.util.zip.GZIPInputStream;

public class Runner {
    // map the gene names to the GoEntry IDs the belong to
    private static HashMap<String, HashSet<String>> geneToGos;
    private static int numGenes;
    private static int diffGenes;

    public static void main(String[] args) {
        Parser parser = new Parser("-obo", "-root", "-mapping", "-mappingtype", "-enrich", "-o", "-minsize", "-maxsize", "-overlapout");
        parser.parseArguments(args);

        String oboFile = parser.getValue("-obo");
        String namespace = parser.getValue("-root");
        String mappingFile = parser.getValue("-mapping");
        String mappingType = parser.getValue("-mappingtype");
        String diffExpFile = parser.getValue("-enrich");
        String outputFile = parser.getValue("-o");
        String min = parser.getValue("-minsize");
        String max = parser.getValue("-maxsize");
        String overlapOutFile = parser.getValue("-overlapout");

        boolean problem = false;
        // check input values
        if (oboFile == null || mappingFile == null || mappingType == null || diffExpFile == null || outputFile == null)
            problem = true;
        if (namespace != null && !namespace.equals("biological_process") && !namespace.equals("cellular_component") &&
                !namespace.equals("molecular_function")) {
            problem = true;
            System.out.println("Please use one of the three GO namespaces: biological_process, cellular_component, molecular_function");
        }
        if (mappingType != null && !mappingType.equals("ensembl") && !mappingType.equals("go")) {
            problem = true;
            System.out.println("Please use [ensembl|go] as the mapping type!");
        }
        int minsize = 0;
        int maxsize = 1;
        if (min == null || max == null) problem = true;
        else {
            try {
                minsize = Integer.parseInt(min);
                maxsize = Integer.parseInt(max);
                if (maxsize < minsize) {
                    problem = true;
                    System.out.println("Please provide minsize < maxsize!");
                }
            } catch (NumberFormatException e) {
                problem = true;
                System.out.println("minsize and maxsize must be integers!");
            }
        }

        // if one of the parameters was missing or not correct, print the usage:
        if (problem) {
            System.out.println("Usage:");
            System.out.println("-obo <obo_file>");
            System.out.println("-root <GO_namespace>\t\t\tmust be in [biological_process|cellular_component|molecular_function]");
            System.out.println("-mapping <gene2go_mapping_file>");
            System.out.println("-mappingtype [ensembl|go]\t\tspecifies format of mapping file");
            System.out.println("-enrich <diffexp_file>");
            System.out.println("-o <output_tsv>");
            System.out.println("-minsize <int>\t\t\t\t\tdefine which GO entries are considered");
            System.out.println("-maxsize <int>\t\t\t\t\tonly consider entries with minsize <= size <= maxzize");
            System.out.println("[-overlapout <overlap_out_tsv>]");
            System.exit(0);
        }


        // initialize static variables:
        geneToGos = new HashMap<>();
        numGenes = 0;
        diffGenes = 0;


        //#### read oboFile and build GoEntry objects ##################################################################
        HashMap<String, GoEntry> goToEntry = new HashMap<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(oboFile));
            String line;
            ArrayList<String> block = new ArrayList<>();
            // iterate over the file and gather the blocks
            while ((line = br.readLine()) != null) {
                if (line.startsWith("[Term]")) {
                    block.add(line);
                } else if (line.equals("") && block.size() > 0) {
                    GoEntry goe = buildGoEntry(block, namespace);
                    if (goe != null) goToEntry.put(goe.getId(), goe);
                    block.clear();
                } else if (block.size() > 0) block.add(line);
            }

            br.close();
        } catch (FileNotFoundException e) {
            System.out.println("The obo-file " + oboFile + " was not found!");
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }


        //#### read the mapping file and add gene names to GoEntry objects #############################################
        // we have two different types of mapping files => one method for each type...
        if (mappingType.equals("ensembl")) parseMappingEnsembl(mappingFile, goToEntry);
        else parseMappingGo(mappingFile, goToEntry);


        //#### read and process enrich file ############################################################################
        ArrayList<Double> allFcs = new ArrayList<>();
        HashSet<String> simulatedTruths = new HashSet<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(diffExpFile));
            String line;

            while ((line = br.readLine()) != null) {
                line = line.replace("\n", "");
                if (line.startsWith("#")) {
                    // ground truth: truly enriched GO classes
                    line = line.substring(1);
                    if (goToEntry.containsKey(line)) {
                        simulatedTruths.add(line);
                        goToEntry.get(line).setSignifTrue();
                    }
                } else if (!line.startsWith("id")) {
                    // all the other lines except from the header line
                    String[] parts = line.split("\t");
                    double fc = 0.0;
                    try {
                        fc = Double.parseDouble(parts[1]);
                    } catch (NumberFormatException e) {
                        System.out.println("WARNING: The enrichment file " + diffExpFile + " contains the log fold change " + parts[1] + " which seems to be no double! Replaced it with 0.0");
                    }
                    // parseBoolean gives us true if the String is some variation of "true", else: false
                    boolean signif = Boolean.parseBoolean(parts[2]);

                    // if the gene is significant => count it for the corresponding GoEntries
                    if (signif && geneToGos.containsKey(parts[0])) {
                        diffGenes++;
                        allFcs.add(fc);
                        for (String go : geneToGos.get(parts[0])) {
                            goToEntry.get(go).incrementNumSignifGenes();
                            goToEntry.get(go).addMeasured(fc);
                        }
                    } else if (geneToGos.containsKey(parts[0])) {
                        allFcs.add(fc);
                        for (String go : geneToGos.get(parts[0])) {
                            goToEntry.get(go).addMeasured(fc);
                        }
                    }
                }
            }

            br.close();
        } catch (FileNotFoundException e) {
            System.out.println("The enrich file " + diffExpFile + " was not found!");
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }
        // for all simulated truth GOEntries: pre-compute all entries on the path to the node (and save with distance)
        for (String go : simulatedTruths) {
            goToEntry.get(go).preComputeAncestors(goToEntry);
        }
        Collections.sort(allFcs);

        //#### generate pvals/do statistical tests #####################################################################
        GoEntry[] allGoEntries = new GoEntry[goToEntry.keySet().size()];
        int k = 0;
        for (String go : goToEntry.keySet()) {
            if (goToEntry.get(go).getSize() >= minsize && goToEntry.get(go).getSize() <= maxsize) {
                allGoEntries[k] = goToEntry.get(go);
                k++;

                goToEntry.get(go).fisherTest(numGenes, diffGenes);
                goToEntry.get(go).hgTest(numGenes, diffGenes);
                goToEntry.get(go).ksTest(allFcs);
            }
        }
        allGoEntries = Arrays.copyOfRange(allGoEntries, 0, k);

        // Perform the multiple testing corrections
        bhHg(allGoEntries);
        bhFischer(allGoEntries);
        bhKs(allGoEntries);


        //#### write output ############################################################################################
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
            bw.write("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true\n");
            // iterate over all the GoEntries and create output for all with the correct size
            for (String go : goToEntry.keySet()) {
                GoEntry goen = goToEntry.get(go);
                int size = goen.getSize();
                if (size >= minsize && size <= maxsize) {
                    String pathInfo = "";
                    if (!goen.getSignifTruth() && simulatedTruths.size() > 0) {
                        // pathInfo is set
                        Integer currentMin = null;
                        String shortestPath = null;
                        for (String truth : simulatedTruths) {
                            String s = getShortestPathTo(goen, goToEntry.get(truth), goToEntry, currentMin);
                            if (s != null) {
                                shortestPath = s;
                                currentMin = shortestPath.split("\\|").length;
                            }
                        }
                        pathInfo = goPathToNamePath(shortestPath, goToEntry);
                    }

                    String line = go + "\t" + goen.getName() + "\t" + size + "\t" + goen.getSignifTruth() + "\t" +
                            goen.getNumSignificantGenes() + "\t" + goen.getHgPval() + "\t" + goen.getHgPadj() + "\t" +
                            goen.getFejPval() + "\t" + goen.getFejPadj() + "\t" + goen.getKsStat() + "\t" +
                            goen.getKsPval() + "\t" + goen.getKsPadj() + "\t" + pathInfo;
                    line += "\n";
                    bw.write(line);
                }
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        //#### Option for overlap information ##########################################################################
        if (overlapOutFile != null) {
            // find all the pairs that have an overlap
            HashSet<String> overlapPairs = new HashSet<>();
            for (String gene : geneToGos.keySet()) {
                String[] gos = new String[geneToGos.get(gene).size()];
                int x = 0;
                for (String g : geneToGos.get(gene)) {
                    gos[x] = g;
                    x++;
                }

                // build all pairs and make sure they are not already in the set of pairs...
                for (int i = 0; i < gos.length; i++) {
                    if (goToEntry.get(gos[i]).getSize() >= minsize && goToEntry.get(gos[i]).getSize() <= maxsize) {
                        for (int j = i + 1; j < gos.length; j++) {
                            if (goToEntry.get(gos[j]).getSize() >= minsize && goToEntry.get(gos[j]).getSize() <= maxsize) {
                                String one = gos[i] + "-" + gos[j];
                                String two = gos[j] + "-" + gos[i];
                                if (!overlapPairs.contains(one) && !overlapPairs.contains(two)) {
                                    overlapPairs.add(one);
                                }
                            }
                        }
                    }
                }
            }

            // compute the overlap between every pair...
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(overlapOutFile));
                bw.write("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent\n");

                for (String pair : overlapPairs) {
                    String go1 = pair.split("-")[0];
                    String go2 = pair.split("-")[1];
                    String s = shortestPathBetween(goToEntry.get(go1), goToEntry.get(go2), goToEntry);
                    boolean isRelative = Boolean.parseBoolean(s.split(":")[1]);
                    int minPathLength = Integer.parseInt(s.split(":")[0]);
                    int numOverlapping = goToEntry.get(go1).numOverlappingGenes(goToEntry.get(go2));
                    double frac = numOverlapping / (double) goToEntry.get(go1).getSize();
                    double frac2 = numOverlapping / (double) goToEntry.get(go2).getSize();
                    if (frac2 > frac) frac = frac2;

                    bw.write(go1 + "\t" + go2 + "\t" + isRelative + "\t" + minPathLength + "\t" + numOverlapping +
                            "\t" + doubleToPercent(frac) + "\n");
                }

                bw.close();
            } catch (FileNotFoundException e) {
                System.out.println("The overlap output file " + overlapOutFile + " was not found!");
                System.exit(0);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    // build GoEntry object from a GO entry
    public static GoEntry buildGoEntry(ArrayList<String> entryLines, String root) {
        // program only deals with one DAG (specified by root) => ignore entry when it does not fit
        String id = null;
        String name = null;
        String namespace = null;
        HashSet<String> parents = new HashSet<>();

        // go through the entry and get all necessary information
        for (String line : entryLines) {
            if (line.startsWith("id:")) {
                id = line.substring(4);
                id = id.replace("\n", "");
            } else if (line.startsWith("name:")) {
                name = line.substring(6);
                name = name.replace("\n", "");
                // Ignore obsolete terms!!!
                if (name.contains("obsolete")) return null;
            } else if (line.startsWith("namespace:")) {
                namespace = line.substring(11);
                namespace = namespace.replace("\n", "");
                if (!namespace.equals(root)) return null;
            } else if (line.startsWith("is_a:")) {
                String parentId = line.substring(6);
                parentId = parentId.substring(0, 10);
                parents.add(parentId);
            }
        }

        return new GoEntry(id, name, parents);
    }

    public static void parseMappingEnsembl(String filePath, HashMap<String, GoEntry> goToEntry) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(filePath));
            String line;

            while ((line = br.readLine()) != null) {
                // header line is not considered...
                if (line.startsWith("ENSG")) {
                    HashSet<String> geneGos = new HashSet<>();
                    line = line.replace("\n", "");
                    String[] parts = line.split("\t");
                    // only handle line if the gene has a name
                    if (!parts[1].equals("")) {
                        String[] gos = parts[2].split("\\|");
                        for (String go : gos) {
                            // check because we only handle entries of one namespace
                            if (goToEntry.containsKey(go)) {
                                goToEntry.get(go).addGene(parts[1], goToEntry, geneGos);
                            }
                        }
                        // now, geneGos should be filled with all the GoEntry IDs where the gene belongs to => save in HashMap
                        if (geneGos.size() > 0) {
                            // only do this if the gene was mapped to GOs that are in our desired namespace
                            if (!geneToGos.containsKey(parts[1])) {
                                numGenes++;
                                geneToGos.put(parts[1], geneGos);
                            } else {
                                geneToGos.get(parts[1]).addAll(geneGos);
                            }
                        }
                    }
                }
            }

            br.close();
        } catch (FileNotFoundException e) {
            System.out.println("The mapping file " + filePath + " was not found!");
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void parseMappingGo(String filePath, HashMap<String, GoEntry> goToEntry) {
        try {
            java.util.zip.GZIPInputStream reader = new GZIPInputStream(new FileInputStream(filePath));
            InputStreamReader isr = new InputStreamReader(reader);
            BufferedReader in = new BufferedReader(isr);

            String line;
            while ((line = in.readLine()) != null) {
                // comment lines start with "!"
                if (!line.startsWith("!")) {
                    line = line.replace("\n", "");
                    String[] parts = line.split("\t");

                    // only use the associations without association qualifier modifier:
                    if (parts[3].equals("")) {
                        HashSet<String> geneGos = new HashSet<>();
                        String geneName = parts[2];
                        String goId = parts[4];
                        if (goToEntry.containsKey(goId)) {
                            goToEntry.get(goId).addGene(geneName, goToEntry, geneGos);
                            // add all the gos, where the gene belongs to to the HashMap
                            if (!geneToGos.containsKey(geneName)) {
                                geneToGos.put(geneName, geneGos);
                                numGenes++;
                            } else geneToGos.get(geneName).addAll(geneGos);
                        }
                    }
                }
            }

            in.close();
            isr.close();
            reader.close();

        } catch (FileNotFoundException e) {
            System.out.println("The mapping file " + filePath + " was not founnd!");
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static HashMap<String, Double> BH(HashMap<String, Double> pvals) {
        String[] gos = new String[pvals.keySet().size()];
        double[] p = new double[gos.length];
        double[] padjusted = new double[gos.length];
        int j = 0;
        for (String go : pvals.keySet()) {
            gos[j] = go;
            p[j] = pvals.get(go);
            j++;
        }
        NaturalRanking nr = new NaturalRanking(TiesStrategy.SEQUENTIAL);
        double[] ranks = nr.rank(p);

        // compute adjusted pvalues
        // they are saved in the order of the ranking!!!
        HashMap<String, Double> padj = new HashMap<>();
        for (int i = 0; i < gos.length; i++) {
            int index = (int) ranks[i] - 1;
            if (index == padjusted.length - 1) padjusted[index] = p[i];
            else padjusted[index] = ((p[i] * p.length) / ranks[i]);
        }

        // go through adjusted pvalues in the order from before and make sure, they are not decreasing
        // actually: take for padj the minimum of itself or all padj that come after it
        for (int i = gos.length - 1; i >= 0; i--) {
            // make sure that padj is not bigger than value after that
            if (i + 1 < gos.length) {
                if (padjusted[i + 1] < padjusted[i]) padjusted[i] = padjusted[i + 1];
            }

            // if padj > 1, set it to 1.0 bc bigger doesn't make sense
            if (padjusted[i] > 1.0) padjusted[i] = 1.0;
        }

        // put the result in the HashMap
        for (int i = 0; i < gos.length; i++) {
            padj.put(gos[i], padjusted[(int) ranks[i] - 1]);
        }

        return padj;
    }

    public static double doubleToPercent(double num) {
        double res = num * 10000;
        res = Math.round(res);
        return res / 100.0;
    }

    public static String shortestPathBetween(GoEntry ge1, GoEntry ge2, HashMap<String, GoEntry> goToEntry) {
        HashMap<String, Integer> ancestorsOfOne = new HashMap<>();
        ancestorsOfOne.put(ge1.getId(), 0);
        boolean ancestor = false;
        Integer foundDist = null;
        int dist = 1;
        HashSet<String> parents = ge1.getParents();
        HashSet<String> nextGen;
        // follow ancestors of ge1 to the root
        while (parents.size() > 0) {
            if (parents.contains(ge2.getId())) {
                if (foundDist == null) foundDist = dist;
                ancestor = true;
                //return dist + ":true";
            }
            nextGen = new HashSet<>();
            for (String p : parents) {
                if (!ancestorsOfOne.containsKey(p)) ancestorsOfOne.put(p, dist);
                nextGen.addAll(goToEntry.get(p).getParents());
            }
            parents = nextGen;
            dist++;
        }
        // follow ancestors of ge2 to the root
        dist = 1;
        parents = ge2.getParents();
        while (parents.size() > 0) {
            if (parents.contains(ge1.getId())) {
                ancestor = true;
                if (foundDist == null) foundDist = dist;
                else foundDist = Math.min(dist, foundDist);
                //return dist + ":true";
            }
            nextGen = new HashSet<>();
            // check if there is overlap with the ancestors of ge1
            for (String p : parents) {
                if (ancestorsOfOne.containsKey(p)) {
                    if (foundDist == null) foundDist = dist + ancestorsOfOne.get(p);
                    else foundDist = Math.min(foundDist, dist + ancestorsOfOne.get(p));
                }
                nextGen.addAll(goToEntry.get(p).getParents());
            }
            //if (foundDist != null) return foundDist + ":false";
            parents = nextGen;
            dist++;
        }
        return foundDist + ":" + ancestor;
    }

    public static String getShortestPathTo(GoEntry goen, GoEntry truth, HashMap<String, GoEntry> goToEntry, Integer currentMin) {
        HashMap<String, Path> trueAncestors = truth.getAncestorMap();
        String goenId = goen.getId();
        String path = null;
        Integer foundLength = null;
        if (trueAncestors.containsKey(goenId)) {
            foundLength = trueAncestors.get(goenId).getLength() + 1;
            path = goenId + "*";
            if (trueAncestors.get(goenId).getPath().length() > 0) path += "|" + trueAncestors.get(goenId).getPath();
        }

        // go through all ancestors of the GoEntry and check if there is an overlap with the keySet of the ancestors of the true GoEntry
        HashSet<String> parents = goen.getParents();
        HashMap<String, String> parentToChild = new HashMap<>();
        for (String p : parents) {
            parentToChild.put(p, goenId);
        }

        int dist = 1;
        HashMap<String, Path> searchPaths = new HashMap<>();
        Path seedPath = new Path();
        searchPaths.put(goenId, seedPath);
        HashSet<String> nextGen;

        while (parents.size() > 0) {
            nextGen = new HashSet<>();
            // check if a new path to the true was found
            for (String p : parents) {
                if (trueAncestors.containsKey(p)) {
                    int len = trueAncestors.get(p).getLength() + dist + 1;
                    if (foundLength == null || len < foundLength) {
                        foundLength = len;
                        path = searchPaths.get(parentToChild.get(p)).getPath() + "|" + parentToChild.get(p) + "|" + p +
                                "*|" + trueAncestors.get(p).getPath();
                    }
                }
                // expand next generation
                if (!searchPaths.containsKey(p)) {
                    Path parentPath = searchPaths.get(parentToChild.get(p)).copy();
                    parentPath.addNodeBack(parentToChild.get(p));
                    searchPaths.put(p, parentPath);
                }
                HashSet<String> grandparents = goToEntry.get(p).getParents();
                nextGen.addAll(grandparents);
                for (String gp : grandparents) {
                    if (!parentToChild.containsKey(gp)) parentToChild.put(gp, p);
                }
            }
            // update parents (go to next generation)
            parents = nextGen;
            dist++;
        }

        if (foundLength != null && (currentMin == null || foundLength < currentMin)) return path;
        else return null;
    }

    private static String goPathToNamePath(String goPath, HashMap<String, GoEntry> goToEntry) {
        String[] nodes = goPath.split("\\|");
        StringBuilder res = new StringBuilder();
        for (String node : nodes) {
            if (!node.equals("")) {
                boolean star = false;
                if (node.endsWith("*")) {
                    star = true;
                    node = node.substring(0, node.length() - 1);
                }
                if (res.length() != 0) res.append("|");
                res.append(goToEntry.get(node).getName());
                if (star) res.append("*");
            }
        }
        return res.toString();
    }

    private static void bhFischer(GoEntry[] goEntries) {
        int m = goEntries.length;
        // sort the goEntries from small to big
        Arrays.sort(goEntries, Comparator.comparingDouble(GoEntry::getFejPval));
        for (int i = m - 1; i >= 0; i--) {
            if (i == m - 1) {
                goEntries[i].setfejPadj(goEntries[i].getFejPval());
            } else {
                // adjusted p-val: (#tests/rank)*unadjusted
                double unadjustedPvalue = goEntries[i].getFejPval();
                int divideByM = i + 1;
                double right = goEntries[i + 1].getFejPadj();
                double adj = (m / (double) divideByM) * unadjustedPvalue;
                // make sure that the adjusted values are still descending
                goEntries[i].setfejPadj(Math.min(right, adj));
            }
        }
    }

    private static void bhHg(GoEntry[] goEntries) {
        int m = goEntries.length;
        // sort the goEntries from small to big
        Arrays.sort(goEntries, Comparator.comparingDouble(GoEntry::getHgPval));
        for (int i = m - 1; i >= 0; i--) {
            if (i == m - 1) {
                goEntries[i].sethgPadj(goEntries[i].getHgPval());
            } else {
                // adjusted p-val: (#tests/rank)*unadjusted
                double unadjustedPvalue = goEntries[i].getHgPval();
                int divideByM = i + 1;
                double right = goEntries[i + 1].getHgPadj();
                double adj = (m / (double) divideByM) * unadjustedPvalue;
                // make sure that the adjusted values are still descending
                goEntries[i].sethgPadj(Math.min(right, adj));
            }
        }
    }

    private static void bhKs(GoEntry[] goEntries) {
        int m = goEntries.length;
        // sort the goEntries from small to big
        Arrays.sort(goEntries, Comparator.comparingDouble(GoEntry::getKsPval));
        for (int i = m - 1; i >= 0; i--) {
            if (i == m - 1) {
                goEntries[i].setksPadj(goEntries[i].getKsPval());
            } else {
                // adjusted p-val: (#tests/rank)*unadjusted
                double unadjustedPvalue = goEntries[i].getKsPval();
                int divideByM = i + 1;
                double right = goEntries[i + 1].getKsPadj();
                double adj = (m / (double) divideByM) * unadjustedPvalue;
                // make sure that the adjusted values are still descending
                goEntries[i].setksPadj(Math.min(right, adj));
            }
        }
    }
}
