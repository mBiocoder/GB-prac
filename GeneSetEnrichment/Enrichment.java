public class Enrichment {
    public String geneId;
    public Double fc;
    public boolean signif;

    public Enrichment(String geneId, Double fc, boolean signif) {
        this.geneId = geneId;
        this.fc = fc;
        this.signif = signif;
    }

}
