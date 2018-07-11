package MZIDparser;

import com.compomics.util.experiment.biology.PTM;
import com.compomics.util.experiment.biology.PTMFactory;
import com.compomics.util.experiment.biology.Peptide;
import com.compomics.util.experiment.identification.matches.ModificationMatch;

import java.util.ArrayList;

/**
 * Created by dengyamei on 10/07/2018.
 */
public class JPSM {
    //member variable

    //spectrum identifier, different from SpectrumIndex
    private String SpectrumID;
    //spectrum index in MGF
    private int SpectrumIndex;
    //peptide sequence
    private String pepSeq;
    //peptide modification, if none, this object equals to null
    private ArrayList<JModification> modificationList = new ArrayList<JModification>();

    private int rank;
    private int charge;
    private double calculatedMassToCharge;
    private double experimentalMassToCharge;
    private double evalue;
    private ArrayList<String> proteins = new ArrayList<String>();
    private boolean isDecoy=false;
    private double rt=-1.0;
    private double massErrorMZ; //precursor feature
    private double IsotopeError=0.0;
    //whether or not to add flank aa into peptide sequence
    public static boolean addFlankAA = false;

    //flank aa
    private String preAA = "";
    private String postAA = "";

    /**
     * record non-rank1 results
     */
    private ArrayList<JPSM> jPsmRank = new ArrayList<JPSM>();

    // engine-specific score
    //tandem
    private double rawHyperScore; // tandem-specific score
    //myrimatch
    private double rawMVH;
    private double rawMzFidelity;
    private int numOfMatchIonsForMyrimatch;
    private int numOfUnmatchIonsForMyrimatch;
    private double rawXCorr;

    /**
     * 1=MS-GF+
     * 2=X!Tandem
     * 3=Myrimatch
     */
    public int searchEngineClass;

    //feature
    private ArrayList< Double> featureValueList = new ArrayList<Double>();
    public static ArrayList<String> featureNameList = new ArrayList<String>();
    //the following parameters are used to calculate deltaScore
    /**
     * -log Evalue
     */
    private double lnEvalue;
    /**
     * -log SpecEValue
     */
    private double lnSpecEValue;
    private double rawScore;
    private double deNovoScore;
    private double qvalue;

    public void setLnEvalue(double e){
        this.lnEvalue = e;
    }
    public double getLnEvalue(){
        return this.lnEvalue;
    }
    public void setLnSpecEValue(double e){
        this.lnSpecEValue = e;
    }
    public double getlnSpecEValue(){
        return this.lnSpecEValue;
    }

    /*constructor*/
    public JPSM(){

    }

    /**
     * get the ArrayList<ModificationMatch> object for Peptide Class
     * @return
     */
    public ArrayList<ModificationMatch> getModificationMatch_will_remove() {
        ArrayList<ModificationMatch> modificationMatches = new ArrayList<ModificationMatch>();
        PTMFactory ptmFactory = PTMFactory.getInstance();
        if(modificationList.size()>=1){
            for(int i=0;i<modificationList.size();i++){

                JModification jModification = modificationList.get(i);
                String ptmName = "CMSPTM:"+jModification.getResidue()+":"+jModification.getModMassDelta();
                if(!ptmFactory.containsPTM(ptmName)){
                    ArrayList<String> residuesArray = new ArrayList<String>();
                    residuesArray.add(jModification.getResidue());
                    @SuppressWarnings("deprecation")
                    PTM ptm = new PTM(1,ptmName,jModification.getModMassDelta(),residuesArray);
                    ptmFactory.addUserPTM(ptm);
                }

                ModificationMatch mm = new ModificationMatch(ptmName, true, jModification.getModLocation());

                modificationMatches.add(mm);
            }

        }
        return(modificationMatches);
    }

    public ArrayList<ModificationMatch> getModificationMatch() {
        ArrayList<ModificationMatch> modificationMatches = new ArrayList<ModificationMatch>();
        PTMFactory ptmFactory = PTMFactory.getInstance();
        if (modificationList.size() >= 1) {
            for (int i = 0; i < modificationList.size(); i++) {
                JModification jModification = modificationList.get(i);
                double modMassDelta = Double.valueOf(String.format("%.4f", jModification.getModMassDelta()));
                String ptmName = "CMSPTM:" + jModification.getResidue() + ":" + modMassDelta;
                if (!ptmFactory.containsPTM(ptmName)) {
                    ArrayList<String> residuesArray = new ArrayList<String>();
                    residuesArray.add(jModification.getResidue());
                    PTM ptm = null;
                    if (jModification.getResidue().equalsIgnoreCase("N-term")) {
                        ptm = new PTM(PTM.MODNP, ptmName, modMassDelta, new ArrayList<String>());
                    } else if (jModification.getResidue().equalsIgnoreCase("C-term")) {
                        ptm = new PTM(PTM.MODCP, ptmName, modMassDelta, new ArrayList<String>());
                    } else {
                        ptm = new PTM(PTM.MODAA, ptmName, modMassDelta, residuesArray);
                    }
                    //System.out.println(ptm.getType() + "," + ptmName + "," + modMassDelta);
                    ptmFactory.addUserPTM(ptm);
                }

                ModificationMatch mm = new ModificationMatch(ptmName, true, jModification.getModLocation());
                //System.out.println(+","+ptmName+","+jModification.getModLocation());

                modificationMatches.add(mm);
            }

        }
        return (modificationMatches);
    }

    /**
     * get peptide mass
     */
    public double getMass() throws InterruptedException {
        Peptide peptideTmp  = new Peptide( getPepSeq(),getModificationMatch());
        return peptideTmp.getMass();
    }


    public String getSpectrumID() {
        return SpectrumID;
    }


    public void setSpectrumID(String spectrumID) {
        SpectrumID = spectrumID;
    }


    public int getSpectrumIndex() {
        return SpectrumIndex;
    }


    public void setSpectrumIndex(int spectrumIndex) {
        SpectrumIndex = spectrumIndex;
    }


    public String getPepSeq() {
        return pepSeq;
    }


    public void setPepSeq(String pepSeq) {
        this.pepSeq = pepSeq;
    }


    public ArrayList<JModification> getModificationList() {
        return modificationList;
    }


    public void addModificationList(JModification modification) {
        this.modificationList.add(modification);
    }


    public int getRank() {
        return rank;
    }


    public void setRank(int rank) {
        this.rank = rank;
    }


    public int getCharge() {
        return charge;
    }


    public void setCharge(int charge) {
        this.charge = charge;
    }


    public double getCalculatedMassToCharge() {
        return calculatedMassToCharge;
    }


    public void setCalculatedMassToCharge(double calculatedMassToCharge) {
        this.calculatedMassToCharge = calculatedMassToCharge;
    }


    public double getExperimentalMassToCharge() {
        return experimentalMassToCharge;
    }


    public void setExperimentalMassToCharge(double experimentalMassToCharge) {
        this.experimentalMassToCharge = experimentalMassToCharge;
    }


    public double getEvalue() {
        return evalue;
    }


    public void setEvalue(double evalue) {
        this.evalue = evalue;
    }

    public ArrayList<String> getProteins() {
        return proteins;
    }

    public void addProteins(String acc) {
        this.proteins.add(acc);
    }


    public double getRawHyperScore() {
        return this.rawHyperScore;
    }

    public void setRawHyperScore(double hyperScore) {
        this.rawHyperScore = hyperScore;
    }


    public double getRawMVH() {
        return rawMVH;
    }

    public void setRawMVH(double rawMVH) {
        this.rawMVH = rawMVH;
    }

    public double getRawMzFidelity() {
        return rawMzFidelity;
    }

    public void setRawMzFidelity(double rawMzFidelity) {
        this.rawMzFidelity = rawMzFidelity;
    }

    public int getNumOfMatchIonsForMyrimatch() {
        return numOfMatchIonsForMyrimatch;
    }

    public void setNumOfMatchIonsForMyrimatch(int numOfMatchIonsForMyrimatch) {
        this.numOfMatchIonsForMyrimatch = numOfMatchIonsForMyrimatch;
    }

    public int getNumOfUnmatchIonsForMyrimatch() {
        return numOfUnmatchIonsForMyrimatch;
    }

    public void setNumOfUnmatchIonsForMyrimatch(int numOfUnmatchIonsForMyrimatch) {
        this.numOfUnmatchIonsForMyrimatch = numOfUnmatchIonsForMyrimatch;
    }

    public ArrayList<Double> getFeatureValueList() {
        return featureValueList;
    }

    public void addFeatureValueList(double va) {
        this.featureValueList.add(va);
    }

    public ArrayList<String> getFeatureNameList() {
        return featureNameList;
    }

    public static void addFeatureNameList(String name) {
        if(!featureNameList.contains(name)){
            featureNameList.add(name);
        }
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public void setDecoy(boolean isDecoy) {
        this.isDecoy = isDecoy;
    }

    public int getLabel(){
        if(this.isDecoy){
            return -1;
        }else{
            return 1;
        }
    }

    public double getRawXCorr() {
        return rawXCorr;
    }

    public void setRawXCorr(double rawXCorr) {
        this.rawXCorr = rawXCorr;
    }

    /**
     * keep non-rank1 result
     * @param jpsm
     */
    public void addPSMRank(JPSM jpsm){
        this.jPsmRank.add(jpsm);
    }

    /**
     * return non-rank1 result
     * @return
     */
    public ArrayList<JPSM> getPSMRank(){
        return this.jPsmRank;
    }
    public double getRawScore() {
        return rawScore;
    }
    public void setRawScore(double rawScore) {
        this.rawScore = rawScore;
    }
    public double getDeNovoScore() {
        return deNovoScore;
    }
    public void setDeNovoScore(double deNovoScore) {
        this.deNovoScore = deNovoScore;
    }

    public JPSM getNextRank(){
        if(this.jPsmRank.size()>=1){
            for(int i=0;i<this.jPsmRank.size();i++){
                if(i==this.jPsmRank.size()-1){
                    return this.jPsmRank.get(i);
                }else{
                    if(!this.jPsmRank.get(i).getPepSeq().equalsIgnoreCase(this.getPepSeq())){
                        return this.jPsmRank.get(i);
                    }
                }
            }
        }else{
            return null;
        }
        return null;
    }

    public double getDeltaLnEvalue(){
        if(this.getNextRank()!=null){
            double delta = this.getLnEvalue() - this.getNextRank().getLnEvalue();
            return delta;
        }else{
            return 0.0;
        }

    }

    public double getDeltaLnSpecEValue(){
        if(this.getNextRank()!=null){
            double delta = this.getlnSpecEValue() - this.getNextRank().getlnSpecEValue();
            return delta;
        }else{
            return 0.0;
        }

    }


    public double getDeltaRawScore(){
        if(this.getNextRank()!=null){
            double delta = this.getRawScore() - this.getNextRank().getRawScore();
            return delta;
        }else{
            return 0.0;
        }

    }

    public double getDeltaDeNovoScore(){
        if(this.getNextRank()!=null){
            double delta = this.getDeNovoScore() - this.getNextRank().getDeNovoScore();
            return delta;
        }else{
            return 0.0;
        }

    }


    /**
     * output feature line
     * @return
     */
    public String getOutFeatureLine(){
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(this.SpectrumIndex + ":" + this.rank+ "\t");
        stringBuilder.append(this.getLabel() + "\t");
        stringBuilder.append(this.SpectrumIndex + "\t");// Scannum
        if(this.rt!=-1){
            stringBuilder.append(this.rt + "\t");// rt
        }
        stringBuilder.append(this.getMassErrorMZ() + "\t");// delta

        for (int i = 0; i < this.getFeatureValueList().size(); i++) {
            stringBuilder.append(this.getFeatureValueList().get(i) + "\t");

        }
        if(addFlankAA){
            stringBuilder.append("X."+this.pepSeq +".X"+ "\t");
        }else{
            stringBuilder.append(this.pepSeq + "\t");
        }
        for (int j = 0; j < this.getProteins().size(); j++) {
            if (j == this.getProteins().size() - 1) {
                stringBuilder.append(this.getProteins().get(j));
            } else {
                stringBuilder.append(this.getProteins().get(j) + "\t");
            }
        }
        return stringBuilder.toString();
    }


    /**
     * output feature line, designed for ScorerTest class
     * @return
     */
    public String getOutFeatureLine4ScorerTest(){
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(this.SpectrumIndex + ":" + this.rank+ "\t");
        stringBuilder.append(this.getLabel() + "\t");
        stringBuilder.append(this.SpectrumIndex + "\t");
        if(this.rt!=-1){
            stringBuilder.append(this.rt + "\t");// rt
        }
        stringBuilder.append(this.getMassErrorMZ() + "\t");// delta

        for (int i = 0; i < this.getFeatureValueList().size(); i++) {
            stringBuilder.append(this.getFeatureValueList().get(i) + "\t");

        }
        if(addFlankAA){
            stringBuilder.append("X."+this.pepSeq +".X"+ "\t");
        }else{
            stringBuilder.append(this.pepSeq + "\t");
        }
        for (int j = 0; j < this.getProteins().size(); j++) {
            if (j == this.getProteins().size() - 1) {
                stringBuilder.append(this.getProteins().get(j));
            } else {
                stringBuilder.append(this.getProteins().get(j) + ";");
            }
        }
        return stringBuilder.toString();
    }


    /**
     * output feature line for specific index and rank
     */
    public String getOutFeatureLine(String hIndex, String hRank, long scannum){
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(hIndex + ":" + hRank + "\t");
        stringBuilder.append(this.getLabel() + "\t");

        stringBuilder.append(scannum + "\t");// Scannum
        if(this.rt!=-1){
            stringBuilder.append(this.rt + "\t");// rt
        }
        stringBuilder.append(this.getMassErrorMZ() + "\t");// delta

        for (int i = 0; i < this.getFeatureValueList().size(); i++) {
            stringBuilder.append(this.getFeatureValueList().get(i) + "\t");

        }
        if(addFlankAA){
            stringBuilder.append("X."+this.pepSeq +".X"+ "\t");
        }else{
            stringBuilder.append(this.pepSeq + "\t");
        }
        for (int j = 0; j < this.getProteins().size(); j++) {
            if (j == this.getProteins().size() - 1) {
                stringBuilder.append(this.getProteins().get(j));
            } else {
                stringBuilder.append(this.getProteins().get(j) + "\t");
            }
        }
        return stringBuilder.toString();
    }

    public String getPreAA() {
        return preAA;
    }
    public void setPreAA(String preAA) {
        this.preAA = preAA;
    }
    public String getPostAA() {
        return postAA;
    }
    public void setPostAA(String postAA) {
        this.postAA = postAA;
    }
    public double getRt() {
        return rt;
    }
    public void setRt(double rt) {
        this.rt = rt;
    }
    public double getIsotopeError() {
        return IsotopeError;
    }
    public void setIsotopeError(double isotopeError) {
        IsotopeError = isotopeError;
    }
    public double getMassErrorMZ() {
        return massErrorMZ;
    }
    public void setMassErrorMZ(double massErrorMZ) {
        this.massErrorMZ = massErrorMZ;
    }
    public double getQvalue() {
        return qvalue;
    }
    public void setQvalue(double qvalue) {
        this.qvalue = qvalue;
    }
}
