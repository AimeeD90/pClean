package MZIDparser;

import Preprocessing.JPeak;
import Preprocessing.JSpectrum;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.NeutralLoss;
import com.compomics.util.experiment.biology.PTMFactory;
import com.compomics.util.experiment.biology.Peptide;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.IonMatch;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.spectrum_annotation.AnnotationSettings;
import com.compomics.util.experiment.identification.spectrum_annotation.SpecificAnnotationSettings;
import com.compomics.util.experiment.identification.spectrum_annotation.spectrum_annotators.PeptideSpectrumAnnotator;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.massspectrometry.Charge;
import com.compomics.util.experiment.massspectrometry.MSnSpectrum;
import org.apache.commons.math.MathException;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Created by dengyamei on 10/07/2018.
 */
public class PeakMatch {
    //member variable
    public boolean bUseDynamicRange=true; // enables using the dynamic range
    public boolean bUseLowestMass=true; // enables the removal of very low m/z peaks from a spectrum
    public boolean bUseMaxPeaks=true; // enables removing low intensity peaks
    public boolean bUseMinMass=true; // sets the minimum parent ion mass allowed
    public boolean bUseMinSize=true; // enables using the minimum number of peaks to exclude spectra
    public boolean bUseNeutralLoss=true;
    public int iMaxPeaks=50; // the maximum number of peaks in a spectrum
    public double dDynamicRange=100; // the normalized intensity of the most intense peak in a spectrum;
    public double dLowestMass = 150.0; // the lowest m/z in a spectrum，低于这个质荷比的峰将会被去掉
    public int iMinSize=5; // the minimum number of peaks in a spectrum;
    public static double ms2tol = 0.5;//MS/MS tol. ±， only da unit is supported
    public double intensityLimit = 0.02;//低于1%的峰将被去掉
    public double ticCutoffPercentage = 0.98; //myrimatch用于谱图预处理过滤的参数

    //the following three parameter must be assigned
    public MSnSpectrum spectrum;
    public ArrayList<IonMatch> matches = new ArrayList<IonMatch>();
    public JSpectrum jSpectrum = new JSpectrum();
    public JPSM jpsm = new JPSM();

    //fragmentation pattern
    public String fragmentMethod = "cid"; //default CID
    public int maxFragmentChargeState = 3;
    //private AnnotationPreferences annotationPreferences = new AnnotationPreferences();
    private boolean lossWaterNH3=false;

    public PeakMatch(){

    }

    /**
     * main constructor
     * @param spectrum
     * @param jSpectrum
     * @param jpsm
     */
    public PeakMatch(MSnSpectrum spectrum, JSpectrum jSpectrum, JPSM jpsm){
        this.spectrum = spectrum;
        this.jSpectrum = jSpectrum;
        this.jpsm = jpsm;
    }

    public void initialize(boolean lossWaterNH3, int maxFragCharge) throws IllegalArgumentException, FileNotFoundException, ClassNotFoundException, IOException, InterruptedException, SQLException, MathException {
        setLossWaterNH3(lossWaterNH3);
        setMaxFragmentChargeState(maxFragCharge);
        getSpectrumAnnotation();

    }


    /**
     * remove low mass peaks in MS/MS spectrum
     * @return
     */
    public int removeLowMasses(){
        int rpn = 0;
        Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();
        while(pListIterator.hasNext()){
            JPeak jPeak = pListIterator.next();
            if( jPeak.getMz() <= this.dLowestMass){
                //System.out.println("remove low mass peak:"+jPeak.getMz());
                rpn++;
                pListIterator.remove();
            }
        }

        return rpn;
    }

    /**
     * select peaks with high intensity, the number of peaks is iMaxPeaks, default 50
     * @return
     */
    public boolean filterByPeakCount(){
        boolean isRemove = false;
        int peaksNumber = this.jSpectrum.getPeaks().size();
        if(peaksNumber>this.iMaxPeaks){
            this.jSpectrum.sortPeaksByIntensity();
            this.jSpectrum.getPeaks().subList(0, peaksNumber-this.iMaxPeaks).clear();
            //System.out.println("remove none top 50 peak!");
            isRemove=true;
        }
        return isRemove;
    }

    /**
     * select peaks with high intensity, the number of peaks is iMaxPeaks, user defined
     * @return
     */
    public boolean filterByPeakCount(int maxPeakCount){
        boolean isRemove = false;
        int peaksNumber = this.jSpectrum.getPeaks().size();
        if(peaksNumber>maxPeakCount){
            this.jSpectrum.sortPeaksByIntensity();
            this.jSpectrum.getPeaks().subList(0, peaksNumber-maxPeakCount).clear();
            //System.out.println("remove none top maxPeakCount peak!");
            isRemove=true;
        }
        return isRemove;
    }

    /**
     * remove peaks whose intensity lower than highest intensity * intensityLimit
     * @return
     */
    public int removeLowIntensityPeak(){
        int rpn = 0;
        double maxIntensity = this.jSpectrum.getMaxIntensityFragmentIonPeak().getIntensity();
        double limitIntensity = 1.0*this.intensityLimit*maxIntensity;
        //System.out.println("limitIntensity:"+limitIntensity);

        Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();
        while(pListIterator.hasNext()){
            JPeak jPeak = pListIterator.next();
            if( jPeak.getIntensity()<limitIntensity){
                //System.out.println("remove peak:"+jPeak.getMz()+"\t"+jPeak.getIntensity());
                //System.out.println("remove low intensity peak:"+jPeak.getMz());
                rpn++;
                pListIterator.remove();
            }
        }

        return rpn;
    }

    /**
     * remove peaks whose intensity lower than highest intensity * intensityLimit,
     * user defined intensityLimit, percent
     * @param interLimit
     * @return
     */
    public int removeLowIntensityPeak( double interLimit){
        int rpn = 0;
        double maxIntensity = this.jSpectrum.getMaxIntensityFragmentIonPeak().getIntensity();
        double limitIntensity = 1.0*interLimit*maxIntensity;
        //System.out.println("limitIntensity:"+limitIntensity);

        Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();
        while(pListIterator.hasNext()){
            JPeak jPeak = pListIterator.next();
            if( jPeak.getIntensity()<limitIntensity){
                rpn++;
                pListIterator.remove();
            }
        }
        return rpn;
    }

    /**
     * remove isotopes removes multiple entries within 0.95 Da of each other, retaining
     * the highest value. this is necessary because of the behavior of some peak
     * finding routines in commercial software
     * @return
     */
    public int removeIsotopes(){
        int rpn = 0;
        if(this.jSpectrum.getPeaks().size()>2){
            this.getJSpectrum().sortPeaksByMZ();// 必须先把峰按mz从小到大排序
            JPeak jPeak1 = this.getJSpectrum().getPeaks().get(0);
            JPeak jPeak2 = this.getJSpectrum().getPeaks().get(1);
            ArrayList< JPeak> tmpJPeaks = new ArrayList<JPeak>();
            double mz1 = jPeak1.getMz();
            for(int i=1;i<this.getJSpectrum().getPeaks().size();i++){
                jPeak2 = this.getJSpectrum().getPeaks().get(i);
                if(jPeak2.getMz() - mz1 >= 0.95 || jPeak2.getMz() < 200){
                    tmpJPeaks.add(jPeak1);
                    jPeak1 = jPeak2;
                    mz1 = jPeak1.getMz();
                }else if (jPeak2.getIntensity() > jPeak1.getIntensity()){
                    jPeak1 = jPeak2;
                    mz1 = jPeak1.getMz();//注意tandem的代码里是没有这个的
                }
            }
            rpn = this.jSpectrum.getPeaks().size() - tmpJPeaks.size();
            tmpJPeaks.add(jPeak1);
            this.jSpectrum.setPeaks(tmpJPeaks);
        }
        return rpn;
    }


    /**
     * clean_isotopes removes peaks that are probably C13 isotopes
     * @return
     */
    public int cleanIsotopes(){
        int rpn = 0;
        if(this.jSpectrum.getPeaks().size()>2){
            this.getJSpectrum().sortPeaksByMZ();// 必须先把峰按mz从小到大排序
            JPeak jPeak1 = this.getJSpectrum().getPeaks().get(0);
            JPeak jPeak2 = this.getJSpectrum().getPeaks().get(1);
            ArrayList< JPeak> tmpJPeaks = new ArrayList<JPeak>();
            double mz1 = jPeak1.getMz();
            for(int i=1;i<this.getJSpectrum().getPeaks().size();i++){
                jPeak2 = this.getJSpectrum().getPeaks().get(i);
                if(jPeak2.getMz() - mz1 >= 1.5 || jPeak2.getMz() < 200){
                    tmpJPeaks.add(jPeak1);
                    jPeak1 = jPeak2;
                    mz1 = jPeak1.getMz();
                }else if (jPeak2.getIntensity() > jPeak1.getIntensity()){
                    jPeak1 = jPeak2;
                    mz1 = jPeak1.getMz();//注意tandem的代码里是没有这个的
                }
            }
            rpn = this.jSpectrum.getPeaks().size() - tmpJPeaks.size();
            tmpJPeaks.add(jPeak1);
            this.jSpectrum.setPeaks(tmpJPeaks);
        }
        return rpn;
    }

    /**
     * normalize spectrum intensity, the maximum intensity set as maxPeakAfterNornalize
     * @param maxPeakAfterNornalize
     */
    public void normalizePeaks(double maxPeakAfterNornalize){
        double maxIntensity = this.jSpectrum.getMaxIntensityFragmentIonPeak().getIntensity();
        for(int i=0;i<this.jSpectrum.getPeaks().size();i++){
            double norIntensity = 1.0*maxPeakAfterNornalize*this.jSpectrum.getPeaks().get(i).getIntensity()/maxIntensity;
            this.jSpectrum.getPeaks().get(i).setIntensity(norIntensity);
        }
    }

    public MSnSpectrum getMSnSpectrum() {
        return spectrum;
    }

    /**
     * set a PSM object using a MSnSpectrum instance, different form Spectrum
     * @param spectrum
     */
    public void setMSnSpectrum(MSnSpectrum spectrum) {
        this.spectrum = spectrum;
    }

    public ArrayList<IonMatch> getIonMatches() {
        return matches;
    }

    public void setIonMatches(ArrayList<IonMatch> matches) {
        this.matches = matches;
    }

    public JSpectrum getJSpectrum() {
        return jSpectrum;
    }

    /**
     * set a JSpectrum object using a JSpectrum instance
     * @param jSpectrum
     */
    public void setJSpectrum(JSpectrum jSpectrum) {
        this.jSpectrum = jSpectrum;
    }

    public ArrayList<IonMatch> getCopyOfIonMatches(){
        ArrayList<IonMatch> newIonMatches = new ArrayList<IonMatch>();
        for(IonMatch ionMatch :this.getIonMatches()){
            newIonMatches.add(ionMatch);

        }
        return newIonMatches;
    }

    public String getFragmentMethod() {
        return fragmentMethod;
    }

    public void setFragmentMethod(String fragmentMethod) {
        this.fragmentMethod = fragmentMethod;
    }

    public JPSM getPSM(){
        return this.jpsm;
    }

    public void setPSM(JPSM j){
        this.jpsm = j;
    }

    //the following functions are mainly used for Myrimatch MS/MS spectrum preprocessing
    /**
     * Filters out the peaks with the lowest intensities until only
     * ticCutoffPercentage of the total ion current remains
     * @return
     */
    public int filterByTIC(){
        int rpn=0;
        // Sort peak list in descending order of intensity while calculating the total ion current in the spectrum.
        // Use a multimap because multiple peaks can have the same intensity.
        if (this.jSpectrum.getPeaks().size() > 7) { //at least more than 7 peaks, because these peaks will be divided into three types, each type has 1, 2, and 4 peaks at least
            double totalIonCurrent = this.jSpectrum.getTotalIonCurrent();
            double relativeIntensity=0.0;
            this.jSpectrum.sortPeaksByIntensity();
            ArrayList< JPeak> tmpJPeaks = new ArrayList<JPeak>();

            for(int i=this.jSpectrum.getPeaks().size()-1;i>=0;i--){
                relativeIntensity+=this.jSpectrum.getPeaks().get(i).getIntensity()/totalIonCurrent;
                if(relativeIntensity<=this.ticCutoffPercentage){
                    tmpJPeaks.add(this.jSpectrum.getPeaks().get(i));
                }else{
                    rpn++;
                }
            }
            this.jSpectrum.setPeaks(tmpJPeaks);
        }
        return rpn;
    }

    public int getMaxFragmentChargeState() {
        return maxFragmentChargeState;
    }

    public void setMaxFragmentChargeState(int maxFragmentChargeState) {
        this.maxFragmentChargeState = maxFragmentChargeState;
    }

    /**
     * set a mz range, then filter fragment ions
     * @param lowBoundmz
     * @param upBoundmz
     * @return
     */
    public int filterPeakByRange(double lowBoundmz,double upBoundmz){
        int rpn = 0;
        Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();
        while(pListIterator.hasNext()){
            JPeak jPeak = pListIterator.next();
            if(jPeak.getMz() < lowBoundmz || jPeak.getMz() > upBoundmz){
                rpn++;
                pListIterator.remove();
            }
        }
        return rpn;
    }


    /**
     * split peaks into fixed-size regions (default 100Da), select topN intensity peaks at each region
     * @param n
     * @param regionMass
     * @return
     * @throws InterruptedException
     */
    public int keepTopNPeaksInRegions(int n, double regionMass) throws InterruptedException{
        int rpn = 0;
        HashMap<Double, Integer> mz2region = new HashMap<Double, Integer>();
        HashMap<Integer, Integer> regionP = new HashMap<Integer, Integer>();
        for(JPeak jPeak:this.jSpectrum.getPeaks()){
            int rg = (int) Math.floor(1.0*jPeak.getMz()/regionMass) - 1;
            if(rg<0){
                continue;
            }
            mz2region.put(jPeak.getMz(), rg);
            regionP.put(rg, 0);
        }
        this.jSpectrum.sortPeaksByIntensity();
        Collections.reverse(this.jSpectrum.getPeaks());
        Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();
        while(pListIterator.hasNext()){
            JPeak jPeak = pListIterator.next();
            double ionMz = jPeak.getMz();
            if(!mz2region.containsKey(ionMz)){
                pListIterator.remove();
                rpn++;
                continue;
            }

            int rg = mz2region.get(ionMz);
            if(regionP.get(rg) >= n){
                pListIterator.remove();
                rpn++;
            }else{
                regionP.put(rg, regionP.get(rg)+1);
            }
        }

        return rpn;
    }


    /**
     * remove mz larger that cutoff
     * @param cutoff
     * @return
     */
    public int removePeakGreaterBy(double cutoff){
        int rpn = 0;
        Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();
        while(pListIterator.hasNext()){
            JPeak jPeak = pListIterator.next();
            if( jPeak.getMz() > cutoff){
                //System.out.println("remove low mass peak:"+jPeak.getMz());
                rpn++;
                pListIterator.remove();
            }
        }

        return rpn;
    }

    /**
     * remove adjecent peaks of precursor
     * @param cutoff +-cutoff
     * @return
     */
    public int removePeakAroundParent(double cutoff){
        int rpn = 0;
        Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();
        double precursor_mz = this.jSpectrum.getParentMassToCharge();
        while(pListIterator.hasNext()){
            JPeak jPeak = pListIterator.next();
            if( jPeak.getMz() > precursor_mz - cutoff && jPeak.getMz() < precursor_mz + cutoff){
                //System.out.println("remove low mass peak:"+jPeak.getMz());
                rpn++;
                pListIterator.remove();
            }
        }

        return rpn;
    }

    /**
     * match spectrum to peptide
     * @return
     * @throws IllegalArgumentException
     * @throws FileNotFoundException
     * @throws ClassNotFoundException
     * @throws IOException
     * @throws InterruptedException
     * @throws SQLException
     */
    public void getSpectrumAnnotation() throws IllegalArgumentException, ClassNotFoundException, IOException, InterruptedException, SQLException, MathException {

        jSpectrum.resetPeaks();
        jSpectrum.sortPeaksByMZ();
        Peptide objPeptide = new Peptide(this.jpsm.getPepSeq(), this.jpsm.getModificationMatch());
        /*
        System.out.println("mass:" + objPeptide.getMass());
        System.out.println("ptms:" + objPeptide.getSequenceWithLowerCasePtms());
        */
        ArrayList<ModificationMatch> mm = objPeptide.getModificationMatches();
        PTMFactory ptmFactory = PTMFactory.getInstance();
        for (int m = 0; m < mm.size(); m++) {
            /*System.out.println();
            System.out.println(mm.get(m).getModificationSite() + "," + mm.get(m).getTheoreticPtm());
            System.out.println(ptmFactory.getPTM(mm.get(m).getTheoreticPtm()).getMass());
            System.out.println(ptmFactory.getPTM(mm.get(m).getTheoreticPtm()).getName());
            System.out.println(ptmFactory.getPTM(mm.get(m).getTheoreticPtm()).getType());*/
        }


        PeptideSpectrumAnnotator peptideSpectrumAnnotator = new PeptideSpectrumAnnotator();

        String spectrumKey = spectrum.getSpectrumKey();
        Charge objCharge = new Charge(Charge.PLUS, this.jpsm.getCharge());
        PeptideAssumption peptideAssumption = new PeptideAssumption(objPeptide, 1, 1, objCharge, 1.0);
        SpecificAnnotationSettings specificAnnotationPreferences = new SpecificAnnotationSettings(spectrumKey, peptideAssumption);

        ArrayList<Integer> charges = new ArrayList<Integer>(4);
        int precursorCharge = peptideAssumption.getIdentificationCharge().value;
        if (precursorCharge == 1) {
            charges.add(precursorCharge);
        } else {
            for (int c = 1; c < precursorCharge; c++) {
                charges.add(c);
            }
        }
        specificAnnotationPreferences.setSelectedCharges(charges);

        specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, PeptideFragmentIon.B_ION);
        specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, PeptideFragmentIon.Y_ION);
        specificAnnotationPreferences.setFragmentIonAccuracy(ms2tol);
        specificAnnotationPreferences.setFragmentIonPpm(false);
        specificAnnotationPreferences.setNeutralLossesAuto(false);
        specificAnnotationPreferences.clearNeutralLosses();
        specificAnnotationPreferences.addNeutralLoss(NeutralLoss.H2O);
        specificAnnotationPreferences.addNeutralLoss(NeutralLoss.NH3);


        AnnotationSettings annotationSettings = new AnnotationSettings();
        annotationSettings.setHighResolutionAnnotation(false);
        annotationSettings.setFragmentIonAccuracy(ms2tol);
        annotationSettings.setFragmentIonPpm(false);
        annotationSettings.setIntensityLimit(intensityLimit);
        ArrayList<IonMatch> matches = peptideSpectrumAnnotator.getSpectrumAnnotation(annotationSettings, specificAnnotationPreferences, spectrum, objPeptide);


        //System.out.println(matches.size()+"fffff");
        //System.out.println(ms2tol);
        if (matches.isEmpty()) {
            System.err.println("No ions matched!");
        }
        this.matches = matches;

    }

	/*
    public Peptide getObjPeptide(){
		if(this.objPeptide==null){
			this.objPeptide = new Peptide( jpsm.getPepSeq(),jpsm.getModificationMatch());
		}
		return this.objPeptide;
	}
	*/

    public boolean isLossWaterNH3() {
        return lossWaterNH3;
    }

    public void setLossWaterNH3(boolean lossWaterNH3) {
        this.lossWaterNH3 = lossWaterNH3;
    }

    public static double getMs2tol() {
        return ms2tol;
    }

    public static void setMs2tol(double va) {
        ms2tol = va;
    }
}
