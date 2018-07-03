import com.compomics.util.experiment.biology.ions.ElementaryIon;

import java.util.*;

/**
 * Object for a MS/MS spectrum
 * provide with various functions to manipulate MS/MS spectrum's properties
 *
 * Created by dengyamei on 30/06/2018.
 */
public class JSpectrum {
    private String spectrumTitle;
    private String scanNumber;
    private double rt;
    private int charge;
    private double parentMassToCharge;
    private double intensity;
    private double parentMass;
    /*store the information of fragment ions for normalization*/
    private ArrayList<JPeak> peaks = new ArrayList<JPeak>();
    /*store the original information of fragment ions*/
    private ArrayList<JPeak> rawPeaks = new ArrayList<JPeak>();

    /*set the list of immonium ions*/
    public static HashMap<Double, String> ImmoniumIons = new HashMap<Double, String>();
    public static void setImmoniumIons(){
        ImmoniumIons.put(30.03, "G");
        ImmoniumIons.put(44.05, "A");
        ImmoniumIons.put(60.04, "S");
        ImmoniumIons.put(70.07, "P/R");
        ImmoniumIons.put(72.08, "V");
        ImmoniumIons.put(74.06, "T");
        ImmoniumIons.put(86.10, "I/L");
        ImmoniumIons.put(87.06, "N");
        ImmoniumIons.put(70.03, "N/D");
        ImmoniumIons.put(88.04, "D");
        ImmoniumIons.put(129.10, "Q");
        ImmoniumIons.put(101.07, "Q");
        ImmoniumIons.put(84.04, "Q/E");
        ImmoniumIons.put(56.05, "Q/K");
        ImmoniumIons.put(129.11, "K/R");
        ImmoniumIons.put(101.11, "K");
        ImmoniumIons.put(84.08, "K");
        ImmoniumIons.put(101.11, "K");
        ImmoniumIons.put(102.05, "E");
        ImmoniumIons.put(104.06, "M");
        ImmoniumIons.put(120.08, "F");
        ImmoniumIons.put(110.07, "H");
        ImmoniumIons.put(115.09, "R");
        ImmoniumIons.put(112.09, "R");
        ImmoniumIons.put(87.09, "R");
        ImmoniumIons.put(60.06, "R");
        ImmoniumIons.put(133.04, "C");
        ImmoniumIons.put(136.08, "Y");
        ImmoniumIons.put(159.09, "W");
        ImmoniumIons.put(132.08, "W");
        ImmoniumIons.put(130.07, "W");
    }

    /*
    * Empty constructor
    * */
    public JSpectrum() {
    }


    /*Getter and Setter*/
    public String getSpectrumTitle() {
        return spectrumTitle;
    }

    public void setSpectrumTitle(String spectrumTitle) {
        this.spectrumTitle = spectrumTitle;
    }

    public String getScanNumber() {
        return scanNumber;
    }

    public void setScanNumber(String scanNumber) {
        this.scanNumber = scanNumber;
    }

    public double getRt() {
        return rt;
    }

    public void setRt(double rt) {
        this.rt = rt;
    }

    public int getCharge() {
        return charge;
    }

    public void setCharge(int charge) {
        this.charge = charge;
    }

    public double getParentMassToCharge() {
        return parentMassToCharge;
    }

    public void setParentMassToCharge(double parentMassToCharge) {
        this.parentMassToCharge = parentMassToCharge;
    }

    public double getIntensity() {
        return intensity;
    }

    public void setIntensity(double intensity) {
        this.intensity = intensity;
    }

    public void addPeak(JPeak peak){
        this.peaks.add(peak);
    }

    public void setPeaks(ArrayList<JPeak> peaks) {
        this.peaks = peaks;
    }

    public double getParentMass() {
        return parentMass;
    }

    public void setParentMass(double parentMass) {
        this.parentMass = parentMass;
    }

    public ArrayList<JPeak> getPeaks() {
        if(this.peaks.size()==0){
            this.resetPeaks();
        }
        return this.peaks;
    }

    /*remove peak from peaks based on the index of this peak*/
    public void removePeak(int index){
        this.peaks.remove(index);
    }

    /*initialize rawPeaks*/
    public void addRawPeak(JPeak peak){
        this.rawPeaks.add(peak);
    }

    /*get the information of peaks from the original MS/MS spectrum*/
    public ArrayList<JPeak> getRawPeaks() {
        return rawPeaks;
    }

    /*remove peak from rawPeaks based on the index of this peak*/
    public void removeRawPeak(int index){
        this.rawPeaks.remove(index);
    }

    /*reset peaks as rawPeaks*/
    public void resetPeaks() {
        this.peaks = new ArrayList<JPeak>();
        for (JPeak jPeak : this.rawPeaks) {
            JPeak tmpJPeak = jPeak.clone();
            this.peaks.add(tmpJPeak);
        }
    }

    /*sort peaks based on m/z in ascending order*/
    public void sortPeaksByMZ(){
        PeakMZComparator pComparator = new PeakMZComparator();
        Collections.sort(this.peaks, pComparator);
    }

    /*sort peaks based on intensity in ascending order*/
    public void sortPeaksByIntensity(){
        PeakIntensityComparator pComparator = new PeakIntensityComparator();
        Collections.sort(this.peaks,pComparator);
    }

    /*print the objects of peaks to screen for checking out*/
    public void printPeaks(){
        this.sortPeaksByMZ();
        int i=1;
        for(JPeak peak: this.peaks){
            System.out.println(i + "\t" + peak.getMz() + "\t" + peak.getIntensity());
        }
    }

    /*calculate the sum intensity of all the fragment peaks*/
    public double getTotalIntensity(){
        double sumInt = 0.0;
        for(JPeak jPeak: this.peaks){
            sumInt += jPeak.getIntensity();
        }
        return sumInt;
    }

    /*get the peak object with the maximum intensity*/
    public JPeak getMaxIntensityFragmentIonPeak(){
        sortPeaksByIntensity();
        return peaks.get(peaks.size()-1);
    }

    /*get the peak object with the maximum m/z*/
    public JPeak getMaxMzFragmentIonPeak(){
        sortPeaksByMZ();
        return peaks.get(peaks.size()-1);
    }

    /*get the peak object with the maximum m/z*/
    public JPeak getMinMzFragmentPeak(){
        sortPeaksByMZ();
        return peaks.get(0);
    }

    /*normalization of intensity*/
    public void normalization(){
        double maxInt = getMaxIntensityFragmentIonPeak().getIntensity();
        for(JPeak jPeak: this.peaks){
            jPeak.setIntensity(jPeak.getIntensity()/maxInt);
        }
    }

    /*format the decimal state of mz value*/
    public void formatMZdeciPoint() {
        for (int i = 0; i < peaks.size(); i++) {
            peaks.get(i).setMz(Double.valueOf(String.format("%.4f", peaks.get(i).getMz())));
        }
    }



    /*
    * Default preprocessing.
    * 1. remove immonium ions
    * */
    public int removeImmoniumIons() {
        Iterator<JPeak> pListIterator = this.getPeaks().iterator();
        int rpn = 0;
        while (pListIterator.hasNext()) {
            JPeak jPeak = pListIterator.next();
            if (jPeak.getMz() >= 170) {
                continue;
            }
            for (Double mion : ImmoniumIons.keySet()) {
                if (Math.abs(jPeak.getMz() - mion) <= Config.ms2tol) {
                    //System.out.println("mion:" + mion + "\tpeak:" + jPeak.getMz());
                    rpn++;
                    pListIterator.remove();
                    break;
                }

            }
        }
        return rpn;
    }



    /*
    * 1. intensity-based preprocessing method
    * a. remove peaks with intensity lower than the threshold: limit * maxIntensity, limit ranges from 0 to 1.
    * */
    public int removeLowIntensityPeak(double intensityLimit){
        int rpn = 0;
        double maxInt = this.getMaxIntensityFragmentIonPeak().getIntensity();
        double limitInt = 1.0 * intensityLimit * maxInt;
        //System.out.println("limitIntensity:"+limitIntensity);

        Iterator<JPeak> pListIterator = this.getPeaks().iterator();
        while(pListIterator.hasNext()){
            JPeak jPeak = pListIterator.next();
            if( jPeak.getIntensity()<limitInt){
                //System.out.println("remove peak:"+jPeak.getMz()+"\t"+jPeak.getIntensity());
                //System.out.println("remove low intensity peak:"+jPeak.getMz());
                rpn++;
                pListIterator.remove();
            }
        }
        return rpn;
    }

    /*
    * 1. intensity-based preprocessing method
    * b. keep TopN peaks in fixed regions, such as keep Top10 in every 50Da scale.
    * */
    public int keepTopNPeaksInRegions(int n, double regionMass) throws InterruptedException{
        int rpn = 0;
        HashMap<Double, Integer> mz2region = new HashMap<Double, Integer>();
        HashMap<Integer, Integer> regionP = new HashMap<Integer, Integer>();
        for (JPeak jPeak : getPeaks()) {
            /*here we use 0Da as the beginning of regionMass, while in the real-life data preprocessing, it is needed to recalculate "rg" with the beginning of regionMass equals to (the minimum m/z of this spectrum + the size of regionMass) */
            int rg = (int) Math.floor(1.0 * jPeak.getMz() / regionMass) - 1;

            /*to avoid wrong setting by user*/
            if(rg<0){
                System.out.println("JSpectrum.java line253: rg is smaller than 0!");
                continue;
            }

            mz2region.put(jPeak.getMz(), rg);
            regionP.put(rg, 0);
        }

        sortPeaksByIntensity(); // sort intensity in ascending order
        Collections.reverse(getPeaks()); // reverse intensity in descending order
        Iterator<JPeak> pListIterator = getPeaks().iterator();
        ArrayList<JPeak> jPeaksToRemove = new ArrayList<JPeak>();
        while(pListIterator.hasNext()){
            JPeak jPeak = pListIterator.next();
            double ionMz = jPeak.getMz();
            if(!mz2region.containsKey(ionMz)){
                pListIterator.remove();
                rpn++;
                continue;
            }
            int rg = mz2region.get(ionMz);
            Integer value = regionP.get(rg) + 1;
            regionP.put(rg, value);
            if (regionP.get(rg) > n) {
                jPeaksToRemove.add(jPeak);
                rpn++;
            }
        }
        getPeaks().removeAll(jPeaksToRemove);
        return rpn;
    }


    /*
    * Module 1. removal of label-associated ions and of b-/y-ions free windows
    *
    * */
    public void module1(String method, Boolean removeRep, Boolean labelFilter, Boolean byFreeWinLow, Boolean byFreeWinHigh) {
        Double isobaric_tag;
        if (method.equals("iTRAQ4plex")) {
            isobaric_tag = 144.102063;
            ArrayList<Double> itraq4Reporters = new ArrayList<Double>();
            itraq4Reporters.add(114.110679);
            itraq4Reporters.add(115.107714);
            itraq4Reporters.add(116.111069);
            itraq4Reporters.add(117.114424);
            if (labelFilter) {
                doFilterLabelAssociatedIons(itraq4Reporters, isobaric_tag, removeRep);
            }
            if (byFreeWinLow || byFreeWinHigh) {
                doFilterBYfreeWins(byFreeWinLow, byFreeWinHigh, isobaric_tag);
            }

        } else if (method.equals("iTRAQ8plex")) {
            isobaric_tag = 304.205360;
            ArrayList<Double> itraq8Reporters = new ArrayList<Double>();
            itraq8Reporters.add(113.107324);
            itraq8Reporters.add(114.110679);
            itraq8Reporters.add(115.107714);
            itraq8Reporters.add(116.111069);
            itraq8Reporters.add(117.114424);
            itraq8Reporters.add(118.111459);
            itraq8Reporters.add(119.114814);
            itraq8Reporters.add(121.121523);
            if (labelFilter) {
                doFilterLabelAssociatedIons(itraq8Reporters, isobaric_tag, removeRep);
            }
            if (byFreeWinLow || byFreeWinHigh) {
                System.out.println(byFreeWinHigh + "-high\tlow-" + byFreeWinLow);
                doFilterBYfreeWins(byFreeWinLow, byFreeWinHigh, isobaric_tag);
            }

        } else if (method.equals("TMT6plex")) {
            isobaric_tag = 229.162932;
            ArrayList<Double> tmt6Reporters = new ArrayList<Double>();
            tmt6Reporters.add(126.127726);
            tmt6Reporters.add(127.131080);
            tmt6Reporters.add(128.134435);
            tmt6Reporters.add(129.137790);
            tmt6Reporters.add(130.141145);
            tmt6Reporters.add(131.138180);
            if (labelFilter) {
                doFilterLabelAssociatedIons(tmt6Reporters, isobaric_tag, removeRep);
            }
            if (byFreeWinLow || byFreeWinHigh) {
                doFilterBYfreeWins(byFreeWinLow, byFreeWinHigh, isobaric_tag);
            }

        } else if (method.equals("TMT10plex")) {
            isobaric_tag = 229.162932;
            ArrayList<Double> tmt10Reporters = new ArrayList<Double>();
            tmt10Reporters.add(126.127726);
            tmt10Reporters.add(127.1247610);
            tmt10Reporters.add(127.1310809);
            tmt10Reporters.add(128.1281158);
            tmt10Reporters.add(128.1344357);
            tmt10Reporters.add(129.1314706);
            tmt10Reporters.add(129.1377905);
            tmt10Reporters.add(130.1348254);
            tmt10Reporters.add(130.1411453);
            tmt10Reporters.add(131.1381802);
            if (labelFilter) {
                doFilterLabelAssociatedIons(tmt10Reporters, isobaric_tag, removeRep);
            }
            if (byFreeWinLow || byFreeWinHigh) {
                doFilterBYfreeWins(byFreeWinLow, byFreeWinHigh, isobaric_tag);
            }

        } else {
            System.out.println("Please check out the proper setting for labeling method.");

        }
    }

    private void doFilterBYfreeWins(Boolean low, Boolean high, Double isobaric_tag) {
        double h = ElementaryIon.proton.getTheoreticMass();
        double gly = 57.021464;
        double arg = 156.101111;
        double lys = 128.094963;
        double h2o = h * 2 + 15.99491463;
        double min_b = isobaric_tag + gly + h;
        double max_b = getParentMass() - arg - h2o;
        double min_y = arg + h2o + h;      // isobaric_tag + lys + h2o + h;
        double max_y = getParentMass() - gly - isobaric_tag - h2o;
        double lowWin = Math.min(min_b, min_y);
        double highWin = Math.max(max_b, max_y);
        //System.out.println(highWin + "-highWin\tlowWin-" + lowWin);
        ArrayList<JPeak> uninformative = new ArrayList<JPeak>();
        for (JPeak jPeak : getPeaks()) {
            double mz = jPeak.getMz();
            if (low) {
                if (mz < lowWin - 0.02) {
                    //System.out.println(jPeak.getMz() + "\t<lowWin");
                    uninformative.add(jPeak);
                }
            }
            if (high) {
                if (mz > highWin + 0.02) {
                    //System.out.println(jPeak.getMz() + "\t>highWin");
                    uninformative.add(jPeak);
                }
            }
        }
        //System.out.println("before-" + getPeaks().size());
        getPeaks().removeAll(uninformative);
        //System.out.println("after-" + getPeaks().size());
    }

    private void doFilterLabelAssociatedIons(ArrayList<Double> reporters, double isobaric_tag, Boolean reporter) {
        ArrayList<JPeak> uninformative = new ArrayList<JPeak>();
        double h = ElementaryIon.proton.getTheoreticMass();
        double labelH = isobaric_tag + h;
        double labelH2 = isobaric_tag + h * 2;
        double o = 15.994915;
        double c = 12.0;
        double precusorMinuslabel = getParentMass() - isobaric_tag;

        for (JPeak jPeak : getPeaks()) {
            double mz = jPeak.getMz();
            /*pClean is designed to preprocess high-resolution MS/MS data, so the */
            /*label-associated ions: label+ */
            if (Config.delta(labelH, mz) <= Config.ms2tol) {
                //System.out.println(jPeak.getMz() + "\t" + "label+");
                uninformative.add(jPeak);
                continue;
            }
            /*label-associated ions: label++ */
            if (Config.delta(labelH2, mz) <= Config.ms2tol) {
                //System.out.println(jPeak.getMz() + "\t" + "label++");
                uninformative.add(jPeak);
                continue;
            }
            /*label-associated ions: precursor-label+ */
            if (Config.delta(precusorMinuslabel, mz) <= Config.ms2tol) {
                //System.out.println(jPeak.getMz() + "\t" + "precursor-label+");
                uninformative.add(jPeak);
                continue;
            }
            for (Double mass : reporters) {
                /*remove reporter ions*/
                if (reporter) {
                    if (Config.delta(mass, mz) <= Config.ms2tol) {
                        //System.out.println(jPeak.getMz() + "\t" + "reporter");
                        uninformative.add(jPeak);
                        break;
                    }
                }

                double repCOH = mass + c + o + h;
                double precusorMinusRepCOH = getParentMass() - repCOH;
                /*label-associated ions: repCO+ */
                if (Config.delta(repCOH, mz) <= Config.ms2tol) {
                    //System.out.println(jPeak.getMz() + "\t" + "repCO+");
                    uninformative.add(jPeak);
                    break;
                }
                /*label-associated ions: precursor-repCO+ */
                if (Config.delta(precusorMinusRepCOH, mz) <= Config.ms2tol) {
                    //System.out.println(jPeak.getMz() + "\t" + "precursor-repCO+");
                    uninformative.add(jPeak);
                    break;
                }
            }
        }
        //System.out.println("before-" + getPeaks().size());
        getPeaks().removeAll(uninformative);
        //System.out.println("after-" + getPeaks().size());
    }


    /*
    * Module 2. Isotopic peak reduction & charge deconvolution
    *
    * */
    public void module2(Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMarge, Boolean filter) {
        if (isoReduction) {
            /*for (JPeak jPeak : getPeaks()) {
                System.out.println("before\t" + jPeak.getMz() + "\t" + jPeak.getCharge() + "\t" + jPeak.getIntensity());
            }*/
            doBackwardSearch();
            /*sortPeaksByMZ();
            System.out.println();
            System.out.println();
            System.out.println();
            for (JPeak jPeak : getPeaks()) {
                System.out.println("after\t" + jPeak.getMz() + "\t" + jPeak.getCharge() + "\t" + jPeak.getIntensity());
            }*/

            doIsotopicPeakReductionChargeDeconvolution(chargeDeconv);

            doIonMarge(ionsMarge);

            doFilterIonsLargerThanPrecursor(filter);
        } else {
            /*System.out.println("before:" + getPeaks().size());
            for (JPeak jPeak : getPeaks()) {
                System.out.println(jPeak.getMz() + "\t" + jPeak.getCharge());
            }*/
            fragmentChargePredictionWithoutDeisotoping();
            /*sortPeaksByMZ();
            System.out.println();
            System.out.println();
            System.out.println();
            for (JPeak jPeak : getPeaks()) {
                System.out.println(jPeak.getMz() + "\t" + jPeak.getCharge() + "\t" + jPeak.getIntensity());
            }
            System.out.println("after:" + getPeaks().size());*/

            doIsotopicPeakReductionChargeDeconvolution(chargeDeconv);
            doIonMarge(ionsMarge);
            doFilterIonsLargerThanPrecursor(filter);
        }
    }

    private void doFindEnvelopCluster() {
        /*
        * Note: Module2 can't handle well with overlap envelops, neither doForwordSearch or doBackwardSearch.
        * The only difference between these two methods is when different charge state is detected for a peak,
        * the finial charge state is set as the same with the larger adjacent peak (doForwordSearch) or as the same with
        * the smaller one (doBackwardSearch). The choice of setting is up to the user.
        * */
        doForwardSearch();
        doBackwardSearch();
    }

    /*predict charge state of fragment ions and reduce isotopic peaks in a forward search manner*/
    private void doForwardSearch() {
        initializeFragmentCharge();
        int maxcharge = getCharge();
        double neutron = 1.0033548;
        sortPeaksByMZ();
        for (int i = 0; i < getPeaks().size(); i++) {
            JPeak jPeak1 = getPeaks().get(i);
            for (int j = i + 1; j < getPeaks().size(); j++) {
                JPeak jPeak2 = getPeaks().get(j);
                int charge = maxcharge;
                while (charge > 0) {
                    if (Math.abs(jPeak2.getMz() - jPeak1.getMz() - neutron / charge) < 0.02 && jPeak1.getIntensity() > (jPeak2.getIntensity() / 2.0)) {
                        jPeak2.setIntensity(0.0);
                        jPeak1.setCharge(charge);
                        jPeak2.setCharge(charge);
                    }
                    charge -= 1;
                }
                if ((jPeak2.getMz() - jPeak1.getMz()) > 2) {
                    break;
                }
            }
        }
    }

    /*predict charge state of fragment ions and reduce isotopic peaks in a backward search manner*/
    private void doBackwardSearch() {
        initializeFragmentCharge();
        int maxcharge = getCharge();
        double neutron = 1.0033548;
        sortPeaksByMZ();
        Collections.reverse(getPeaks());
        for (int i = 0; i < getPeaks().size(); i++) {
            JPeak jPeak1 = getPeaks().get(i);
            for (int j = i + 1; j < getPeaks().size(); j++) {
                JPeak jPeak2 = getPeaks().get(j);
                int charge = maxcharge;
                while (charge > 0) {
                    if (Math.abs(jPeak1.getMz() - jPeak2.getMz() - neutron / charge) < 0.02 && jPeak2.getIntensity() > (jPeak1.getIntensity() / 2.0)) {
                        jPeak1.setIntensity(0.0);
                        jPeak1.setCharge(charge);
                        jPeak2.setCharge(charge);
                    }
                    charge -= 1;
                }
                if ((jPeak1.getMz() - jPeak2.getMz()) > 2) {
                    break;
                }
            }
        }
    }

    /*the function is limited to predict charge state of fragment ions, and has done nothing with peak reduction*/
    private void fragmentChargePredictionWithoutDeisotoping() {
        initializeFragmentCharge();
        int maxcharge = getCharge();
        double neutron = 1.0033548;
        sortPeaksByMZ();
        Collections.reverse(getPeaks());
        for (int i = 0; i < getPeaks().size(); i++) {
            JPeak jPeak1 = getPeaks().get(i);
            for (int j = i + 1; j < getPeaks().size(); j++) {
                JPeak jPeak2 = getPeaks().get(j);
                int charge = maxcharge;
                while (charge > 0) {
                    if (Math.abs(jPeak1.getMz() - jPeak2.getMz() - neutron / charge) < Config.ms2tol && jPeak2.getIntensity() > (jPeak1.getIntensity() / 2.0)) {
                        jPeak1.setCharge(charge);
                        jPeak2.setCharge(charge);
                    }
                    charge -= 1;
                }
                if ((jPeak1.getMz() - jPeak2.getMz()) > 2) {
                    break;
                }
            }
        }

        /*peaks with charge state 0 are reset to 1*/
        for (JPeak jPeak : getPeaks()) {
            if (jPeak.getCharge() == 0) {
                jPeak.setCharge(1);
            }
        }
    }

    private void initializeFragmentCharge() {
        for (JPeak jPeak : getPeaks()) {
            jPeak.setCharge(0);
        }
    }

    /*remove peaks with 0.0 intensity (heavy isotopic peaks) and deconvolute higher charges*/
    private void doIsotopicPeakReductionChargeDeconvolution(Boolean chargeDeconv) {
        ArrayList<JPeak> filter = new ArrayList<JPeak>();
        for (JPeak jPeak : getPeaks()) {
            if (jPeak.getIntensity() > 0.0) {
                if (chargeDeconv) {
                    if (jPeak.getCharge() > 1) {
                        double mass = jPeak.getMz() * jPeak.getCharge() - (jPeak.getCharge() - 1) * ElementaryIon.proton.getTheoreticMass();
                        //System.out.println("before-" + jPeak.getMz() + "\tafter-" + mass + "\t" + jPeak.getCharge());
                        jPeak.setMz(mass);
                        jPeak.setCharge(1);
                    } else if (jPeak.getCharge() == 0) {
                        jPeak.setCharge(1);
                    }
                }
            } else {
                filter.add(jPeak);
            }
        }
        //System.out.println(getPeaks().size());
        getPeaks().removeAll(filter);
        //System.out.println(getPeaks().size());
        sortPeaksByMZ();
    }

    /*marge similar ions in a direct way*/
    private void doIonMarge(Boolean ionsMarge) {
        ArrayList<JPeak> jPeaks = new ArrayList<JPeak>();
        double previousMZ = getPeaks().get(0).getMz();
        double previousIntensity = getPeaks().get(0).getIntensity();
        for (int i = 0; i < getPeaks().size(); i++) {
            double mz = getPeaks().get(i).getMz();
            double intensity = getPeaks().get(i).getIntensity();

            if (Math.abs(previousMZ - mz) <= 0.05) { /*0.01 or 0.02 or 0.05???*/
                //System.out.println("pre-" + previousMZ + "\taft-" + mz);
                double mzMarge = (mz + previousMZ) / 2;
                double intensityMarge = (intensity + previousIntensity) / 2;
                previousMZ = mzMarge;
                previousIntensity = intensityMarge;
            } else {
                JPeak jPeak = new JPeak(previousMZ, previousIntensity);
                jPeaks.add(jPeak);
                previousMZ = mz;
                previousIntensity = intensity;

                if (i == (getPeaks().size() - 1)) {
                    JPeak jPeakLastOne = new JPeak(previousMZ, previousIntensity);
                    jPeaks.add(jPeakLastOne);
                }
            }
        }
        setPeaks(jPeaks);
    }

    /*filter out ions with masses larger than precursor*/
    private void doFilterIonsLargerThanPrecursor(Boolean filter) {
        ArrayList<JPeak> filer = new ArrayList<JPeak>();
        for (JPeak jPeak : getPeaks()) {
            if (jPeak.getMz() > getParentMass()) {
                //System.out.println(jPeak.getMz() + "\tlargeThanPrecursor");
                filer.add(jPeak);
            }
        }
        getPeaks().removeAll(filer);
    }


    /*marge similar ions using charge information*/
    /*private void doIonMarge(Boolean ionsMarge) {
        ArrayList<JPeak> jPeaks = new ArrayList<JPeak>();
        double previousMZ = getPeaks().get(0).getMz();
        double previousIntensity = getPeaks().get(0).getIntensity();
        int previousCharge = getPeaks().get(0).getCharge();
        for (int i = 1; i < getPeaks().size(); i++) {
            double mz = getPeaks().get(i).getMz();
            double intensity = getPeaks().get(i).getIntensity();
            int charge = getPeaks().get(i).getCharge();

            if (Math.abs(previousMZ - mz) < 0.01) {
                if (previousCharge > 0 && charge > 0) {
                    double mzMarge = (mz + previousMZ) / 2;
                    double intensityMarge = (intensity + previousIntensity);

                    previousMZ = mzMarge;
                    previousIntensity = intensityMarge;
                    previousCharge = charge;
                } else {
                    //System.out.println("Pay Attention!");
                    //System.out.println(previousMZ + "\t" + mz + "\t" + previousIntensity + "\t" + intensity + "\t" + previousCharge + "\t" + charge);
                    *//*if two adjacent ions have similar mz and the charge of both of them is 0, they should be kept the original state.*//*
                    if (previousCharge == 0 && charge == 0) {
                        if (intensity > previousIntensity) {
                            previousMZ = mz;
                            previousIntensity = intensity;
                            previousCharge = charge;
                        }
                    } else {
                        double mzMarge = (mz + previousMZ) / 2;
                        double intensityMarge = (intensity + previousIntensity);
                        previousMZ = mzMarge;
                        previousIntensity = intensityMarge;
                        previousCharge = charge;
                    }
                }
            } else {
                JPeak jPeak = new JPeak(previousMZ, previousIntensity, previousCharge);
                jPeaks.add(jPeak);
                previousMZ = mz;
                previousIntensity = intensity;
                previousCharge = charge;

                if (i == (getPeaks().size() - 1)) {
                    JPeak jPeakLastOne = new JPeak(previousMZ, previousIntensity, previousCharge);
                    jPeaks.add(jPeakLastOne);
                }
            }
        }
        setPeaks(jPeaks);
    }*/

}