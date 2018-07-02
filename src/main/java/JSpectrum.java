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
    private void resetPeaks() {
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
    public void module1(String method, Boolean labelAssociated, Boolean byFreeWinLow, Boolean byFreeWinHigh) {
        Double isobaric_tag;
        if (method.equals("iTRAQ4plex")) {
            isobaric_tag = 144.102063;
            ArrayList<Double> itraq4Reporters = new ArrayList<Double>();
            itraq4Reporters.add(114.110679);
            itraq4Reporters.add(115.107714);
            itraq4Reporters.add(116.111069);
            itraq4Reporters.add(117.114424);
            if (!labelAssociated) {
                filterLabelAssociatedIons(itraq4Reporters, isobaric_tag);
            }
            if (!byFreeWinLow || !byFreeWinHigh) {
                filterBYfreeWins(byFreeWinLow, byFreeWinHigh, isobaric_tag);
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
            if (!labelAssociated) {
                filterLabelAssociatedIons(itraq8Reporters, isobaric_tag);
            }
            if (!byFreeWinLow || !byFreeWinHigh) {
                filterBYfreeWins(byFreeWinLow, byFreeWinHigh, isobaric_tag);
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
            if (!labelAssociated) {
                filterLabelAssociatedIons(tmt6Reporters, isobaric_tag);
            }
            if (!byFreeWinLow || !byFreeWinHigh) {
                filterBYfreeWins(byFreeWinLow, byFreeWinHigh, isobaric_tag);
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
            if (!labelAssociated) {
                filterLabelAssociatedIons(tmt10Reporters, isobaric_tag);
            }
            if (!byFreeWinLow || !byFreeWinHigh) {
                filterBYfreeWins(byFreeWinLow, byFreeWinHigh, isobaric_tag);
            }

        } else {
            System.out.println("Please check out the proper setting for labeling method.");

        }
    }

    private void filterBYfreeWins(Boolean low, Boolean high, Double isobaric_tag) {
        Double h = ElementaryIon.proton.getTheoreticMass();
        Double gly = 57.021464;
        Double arg = 156.101111;
        Double lys = 128.094963;
        Double h2o = h * 2 + 15.99491463;
        Double min_b = isobaric_tag + gly + h;
        Double max_b = getParentMass() - arg - h2o;
        Double min_y = arg + h2o + h;      // isobaric_tag + lys + h2o + h;
        Double max_y = getParentMass() - gly - isobaric_tag - h2o;
        Double lowWin = Math.min(min_b, min_b);
        Double highWin = Math.max(max_b, max_y);
        ArrayList<JPeak> uninformative = new ArrayList<JPeak>();
        for (JPeak jPeak : getPeaks()) {
            Double mz = jPeak.getMz();
            if (low) {
                if (mz < lowWin - 0.02) {
                    uninformative.add(jPeak);
                }
            }
            if (high) {
                if (mz > highWin + 0.02) {
                    uninformative.add(jPeak);
                }
            }
        }
        getPeaks().removeAll(uninformative);
    }

    private void filterLabelAssociatedIons(ArrayList<Double> reporters, Double isobaric_tag) {
        ArrayList<JPeak> uninformative = new ArrayList<JPeak>();
        Double labelH = isobaric_tag + ElementaryIon.proton.getTheoreticMass();
        Double precusorMinuslabel = getParentMass() - isobaric_tag;
        for (JPeak jPeak : getPeaks()) {
            Double mz = jPeak.getMz();
            if (Config.deltaPPM(labelH, mz) <= 20) {
                uninformative.add(jPeak);
                continue;
            }
            if (Config.deltaPPM(precusorMinuslabel, mz) <= 20) {
                uninformative.add(jPeak);
                continue;
            }
            for (Double mass : reporters) {
                if (Config.deltaPPM(mass, mz) <= 20) {
                    uninformative.add(jPeak);
                    break;
                }
            }
        }
        getPeaks().removeAll(uninformative);
    }



    public double getParentMass() {
        return parentMass;
    }

    public void setParentMass(double parentMass) {
        this.parentMass = parentMass;
    }

}



/*
*
* to sort peaks in descending order of m/z
* sortPeaksByIntensity(); Collections.reverse(getPeaks());
* */