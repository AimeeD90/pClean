package MZIDparser;

import Preprocessing.Config;
import Preprocessing.JPeak;
import Preprocessing.JSpectrum;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.IonMatch;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by dengyamei on 10/07/2018.
 */
public class CalcFeature {
    private double nearTol = 0.3;
    private JSpectrum jSpectrum;
    public HashMap<Double, ArrayList<Double>> featureVal = new HashMap<Double, ArrayList<Double>>();
    public HashMap<Double, ArrayList<String>> featureTag = new HashMap<Double, ArrayList<String>>();
    public CalcFeature(JSpectrum js){
        this.jSpectrum = js;

    }

    public JSpectrum getSpectrum(){
        return this.jSpectrum;
    }

    public boolean classifyPeakIntensities(){
        boolean yes =false;
        int totalPeaks = jSpectrum.getPeaks().size();
        if(totalPeaks>=7){
            yes = true;
            //at present, three types are set by default, the ratio is 1:2:4
            int npeak = (int) Math.floor(1.0*totalPeaks/7.0);
            jSpectrum.sortPeaksByIntensity();
            int nclass1 = npeak;
            int nclass2 = 2*npeak;
            for(int i=totalPeaks-1;i>=totalPeaks-nclass1;i--){
                //System.out.println("i="+i+","+totalPeaks+","+nclass1);
                jSpectrum.getPeaks().get(i).setIntenClass(0);
            }
            //class two
            for(int i=totalPeaks-1-nclass1;i>=totalPeaks-nclass1-nclass2;i--){
                jSpectrum.getPeaks().get(i).setIntenClass(1);
            }
            //class three
            for(int i=totalPeaks-1-nclass1-nclass2;i>=0;i--){
                jSpectrum.getPeaks().get(i).setIntenClass(2);
            }
        }else{
            for(int i=0;i<jSpectrum.getPeaks().size();i++){
                jSpectrum.getPeaks().get(i).setIntenClass(-1);
            }
        }
        return yes;
    }



    public void run(){
        classifyPeakIntensities();
        jSpectrum.sortPeaksByMZ();
        int npeaks = jSpectrum.getPeaks().size();
        double maxInt = jSpectrum.getMaxIntensityFragmentIonPeak().getIntensity();
        double sumInt = jSpectrum.getTotalIntensity();

        PAminoAcid pAminoAcid = new PAminoAcid();
        for(int i=0;i<npeaks;i++){
            JPeak jPeak = jSpectrum.getPeaks().get(i);
            double mz = jPeak.getMz();
            double intensity = jPeak.getIntensity();
            featureVal.put(mz, new ArrayList<Double>());
            featureTag.put(mz, new ArrayList<String>());

            featureVal.get(mz).add(mz);
            featureTag.get(mz).add("mz");
            featureVal.get(mz).add(Math.log(intensity));
            featureTag.get(mz).add("log(rawInt)");
            featureVal.get(mz).add(intensity/maxInt);
            featureTag.get(mz).add("ratio2MaxInt");

            featureVal.get(mz).add(Math.log(intensity/maxInt));
            featureTag.get(mz).add("log(ratio2MaxInt)");

            featureVal.get(mz).add((double) jPeak.getIntenClass());
            featureTag.get(mz).add("peakClass");

            double mzDelta2right=0;
            double mzDelta2left=0;
            if(i==0){
                mzDelta2left = 0;
                mzDelta2right = jSpectrum.getPeaks().get(1).getMz() - mz;
            }else if(i==npeaks-1){
                mzDelta2right = 0;
                mzDelta2left = mz - jSpectrum.getPeaks().get(i-1).getMz();
            }else{
                mzDelta2right = jSpectrum.getPeaks().get(i+1).getMz() - mz;
                mzDelta2left = mz - jSpectrum.getPeaks().get(i-1).getMz();
            }

            featureVal.get(mz).add(mzDelta2left);
            featureTag.get(mz).add("mzDelta2left");

            featureVal.get(mz).add(mzDelta2right);
            featureTag.get(mz).add("mzDelta2right");

            double mzLeftBorder = (mz-50)<0?0:mz-50;
            double mzRightBorder = mz+50;
            double nRegionPeak = 1;
            double intSumRegionPeak = intensity;
            //int rankInRegion = 1;
            double maxIntRegion = intensity;

            double nearPeaks = 0;
            double greatPeaks = 0;
            double match2aa = 0;
            double match2aaInt =0 ;
            double match2aaIntGreat50 = 0;
            double match2aaIntGreat30 = 0;
            double match2aaIntGreat20 = 0;
            double match2aaIntGreat10 = 0;
            double match2aaIntGreat5 = 0;
            double maxInt50 = 0.5*maxInt;
            double maxInt30 = 0.3*maxInt;
            double maxInt20 = 0.2*maxInt;
            double maxInt10 = 0.1*maxInt;
            double maxInt5 = 0.05*maxInt;
            double minms2Delta = 100;

            double pairPeaks = 0;
            for(int j=0;j<npeaks;j++){
                if(i==j){
                    continue;
                }
                double mzNext = jSpectrum.getPeaks().get(j).getMz();
                double intNext = jSpectrum.getPeaks().get(j).getIntensity();

                long mzDelta = Math.round(Math.abs(mz-mzNext));

                if(pAminoAcid.aa.containsKey((mzDelta)) && Math.abs(Math.abs(mz-mzNext)-pAminoAcid.aa.get(mzDelta)) <= Config.ms2tol){
                    double tol = Math.abs(Math.abs(mz-mzNext)-pAminoAcid.aa.get(mzDelta));
                    match2aa++;
                    match2aaInt+=intNext;
                    minms2Delta = tol < minms2Delta?tol:minms2Delta;
                    if(intNext >= maxInt50){
                        match2aaIntGreat50++;
                    }
                    if(intNext >= maxInt30){
                        match2aaIntGreat30++;
                    }
                    if(intNext >= maxInt20){
                        match2aaIntGreat20++;
                    }
                    if(intNext >= maxInt10){
                        match2aaIntGreat10++;
                    }
                    if(intNext >= maxInt5){
                        match2aaIntGreat5++;
                    }
                    //System.out.println(jSpectrum.getSpectrumTitle());
                    //if(jSpectrum.getSpectrumTitle().equals("40")){
                    System.err.println(mz+"="+mzNext+","+Math.abs(mz-mzNext)+","+pAminoAcid.aa.get(mzDelta));


                }

                if((mz+mzNext-jSpectrum.getParentMassToCharge()) <=2){
                    pairPeaks++;
                }

                if(Math.abs(intensity-intNext) <this.nearTol*intensity){
                    nearPeaks++;
                }
                if(intNext>=intensity){
                    greatPeaks++;
                }
                //double mzDelta = mzNext - mz;
                if(mzNext >=mzLeftBorder && mzNext <= mzRightBorder ){
                    nRegionPeak++;
                    intSumRegionPeak+=intNext;
                    maxIntRegion = intNext>maxIntRegion?intNext:maxIntRegion;
                }
            }

            featureVal.get(mz).add(match2aa);
            featureTag.get(mz).add("match2aa");

            featureVal.get(mz).add(match2aa/npeaks);
            featureTag.get(mz).add("match2aa2npeaks");

            featureVal.get(mz).add(minms2Delta);
            featureTag.get(mz).add("minms2Delta");

            featureVal.get(mz).add(minms2Delta*1000000/mz);
            featureTag.get(mz).add("minms2DeltaPPM");

            featureVal.get(mz).add(match2aaInt);
            featureTag.get(mz).add("match2aaInt");

            featureVal.get(mz).add(Math.log(match2aaInt+1));
            featureTag.get(mz).add("log(match2aaInt)");

            featureVal.get(mz).add(match2aaInt/maxInt);
            featureTag.get(mz).add("match2aaInt2maxInt");

            featureVal.get(mz).add(match2aaInt/sumInt);
            featureTag.get(mz).add("match2aaInt2sumInt");

            featureVal.get(mz).add(match2aaIntGreat50);
            featureTag.get(mz).add("match2aaIntGreat50");
            featureVal.get(mz).add(match2aaIntGreat30);
            featureTag.get(mz).add("match2aaIntGreat30");
            featureVal.get(mz).add(match2aaIntGreat20);
            featureTag.get(mz).add("match2aaIntGreat20");
            featureVal.get(mz).add(match2aaIntGreat10);
            featureTag.get(mz).add("match2aaIntGreat10");
            featureVal.get(mz).add(match2aaIntGreat5);
            featureTag.get(mz).add("match2aaIntGreat5");

            featureVal.get(mz).add(match2aaIntGreat50/npeaks);
            featureTag.get(mz).add("match2aaIntGreat502npeaks");
            featureVal.get(mz).add(match2aaIntGreat30/npeaks);
            featureTag.get(mz).add("match2aaIntGreat302npeaks");
            featureVal.get(mz).add(match2aaIntGreat20/npeaks);
            featureTag.get(mz).add("match2aaIntGreat202npeaks");
            featureVal.get(mz).add(match2aaIntGreat10/npeaks);
            featureTag.get(mz).add("match2aaIntGreat102npeaks");
            featureVal.get(mz).add(match2aaIntGreat5/npeaks);
            featureTag.get(mz).add("match2aaIntGreat52npeaks");

            featureVal.get(mz).add(pairPeaks);
            featureTag.get(mz).add("pairPeaks");

            featureVal.get(mz).add(nearPeaks);
            featureTag.get(mz).add("nearPeaks");

            featureVal.get(mz).add(nearPeaks/npeaks);
            featureTag.get(mz).add("nearPeaks2npeaks");

            featureVal.get(mz).add(greatPeaks);
            featureTag.get(mz).add("greatPeaks");

            featureVal.get(mz).add(greatPeaks/npeaks);
            featureTag.get(mz).add("greatPeaks2npeaks");


            featureVal.get(mz).add(intensity/intSumRegionPeak);
            featureTag.get(mz).add("intRatio2regionSum");

            featureVal.get(mz).add(nRegionPeak);
            featureTag.get(mz).add("nRegionPeak");

            featureVal.get(mz).add(intensity/maxIntRegion);
            featureTag.get(mz).add("intRatio2regionMax");

            featureVal.get(mz).add(intSumRegionPeak/sumInt);
            featureTag.get(mz).add("intRegionRatio2totalSum");

            featureVal.get(mz).add(nRegionPeak/npeaks);
            featureTag.get(mz).add("nRegionPeak2npeaks");

        }

    }



    public void addMatch(ArrayList<IonMatch> ionMatchs){
        HashMap<Double,IonMatch> mzMatch = new HashMap<Double,IonMatch>();
        for(IonMatch im:ionMatchs){

            mzMatch.put(im.peak.getMz(),im);
        }
        for(int i=0;i<jSpectrum.getPeaks().size();i++){
            if(mzMatch.containsKey(jSpectrum.getPeaks().get(i).getMz())){
                jSpectrum.getPeaks().get(i).isMatch = 1;
                IonMatch ionMatch = mzMatch.get(jSpectrum.getPeaks().get(i).getMz());
                PeptideFragmentIon fragmentIon = ((PeptideFragmentIon) ionMatch.ion);
                String label = ionMatch.ion.getSubTypeAsString() + fragmentIon.getNumber() + ionMatch.ion.getNeutralLossesAsString() + ionMatch.charge.intValue();

                jSpectrum.getPeaks().get(i).setIonType(label);
            }
        }
    }
}
