import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.massspectrometry.MSnSpectrum;
import com.compomics.util.experiment.massspectrometry.Peak;
import com.compomics.util.experiment.massspectrometry.SpectrumFactory;
import org.apache.commons.cli.*;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshallerException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 *
 * Created by dengyamei on 30/06/2018.
 */


public class pCleanMainClass {
    public static void main(String[] args) throws ParseException, IOException, MzMLUnmarshallerException {
        Options options = new Options();
        options.addOption("i", true, "MS/MS data in MGF format");
        options.addOption("o", true, "Output directory");
        options.addOption("itol", true, "Fragment tolerance");
        options.addOption("mionFilter", false, "immonium ions filter");
        options.addOption("labelMethod", true, "Peptide labeling method, iTRAQ4plex, iTRAQ8plex or TMT6, TMT12");
        options.addOption("repFilter", false, "remove reporter ions");
        options.addOption("labelFilter", false, "remove label-associated ions");
        options.addOption("low", false, "removal of low b-/y-free window");
        options.addOption("high", false, "removal of low b-/y-free window");
        options.addOption("isoReduction", false, "reduction of heavy isotopic peaks");
        options.addOption("chargeDeconv", false, "high charge deconvolution");
        options.addOption("ionsMerge", false, "merge two adjacent peaks within a mass tolerance of 20ppm");
        options.addOption("largerThanPrecursor", false, "remove peaks larger than precursor");
        options.addOption("a2", false, "Consider gap masses of two amino acids");
        options.addOption("m", true, "mzIdentML file");
        options.addOption("h", false, "Help info");

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();
            System.out.println("java -Xmx2G pClean.jar");
            f.printHelp("Options", options);
            System.exit(0);
        }
        String mgf = cmd.getOptionValue("i");
        String outdir = cmd.getOptionValue("o");
        File OD = new File(outdir);
        if (OD.isDirectory() && OD.exists()) {

        } else {
            OD.mkdirs();
        }

        if (cmd.hasOption("itol")) {
            Config.ms2tol = Double.valueOf(cmd.getOptionValue("itol"));
        } else {
            Config.ms2tol = 0.05;
        }

        /*Default treatments for MS/MS spectra before applying a module*/
        JSpectrum.setImmoniumIons();
        Boolean imonFilter = cmd.hasOption("mionFilter");

        String labelMethod = null;
        if (cmd.hasOption("labelMethod")) {
            labelMethod = cmd.getOptionValue("labelMethod");
        }
        Boolean repFilter = cmd.hasOption("repFilter");
        Boolean labelFilter = cmd.hasOption("labelFilter");
        Boolean lowWinFilter = cmd.hasOption("low");
        Boolean highWinFilter = cmd.hasOption("high");
        Boolean isoReduction = cmd.hasOption("isoReduction");
        Boolean chargeDeconv = cmd.hasOption("chargeDeconv");
        Boolean ionsMerge = cmd.hasOption("ionsMerge");
        Boolean largerThanPrecursor = cmd.hasOption("largerThanPrecursor");

        if (cmd.hasOption("a2")) {
            DeltaMassDB.consider2aa = true;
        } else {
            DeltaMassDB.consider2aa = false;
        }

        String outlog = outdir+"/spectrumInfor.txt";
        if (cmd.hasOption("m")) {
            String mzID = cmd.getOptionValue("m");

        } else {
            if (labelMethod != null) {
                //System.out.println(mgf + "\n" + outdir + "\n" + outlog + "\n" + imonFilter + "\n" + labelMethod + "\n" + repFilter + "\n" + labelFilter + "\n" + lowWinFilter + "\n" + highWinFilter + "\n" + isoReduction + "\n" + chargeDeconv + "\n" + ionsMarge + "\n" + largerThanPrecursor);
                //doModul1Test(mgf, outdir, outlog, imonFilter, labelMethod, repFilter, labelFilter, lowWinFilter, highWinFilter); /*finish test, pass*/
                //doModule123Test(mgf, outdir, outlog, imonFilter, labelMethod, repFilter, labelFilter, lowWinFilter, highWinFilter, isoReduction, chargeDeconv, ionsMarge, largerThanPrecursor); /*finish test, pass*/

                doPreprocessing123(mgf, outdir, outlog, imonFilter, labelMethod, repFilter, labelFilter, lowWinFilter, highWinFilter, isoReduction, chargeDeconv, ionsMerge, largerThanPrecursor);
            } else {
                //doImoniumIonsFilterTest(mgf, outdir, outlog, imonFilter); /*finish test, pass*/
                //doModule2Test(mgf, outdir, outlog, imonFilter, isoReduction, chargeDeconv, ionsMarge, largerThanPrecursor); /*finish test, pass*/
                //doPreprocessing23Test(mgf, outdir, outlog, imonFilter, isoReduction, chargeDeconv, ionsMarge, largerThanPrecursor); /*finish test, pass*/

                doPreprocessing23(mgf, outdir, outlog, imonFilter, isoReduction, chargeDeconv, ionsMerge, largerThanPrecursor);
            }
        }
    }


    /*For high-resolution label-based MS/MS data: implemented module 1, 2, and 3*/
    private static void doPreprocessing123(String mgf, String outdir, String outlog, Boolean imonFilter, String labelMethod, Boolean repFilter, Boolean labelFilter, Boolean lowWinFilter, Boolean highWinFilter, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMerge, Boolean largerThanPrecursor) throws IOException, MzMLUnmarshallerException {
        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        spectrumFactory.addSpectra(mgfFile, null);

        /*spectrumInfor.txt*/
        StringBuilder logStrBuilder = new StringBuilder();
        logStrBuilder.append("index\ttitle\tedge\tvertex\tmz\tintensity\tcharge\n");

        /*parsing MS/MS spectrum one by one*/
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());
        for (int k=0;k<tList.size();k++) {
            String outprefix = "spectrum" + k;
            /*peak annotation file*/
            String peakfile = outdir + "/" + outprefix + "_peak.txt";
            BufferedWriter bWriter = new BufferedWriter(new FileWriter(new File(peakfile)));
            bWriter.write("name\ttype\tintensity\n");

            MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(mgfFile.getName(),tList.get(k));
            /*construct a JSpectrum object*/
            JSpectrum jSpectrum = new JSpectrum();
            Integer ch = spectrum.getPrecursor().getPossibleCharges().get(0).value;
            jSpectrum.setParentMass(spectrum.getPrecursor().getMassPlusProton(ch));
            jSpectrum.setParentMassToCharge(spectrum.getPrecursor().getMz());
            jSpectrum.setCharge(ch);
            jSpectrum.setSpectrumTitle(spectrum.getSpectrumTitle());
            jSpectrum.setIntensity(spectrum.getPrecursor().getIntensity());
            jSpectrum.setRt(spectrum.getPrecursor().getRt());
            for (Peak p : spectrum.getPeakList()) {
                JPeak jPeak = new JPeak(p.getMz(), p.getIntensity());
                jSpectrum.addRawPeak(jPeak);
            }
            jSpectrum.resetPeaks();
            if (imonFilter) {
                jSpectrum.removeImmoniumIons();
            }

            /*module1 treatment*/
            jSpectrum.sortPeaksByMZ();
            jSpectrum.module1(labelMethod, repFilter, labelFilter, lowWinFilter, highWinFilter);

            /*moudle2 treatment*/
            jSpectrum.sortPeaksByMZ();
            jSpectrum.module2(isoReduction, chargeDeconv, ionsMerge, largerThanPrecursor);

            /*module3 treatment*/
            jSpectrum.sortPeaksByMZ();
            StringBuilder peakBuilder = new StringBuilder();
            for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
                JPeak jPeak = jSpectrum.getPeaks().get(i);
                peakBuilder.append(jPeak.getMz());
                peakBuilder.append("\tNA\t");
                peakBuilder.append(jPeak.getIntensity());
                peakBuilder.append("\n");
            }
            bWriter.write(peakBuilder.toString() + "\n");
            bWriter.close();

            /*individual and original MS/MS spectrum*/
            /*String outMS1 = outdir + "/" + outprefix + "_ms2-raw.mgf";
            BufferedWriter msWriter1 = new BufferedWriter(new FileWriter(new File(outMS1)));
            msWriter1.write(spectrum.asMgf() + "\n");
            msWriter1.close();*/

            for(int i=0;i<jSpectrum.getPeaks().size();i++){
                JPeak jPeak = jSpectrum.getPeaks().get(i);
                jPeak.setID(String.valueOf(jPeak.getMz()));
                /*
                * if user did not set charge deconvolution, charge deconvolution will be temporarily transformed to charge 1,
                * this is necessary for graph-based network filtration.
                * */
                if (!chargeDeconv) {
                    jPeak.setMz(jPeak.getMz() * jPeak.getCharge() - (jPeak.getCharge() - 1) * ElementaryIon.proton.getTheoreticMass());
                }
            }

            ArrayList<JPeakPair> jPeakPairs = findPeakPair(jSpectrum, labelMethod);
            /*output edge.txt*/
            StringBuilder edgeStrBuilder = new StringBuilder();
            edgeStrBuilder.append("From\tTo\tmass1\tmass2\tcharge1\tcharge2\tmztol\tdelta\tdeltaName\tintensity\n");
            for(JPeakPair jpp:jPeakPairs){
                edgeStrBuilder.append(jpp.print());
                edgeStrBuilder.append("\n");
            }
            String edgefile = outdir + "/" + outprefix + "_edge.txt";
            BufferedWriter edgeWriter = new BufferedWriter(new FileWriter(new File(edgefile)));
            edgeWriter.write(edgeStrBuilder.toString());
            edgeWriter.close();

            /*for spectrumInfor.txt*/
            logStrBuilder.append(outprefix + "\t" + spectrum.getSpectrumTitle() + "\t" + edgefile + "\t" + peakfile + "\t" + spectrum.getPrecursor().getMz() + "\t" + spectrum.getPrecursor().getIntensity() + "\t" + ch + "\n");
        }
        BufferedWriter logWriter = new BufferedWriter(new FileWriter(new File(outlog)));
        logWriter.write(logStrBuilder.toString());
        logWriter.close();
    }

    /*For high-resolution label-free MS/MS data: implemented module 2 and 3*/
    private static void doPreprocessing23(String mgf, String outdir, String outlog, Boolean imon, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMerge, Boolean largerThanPrecursor) throws IOException, MzMLUnmarshallerException {
        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        spectrumFactory.addSpectra(mgfFile, null);

        /*spectrumInfor.txt*/
        StringBuilder logStrBuilder = new StringBuilder();
        logStrBuilder.append("index\ttitle\tedge\tvertex\tmz\tintensity\tcharge\n");

        /*parsing MS/MS spectrum one by one*/
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());
        for (int k=0;k<tList.size();k++) {
            String outprefix = "spectrum" + k;
            /*peak annotation file*/
            String peakfile = outdir + "/" + outprefix + "_peak.txt";
            BufferedWriter bWriter = new BufferedWriter(new FileWriter(new File(peakfile)));
            bWriter.write("name\ttype\tintensity\n");

            MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(mgfFile.getName(),tList.get(k));
            /*construct a JSpectrum object*/
            JSpectrum jSpectrum = new JSpectrum();
            Integer ch = spectrum.getPrecursor().getPossibleCharges().get(0).value;
            jSpectrum.setParentMass(spectrum.getPrecursor().getMassPlusProton(ch));
            jSpectrum.setParentMassToCharge(spectrum.getPrecursor().getMz());
            jSpectrum.setCharge(ch);
            jSpectrum.setSpectrumTitle(spectrum.getSpectrumTitle());
            jSpectrum.setIntensity(spectrum.getPrecursor().getIntensity());
            jSpectrum.setRt(spectrum.getPrecursor().getRt());
            for (Peak p : spectrum.getPeakList()) {
                JPeak jPeak = new JPeak(p.getMz(), p.getIntensity());
                jSpectrum.addRawPeak(jPeak);
            }
            jSpectrum.resetPeaks();
            if (imon) {
                jSpectrum.removeImmoniumIons();
            }

            /*moudle2 treatment*/
            jSpectrum.sortPeaksByMZ();
            jSpectrum.module2(isoReduction, chargeDeconv, ionsMerge, largerThanPrecursor);

            /*module3 treatment*/
            jSpectrum.sortPeaksByMZ();
            StringBuilder peakBuilder = new StringBuilder();
            for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
                JPeak jPeak = jSpectrum.getPeaks().get(i);
                peakBuilder.append(jPeak.getMz());
                peakBuilder.append("\tNA\t");
                peakBuilder.append(jPeak.getIntensity());
                peakBuilder.append("\n");
            }
            bWriter.write(peakBuilder.toString() + "\n");
            bWriter.close();

            /*individual and original MS/MS spectrum*/
            /*String outMS1 = outdir + "/" + outprefix + "_ms2-raw.mgf";
            BufferedWriter msWriter1 = new BufferedWriter(new FileWriter(new File(outMS1)));
            msWriter1.write(spectrum.asMgf() + "\n");
            msWriter1.close();*/

            for(int i=0;i<jSpectrum.getPeaks().size();i++){
                JPeak jPeak = jSpectrum.getPeaks().get(i);
                jPeak.setID(String.valueOf(jPeak.getMz()));
                /*
                * if user did not set charge deconvolution, charge deconvolution will be temporarily transformed to charge 1,
                * this is necessary for graph-based network filtration.
                * */
                if (!chargeDeconv) {
                    jPeak.setMz(jPeak.getMz() * jPeak.getCharge() - (jPeak.getCharge() - 1) * ElementaryIon.proton.getTheoreticMass());
                }
            }

            ArrayList<JPeakPair> jPeakPairs = findPeakPair(jSpectrum, null);
            /*output edge.txt*/
            StringBuilder edgeStrBuilder = new StringBuilder();
            edgeStrBuilder.append("From\tTo\tmass1\tmass2\tcharge1\tcharge2\tmztol\tdelta\tdeltaName\tintensity\n");
            for(JPeakPair jpp:jPeakPairs){
                edgeStrBuilder.append(jpp.print());
                edgeStrBuilder.append("\n");
            }
            String edgefile = outdir + "/" + outprefix + "_edge.txt";
            BufferedWriter edgeWriter = new BufferedWriter(new FileWriter(new File(edgefile)));
            edgeWriter.write(edgeStrBuilder.toString());
            edgeWriter.close();

            /*for spectrumInfor.txt*/
            logStrBuilder.append(outprefix + "\t" + spectrum.getSpectrumTitle() + "\t" + edgefile + "\t" + peakfile + "\t" + spectrum.getPrecursor().getMz() + "\t" + spectrum.getPrecursor().getIntensity() + "\t" + ch + "\n");
        }
        BufferedWriter logWriter = new BufferedWriter(new FileWriter(new File(outlog)));
        logWriter.write(logStrBuilder.toString());
        logWriter.close();
    }

    private static ArrayList<JPeakPair> findPeakPair(JSpectrum jSpectrum, String labelMethod) {
        ArrayList<JPeakPair> jPeakPairList = new ArrayList<JPeakPair>();
        DeltaMassDB deltaMassDB = new DeltaMassDB(labelMethod);
        deltaMassDB.init();
        jSpectrum.sortPeaksByMZ();
        int npeaks = jSpectrum.getPeaks().size();
        HashSet<String> usedPeak = new HashSet<String>();
        HashSet<String> findPairs = new HashSet<String>();
        for(int i=0;i<npeaks;i++){
            JPeak jPeak = jSpectrum.getPeaks().get(i);
            double mz = jPeak.getMz();
            String pID = jPeak.getID();
            if(usedPeak.contains(pID)){
                continue;
            }
            for (int j = i; j < npeaks; j++) {
                if (i == j) {
                    continue;
                }
                double mzNext = jSpectrum.getPeaks().get(j).getMz();
                String pID2 = jSpectrum.getPeaks().get(j).getID();
                if (pID.equals(pID2)) {
                    continue;
                }

                /*formula 1: Mass(bi)+Mass(yn-i) = Mass(protonated precursor)+1*/
                double delta2parent = Math.abs((mz + mzNext - 2) - jSpectrum.getParentMass());
                if (delta2parent <= 2 * Config.ms2tol) {
                    double sumInt = jPeak.getIntensity() + jSpectrum.getPeaks().get(j).getIntensity();
                    /*make sure the peak-peak-match is unique*/
                    if (Double.valueOf(pID) >= Double.valueOf(pID2)) {
                        if (!findPairs.contains(pID2 + ";" + pID)) {
                            findPairs.add(pID2 + ";" + pID);
                        } else {
                            continue;
                        }
                    } else {
                        if (!findPairs.contains(pID + ";" + pID2)) {
                            findPairs.add(pID + ";" + pID2);
                        } else {
                            continue;
                        }
                    }
                    JPeakPair jPeakPair = new JPeakPair();
                    jPeakPair.setFrom(pID);
                    jPeakPair.setTo(pID2);
                    jPeakPair.mass1 = mz;
                    jPeakPair.mass2 = mzNext;
                    jPeakPair.charge1 = jPeak.getCharge();
                    jPeakPair.charge2 = jSpectrum.getPeaks().get(j).getCharge();
                    jPeakPair.delta = delta2parent;
                    jPeakPair.mzdelta = delta2parent;
                    jPeakPair.deltaName = "Parent";
                    jPeakPair.intensity = sumInt;
                    jPeakPairList.add(jPeakPair);

                    usedPeak.add(pID);
                } else {
                    /*
                    * formula 2: |Mass(bi or yi) - Mass(bi-1 or yi-1)| ∈ S{AMnative, AMvariant}
                    * formula 4: |Mass(bi or yi) - Mass(bi-1 or yi-1)| ∈ S{AAMnative, AAMvariant}
                    * */
                    ArrayList<PMass> matchPMass = deltaMassDB.searchDB(Math.abs(mz - mzNext), Config.ms2tol, 2);
                    if (matchPMass.size() == 1) {
                        double sumInt = jPeak.getIntensity() + jSpectrum.getPeaks().get(j).getIntensity();
                        if (Double.valueOf(pID) >= Double.valueOf(pID2)) {
                            if (!findPairs.contains(pID2 + ";" + pID)) {
                                findPairs.add(pID2 + ";" + pID);
                            } else {
                                continue;
                            }
                        } else {
                            if (!findPairs.contains(pID + ";" + pID2)) {
                                findPairs.add(pID + ";" + pID2);
                            } else {
                                continue;
                            }
                        }
                        JPeakPair jPeakPair = new JPeakPair();
                        jPeakPair.setFrom(pID);
                        jPeakPair.setTo(pID2);
                        jPeakPair.mass1 = mz;
                        jPeakPair.mass2 = mzNext;
                        jPeakPair.charge1 = jPeak.getCharge();
                        jPeakPair.charge2 = jSpectrum.getPeaks().get(j).getCharge();
                        jPeakPair.delta = matchPMass.get(0).deltaMZ;
                        jPeakPair.mzdelta = matchPMass.get(0).mzTol;
                        jPeakPair.deltaName = matchPMass.get(0).name;
                        jPeakPair.intensity = sumInt;
                        jPeakPairList.add(jPeakPair);

                        usedPeak.add(pID);
                    } else if (!matchPMass.isEmpty() && matchPMass.size() > 1) {
                        double sumInt = jPeak.getIntensity() + jSpectrum.getPeaks().get(j).getIntensity();
                        if (Double.valueOf(pID) >= Double.valueOf(pID2)) {
                            if (!findPairs.contains(pID2 + ";" + pID)) {
                                findPairs.add(pID2 + ";" + pID);
                            } else {
                                continue;
                            }
                        } else {
                            if (!findPairs.contains(pID + ";" + pID2)) {
                                findPairs.add(pID + ";" + pID2);
                            } else {
                                continue;
                            }
                        }

                        /*
                        * check out the returned matchPMass by searchDB when searching deltaMassDB library is corresponding to an amino acid's mass,
                        * or the sum mass of aa.
                        * The edge with an amino-acid relationship is preferable.
                        * */
                        double deltaTmp = matchPMass.get(0).deltaMZ;
                        double mzdeltaTmp = matchPMass.get(0).mzTol;
                        String deltaNameTmp = matchPMass.get(0).name;
                        for (PMass pm : matchPMass) {
                            if (pm.name.length() == 1) {
                                deltaTmp = pm.deltaMZ;
                                mzdeltaTmp = pm.mzTol;
                                deltaNameTmp = pm.name;
                            }
                        }
                        JPeakPair jPeakPair = new JPeakPair();
                        jPeakPair.setFrom(pID);
                        jPeakPair.setTo(pID2);
                        jPeakPair.mass1 = mz;
                        jPeakPair.mass2 = mzNext;
                        jPeakPair.charge1 = jPeak.getCharge();
                        jPeakPair.charge2 = jSpectrum.getPeaks().get(j).getCharge();
                        jPeakPair.delta = deltaTmp;
                        jPeakPair.mzdelta = mzdeltaTmp;
                        jPeakPair.deltaName = deltaNameTmp;
                        jPeakPair.intensity = sumInt;
                        jPeakPairList.add(jPeakPair);

                        usedPeak.add(pID);
                    } else {
                        /*
                        * if none matchPMass is returned, check out amino acid mass complementation relationships
                        *
                        * formula 3: |Mass(protonated precursor) - (Mass(bi or yi) - Mass(yn-i-1 or bn-i-1))| ∈ S{AMnative, AMvariant}
                        * formula 5: |Mass(protonated precursor) - (Mass(bi or yi) - Mass(yn-i-1 or bn-i-1))| ∈ S{AAMnative, AAMvariant}
                        *
                        * */
                        ArrayList<PMass> matchPMass2 = deltaMassDB.searchDB(delta2parent, Config.ms2tol, 2);
                        if (!matchPMass2.isEmpty() && matchPMass2.size() >= 1) {
                            double sumInt = jPeak.getIntensity() + jSpectrum.getPeaks().get(j).getIntensity();
                            if (Double.valueOf(pID) >= Double.valueOf(pID2)) {
                                if (!findPairs.contains(pID2 + ";" + pID)) {
                                    findPairs.add(pID2 + ";" + pID);
                                } else {
                                    continue;
                                }
                            } else {
                                if (!findPairs.contains(pID + ";" + pID2)) {
                                    findPairs.add(pID + ";" + pID2);
                                } else {
                                    continue;
                                }
                            }
                            String deltaName = "Parent-" + matchPMass2.get(0).name;
                            JPeakPair jPeakPair = new JPeakPair();
                            jPeakPair.setFrom(pID);
                            jPeakPair.setTo(pID2);
                            jPeakPair.mass1 = mz;
                            jPeakPair.mass2 = mzNext;
                            jPeakPair.charge1 = jPeak.getCharge();
                            jPeakPair.charge2 = jSpectrum.getPeaks().get(j).getCharge();
                            jPeakPair.delta = matchPMass2.get(0).deltaMZ;
                            jPeakPair.mzdelta = matchPMass2.get(0).mzTol;
                            jPeakPair.deltaName = deltaName;
                            jPeakPair.intensity = sumInt;
                            jPeakPairList.add(jPeakPair);

                            usedPeak.add(pID);
                        }
                    }
                }
            }
        }
        return jPeakPairList;
    }

}
