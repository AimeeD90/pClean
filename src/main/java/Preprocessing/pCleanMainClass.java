package Preprocessing;

import DatParser.MascotDatparser;
import DatParser.SpectrumAnnotationContainer;
import MZIDparser.CalcFeature;
import MZIDparser.JPSM;
import MZIDparser.PeakMatch;
import MZIDparser.ReadmzID;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.massspectrometry.MSnSpectrum;
import com.compomics.util.experiment.massspectrometry.Peak;
import com.compomics.util.experiment.massspectrometry.SpectrumFactory;
import org.apache.commons.cli.*;
import org.apache.commons.math.MathException;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshallerException;

import javax.xml.bind.JAXBException;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * Created by dengyamei on 30/06/2018.
 */


public class pCleanMainClass {
    public static final String TOOL = "pClean";
    public static final String VERSION = "Release (v2018.08.06)";
    public static final String RELEASE_DATE = "6 August 2018";

    public static void main(String[] args) throws ParseException, IOException, MzMLUnmarshallerException, JAXBException, ClassNotFoundException, SQLException, MathException, InterruptedException {
        //long time = System.currentTimeMillis();

        ParameterManager params = new ParameterManager(pCleanMainClass.TOOL, pCleanMainClass.VERSION, pCleanMainClass.RELEASE_DATE);
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(params.getOptions(), args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            params.printUsageInfo();
            System.exit(0);
        }
        params.InitializeParameters(cmd);

        /*Default treatments for MS/MS spectra before applying a module*/
        JSpectrum.setImmoniumIons();

        String outlog = params.getOutdir() + "/spectrumInfor.txt";
        if (params.getIdres() == null) { // don't plot png
            if (params.doNetworkFilter()) { // do graph-based network filtration
                if (params.getLabelMethod() != null) { // label-based MS/MS data preprocessing
                    String labelmethod = params.getLabelMethod();
                    doPreprocessing123(params.getInput(), params.getOutdir(), outlog, params.doImionFilter(), params.getLabelMethod(), params.doRepFilter(), params.doLabelFilter(), params.doLowWindowFilter(), params.doHighWindowFilter(), params.doIsoReduction(), params.doChargeDeconv(), params.doIonsMerge(), params.doLargerThanPrecursor());

                } else { // label-free MS/MS data preprocessing
                    doPreprocessing23(params.getInput(), params.getOutdir(), outlog, params.doImionFilter(), params.doIsoReduction(), params.doChargeDeconv(), params.doIonsMerge(), params.doLargerThanPrecursor());

                }

            } else { // don't do graph-based network filtration, and export MGF directly
                if (params.getLabelMethod() != null) {
                    // do module1 and module 2
                    if (params.doIsoReduction() || params.doChargeDeconv() || params.doIonsMerge() || params.doLargerThanPrecursor()) { // true: do module1, module2 and module3
                        doModule12(params.getInput(), params.getOutdir(), params.doImionFilter(), params.getLabelMethod(), params.doRepFilter(), params.doLabelFilter(), params.doLowWindowFilter(), params.doHighWindowFilter(), params.doIsoReduction(), params.doChargeDeconv(), params.doIonsMerge(), params.doLargerThanPrecursor());

                    } else { // do module 1
                        doModule1(params.getInput(), params.getOutdir(), params.doImionFilter(), params.getLabelMethod(), params.doRepFilter(), params.doLabelFilter(), params.doLowWindowFilter(), params.doHighWindowFilter());

                    }

                } else {
                    doModule2(params.getInput(), params.getOutdir(), params.doImionFilter(), params.doIsoReduction(), params.doChargeDeconv(), params.doIonsMerge(), params.doLargerThanPrecursor());
                }


            }
        } else {
            String idres = params.getIdres();
            if (idres.endsWith(".mzid")) {
                if (params.getLabelMethod() != null) {
                    doParseMzid123(params.getInput(), params.getOutdir(), outlog, params.doImionFilter(), params.getLabelMethod(), params.doRepFilter(), params.doLabelFilter(), params.doLowWindowFilter(), params.doHighWindowFilter(), params.doIsoReduction(), params.doChargeDeconv(), params.doIonsMerge(), params.doLargerThanPrecursor(), idres);

                } else {
                    doParseMzid23(params.getInput(), params.getOutdir(), outlog, params.doImionFilter(), params.doIsoReduction(), params.doChargeDeconv(), params.doIonsMerge(), params.doLargerThanPrecursor(), idres);

                }
            } else if (idres.endsWith(".dat")) {
                if (params.getLabelMethod() != null) {
                    doParseDat123(params.getInput(), params.getOutdir(), outlog, params.doImionFilter(), params.getLabelMethod(), params.doRepFilter(), params.doLabelFilter(), params.doLowWindowFilter(), params.doHighWindowFilter(), params.doIsoReduction(), params.doChargeDeconv(), params.doIonsMerge(), params.doLargerThanPrecursor(), idres);

                } else {
                    doParseDat23(params.getInput(), params.getOutdir(), outlog, params.doImionFilter(), params.doIsoReduction(), params.doChargeDeconv(), params.doIonsMerge(), params.doLargerThanPrecursor(), idres);

                }

            }
        }

        //System.out.format("pClean complete (total elapsed time: %.2f sec)\n", (System.currentTimeMillis() - time) / (float) 1000);
    }


    /*
    * perform module 1 preprocessing, at the same time generate a resultant mgf.
    * */
    private static void doModule1(String mgf, String outdir, Boolean imonFilter, String labelMethod, Boolean repFilter, Boolean labelFilter, Boolean lowWinFilter, Boolean highWinFilter) throws IOException, MzMLUnmarshallerException {
        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        String outmgf = outdir + System.getProperty("file.separator") + mgfFile.getName().replace(".mgf", ".pClean_M1.mgf");
        BufferedWriter bwmgf = new BufferedWriter(new FileWriter(new File(outmgf)));

        spectrumFactory.addSpectra(mgfFile, null);
        /*parsing MS/MS spectrum one by one*/
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());
        for (int k = 0; k < tList.size(); k++) {
            /*peak annotation file*/
            MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(mgfFile.getName(), tList.get(k));
            /*construct a Preprocessing.JSpectrum object*/
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

            StringBuilder sb = jSpectrum.toMgf();
            bwmgf.write(sb.toString());
        }
        bwmgf.close();
    }


    /*
    * perform module 2 preprocessing, at the same time generate a resultant mgf.
    * */
    private static void doModule2(String mgf, String outdir, Boolean imonFilter, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMerge, Boolean largerThanPrecursor) throws IOException, MzMLUnmarshallerException {
        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        String outmgf = outdir + System.getProperty("file.separator") + mgfFile.getName().replace(".mgf", ".pClean_M2.mgf");
        BufferedWriter bwmgf = new BufferedWriter(new FileWriter(new File(outmgf)));

        spectrumFactory.addSpectra(mgfFile, null);
        /*parsing MS/MS spectrum one by one*/
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());
        for (int k=0;k<tList.size();k++) {
            MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(mgfFile.getName(),tList.get(k));
            /*construct a Preprocessing.JSpectrum object*/
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

            /*moudle2 treatment*/
            jSpectrum.sortPeaksByMZ();
            jSpectrum.module2(isoReduction, chargeDeconv, ionsMerge, largerThanPrecursor);

            StringBuilder sb = jSpectrum.toMgf();
            bwmgf.write(sb.toString());
        }
        bwmgf.close();
    }

    /*
    * perform module 1 and module 2 preprocessing, at the same time generate a resultant mgf.
    * */
    private static void doModule12(String mgf, String outdir, Boolean imonFilter, String labelMethod, Boolean repFilter, Boolean labelFilter, Boolean lowWinFilter, Boolean highWinFilter, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMerge, Boolean largerThanPrecursor) throws IOException, MzMLUnmarshallerException {
        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        String outmgf = outdir + System.getProperty("file.separator") + mgfFile.getName().replace(".mgf", ".pClean_M12.mgf");
        BufferedWriter bwmgf = new BufferedWriter(new FileWriter(new File(outmgf)));

        spectrumFactory.addSpectra(mgfFile, null);
        /*parsing MS/MS spectrum one by one*/
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());
        for (int k = 0; k < tList.size(); k++) {
            /*peak annotation file*/
            MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(mgfFile.getName(), tList.get(k));
            /*construct a Preprocessing.JSpectrum object*/
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

            StringBuilder sb = jSpectrum.toMgf();
            bwmgf.write(sb.toString());
        }
        bwmgf.close();
    }

    /*
    * perform module 1 and module 2 preprocessing, and not generate a resultant mgf but export _peak.txt and _edge.txt file.
    * */
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
            /*construct a Preprocessing.JSpectrum object*/
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


    /*
    * perform module 2 preprocessing, and not generate a resultant mgf but export _peak.txt and _edge.txt file.
    * */
    private static void doPreprocessing23(String mgf, String outdir, String outlog, Boolean imonFilter, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMerge, Boolean largerThanPrecursor) throws IOException, MzMLUnmarshallerException {
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
            /*construct a Preprocessing.JSpectrum object*/
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


    private static void doParseMzid123(String mgf, String outdir, String outlog, Boolean imonFilter, String labelMethod, Boolean repFilter, Boolean labelFilter, Boolean lowWinFilter, Boolean highWinFilter, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMerge, Boolean largerThanPrecursor, String mzid) throws IOException, MzMLUnmarshallerException, JAXBException, ClassNotFoundException, SQLException, MathException, InterruptedException {
        /*read mzid*/
        ReadmzID readmzID = new ReadmzID(100);
        ArrayList<JPSM> psmArrayList = readmzID.readmzIdentML(mzid);
        HashMap<Integer, JPSM> psmMap = new HashMap<Integer, JPSM>();
        for (JPSM psm : psmArrayList) {
            if (psm.getRank() >= 2) {
                continue;
            }
            psmMap.put(psm.getSpectrumIndex(), psm);
        }

        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        spectrumFactory.addSpectra(mgfFile, null);

        /*spectrumInfor.txt*/
        StringBuilder logStrBuilder = new StringBuilder();
        logStrBuilder.append("index\ttitle\tedge\tvertex\tmz\tintensity\tcharge\n");

        /*parsing MS/MS spectrum one by one*/
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());
        for (int k = 0; k < tList.size(); k++) {
            String outprefix = "spectrum" + k;
            /*peak annotation file*/
            String peakfile = outdir + "/" + outprefix + "_peak.txt";
            BufferedWriter bWriter = new BufferedWriter(new FileWriter(new File(peakfile)));
            bWriter.write("name\ttype\tintensity\n");

            MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(mgfFile.getName(), tList.get(k));
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
            if (psmMap.containsKey(k)) {
                JPSM psm = psmMap.get(k);
                PeakMatch peakMatch = new PeakMatch(spectrum, jSpectrum, psm);
                int maxFragCharge = 2;
                PeakMatch.setMs2tol(Config.ms2tol);
                peakMatch.initialize(true, maxFragCharge);
                CalcFeature calcFeature = new CalcFeature(jSpectrum);
                calcFeature.addMatch(peakMatch.getIonMatches());
                JSpectrum js = calcFeature.getSpectrum();

                /*module1 treatment*/
                jSpectrum.sortPeaksByMZ();
                jSpectrum.module1(labelMethod, repFilter, labelFilter, lowWinFilter, highWinFilter);

                /*moudle2 treatment*/
                js.sortPeaksByMZ();
                js.module2(isoReduction, chargeDeconv, ionsMerge, largerThanPrecursor);
                js.sortPeaksByMZ();
                for (int i = 0; i < js.getPeaks().size(); i++) {
                    JPeak jPeak = js.getPeaks().get(i);
                    if (jPeak.isMatch == 1) {
                        peakBuilder.append(jPeak.getMz());
                        peakBuilder.append("\t" + jPeak.getIonType() + "\t");
                        peakBuilder.append(jPeak.getIntensity());
                        peakBuilder.append("\n");
                    } else {
                        peakBuilder.append(jPeak.getMz());
                        peakBuilder.append("\tNA\t");
                        peakBuilder.append(jPeak.getIntensity());
                        peakBuilder.append("\n");
                    }
                }
            } else {
                for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
                    JPeak jPeak = jSpectrum.getPeaks().get(i);
                    peakBuilder.append(jPeak.getMz());
                    peakBuilder.append("\tNA\t");
                    peakBuilder.append(jPeak.getIntensity());
                    peakBuilder.append("\n");
                }
            }
            bWriter.write(peakBuilder.toString() + "\n");
            bWriter.close();

            /*individual and original MS/MS spectrum*/
            /*String outMS1 = outdir + "/" + outprefix + "_ms2-raw.mgf";
            BufferedWriter msWriter1 = new BufferedWriter(new FileWriter(new File(outMS1)));
            msWriter1.write(spectrum.asMgf() + "\n");
            msWriter1.close();*/

            for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
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
            for (JPeakPair jpp : jPeakPairs) {
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

    private static void doParseMzid23(String mgf, String outdir, String outlog, Boolean imonFilter, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMerge, Boolean largerThanPrecursor, String mzid) throws JAXBException, IOException, MzMLUnmarshallerException, ClassNotFoundException, SQLException, MathException, InterruptedException {
        /*read mzid*/
        ReadmzID readmzID = new ReadmzID(100);
        ArrayList<JPSM> psmArrayList = readmzID.readmzIdentML(mzid);
        HashMap<Integer, JPSM> psmMap = new HashMap<Integer, JPSM>();
        for (JPSM psm : psmArrayList) {
            if (psm.getRank() >= 2) {
                continue;
            }
            psmMap.put(psm.getSpectrumIndex(), psm);
        }

        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        spectrumFactory.addSpectra(mgfFile, null);

        /*spectrumInfor.txt*/
        StringBuilder logStrBuilder = new StringBuilder();
        logStrBuilder.append("index\ttitle\tedge\tvertex\tmz\tintensity\tcharge\n");

        /*parsing MS/MS spectrum one by one*/
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());
        for (int k = 0; k < tList.size(); k++) {
            String outprefix = "spectrum" + k;
            /*peak annotation file*/
            String peakfile = outdir + "/" + outprefix + "_peak.txt";
            BufferedWriter bWriter = new BufferedWriter(new FileWriter(new File(peakfile)));
            bWriter.write("name\ttype\tintensity\n");

            MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(mgfFile.getName(), tList.get(k));
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

            /*moudle2 treatment*/
            jSpectrum.sortPeaksByMZ();
            jSpectrum.module2(isoReduction, chargeDeconv, ionsMerge, largerThanPrecursor);

            /*module3 treatment*/
            jSpectrum.sortPeaksByMZ();
            StringBuilder peakBuilder = new StringBuilder();
            if (psmMap.containsKey(k)) {
                JPSM psm = psmMap.get(k);
                PeakMatch peakMatch = new PeakMatch(spectrum, jSpectrum, psm);
                int maxFragCharge = 2;
                PeakMatch.setMs2tol(Config.ms2tol);
                peakMatch.initialize(true, maxFragCharge);
                CalcFeature calcFeature = new CalcFeature(jSpectrum);
                calcFeature.addMatch(peakMatch.getIonMatches());
                JSpectrum js = calcFeature.getSpectrum();

                /*moudle2 treatment*/
                js.sortPeaksByMZ();
                js.module2(isoReduction, chargeDeconv, ionsMerge, largerThanPrecursor);
                js.sortPeaksByMZ();
                for (int i = 0; i < js.getPeaks().size(); i++) {
                    JPeak jPeak = js.getPeaks().get(i);
                    if (jPeak.isMatch == 1) {
                        peakBuilder.append(jPeak.getMz());
                        peakBuilder.append("\t" + jPeak.getIonType() + "\t");
                        peakBuilder.append(jPeak.getIntensity());
                        peakBuilder.append("\n");
                    } else {
                        peakBuilder.append(jPeak.getMz());
                        peakBuilder.append("\tNA\t");
                        peakBuilder.append(jPeak.getIntensity());
                        peakBuilder.append("\n");
                    }
                }
            } else {
                for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
                    JPeak jPeak = jSpectrum.getPeaks().get(i);
                    peakBuilder.append(jPeak.getMz());
                    peakBuilder.append("\tNA\t");
                    peakBuilder.append(jPeak.getIntensity());
                    peakBuilder.append("\n");
                }
            }
            bWriter.write(peakBuilder.toString() + "\n");
            bWriter.close();

            /*individual and original MS/MS spectrum*/
            /*String outMS1 = outdir + "/" + outprefix + "_ms2-raw.mgf";
            BufferedWriter msWriter1 = new BufferedWriter(new FileWriter(new File(outMS1)));
            msWriter1.write(spectrum.asMgf() + "\n");
            msWriter1.close();*/

            for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
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
            for (JPeakPair jpp : jPeakPairs) {
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

    private static void doParseDat123(String mgf, String outdir, String outlog, Boolean imonFilter, String labelMethod, Boolean repFilter, Boolean labelFilter, Boolean lowWinFilter, Boolean highWinFilter, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMerge, Boolean largerThanPrecursor, String dat) throws IOException, MzMLUnmarshallerException {
        MascotDatparser mascotDatparser = new MascotDatparser();
        ArrayList<SpectrumAnnotationContainer> dpsms = mascotDatparser.fragAnnotation(dat);
        HashMap<String, SpectrumAnnotationContainer> psmMap = new HashMap<String, SpectrumAnnotationContainer>();
        for (SpectrumAnnotationContainer psm : dpsms) {
            psmMap.put(psm.getSpectrumID(), psm); /*spectrumID is spectrumtitle*/
        }

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
            /*construct a Preprocessing.JSpectrum object*/
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

            if (psmMap.containsKey(jSpectrum.getSpectrumTitle())) {
                SpectrumAnnotationContainer psm = psmMap.get(jSpectrum.getSpectrumTitle());
                HashMap<Double, String> frags = psm.getFragannotation();
                for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
                    JPeak jPeak = jSpectrum.getPeaks().get(i);
                    if (frags.containsKey(jPeak.getMz())) {
                        peakBuilder.append(jPeak.getMz());
                        peakBuilder.append("\t" + frags.get(jPeak.getMz()) + "\t");
                        peakBuilder.append(jPeak.getIntensity());
                        peakBuilder.append("\n");
                        System.out.println(jSpectrum.getSpectrumTitle() + "\t" + frags.get(jPeak.getMz()) + "\t" + jPeak.getMz());
                    } else {
                        peakBuilder.append(jPeak.getMz());
                        peakBuilder.append("\tNA\t");
                        peakBuilder.append(jPeak.getIntensity());
                        peakBuilder.append("\n");
                    }
                }

            } else {
                for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
                    JPeak jPeak = jSpectrum.getPeaks().get(i);
                    peakBuilder.append(jPeak.getMz());
                    peakBuilder.append("\tNA\t");
                    peakBuilder.append(jPeak.getIntensity());
                    peakBuilder.append("\n");
                }
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

    private static void doParseDat23(String mgf, String outdir, String outlog, Boolean imonFilter, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMerge, Boolean largerThanPrecursor, String dat) throws IOException, MzMLUnmarshallerException {
        MascotDatparser mascotDatparser = new MascotDatparser();
        ArrayList<SpectrumAnnotationContainer> dpsms = mascotDatparser.fragAnnotation(dat);
        HashMap<String, SpectrumAnnotationContainer> psmMap = new HashMap<String, SpectrumAnnotationContainer>();
        for (SpectrumAnnotationContainer psm : dpsms) {
            psmMap.put(psm.getSpectrumID(), psm); /*spectrumID is spectrumtitle*/
        }

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
            /*construct a Preprocessing.JSpectrum object*/
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

            /*moudle2 treatment*/
            jSpectrum.sortPeaksByMZ();
            jSpectrum.module2(isoReduction, chargeDeconv, ionsMerge, largerThanPrecursor);

            /*module3 treatment*/
            jSpectrum.sortPeaksByMZ();
            StringBuilder peakBuilder = new StringBuilder();

            if (psmMap.containsKey(jSpectrum.getSpectrumTitle())) {
                SpectrumAnnotationContainer psm = psmMap.get(jSpectrum.getSpectrumTitle());
                HashMap<Double, String> frags = psm.getFragannotation();
                for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
                    JPeak jPeak = jSpectrum.getPeaks().get(i);
                    if (frags.containsKey(jPeak.getMz())) {
                        peakBuilder.append(jPeak.getMz());
                        peakBuilder.append("\t" + frags.get(jPeak.getMz()) + "\t");
                        peakBuilder.append(jPeak.getIntensity());
                        peakBuilder.append("\n");
                        System.out.println(jSpectrum.getSpectrumTitle() + "\t" + frags.get(jPeak.getMz()) + "\t" + jPeak.getMz());
                    } else {
                        peakBuilder.append(jPeak.getMz());
                        peakBuilder.append("\tNA\t");
                        peakBuilder.append(jPeak.getIntensity());
                        peakBuilder.append("\n");
                    }
                }

            } else {
                for (int i = 0; i < jSpectrum.getPeaks().size(); i++) {
                    JPeak jPeak = jSpectrum.getPeaks().get(i);
                    peakBuilder.append(jPeak.getMz());
                    peakBuilder.append("\tNA\t");
                    peakBuilder.append(jPeak.getIntensity());
                    peakBuilder.append("\n");
                }
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


    private static void doPrecursorIsotopesTest(String mgf, String outdir, Boolean b) throws IOException, MzMLUnmarshallerException {
        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        spectrumFactory.addSpectra(mgfFile, null);


        String outmgf = outdir + "/precursorIsotope_" + mgfFile.getName();
        BufferedWriter bwmgf = new BufferedWriter(new FileWriter(new File(outmgf)));


        /*parsing MS/MS spectrum one by one*/
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());
        for (int k=0;k<tList.size();k++) {
            String outprefix = "spectrum" + k;
            MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(mgfFile.getName(),tList.get(k));
            /*construct a Preprocessing.JSpectrum object*/
            JSpectrum jSpectrum = new JSpectrum();
            Integer ch = spectrum.getPrecursor().getPossibleCharges().get(0).value;
            jSpectrum.setParentMass(spectrum.getPrecursor().getMassPlusProton(ch));
            jSpectrum.setParentMassToCharge(spectrum.getPrecursor().getMz());
            jSpectrum.setCharge(ch);
            jSpectrum.setSpectrumTitle(spectrum.getSpectrumTitle());
            jSpectrum.setIntensity(spectrum.getPrecursor().getIntensity());
            jSpectrum.setRt(spectrum.getPrecursor().getRt());
            for (Peak p : spectrum.getPeakList()) {
                if (b) {
                    Boolean q = true;
                    for (double preciso : jSpectrum.predictPrecEnvelope()) {
                        if (Config.delta(p.getMz(), preciso) < Config.ms2tol) {
                            System.out.println(jSpectrum.getSpectrumTitle() + "\t" + preciso);
                            q = false;
                            break;
                        }
                    }
                    if (q) {
                        JPeak jPeak = new JPeak(p.getMz(), p.getIntensity());
                        jSpectrum.addRawPeak(jPeak);
                    }

                } else {
                    JPeak jPeak = new JPeak(p.getMz(), p.getIntensity());
                    jSpectrum.addRawPeak(jPeak);
                }
            }
            jSpectrum.resetPeaks();

            /*module1 treatment*/
            jSpectrum.sortPeaksByMZ();
            StringBuilder sb = jSpectrum.toMgf();
            bwmgf.write(sb.toString());
        }
        bwmgf.close();
    }



}
