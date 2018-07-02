import com.compomics.util.experiment.biology.aminoacids.B;
import com.compomics.util.experiment.massspectrometry.SpectrumFactory;
import com.sun.org.apache.xpath.internal.operations.Bool;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Copyright (C) Yamei Deng and Bo Wen
 *
 * Created by dengyamei on 30/06/2018.
 */


public class pCleanMainClass {
    public static void main(String[] args) throws ParseException, IOException {
        Options options = new Options();
        options.addOption("i", true, "MS/MS data in MGF format");
        options.addOption("o", true, "Output directory");
        options.addOption("itol", true, "Fragment tolerance");
        options.addOption("mion", true, "True: remove immonium ions; False: undo");
        options.addOption("labelMethod", true, "Peptide labeling method, iTRAQ4plex, iTRAQ8plex or TMT6, TMT12");
        options.addOption("rep+", true, "True: remove reporter ions; False: undo");
        options.addOption("label+", true, "True: remove label-associated ions; False: undo");
        options.addOption("low", true, "True: removal of low b-/y-free window; False: undo");
        options.addOption("high", true, "True: removal of low b-/y-free window; False: undo");
        options.addOption("isoReduction", true, "True: reduction of heavy isotopic peaks; False: undo");
        options.addOption("chargeDeconv", true, "True: high charge deconvolution; False: undo");
        options.addOption("ionsMarge", true, "True: marge two adjacent peaks within a mass tolerance of 20ppm; False: undo");
        options.addOption("largerThanPrecursor", true, "True: remove peaks larger than precursor; False: undo");
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
        Boolean imon = cmd.hasOption("mion");

        String labelMethod = null;
        if (cmd.hasOption("label")) {
            labelMethod = cmd.getOptionValue("labelMethod");
        }

        Boolean repFilter = cmd.hasOption("rep+");
        Boolean labelFilter = cmd.hasOption("label+");
        Boolean lowWinFilter = cmd.hasOption("low");
        Boolean highWinFilter = cmd.hasOption("high");
        Boolean isoReduction = cmd.hasOption("isoReduction");
        Boolean chargeDeconv = cmd.hasOption("chargeDeconv");
        Boolean ionsMarge = cmd.hasOption("ionsMarge");
        Boolean largerThanPrecursor = cmd.hasOption("largerThanPrecursor");


        String outlog = outdir+"/spectrumInfor.txt";
        if (cmd.hasOption("m")) {
            String mzID = cmd.getOptionValue("m");

        } else {
            if (labelMethod != null) {
                doPreprocessing123(mgf, outdir, outlog, imon, labelMethod, repFilter, labelFilter, lowWinFilter, highWinFilter, isoReduction, chargeDeconv, ionsMarge, largerThanPrecursor);
            } else {
                doPreprocessing23(mgf, outdir, outlog, imon, isoReduction, chargeDeconv, ionsMarge, largerThanPrecursor);
            }
        }
    }

    /*For high-resolution label-based MS/MS data: implemented module 1, 2, and 3*/
    private static void doPreprocessing123(String mgf, String outdir, String outlog, Boolean imon, String labelMethod, Boolean repFilter, Boolean labelFilter, Boolean lowWinFilter, Boolean highWinFilter, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMarge, Boolean largerThanPrecursor) throws IOException {
        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        spectrumFactory.addSpectra(mgfFile, null);
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());


    }

    /*For high-resolution label-free MS/MS data: implemented module 2 and 3*/
    private static void doPreprocessing23(String mgf, String outdir, String outlog, Boolean imon, Boolean isoReduction, Boolean chargeDeconv, Boolean ionsMarge, Boolean largerThanPrecursor) {

    }

}
