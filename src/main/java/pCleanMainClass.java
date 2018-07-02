import org.apache.commons.cli.*;

import java.io.File;

/**
 * Copyright (C) Yamei Deng and Bo Wen
 *
 * Created by dengyamei on 30/06/2018.
 */


public class pCleanMainClass {
    public static void main(String[] args) throws ParseException {
        Options options = new Options();
        options.addOption("i", true, "MS/MS data in MGF format");
        options.addOption("o", true, "Output directory");
        options.addOption("itol", true, "Fragment tolerance");
        options.addOption("m", true, "mzIdentML file");
        options.addOption("label", true, "Peptide labeling method, iTRAQ4plex, iTRAQ8plex or TMT6, TMT12");
        options.addOption("rep+", true, "True: remove reporter ions; False: keep reporter ions");
        options.addOption("low", true, "True: removal of low b-/y-free window");
        options.addOption("high", true, "True: removal of low b-/y-free window");
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

        String label = null;
        if (cmd.hasOption("label")) {
            label = cmd.getOptionValue("label");
        }


        /*Default treatments for MS/MS spectra before applying a module*/
        JSpectrum.setImmoniumIons();

        String outlog = outdir+"/spectrumInfor.txt";
        if (cmd.hasOption("m")) {
            String mzID = cmd.getOptionValue("m");

        } else {
//            runFindPeakPair(mgf, outdir, outlog, label);
        }







    }

}
