package Preprocessing;

import org.apache.commons.cli.*;

import java.io.File;

/**
 * Created by dengyamei on 06/08/2018.
 */
public class ParameterManager {
    private String tool;
    private String version;
    private String releasedate;
    private Options options;

    public ParameterManager() {
    }

    /*
    * main constructor
    * */
    public ParameterManager(String tool, String version, String releasedate) {
        this.tool = tool;
        this.version = version;
        this.releasedate = releasedate;
        ParameterOptions();
    }

    public void printToolInfo() {
        System.out.println(this.tool + " " + this.version + " (" + this.releasedate + ")");
    }

    public void printUsageInfo() {
        printToolInfo();
        System.out.println();
        HelpFormatter help = new HelpFormatter();
        System.out.println("java -Xmx2G pClean.jar");
        help.printHelp("Options", this.options);
        System.out.println();
    }

    public void ParameterOptions() {
        Options options = new Options();
        options.addOption("h", false, "Help info");
        options.addOption("i", true, "MS/MS data in MGF format");
        options.addOption("idres", true, "Identification result");
        options.addOption("o", true, "Output directory");
        options.addOption("itol", true, "Fragment tolerance");
        options.addOption("a2", false, "Consider gap masses of two amino acids");
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
        options.addOption("network", false, "Graph-based network filtration");
        this.options = options;
    }

    public Options getOptions() {
        return options;
    }

    /*
    * initialize parameters
    * */
    private boolean a2 = true;
    private String idres = null;
    private String input = null;
    private String outdir = null;
    private double itol;
    private boolean mionFilter;
    private String labelMethod = null;
    private boolean repFilter;
    private boolean labelFilter;
    private boolean low;
    private boolean high;
    private boolean isoReduction;
    private boolean chargeDeconv;
    private boolean ionsMerge;
    private boolean largerThanPrecursor;
    private boolean network;

    public void InitializeParameters(CommandLine cmd) {
        if (cmd.hasOption("a2")) {
            DeltaMassDB.consider2aa = true;
        } else {
            DeltaMassDB.consider2aa = false;
        }
        this.a2 = DeltaMassDB.consider2aa;

        if (cmd.hasOption("idres")) {
            this.idres = cmd.getOptionValue("idres");
        }

        if (cmd.hasOption("i")) {
            this.input = cmd.getOptionValue("i");
        }

        if (cmd.hasOption("o")) {
            this.outdir = cmd.getOptionValue("o");
            File OD = new File(outdir);
            if (OD.isDirectory() && OD.exists()) {

            } else {
                OD.mkdirs();
            }
        }

        if (cmd.hasOption("itol")) {
            Config.ms2tol = Double.valueOf(cmd.getOptionValue("itol"));
        } else {
            Config.ms2tol = 0.05;
        }
        this.itol = Config.ms2tol;
        this.mionFilter = cmd.hasOption("mionFilter");

        if (cmd.hasOption("labelMethod")) {
            this.labelMethod = cmd.getOptionValue("labelMethod");
        }
        this.repFilter = cmd.hasOption("repFilter");
        this.labelFilter = cmd.hasOption("labelFilter");
        this.low = cmd.hasOption("low");
        this.high = cmd.hasOption("high");
        this.isoReduction = cmd.hasOption("isoReduction");
        this.chargeDeconv = cmd.hasOption("chargeDeconv");
        this.ionsMerge = cmd.hasOption("ionsMerge");
        this.largerThanPrecursor = cmd.hasOption("largerThanPrecursor");
        this.network = cmd.hasOption("network");
    }

    public boolean considerA2() {
        return a2;
    }

    public String getIdres() {
        return idres;
    }

    public String getInput() {
        return input;
    }

    public String getOutdir() {
        return outdir;
    }

    public double getItol() {
        return itol;
    }

    public boolean doImionFilter() {
        return mionFilter;
    }

    public String getLabelMethod() {
        return labelMethod;
    }

    public boolean doRepFilter() {
        return repFilter;
    }

    public boolean doLabelFilter() {
        return labelFilter;
    }

    public boolean doLowWindowFilter() {
        return low;
    }

    public boolean doHighWindowFilter() {
        return high;
    }

    public boolean doIsoReduction() {
        return isoReduction;
    }

    public boolean doChargeDeconv() {
        return chargeDeconv;
    }

    public boolean doIonsMerge() {
        return ionsMerge;
    }

    public boolean doLargerThanPrecursor() {
        return largerThanPrecursor;
    }

    public boolean doNetworkFilter() {
        return network;
    }
}
