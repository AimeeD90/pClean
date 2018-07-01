/**
 * Object for each fragment ion
 * <p>
 * Created by dengyamei on 30/06/2018.
 */
public class JPeak implements Cloneable {
    private double mz;
    private double intensity;
    private int charge = 0;
    private String ionType = "-";
    private String ID = ""; // the unique identifier of this peak
    private boolean isotopePrediction = false;
    public double isMatch = -1;
    public boolean remove = false;

    public boolean isIsotopePrediction() {
        return isotopePrediction;
    }

    public void setIsotopePrediction(boolean isotopePrediction) {
        this.isotopePrediction = isotopePrediction;
    }

    public String getID() {
        return ID;
    }

    public void setID(String ID) {
        this.ID = ID;
    }

    public double getMz() {
        return mz;
    }

    public void setMz(double mz) {
        this.mz = mz;
    }

    public double getIntensity() {
        return intensity;
    }

    public void setIntensity(double intensity) {
        this.intensity = intensity;
    }

    public int getCharge() {
        return charge;
    }

    public void setCharge(int charge) {
        this.charge = charge;
    }

    public String getIonType() {
        return ionType;
    }

    public void setIonType(String ionType) {
        this.ionType = ionType;
    }

    public JPeak clone() {
        JPeak jPeak = null;
        try {
            jPeak = (JPeak) super.clone();
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();
        }
        return jPeak;
    }

    /*
    * Empty constructor
    * */
    public JPeak() {
    }

    /*
    * Main constructor:
    * @param mz
    * @param intensity
    * */
    public JPeak(double mz, double intensity) {
        this.mz = mz;
        this.intensity = intensity;
    }

    /*
    * Main constructor:
    * @param JPeak
    * */
    public JPeak(JPeak jPeak) {
        this.mz = jPeak.getMz();
        this.intensity = jPeak.getIntensity();
    }

}

/*
* The pending problem:
* public double isMatch = -1; public boolean remove = false;
* these two fields have not be set with getter or setter.
* */
