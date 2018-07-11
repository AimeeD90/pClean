package MZIDparser;

/**
 * Created by dengyamei on 10/07/2018.
 */
public class JModification {
    //member variable
    private int modLocation; //0-based
    private String residue;
    private double modMassDelta;

    // constructor
    public JModification(){

    }

    // setter and getter
    public int getModLocation() {
        return modLocation;
    }
    public void setModLocation(int modLocation) {
        this.modLocation = modLocation;
    }
    public String getResidue() {
        return residue;
    }
    public void setResidue(String residue) {
        this.residue = residue;
    }
    public double getModMassDelta() {
        return modMassDelta;
    }
    public void setModMassDelta(double modMassDelta) {
        this.modMassDelta = modMassDelta;
    }
}
