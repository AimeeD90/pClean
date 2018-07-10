package Preprocessing;

/**
 * Created by dengyamei on 03/07/2018.
 */
public class PMass {
    //public double deltaMass = 0;
    public double deltaMZ = 0;
    public String name = "";

    /*store match error*/
    public double mzTol = 0.0;

    public PMass(){

    }

    public PMass(double m1) {
        this.deltaMZ = m1;
    }
}
