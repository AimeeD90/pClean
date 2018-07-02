/**
 * Config.java provides with three transformation functions for mass calculation:
 * 1) from Da to ppm,
 * 2) from ppm to Da,
 * 3) to calculate the delta mass of two peaks
 *
 * Created by dengyamei on 30/06/2018.
 */
public class Config {
    /*the fragment tolerance value*/
    public static double ms2tol;

    /*convert the fragment tolerance value for a specific peak from ppm to Da*/
    public static double ppm2mz(double mz, double ppmTolerence) {
        return mz * ppmTolerence / 1000000;
    }

    /*convert the fragment tolerance value for a specific peak from Da to ppm*/
    public static double mz2ppm(double mz, double mzTolerence) {
        return mzTolerence * 1000000 / mz;
    }

    /*calculate the mass difference in ppm for two peaks, and the reference peak is jPeak 1*/
    public static double deltaPPM(JPeak jPeak1, JPeak jPeak2) {
        double delta = 0;
        try {
            delta = Math.abs((jPeak1.getMz() - jPeak2.getMz()) / jPeak1.getMz() * 1000000);
        } catch (NullPointerException b) {
            System.out.println("Exception");
        } catch (Exception b) {
            System.out.println("capture exception");
        }
        return delta;
    }

    /*calculate the mass difference in ppm for two mass, and the reference mass is mass 1*/
    public static double deltaPPM(Double mass1, Double mass2) {
        double delta = 0;
        try {
            delta = Math.abs((mass1 - mass2) / mass1 * 1000000);
        } catch (NullPointerException b) {
            System.out.println("Exception");
        } catch (Exception b) {
            System.out.println("capture exception");
        }
        return delta;
    }

}
