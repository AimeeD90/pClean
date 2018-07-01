import java.util.Comparator;

/**
 * Sort a JPeak-objects arrayList according to the peak intensity, from the minimum to the maximum
 *
 * Created by dengyamei on 30/06/2018.
 */
public class PeakIntensityComparator implements Comparator<JPeak> {
    public int compare(JPeak p1, JPeak p2) {
        double  i1 = p1.getIntensity();
        double  i2 = p2.getIntensity();
        if(i1>i2){
            return 1;
        }else if(i1==i2){
            return 0;
        }else{
            return -1;
        }
    }

}
