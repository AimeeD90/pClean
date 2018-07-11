package MZIDparser;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by dengyamei on 10/07/2018.
 */
public class PAminoAcid {
    //public HashSet<Long> aa = new HashSet<Long>();
    public HashMap<Long, Double> aa =  new HashMap<Long, Double>();

    public PAminoAcid(){
        //aa.add(Integer.valueOf("71.037114"));
        aa.put(Math.round(71.037114),71.037114); //A
        aa.put(Math.round(156.101111),156.101111); //R
        aa.put(Math.round(114.042927),114.042927); //N
        aa.put(Math.round(115.026943),115.026943); //D
        aa.put(Math.round(103.009185),103.009185); //C
        aa.put(Math.round(129.042593),129.042593); //E
        aa.put(Math.round(128.058578),128.058578); //Q
        aa.put(Math.round(57.021464),57.021464); //G
        aa.put(Math.round(137.058912),137.058912); //H
        aa.put(Math.round(113.084064),113.084064); //I
        aa.put(Math.round(113.084064),113.084064); //L
        aa.put(Math.round(128.094963),128.094963); //K
        aa.put(Math.round(131.040485),131.040485); //M
        aa.put(Math.round(147.068414),147.068414); //F
        aa.put(Math.round(97.052764),97.052764); //P
        aa.put(Math.round(87.032028),87.032028); //S
        aa.put(Math.round(101.047679),101.047679); //T
        aa.put(Math.round(186.079313),186.079313); //W
        aa.put(Math.round(163.063329),163.063329); //Y
        aa.put(Math.round(99.068414),99.068414); //V

        //HashMap<Long, Double> aatmp =  new HashMap<Long, Double>();

        ArrayList<Long> masses  = new ArrayList<Long>();
        for(long mass :aa.keySet()){
            masses.add(mass);
        }

        for(long mass:masses){

            aa.put(mass-17,aa.get(mass)-17);
            aa.put(mass-18,aa.get(mass)-18);
        }


    }
}
