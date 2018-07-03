import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

/**
 * the mass gap database used in pClean
 * we enlarged this database for interaction-calculation.
 * Created by dengyamei on 03/07/2018.
 */
public class DeltaMassDB {
    public String labelMethod;
    public static boolean consider2aa = true;
    public ArrayList<PMass> massDB = new ArrayList<PMass>();
    public HashMap<Integer, Integer> mass2index = new HashMap<Integer, Integer>();
    public HashMap<Integer, Integer> mass2indexUp = new HashMap<Integer, Integer>();
    public HashMap<Double, String> deltaMass2nameMap = new HashMap<Double, String>();

    /*empty constructor*/
    public DeltaMassDB() {
    }

    /*constructor corresponding to specific peptide labeling method*/
    public DeltaMassDB(String labelMethod) {
        this.labelMethod = labelMethod;
    }

    public int IntegerMass(double mass) {
        return (int) (Math.ceil(mass) + 2);
    }

    public int IntegerMassUp(double mass) {
        return (int) (Math.floor(mass) - 2);
    }

    public int IntegerMassForSearch(double mass) {
        return (int) (Math.floor(mass));
    }

    public void init(){
        deltaMass2nameMap.put(71.037114, "A"); //A
        deltaMass2nameMap.put(156.101111, "R"); //R
        deltaMass2nameMap.put(114.042927, "N"); //N
        deltaMass2nameMap.put(114.042927 + 0.984016, "N"); //N Deamidated
        deltaMass2nameMap.put(115.026943, "D"); //D
        deltaMass2nameMap.put(115.026943 + 0.984016, "D"); //D Deamidated
        deltaMass2nameMap.put(103.009185, "C"); //C
        deltaMass2nameMap.put(103.009185 + 57.021464, "C"); //C
        deltaMass2nameMap.put(129.042593, "E"); //E
        deltaMass2nameMap.put(128.058578, "Q"); //Q
        deltaMass2nameMap.put(57.021464, "G"); //G
        deltaMass2nameMap.put(137.058912, "H"); //H
        //deltaMass2nameMap.put( 113.084064,"I"); //I
        deltaMass2nameMap.put(113.084064, "L"); //L
        deltaMass2nameMap.put(128.094963, "K"); //K
        deltaMass2nameMap.put(131.040485, "M"); //M
        deltaMass2nameMap.put(131.040485 + 15.99492, "M"); //M ????????
        deltaMass2nameMap.put(147.068414, "F"); //F
        deltaMass2nameMap.put(97.052764, "P"); //P
        deltaMass2nameMap.put(87.032028, "S"); //S
        deltaMass2nameMap.put(101.047679, "T"); //T
        deltaMass2nameMap.put(186.079313, "W"); //W
        deltaMass2nameMap.put(163.063329, "Y"); //Y
        deltaMass2nameMap.put(99.068414, "V"); //V

        deltaMass2nameMap.put(17.0265, "NH3"); //
        deltaMass2nameMap.put(18.0106, "H2O"); //

        if (labelMethod != null) {
            if (labelMethod.equals("iTRAQ4plex")) {
                deltaMass2nameMap.put(128.094963 + 144.1021, "K"); //K
                //deltaMass2nameMap.put(163.063329 + 144.1021, "Y"); //Y
            }
            if (labelMethod.equals("iTRAQ8plex")) {
                deltaMass2nameMap.put(128.094963 + 304.2054, "K"); //K
                //deltaMass2nameMap.put(163.063329 + 304.2054, "Y"); //Y
            }
            if (labelMethod.equals("TMT6plex") || labelMethod.equals("TMT10plex")) {
                deltaMass2nameMap.put(128.094963 + 229.162932, "k");
            }
        }

        if(consider2aa){
            ArrayList<Double> aaMass = new ArrayList<Double>();
            for (Double mass : deltaMass2nameMap.keySet()) {
                aaMass.add(mass);
            }
            for (int i = 0; i < aaMass.size(); i++) {
                double mass1 = aaMass.get(i);
                for (int j = i; j < aaMass.size(); j++) {
                    double mass2 = aaMass.get(j);
                    double aa2mass = mass1 + mass2;
                    String aa2 = deltaMass2nameMap.get(mass1) + deltaMass2nameMap.get(mass2);
                    if (!deltaMass2nameMap.containsKey(aa2mass)) {
                        deltaMass2nameMap.put(aa2mass, aa2);
                    }
                }
            }
        }

        boolean lossHN = true;
        if (lossHN) {

        }

        for(Double mass:deltaMass2nameMap.keySet()){
            PMass pMass = new PMass(mass);
            pMass.name = deltaMass2nameMap.get(mass);
            massDB.add(pMass);
        }

        Collections.sort(massDB, comparator);
        int maxMass = (int) Math.ceil(massDB.get(massDB.size() - 1).deltaMZ + 1);
        for(int i=0;i<massDB.size();i++){
            int minMass = IntegerMass(massDB.get(i).deltaMZ);
            if(!mass2index.containsKey(minMass)){
                mass2index.put(minMass, i);
            }
        }

        int lastIndex=0;
        for(int i=1;i<=maxMass;i++){
            if (!mass2index.containsKey(i)) {
                mass2index.put(i, lastIndex);
            } else {
                lastIndex = mass2index.get(i);
            }
        }

        for(int i=massDB.size()-1;i>=0;i--){
            int minMass = IntegerMassUp(massDB.get(i).deltaMZ);
            if(!mass2indexUp.containsKey(minMass)){
                mass2indexUp.put(minMass, i);
            }
        }

        lastIndex = massDB.size() - 1;
        for(int i=maxMass;i>=1;i--){
            if (!mass2indexUp.containsKey(i)) {
                mass2indexUp.put(i, lastIndex);
            } else {
                lastIndex = mass2indexUp.get(i);
            }
        }
    }

    public Comparator<PMass> comparator = new Comparator<PMass>() {
        public int compare(PMass s1, PMass s2) {
            if(s2.deltaMZ > s1.deltaMZ){
                return -1;
            } else if (s2.deltaMZ == s1.deltaMZ) {
                return 0;
            } else {
                return 1;
            }

        }
    };

    public ArrayList<PMass> searchDB(double mass, double delta, int unit){
        ArrayList<PMass> result = new ArrayList<PMass>();
        int index = IntegerMassForSearch(mass);
        if(!mass2index.containsKey(index)){
            return result;
        }
        int ind = mass2index.get(index);
        if(!mass2indexUp.containsKey(index)){
            return result;
        }
        int indMax = mass2indexUp.get(index);
        for(int i=ind;i<=indMax;i++){
            if(indMax>(massDB.size()-1)){
                continue;
            }
            PMass pMass = massDB.get(i);
            double mdel = mass - pMass.deltaMZ;
            if(unit ==1){
                mdel = mdel / pMass.deltaMZ * 1000000;
            }
            if(Math.abs(mdel) <= delta){
                /*store mass match error*/
                pMass.mzTol = mdel;
                result.add(pMass);
            }
        }
        return result;
    }
}
