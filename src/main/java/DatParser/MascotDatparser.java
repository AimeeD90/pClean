package DatParser;

import com.compomics.mascotdatfile.util.interfaces.FragmentIon;
import com.compomics.mascotdatfile.util.mascot.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

/**
 * Created by dengyamei on 10/07/2018.
 */
public class MascotDatparser {
    public MascotDatparser() {
    }

    public ArrayList<dPSM> fragAnnotation(String dat) {

        MascotDatfile mdf = new MascotDatfile(dat);
        PeptideToQueryMap lPeptide2Q = mdf.getPeptideToQueryMap();
        QueryToPeptideMap lQuery2P = mdf.getQueryToPeptideMap();

        ArrayList<dPSM> psms = new ArrayList<dPSM>();
        for (int j = 1; j <= lQuery2P.getNumberOfQueries(); j++) {
            Query lQuery = mdf.getQuery(j);
            PeptideHit ph = lQuery2P.getPeptideHitOfOneQuery(j, 1);
            if (ph != null) {
                PeptideHitAnnotation lPha = new PeptideHitAnnotation(ph.getSequence(), ph.getModifications(), mdf.getMasses(), mdf.getParametersSection(), ph.getIonSeriesFound());
                Vector lMascotMatchedIons = lPha.getMatchedIonsByMascot(lQuery.getPeakList(), ph.getPeaksUsedFromIons1());

                HashMap<Double, String> fragannotation = new HashMap<>();
                for (int k = 0; k < lMascotMatchedIons.size(); k++) {
                    FragmentIon fm = (FragmentIon) lMascotMatchedIons.get(k);
                    String fraglabel = fm.getType() + fm.getNumber();
                    double fragmz = fm.getMZ() + fm.getTheoreticalExperimantalMassError();
                    fragannotation.put(fragmz, fraglabel);
                    //System.out.println(lQuery.getFilename() + fm.getType() + fm.getNumber() + "\t" + fragmz);
                }

                dPSM psm = new dPSM(lQuery.getFilename(),fragannotation);
                psms.add(psm);
            }
        }
        return psms;
    }
}
