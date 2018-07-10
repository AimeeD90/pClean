import Preprocessing.JPeak;
import Preprocessing.JSpectrum;
import com.compomics.util.experiment.massspectrometry.MSnSpectrum;
import com.compomics.util.experiment.massspectrometry.Peak;
import com.compomics.util.experiment.massspectrometry.SpectrumFactory;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshallerException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by dengyamei on 01/07/2018.
 */
public class MGFparser {
    public static void main(String[] args) throws IOException, MzMLUnmarshallerException {
        String mgf = args[0];
        SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(mgf);
        spectrumFactory.addSpectra(mgfFile, null);
        ArrayList<String> tList = spectrumFactory.getSpectrumTitles(mgfFile.getName());
        for (int k=0;k<tList.size();k++) {
            String outprefix = "spectrum" + k;
            MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(mgfFile.getName(), tList.get(k));
            JSpectrum jSpectrum = new JSpectrum();
            /*set values for a Preprocessing.JSpectrum object*/
            Integer ch = spectrum.getPrecursor().getPossibleCharges().get(0).value;
            jSpectrum.setParentMass(spectrum.getPrecursor().getMassPlusProton(ch));
            jSpectrum.setParentMassToCharge(spectrum.getPrecursor().getMz());
            jSpectrum.setCharge(ch);
            jSpectrum.setSpectrumTitle(spectrum.getSpectrumTitle());
            jSpectrum.setIntensity(spectrum.getPrecursor().getIntensity());
            jSpectrum.setRt(spectrum.getPrecursor().getRt());
            for (Peak p : spectrum.getPeakList()) {
                JPeak jPeak = new JPeak(p.getMz(), p.getIntensity());
                jSpectrum.addRawPeak(jPeak);
            }
            jSpectrum.setSpectrumTitle("" + k);
            try {
                jSpectrum.keepTopNPeaksInRegions(5, 100);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }
}
