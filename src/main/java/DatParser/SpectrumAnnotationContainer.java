package DatParser;

import java.util.HashMap;

/**
 * Created by dengyamei on 10/07/2018.
 */
public class SpectrumAnnotationContainer {
    private String spectrumID;
    private HashMap<Double, String> fragannotation = new HashMap<>();

    public SpectrumAnnotationContainer(String spectrumID, HashMap<Double, String> fragannotation) {
        this.spectrumID = spectrumID;
        this.fragannotation = fragannotation;
    }

    public SpectrumAnnotationContainer() {
    }

    public String getSpectrumID() {
        return spectrumID;
    }

    public void setSpectrumID(String spectrumID) {
        this.spectrumID = spectrumID;
    }

    public HashMap<Double, String> getFragannotation() {
        return fragannotation;
    }

    public void setFragannotation(HashMap<Double, String> fragannotation) {
        this.fragannotation = fragannotation;
    }
}
