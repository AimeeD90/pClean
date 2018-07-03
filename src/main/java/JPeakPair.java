/**
 * Format output information in _edge.txt files.
 *
 * Created by dengyamei on 03/07/2018.
 */
public class JPeakPair {
    public double from;
    public double to;
    public int charge1;
    public int charge2;
    public double mass1;
    public double mass2;
    public double delta;
    public double mzdelta;
    public String deltaName;
    public double intensity;

    public JPeakPair() {
    }

    public void setFrom(String value) {
        this.from = Double.valueOf(value);
    }

    public void setFrom(double value) {
        this.from = value;
    }

    public void setTo(String value) {
        this.to = Double.valueOf(value);
    }

    public void setTo(double value) {
        this.to = value;
    }

    public String print(){
        String line = from + "\t" + to + "\t" + mass1 + "\t" + mass2 + "\t" + charge1 + "\t" + charge2 + "\t" + mzdelta + "\t" + delta + "\t" + deltaName + "\t" + intensity;
        return (line);
    }
}
