package MZIDparser;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

import javax.xml.bind.JAXBException;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by dengyamei on 10/07/2018.
 */
public class ReadmzID {

    public double qvalue = 0.001;

    public ReadmzID(double qvalue){
        this.qvalue = qvalue;
    }

    public ArrayList<JPSM> readmzIdentML(String mzidFileString) throws JAXBException {
        File mzidFile = new File(mzidFileString);
        if(!mzidFile.isFile()){
            System.err.println("mzIdentML file: "+mzidFileString+" is not exist!");
            System.exit(0);
        }

        ArrayList<JPSM> psmArrayList = new ArrayList<JPSM>();

        MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(mzidFile);
        DataCollection dataCollection=unmarshaller.unmarshal(MzIdentMLElement.DataCollection);
        if(dataCollection==null){
            System.err.println("Failed! DataCollection is null!");
            System.exit(0);
        }

        //check which type of search, auto-decoy or separate decoy
        boolean autoTargetDecoySearch = false;
        //check which type of enzyme
        ArrayList<String> enzyme = new ArrayList<String>();

        AnalysisProtocolCollection apc=unmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
        for(SpectrumIdentificationProtocol sProtocol:apc.getSpectrumIdentificationProtocol()){
            for(UserParam userParam:sProtocol.getAdditionalSearchParams().getUserParam()){
                if(userParam.getName().contentEquals("TargetDecoyApproach")){
                    autoTargetDecoySearch = Boolean.valueOf(userParam.getValue());
                    break;
                }
            }
            for(Enzyme en:sProtocol.getEnzymes().getEnzyme()){
                for(CvParam cvParam:en.getEnzymeName().getCvParam()){
                    enzyme.add(cvParam.getAccession());
                }
            }
        }

        if(autoTargetDecoySearch){
            System.out.println("Auto target-decoy search!");
        }else{
            System.out.println("Must searching with target-decoy mode!");
            System.exit(0);
        }

        // Get the list of SpectrumIdentification elements
        List<SpectrumIdentificationList> sil = dataCollection.getAnalysisData().getSpectrumIdentificationList();
        for (SpectrumIdentificationList sIdentList : sil) {
            for (SpectrumIdentificationResult spectrumIdentResult: sIdentList.getSpectrumIdentificationResult()) {
                // Get the name of SpectrumIdentificationResult
                String spectrumID =  spectrumIdentResult.getSpectrumID();
                //System.out.println("ReadmzID.java spectrumID="+spectrumID); // the iterm stored in mzid are named as index=2093
                //read each rank of a spectrum
                for (SpectrumIdentificationItem spectrumIdentItem: spectrumIdentResult.getSpectrumIdentificationItem()) {
                    // Get the following information for SpectrumIdentificationItem element
                    // String spectrumIdItem = spectrumIdentItem.getId();
                    Double calculatedMassToCharge =  spectrumIdentItem.getCalculatedMassToCharge();
                    Double experimentalMassToCharge = spectrumIdentItem.getExperimentalMassToCharge();
                    int rank = spectrumIdentItem.getRank(); //select rank1 result
                    if(rank>=2){
                        continue;
                    }
                    int charge = spectrumIdentItem.getChargeState();
                    Peptide peptide = unmarshaller.unmarshal(Peptide.class, spectrumIdentItem.getPeptideRef());
                    String peptideSequence = peptide.getPeptideSequence();

                    //construct a JPSM object
                    JPSM psm = new JPSM();
                    //psm.setSpectrumIndex(spectrumIdentItem.getHid());
                    psm.setCalculatedMassToCharge(calculatedMassToCharge);
                    psm.setCharge(charge);
                    psm.setExperimentalMassToCharge(experimentalMassToCharge);
                    psm.setRank(rank);
                    psm.setSpectrumID(spectrumID);

                    String[] tmpIndexArrayList = spectrumID.split("=");
                    psm.setSpectrumIndex(Integer.valueOf(tmpIndexArrayList[tmpIndexArrayList.length-1]));

                    //psm.setSpectrumIndex(Integer.valueOf(spectrumID.split("=")[1]));
                    psm.setPepSeq(peptideSequence);

                    for(PeptideEvidenceRef peptideEvidenceRef:spectrumIdentItem.getPeptideEvidenceRef()){
                        PeptideEvidence pEvidence = unmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                        //System.out.println(peptideEvidenceRef.getPeptideEvidence().getPre());
                        DBSequence dbSequence = unmarshaller.unmarshal(DBSequence.class, pEvidence.getDBSequenceRef());
                        psm.addProteins(dbSequence.getAccession());
                        if(pEvidence.isIsDecoy()){
                            psm.setDecoy(true);
                        }
                    }
                    if(psm.isDecoy()){
                        continue;
                    }

                    for( CvParam cvParam: spectrumIdentItem.getCvParam()){
                        //System.out.println(cvParam.getName()+"\t"+cvParam.getValue());
                        if(cvParam.getAccession().contentEquals("MS:1002054")){ //MS-GF:RawScore false；actually is MS-GF:QValue
                            psm.setQvalue(Double.valueOf(cvParam.getValue()));

                            //System.out.println("ReadID.java MS:1002054 QValue = "+cvParam.getValue());
                        }
                    }

                    if(psm.getQvalue() > this.qvalue){
                        continue;
                    }

                    List <Modification> modificationList = peptide.getModification();
                    for(Modification modObj: modificationList){

                        JModification modification = new JModification();
                        modification.setModLocation(modObj.getLocation());
                        String modaa = "";
                        if(modObj.getLocation()==0){
                            modaa = "N-term";
                        }else if(modObj.getLocation()==(peptideSequence.length()+1)){
                            modaa = "C-term";
                        }else{
                            modaa = String.valueOf(peptideSequence.charAt(modObj.getLocation()-1));
                        }
                        modification.setResidue(modaa);
                        modification.setModMassDelta(modObj.getMonoisotopicMassDelta());
                        psm.addModificationList(modification);

                        //System.out.println(modObj.getAvgMassDelta()+"\t"+modObj.getHid()+"\t"+modObj.getLocation()+"\t"+modObj.getMonoisotopicMassDelta()+"\t"+modObj.getResidues().size()+"\t"+modaa);

                    }
                    if(rank==1){
                        psmArrayList.add(psm);
                    }else{
                        psmArrayList.get(psmArrayList.size()-1).addPSMRank(psm);
                    }

                } // end spectrum identification item
            } // end spectrum identification results
        }

        return(psmArrayList);
    }
}
