/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class Variant extends Base {//implements Comparable<Object> {

//    private int position;
//    private char base;
    private char referenceBase;
    private static final char FIELDS_SEPERATOR = '\t';
    private static final int DIPLOID_COVERAGE_THRESHOLD = 3;
    private static final int POLYPLOID_BASE_COVERAGE_THRESHOLD = 2;
    private static final int DIPLOID_BASE_COVERAGE_THRESHOLD = 3;
    private static final double POLYPLOID_BASE_PROPORTION_THRESHOLD = 0.05;
    private static final double DIPLOID_BASE_PROPORTION_THRESHOLD = 0.30;

    public Variant() {
        super();
        referenceBase = DNA.BASE_BLANK;
    }

    public Variant(int position, char base, char refBase) {
        super(position, base);
        this.referenceBase = refBase;
    }

    public char getReferenceBase() {
        return referenceBase;
    }

    public void setReferenceBase(char referenceBase) {
        this.referenceBase = referenceBase;
    }

    public static HashMap initialiseVariantList(SAMSequenceDictionary sequenceDictionary, boolean orderingRequired) {
        HashMap variantList = new HashMap();

        Iterator<SAMSequenceRecord> sequenceIterator = sequenceDictionary.getSequences().iterator();
        while (sequenceIterator.hasNext()) {
            // get the next sequence
            SAMSequenceRecord sequenceEntry = sequenceIterator.next();
            // Initialise a new variant list
            variantList.put(sequenceEntry.getSequenceName(), (orderingRequired ? new LinkedHashMap<Integer, Variant>() : new HashMap<Integer, Variant>()));
        } // end while

        return variantList;
    }

    public static int readVariantList(String vcfFile, HashMap<String, HashMap<Integer, Variant>> variantList, HashMap<String, LinkedHashMap<Integer, Variant>> referenceBaseList) {

        if (vcfFile == null || vcfFile.isEmpty()) {
            return 0;
        }

        int count = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(vcfFile));
            String previousSequenceName = "";
            HashMap<Integer, Variant> seqVariantList = null;
            HashMap<Integer, Variant> seqReferenceBaseList = null;

            String strFieldSeperator = Character.toString(FIELDS_SEPERATOR);
            String line;
            while ((line = br.readLine()) != null) {
                // ignore empty lines
                if (line.isEmpty()) {
                    continue;
                }
                if (line.startsWith("#")) {
                    continue;
                }

                String[] entries = line.split(strFieldSeperator);

                if (entries.length < 5) {
                    System.err.println("Incorrect VCF format");
                    return -1;
                }

                // ignore non variant positions
                if (entries[4].equals(".")) {
                    continue;
                }

                // ignore indels
                if ((entries.length > 7 && entries[7].startsWith("INDEL")) || (entries[3].length() != entries[4].length() && !entries[4].contains(","))) {
                    continue;
                }
                // also ignore rows with more than one base in the reference
                if (entries[3].length() > 1) {
                    continue;
                }

                String sequenceName = entries[0];
                int position = Integer.parseInt(entries[1]);
                String consensusBase;
                if (entries.length > 8 && entries[8].startsWith("GT")) {

                    // accept the alternate base if "N" in the reference, or reference base is not part of the consensus                    
                    TreeSet bases = new TreeSet();
                    if (entries[3].equals("N") || entries[9].startsWith("1")) {
                        // don't do anything extra. only need to add alternate base(s)
                    } else if (entries[9].startsWith("0")) { // reference base is to included in the consensus
                        // insert the reference base
                        bases.add(entries[3].charAt(0));
                    }

                    // get the bases corresponding to the alternate base (alternate base could be multiple characters seperated by commas)
                    String[] altBases = entries[4].split(",");
                    for (String altBase : altBases) {
                        bases.add(altBase.charAt(0));
                    }

                    // get the consensus base
                    consensusBase = ((Character) DNA.getReverseNucleotideMap().get(bases.toString())).toString();
                } else {
                    consensusBase = entries[4];
                }

                if (!sequenceName.equals(previousSequenceName)) {
                    if (!previousSequenceName.isEmpty()) {
                        variantList.put(previousSequenceName, seqVariantList);
                        //referenceBaseList.put(previousSequenceName, seqReferenceBaseList);
                    }
                    seqVariantList = (HashMap<Integer, Variant>) variantList.get(sequenceName);
                    if (referenceBaseList != null) {
                        seqReferenceBaseList = (HashMap<Integer, Variant>) referenceBaseList.get(sequenceName);
                    }
                    if (seqVariantList == null || (referenceBaseList != null && seqReferenceBaseList == null)) {
                        System.err.println("Invalid reference name in the VCF File");
                        return -1;
                    }
                    previousSequenceName = sequenceName;
                }

                // set the variant base
                seqVariantList.put(position, new Variant(position, consensusBase.charAt(0), entries[3].charAt(0)));
                // also set the reference base (if provided) {
                if (seqReferenceBaseList != null) {
                    seqReferenceBaseList.put(position, new Variant(position, DNA.BASE_BLANK, entries[3].charAt(0)));
                }

                count++;

            } // end while

            variantList.put(previousSequenceName, seqVariantList);

            br.close();

        } catch (FileNotFoundException ex) {
            System.err.println("File: " + vcfFile + " not found!");
            return -1;
            //Logger.getLogger(SNP.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Variant.class.getName()).log(Level.SEVERE, null, ex);
            return -1;
        }

        return count;
    }

    public static int checkPolyploidHSPList(HashMap<String, HashMap<Integer, Variant>> polyploidHSPList,
            HashMap<String, BaseCoverage[]> polyploidBaseCoverageList, HashMap<String, LinkedHashMap<Integer, Variant>> positionList) {

        int count = 0;
        Iterator<Map.Entry<String, LinkedHashMap<Integer, Variant>>> itPositionList = positionList.entrySet().iterator();
        while (itPositionList.hasNext()) {
            // get the next entry
            Map.Entry<String, LinkedHashMap<Integer, Variant>> entry = itPositionList.next();

            // get sequence name
            String sequenceName = entry.getKey();
            // get snp list for this sequence
            HashMap<Integer, Variant> seqPositionList = entry.getValue();

            // and the reference base list for this sequence
            HashMap<Integer, Variant> seqPolyPloidHSPList = polyploidHSPList.get(sequenceName);
            // and the base distribution list
            BaseCoverage[] seqBaseCoverage = polyploidBaseCoverageList.get(sequenceName);

            Iterator<Map.Entry<Integer, Variant>> itPosition = seqPositionList.entrySet().iterator();
            while (itPosition.hasNext()) {
                // get the next entry
                Map.Entry<Integer, Variant> positionEntry = itPosition.next();
                // get the position and corresponding referencey base from the position list
                int position = positionEntry.getKey();
                Variant referenceBase = positionEntry.getValue();

                // get the HSP base at this position
                Variant hsp = seqPolyPloidHSPList.get(position);
                if (hsp == null) { // hsp was not called at this position. create a new hsp 
                    hsp = new Variant(position, DNA.BASE_BLANK, referenceBase.getReferenceBase());
                }

                // get the base distribution for this position
                BaseCoverage baseCoverage = seqBaseCoverage[position - 1]; //AMT-CHECK: 0-based indexing in base coverage list
                if (baseCoverage == null) { // base-coverage now found for this position. ignore
                    continue;
                }
                int coverage = baseCoverage.getCoverage();

                TreeSet newBases = new TreeSet();
                if (checkPolyploidBase(baseCoverage.getCoverageA(), coverage)) {
                    newBases.add('A');
                }
                if (checkPolyploidBase(baseCoverage.getCoverageC(), coverage)) {
                    newBases.add('C');
                }
                if (checkPolyploidBase(baseCoverage.getCoverageG(), coverage)) {
                    newBases.add('G');
                }
                if (checkPolyploidBase(baseCoverage.getCoverageT(), coverage)) {
                    newBases.add('T');
                }

                if (!newBases.isEmpty()) {
                    char newBase = (Character) DNA.getReverseNucleotideMap().get(newBases.toString());
                    if (hsp.getNucleotide() != newBase) {
                        hsp.setNucleotide(newBase);
                        seqPolyPloidHSPList.put(position, hsp); // add to the polyploid hsp list
                    }
                    if (DNA.isExtendedNucleotide(newBase)) {
                        count++;
                    }
                } else { // no valid base found after checking so remove the HSP
                    seqPolyPloidHSPList.remove(position);
                }

            } // end for each position

        }

        return count;
    }

    public static void checkDiploidVariantList(HashMap<String, HashMap<Integer, Variant>> diploidSNPList,
            HashMap<String, BaseCoverage[]> diploidBaseDistributionList, HashMap<String, HashMap<Integer, Variant>> polyploidHSPList,
            HashMap<String, LinkedHashMap<Integer, Variant>> positionList) {

        Iterator<Map.Entry<String, LinkedHashMap<Integer, Variant>>> itPositionList = positionList.entrySet().iterator();
        while (itPositionList.hasNext()) {
            // get the next entry
            Map.Entry<String, LinkedHashMap<Integer, Variant>> entry = itPositionList.next();

            // get sequence name
            String sequenceName = entry.getKey();
            // get snp list for this sequence
            HashMap<Integer, Variant> seqPositionList = entry.getValue();

            // also get diploid snp list for this sequence
            HashMap<Integer, Variant> seqDiploidSNPList = diploidSNPList.get(sequenceName);
            // and the reference base list for this sequence
            HashMap<Integer, Variant> seqPolyploidHSPList = polyploidHSPList.get(sequenceName);
            // and the base distribution list
            BaseCoverage[] seqBaseCoverage = diploidBaseDistributionList.get(sequenceName);

            Iterator<Map.Entry<Integer, Variant>> itPosition = seqPositionList.entrySet().iterator();
            while (itPosition.hasNext()) {
                // get the next entry
                Map.Entry<Integer, Variant> positionEntry = itPosition.next();
                // get the position and corresponding referencey base from the position list
                int position = positionEntry.getKey();
                Variant referenceBase = positionEntry.getValue();

                // get the HSP at this position
                Variant hsp = seqPolyploidHSPList.get(position);
                if (hsp == null) {// no HSP in polyploid at this position, move to the next position
                    continue;
                }

                // get the base in diploid
                Variant diploidVariant = seqDiploidSNPList.get(position);

                if (diploidVariant != null) { // SNP has been called at this position
                    // check if its an ambiguous base
                    if (DNA.isExtendedNucleotide(diploidVariant.getNucleotide())) {
                        diploidVariant.setNucleotide(DNA.BASE_AMBIGUOUS);
                    }
                    //move to the next position
                    continue;
                }

                // create a new diploid variant
                diploidVariant = new Variant(position, DNA.BASE_BLANK, referenceBase.getReferenceBase());
                // we are now checking the position which is not called a SNP in the diploid.
                // check if multiple bases are present at this position (possibly, heterozygous)
                // get the base distribution for this position
                BaseCoverage baseCoverage = seqBaseCoverage[position - 1]; //AMT-CHECK: 0-based indexing in base coverage list
                if (baseCoverage == null) { // base-coverage now found for this position. ignore
                    continue;
                }
                int coverage = baseCoverage.getCoverage();

                if (coverage == 0) {
                    diploidVariant.setNucleotide(DNA.BASE_ZERO_COVERAGE);
                } else if (coverage < DIPLOID_COVERAGE_THRESHOLD) {
                    diploidVariant.setNucleotide(DNA.BASE_LOW_COVERAGE);
                } else {
                    TreeSet newBases = new TreeSet();
                    if (checkDiploidBase(baseCoverage.getCoverageA(), coverage)) {
                        newBases.add('A');
                    }
                    if (checkDiploidBase(baseCoverage.getCoverageC(), coverage)) {
                        newBases.add('C');
                    }
                    if (checkDiploidBase(baseCoverage.getCoverageG(), coverage)) {
                        newBases.add('G');
                    }
                    if (checkDiploidBase(baseCoverage.getCoverageT(), coverage)) {
                        newBases.add('T');
                    }

                    if (newBases.isEmpty()) {// ideally, this should not happen 
                        diploidVariant.setNucleotide(referenceBase.getReferenceBase());
                    } else if (newBases.size() == 1) {
                        // single base found. set the variant base
                        diploidVariant.setNucleotide((Character) newBases.first());
                    } else if (newBases.size() > 1) { // multiple bases found at this position, mark this position as '*'
                        diploidVariant.setNucleotide(DNA.BASE_AMBIGUOUS);
                    }
                }
                // add the base to the list
                seqDiploidSNPList.put(position, diploidVariant);
            } // end for each position
        }
    }

    private static boolean checkPolyploidBase(int baseCoverage, int totalCoverage) {
        return baseCoverage >= POLYPLOID_BASE_COVERAGE_THRESHOLD && (double) baseCoverage / (double) totalCoverage >= POLYPLOID_BASE_PROPORTION_THRESHOLD;
    }

    private static boolean checkDiploidBase(int baseCoverage, int totalCoverage) {
        return baseCoverage >= DIPLOID_BASE_COVERAGE_THRESHOLD && (double) baseCoverage / (double) totalCoverage >= DIPLOID_BASE_PROPORTION_THRESHOLD;
    }

    public static int readPositionList(String positionFile, HashMap<String, LinkedHashMap<Integer, Variant>> positionList) {
        int count = 0; //number of positions

        try {
            BufferedReader br = new BufferedReader(new FileReader(positionFile));
            String previousSequenceName = "";
            HashMap<Integer, Variant> seqRefBaseList = null;
            String line;
            while ((line = br.readLine()) != null) {
                // ignore empty lines
                if (line.isEmpty()) {
                    continue;
                }

                int pos = 0, end;
                // sequence name 
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                String sequenceName = line.substring(pos, end);
                // position
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                int position = Integer.parseInt(line.substring(pos, end));
                // reference base
                pos = end + 1;
                end = line.indexOf(FIELDS_SEPERATOR, pos);
                String referenceBase;
                if (end > -1) {
                    referenceBase = line.substring(pos, end);
                } else {
                    referenceBase = line.substring(pos);
                }

                if (!sequenceName.equals(previousSequenceName)) {
                    seqRefBaseList = (HashMap<Integer, Variant>) positionList.get(sequenceName);
                    previousSequenceName = sequenceName;
                }

                // set the base
                seqRefBaseList.put(position, new Variant(position, DNA.BASE_BLANK, referenceBase.charAt(0)));
                count++;
            } // end while

        } catch (FileNotFoundException ex) {
            Logger.getLogger(Variant.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Variant.class.getName()).log(Level.SEVERE, null, ex);
        }

        return count;

    }

}
