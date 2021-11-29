/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class BaseCoverage {

    // constants
    private static final char FIELDS_SEPERATOR = '\t';
    // variables
    private int coverageA;
    private int coverageC;
    private int coverageG;
    private int coverageT;
    private int coverageN;

    public BaseCoverage() {
        this.coverageA = 0;
        this.coverageC = 0;
        this.coverageG = 0;
        this.coverageT = 0;
        this.coverageN = 0;
    }

    public BaseCoverage(int coverageA, int coverageT, int coverageG, int coverageC, int coverageN) {
        this.coverageA = coverageA;
        this.coverageC = coverageC;
        this.coverageG = coverageG;
        this.coverageT = coverageT;
        this.coverageN = coverageN;
    }

    public int getCoverage() {
        return getCoverageA() + getCoverageC() + getCoverageG() + getCoverageT() + getCoverageN();
    }

    public void setCoverage(int coverageA, int coverageT, int coverageG, int coverageC, int coverageN) {
        this.coverageA = coverageA;
        this.coverageC = coverageC;
        this.coverageG = coverageG;
        this.coverageT = coverageT;
        this.coverageN = coverageN;
    }

    @Override
    public String toString() {
        return Integer.toString(coverageA) + "\t" + Integer.toString(coverageT) + "\t" + Integer.toString(coverageG) + "\t" + Integer.toString(coverageC) + "\t" + Integer.toString(coverageN);
    }

    public static LinkedHashMap<String, BaseCoverage[]> read(String baseCoverageFile, SAMSequenceDictionary sequenceDictionary) {

        LinkedHashMap<String, BaseCoverage[]> baseCoverageList = initialise(sequenceDictionary);

        try {
            BufferedReader br = new BufferedReader(new FileReader(baseCoverageFile));

            String previousSeqName = "";
            BaseCoverage[] seqBaseCoverageList = null;
            String line;
            while ((line = br.readLine()) != null) {
                // ignore empty lines
                if (line.isEmpty()) {
                } else if (line.startsWith("@")) {
                    // ignore header
                } else {

                    int pos = 0, end;
                    // sequence name 
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    String sequenceName = line.substring(pos, end);
                    // position
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer position = Integer.parseInt(line.substring(pos, end));
                    // coverage A
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer coverageA = Integer.parseInt(line.substring(pos, end));
                    // coverage T
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer coverageT = Integer.parseInt(line.substring(pos, end));
                    // coverage G
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer coverageG = Integer.parseInt(line.substring(pos, end));
                    // coverage C
                    pos = end + 1;
                    end = line.indexOf(FIELDS_SEPERATOR, pos);
                    Integer coverageC = Integer.parseInt(line.substring(pos, end));
                    // coverage N
                    pos = end + 1;
                    Integer coverageN = Integer.parseInt(line.substring(pos));

                    if (!sequenceName.equals(previousSeqName)) {
                        if (!previousSeqName.isEmpty()) {
                            baseCoverageList.put(previousSeqName, seqBaseCoverageList);
                        }
                        seqBaseCoverageList = (BaseCoverage[]) baseCoverageList.get(sequenceName);

                        // set the sequence name as previous sequence name
                        previousSeqName = sequenceName;
                    }
                    seqBaseCoverageList[position - 1] = new BaseCoverage(coverageA, coverageT, coverageG, coverageC, coverageN);
                }
            } // end while
            baseCoverageList.put(previousSeqName, seqBaseCoverageList);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(BaseCoverage.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(BaseCoverage.class.getName()).log(Level.SEVERE, null, ex);
        }

        return baseCoverageList;

    }

    public static LinkedHashMap<String, BaseCoverage[]> initialise(SAMSequenceDictionary sequenceDictionary) {

        LinkedHashMap<String, BaseCoverage[]> baseCoverageList = new LinkedHashMap();
        if (sequenceDictionary == null || sequenceDictionary.isEmpty()) {
            return baseCoverageList;
        }

        Iterator<SAMSequenceRecord> sequenceIterator = sequenceDictionary.getSequences().iterator();
        while (sequenceIterator.hasNext()) {
            // get the next sequence
            SAMSequenceRecord sequenceEntry = sequenceIterator.next();
            // Get sequence details
            String sequenceName = sequenceEntry.getSequenceName();
            int sequenceLength = sequenceEntry.getSequenceLength();

            // initialise the base coverage list for this sequence
            BaseCoverage[] seqBaseCoverageList = new BaseCoverage[sequenceLength];
            // add it to the list
            baseCoverageList.put(sequenceName, seqBaseCoverageList);
        }

        return baseCoverageList;
    }

    public void incrementCoverage(char base) {
        switch (base) {
            case 'A':
                coverageA++;
                break;
            case 'T':
                coverageT++;
                break;
            case 'G':
                coverageG++;
                break;
            case 'C':
                coverageC++;
                break;
            case 'N':
                coverageN++;
                break;
        }
    }

    public static void write(LinkedHashMap<String, BaseCoverage[]> baseCoverage, String outFile) {

        try {
            // open the output writer
            BufferedWriter out = new BufferedWriter(new FileWriter(outFile));

            // write the header
            Iterator<Map.Entry<String, BaseCoverage[]>> itBaseCoverage = baseCoverage.entrySet().iterator();
            while (itBaseCoverage.hasNext()) {
                // get the next entry
                Map.Entry<String, BaseCoverage[]> entry = itBaseCoverage.next();

                String sequenceName = entry.getKey();
                int sequenceLength = entry.getValue().length;

                out.write("@SQ\tSN:" + sequenceName + "\tLN:" + Integer.toString(sequenceLength));
                out.newLine();
            }

            // write the base distribution
            itBaseCoverage = baseCoverage.entrySet().iterator();
            while (itBaseCoverage.hasNext()) {
                // get the next entry
                Map.Entry<String, BaseCoverage[]> entry = itBaseCoverage.next();

                String sequence_name = entry.getKey();
                BaseCoverage[] seqBaseCoverage = entry.getValue();

                String emptyPosition = new BaseCoverage().toString();
                for (int i = 0; i < seqBaseCoverage.length; i++) {
                    if (seqBaseCoverage[i] != null) {
                        out.write(sequence_name + "\t" + Integer.toString(i + 1) + "\t" + seqBaseCoverage[i].toString());
                        out.newLine();
                    } else {
                        out.write(sequence_name + "\t" + Integer.toString(i + 1) + "\t" + emptyPosition);
                        out.newLine();
                    }
                }
            }

            // close the output buffer
            out.close();
        } catch (IOException ex) {
            // do nothing
        }
    }

    /**
     * @return the coverageA
     */
    public int getCoverageA() {
        return coverageA;
    }

    /**
     * @return the coverageC
     */
    public int getCoverageC() {
        return coverageC;
    }

    /**
     * @return the coverageG
     */
    public int getCoverageG() {
        return coverageG;
    }

    /**
     * @return the coverageT
     */
    public int getCoverageT() {
        return coverageT;
    }

    /**
     * @return the coverageN
     */
    public int getCoverageN() {
        return coverageN;
    }

    public static void calculateBaseCoverages(String samFile, String outBaseCoverageFile, int baseQualityThreshold) {
        /* Calculate base coverages */

        // Read the SAM/BAM File
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader samReader = samReaderFactory.open(new File(samFile));
        SAMFileHeader samFileHeader = samReader.getFileHeader();
        SAMSequenceDictionary sequenceDictionary = samFileHeader.getSequenceDictionary();

        //initialise base coverage list
        LinkedHashMap<String, BaseCoverage[]> baseDistributionList = BaseCoverage.initialise(sequenceDictionary);

        // read the same file
        Iterator<SAMRecord> samRecordIterator = samReader.iterator();
        while (samRecordIterator.hasNext()) {
            // get the next alignment
            SAMRecord theSAMAlignment = samRecordIterator.next();

            String readSequence = theSAMAlignment.getReadString();
            String baseQualities = theSAMAlignment.getBaseQualityString();

            // find the occurrence of "D","I" or "S" in the cigar
            String cigar = theSAMAlignment.getCigarString();
            if (cigar.contains("D") || cigar.contains("I") || cigar.contains("S")) {
                // process this read

                String newReadSequence = "";
                String newBaseQualities = "";
                String length = "";
                int lastNonDigitPosition = -1;
                int sum = 0;
                int i = 0;
                int start = 0;
                // int total_deleted = 0;
                while (i < cigar.length()) {
                    for (i = start; i < cigar.length(); i++) {
                        if (Character.isDigit(cigar.charAt(i))) {
                            length += cigar.charAt(i);
                        } else if (cigar.charAt(i) == 'S') {
                            if (i == cigar.length() - 1) {
                                for (int j = readSequence.length() - 1; j >= readSequence.length() - Integer.parseInt(length); j--) {
                                    readSequence = replace(readSequence, j, '-');
                                    baseQualities = replace(baseQualities, j, '-');
                                }
                            } else {
                                readSequence = readSequence.substring(Integer.parseInt(length)); //,read.length());
                                baseQualities = baseQualities.substring(Integer.parseInt(length));
                            }
                            length = "";
                            lastNonDigitPosition = i;
                        } else if (cigar.charAt(i) == 'D') {
                            int deletionSize = Integer.parseInt(cigar.substring(lastNonDigitPosition + 1, i));
                            lastNonDigitPosition = i;
                            int deletionStart = sum;
                            start = i + 1;
                            length = "";

                            char[] gaps = new char[deletionSize];
                            Arrays.fill(gaps, '.');
                            String gap = new String(gaps);

                            newReadSequence += readSequence.substring(0, deletionStart) + gap;
                            readSequence = readSequence.substring(deletionStart);

                            newBaseQualities += baseQualities.substring(0, deletionStart) + gap;
                            baseQualities = baseQualities.substring(deletionStart, baseQualities.length());

                            //total_deleted += deletion_start;
                            sum = 0;
                            break;
                        } else if (cigar.charAt(i) == 'I') {
                            int insertionSize = Integer.parseInt(cigar.substring(lastNonDigitPosition + 1, i));
                            lastNonDigitPosition = i;
                            int insertionStart = sum;
                            start = i + 1;
                            sum = 0;
                            length = "";

                            newReadSequence += readSequence.substring(0, insertionStart);
                            readSequence = readSequence.substring(insertionStart + insertionSize);

                            newBaseQualities += baseQualities.substring(0, insertionStart);
                            baseQualities = baseQualities.substring(insertionStart + insertionSize);

                            break;
                        } else if (cigar.charAt(i) == 'M') {
                            sum += Integer.parseInt(length);
                            length = "";
                            lastNonDigitPosition = i;
                        } else if (cigar.charAt(i) == 'N') { // introns -- ignore
                            length = "";
                            lastNonDigitPosition = i;
                        }
                    }
                } // end while
                newReadSequence += readSequence;
                newBaseQualities += baseQualities;

                readSequence = newReadSequence;
                baseQualities = newBaseQualities;
            }

            // get the base distribution for the sequence to which the read is mapped
            BaseCoverage[] seqBaseDistributionList = (BaseCoverage[]) baseDistributionList.get(theSAMAlignment.getReferenceName());
            for (int i = 0; i < readSequence.length(); i++) {
                if (readSequence.charAt(i) == '-') { // soft-clip .. ignore this position in this read
                    continue;
                } else if (readSequence.charAt(i) == '.') { // gap due to deletion .. ignore
                    continue;
                }

                if (theSAMAlignment.getAlignmentStart() + i > seqBaseDistributionList.length) {
                    break;
                }

                // Ignore this base if the sequencing quality is low
                int baseQuality = Integer.valueOf(baseQualities.charAt(i)) - 33; // base qualities are ascii - 33
                if (baseQuality < baseQualityThreshold) {
                    continue;
                }

                if (seqBaseDistributionList[theSAMAlignment.getAlignmentStart() + i - 1] == null) {
                    seqBaseDistributionList[theSAMAlignment.getAlignmentStart() + i - 1] = new BaseCoverage();
                }
                seqBaseDistributionList[theSAMAlignment.getAlignmentStart() + i - 1].incrementCoverage(readSequence.charAt(i));
                //the_base_distribution->set_base_distribution_list(sequence_name, seq_base_distribution_ptr);
            } // end for (each character in the read) 

        } // end while (ifs_sam.good())

        try {
            samReader.close();
        } catch (IOException ex) {
            Logger.getLogger(BaseCoverage.class.getName()).log(Level.SEVERE, null, ex);
        }

        // write to file, if a file is provided
        if (!outBaseCoverageFile.isEmpty()) {
            BaseCoverage.write(baseDistributionList, outBaseCoverageFile);
        }
    }

    private static String replace(String str, int index, char replace) {
        if (str == null) {
            return str;
        } else if (index < 0 || index >= str.length()) {
            return str;
        }
        char[] chars = str.toCharArray();
        chars[index] = replace;
        return String.valueOf(chars);
    }

}
