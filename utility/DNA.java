/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import static utility.GFFEntry.FASTA_TAG;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class DNA {

    public static final char BASE_BLANK = ' ';
    public static final char BASE_ZERO_COVERAGE = '0';
    public static final char BASE_LOW_COVERAGE = '<';
    public static final char BASE_AMBIGUOUS = '*';
    public static final char BASE_N = 'N';
    public static final char BASE_MISSING = '?';
    private static final int BLOCK_SIZE = 80;
    private static HashMap nucleotideMap;
    private static HashMap reverseNucleotideMap;

    private static void initialiseNucleotideMap() {

        TreeSet A = new TreeSet();
        A.add('A');
        TreeSet C = new TreeSet();
        C.add('C');
        TreeSet G = new TreeSet();
        G.add('G');
        TreeSet T = new TreeSet();
        T.add('T');
        TreeSet R = new TreeSet();
        R.add('A');
        R.add('G');
        TreeSet Y = new TreeSet();
        Y.add('C');
        Y.add('T');
        TreeSet M = new TreeSet();
        M.add('A');
        M.add('C');
        TreeSet K = new TreeSet();
        K.add('G');
        K.add('T');
        TreeSet S = new TreeSet();
        S.add('C');
        S.add('G');
        TreeSet W = new TreeSet();
        W.add('A');
        W.add('T');
        TreeSet H = new TreeSet();
        H.add('A');
        H.add('C');
        H.add('T');
        TreeSet B = new TreeSet();
        B.add('C');
        B.add('G');
        B.add('T');
        TreeSet D = new TreeSet();
        D.add('A');
        D.add('G');
        D.add('T');
        TreeSet V = new TreeSet();
        V.add('A');
        V.add('C');
        V.add('G');
        TreeSet N = new TreeSet();
        N.add('N');

        if (nucleotideMap == null) {
            nucleotideMap = new HashMap();
        }

        nucleotideMap.put('A', A);
        nucleotideMap.put('C', C);
        nucleotideMap.put('G', G);
        nucleotideMap.put('T', T);
        nucleotideMap.put('R', R);
        nucleotideMap.put('Y', Y);
        nucleotideMap.put('M', M);
        nucleotideMap.put('K', K);
        nucleotideMap.put('S', S);
        nucleotideMap.put('W', W);
        nucleotideMap.put('H', H);
        nucleotideMap.put('B', B);
        nucleotideMap.put('D', D);
        nucleotideMap.put('V', V);
        nucleotideMap.put('N', N);

    }

    private static void initialiseReverseNucleotideMap() {

        initialiseNucleotideMap();

        if (reverseNucleotideMap == null) {
            reverseNucleotideMap = new HashMap();
        }

        Iterator it = nucleotideMap.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry entry = (Map.Entry) it.next();
            // add the reverse pair
            reverseNucleotideMap.put(entry.getValue().toString(), entry.getKey());
        }

        // also add all four bases as 'N'
        TreeSet N = new TreeSet();
        N.add('A');
        N.add('C');
        N.add('G');
        N.add('T');
        reverseNucleotideMap.put(N.toString(), 'N');

    }

    public static HashMap getNucleotideMap() {

        if (nucleotideMap == null || nucleotideMap.isEmpty()) {
            initialiseNucleotideMap();
        }

        return nucleotideMap;
    }

    public static HashMap getReverseNucleotideMap() {

        if (reverseNucleotideMap == null || reverseNucleotideMap.isEmpty()) {
            initialiseReverseNucleotideMap();
        }

        return reverseNucleotideMap;
    }

    public static boolean isNucleotide(char base) {
        if (base == 'A' || base == 'T' || base == 'G' || base == 'C'
                || base == 'a' || base == 't' || base == 'g' || base == 'c') {
            return true;
        } else {
            return false;
        }
    }

    public static boolean isExtendedNucleotide(char base) {

        char upperBase = Character.toUpperCase(base);
        if (upperBase == 'R' || upperBase == 'Y' || upperBase == 'M' || upperBase == 'K' || upperBase == 'S' || upperBase == 'W' || upperBase == 'H' || upperBase == 'B' || upperBase == 'D' || upperBase == 'V' || upperBase == 'N') {
            return true;
        } else {
            return false;
        }
    }

    public static int multiFastaToFasta(File multiFastaFile, File outFastaFile, String fastaHeader, int gapSize) {
        BufferedReader br = null;
        BufferedWriter out = null;
        int sequenceCount = 0;
        try {
            // create the separator
            char[] gaps = new char[gapSize];
            Arrays.fill(gaps, BASE_N);
            String separator = new String(gaps);

            // get the fasta fame
            String fastaName = extractSequenceName(fastaHeader);

            GFF theGFF = new GFF();
            ArrayList<GFFEntry> seqGFFList = new ArrayList<GFFEntry>();
            theGFF.setGFFList(fastaName, seqGFFList);

            String line;
            String sequenceName = "";
            String sequenceHeader = "";
            StringBuilder sequence = new StringBuilder(separator);
            int sequenceStart = separator.length() + 1;
            int sequenceLength = 0;
            // read the input file
            br = new BufferedReader(new FileReader(multiFastaFile));
            // open the output writer
            out = new BufferedWriter(new FileWriter(outFastaFile));

            while ((line = br.readLine()) != null) {
                // ignore empty lines
                if (line.isEmpty()) {
                    continue;
                }

                // match the reference
                if (line.charAt(0) == '>') { // we are reading a new sequence 
                    if (sequenceLength > 0) {

                        int sequenceEnd = sequenceStart + sequenceLength - 1;
                        // create the gff entry object and add it to the list of gff entries
                        seqGFFList.add(new GFFEntry(fastaName, sequenceStart, sequenceEnd, GFFEntry.addIDAttribute(sequenceName, sequenceHeader)));

        
                        sequence.append(separator);
                        sequenceStart += sequenceLength + separator.length();
                        sequenceLength = 0;
                    }

                    // get the details for this sequence
                    sequenceName = extractSequenceName(line);
                    sequenceHeader = line.substring(1);

                    // increase the sequence count
                    sequenceCount++;
                } else {
                    String str = line.trim().toUpperCase();
                    sequence.append(str);
                    sequenceLength += str.length();
                }
            } // end while ifs_fasta.good()
            // process the last sequence
            int sequenceEnd = sequenceStart + sequenceLength - 1;
            // create the gff entry object and add it to the list of gff entries
            seqGFFList.add(new GFFEntry(fastaName, sequenceStart, sequenceEnd, GFFEntry.addIDAttribute(sequenceName, sequenceHeader)));
            sequence.append(separator);

            // close the input file
            br.close();

            //AMT System.out.println("\tWriting the fasta file");
            // write the fasta file
            // open the output writer
            out = new BufferedWriter(new FileWriter(outFastaFile));
            // write to the output
            write(fastaHeader, sequence.toString(), out);
            // close the output file
            out.close();

            //AMT System.out.println("\tWriting the GFF file");
            // write the GFF file
            theGFF.write(outFastaFile + ".gff");

        } catch (FileNotFoundException ex) {
            Logger.getLogger(DNA.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(DNA.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(DNA.class.getName()).log(Level.SEVERE, null, ex);
            }
            try {
                out.close();
            } catch (IOException ex) {
                Logger.getLogger(DNA.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return sequenceCount;


    }
    
    private static String extractSequenceName(String string) {
        int start = (string.charAt(0) == '>' ? 1 : 0);
        int position;
        if ((position = string.indexOf(' ')) >= 0) {
            return string.substring(start, position).trim();
        } else {
            return string.substring(start).trim();
        }

    }

    private static void write(String sequenceHeader, String sequence, BufferedWriter out) throws IOException {

        // write the header if supplied
        if (!sequenceHeader.isEmpty()) {
            if (sequenceHeader.charAt(0) != '>') {
                out.write(">");
            }
            out.write(sequenceHeader);
            out.newLine();
        }
        // write the sequence
        int position = 0;
        while (position < sequence.length()) {
            out.write(sequence.substring(position, Math.min(position + BLOCK_SIZE, sequence.length())));
            out.newLine();
            position += BLOCK_SIZE;
        }
    }


}
