/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hands2;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import utility.Base;
import utility.BaseCoverage;
import utility.DNA;
import utility.GFF;
import utility.GFFEntry;
import utility.Variant;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class HANDS2 {

    // constant
    static final String HANDS2_VERSION = "v1.1.1";
    static final int MAX_DIPLOIDS = 10;
    static final char BASE_PATTERN_PROPORTION_MODE_MAXIMUM = 'M';
    static final char BASE_PATTERN_PROPORTION_MODE_ADDITIVE = 'A';
    static final int BASE_QUALITY_THRESHOLD = 20;
    // Global variables
    static double SNP_PAIR_PROPORTION_THRESHOLD = 0.05;
    static double SNP_PATTERN_MATCHING_THRESHOLD = 0.5;
    static double SINGLE_POSITION_SNP_PATTERN_COUNT_THRESHOLD = 0;
    static boolean RECTIFY_USING_REFERENCE = false;
    static boolean MERGE_BASE_PATTERNS = false;
    static boolean OUTPUT_VCF = true;
    static boolean TRY_BASE_SET_CLOSURE = true;
    static int DISTANT_GENOME = 0;
    static int DIPLOID_COUNT = 0;
    static int MISSING_DIPLOID = 0;
    static char BASE_PATTERN_PROPORTION_MODE = BASE_PATTERN_PROPORTION_MODE_MAXIMUM;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        //       args = "assign -s ../example/wheat/polyploid.sam -g ../example/wheat/reference.fa.gff -hsp ../example/wheat/polyploid.hsp -snp1 ../example/wheat/diploid1.snp -snp2 ../example/wheat/diploid2.snp -snp3 ../example/wheat/diploid3.snp -bc ../example/wheat/polyploid.bc -bc1 ../example/wheat/diploid1.bc -bc2 ../example/wheat/diploid2.bc -bc3 ../example/wheat/diploid3.bc -out1 ../example/out1.txt -out2 ../example/out2.txt -out3 ../example/out3.txt".split(" +");
//        args = "-s /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR833617-20.sam -g /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-20.fa.gff -hsp /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR833617-20.snp.accepted.con -snp2 /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR1213099-20.snp.accepted -snp1 \"\" -bd /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR833617-20.base.dist -bd2 /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR1213099-20.base.dist -bd1 \"\" -out1 /Volumes/Users/data/brassica/projects/hands/carinata-test/out1.txt -out2 /Volumes/Users/data/brassica/projects/hands/carinata-test/out2.txt -m true".split(" +");
//        args = "-s /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR833617-20.sam -g /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-20.fa.gff -hsp /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR833617-20.snp.accepted.con -snp2 /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR1213099-20.snp.accepted -snp1 \"\" -bd /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR833617-20.base.dist -bd2 /Volumes/Users/data/brassica/projects/hands/carinata-test/Bra1.27-SRR1213099-20.base.dist -bd1 \"\" -out1 /Volumes/Users/data/brassica/projects/hands/carinata-test/out1.txt -out2 /Volumes/Users/data/brassica/projects/hands/carinata-test/out2.txt -m true -p /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.snp.accepted.con".split(" +");
        //args = "-s /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.sam -g /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-20.fa.gff -hsp /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.snp.accepted.con -snp1 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR6436235-20.snp.accepted -snp2 \"\" -bd /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.base.dist -bd1 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR6436235-20.base.dist -bd2 \"\" -out1 /Volumes/Users/data/brassica/projects/hands/napus-test/out1.txt -out2 /Volumes/Users/data/brassica/projects/hands/napus-test/out2.txt".split(" +");
//        args = "assign -s /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.sam -g /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-20.fa.gff -hsp /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.snp.new -snp1 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR6436235-20.snp.new -snp2 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR1213099-20.snp.new -bc /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.base.dist -bc1 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR6436235-20.base.dist -bc2 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR1213099-20.base.dist -out1 /Volumes/Users/data/brassica/projects/hands/napus-test/out1.txt -out2 /Volumes/Users/data/brassica/projects/hands/napus-test/out2.txt -m true -vcf true".split(" +");
//        args = "assign -s /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.sam -g /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-20.fa.gff -hsp /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.snp.new -snp1 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR6436235-20.snp.new -snp2 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR1213099-20.snp.new -bc /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR948826-20.base.dist -bc1 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR6436235-20.base.dist -bc2 /Volumes/Users/data/brassica/projects/hands/napus-test/Bra1.27-SRR1213099-20.base.dist -out1 /Volumes/Users/data/brassica/projects/hands/napus-test/out1.txt -out2 /Volumes/Users/data/brassica/projects/hands/napus-test/out2.txt -m true".split(" +");
        if (args.length == 0 || args[0].equals("help")) {
            printHelp();
        } else if (args[0].equals("seq2ref")) {
            // multifasta to fasta
            multiFastaToFasta(args);
        } else if (args[0].equals("coverage")) {
            // base coverage
            calculateBaseCoverages(args);
        } else if (args[0].equals("assign")) {
            // assign homoeallelic base-identities
            assignBases(args);
        } else {
            System.err.println("Invalid command: " + args[0]);
            printHelp();
        }

    }

    private static void assignBases(String[] args) {
        String alignmentFile = "";
        String GFFFile = "";
        String polyploidHSPFile = "";
        ArrayList<String> diploidSNPFiles = new ArrayList<String>(MAX_DIPLOIDS);
        String polyploidBaseCoverageFile = "";
        ArrayList<String> diploidBaseCoverageFiles = new ArrayList<String>(MAX_DIPLOIDS);
        ArrayList<String> outFiles = new ArrayList<String>(MAX_DIPLOIDS);
        String positionFile = ""; // to specify extra positions to be checked in pre-processing
        char merge = ' '; // input for merge
        char vcf = ' ';

        for (int i = 0; i < MAX_DIPLOIDS; i++) {
            diploidSNPFiles.add(null);
            diploidBaseCoverageFiles.add(null);
            outFiles.add(null);
        }

        // read input arguments
        int idx = 1;
        int nArgs = args.length;
        while (idx < nArgs) {
            String parameter = args[idx++];
            if (parameter.equals("-h") || parameter.equals("-help")) { // Help
                if (nArgs == 2) {
                    printHelpAssign();
                } else {
                    System.err.println("Invalid option used with Help!");
                }
                return;
            } else if (parameter.equals("-i")) { // SAM/BAM file
                alignmentFile = args[idx++];
            } else if (parameter.equals("-g")) { // GFF file
                GFFFile = args[idx++];
            } else if (parameter.equals("-hsp")) { // polyploid HSP File
                polyploidHSPFile = args[idx++];
            } else if (parameter.startsWith("-snp")) { // Diploid Variant Files
                int diploidNo = Integer.parseInt(parameter.substring(4)) - 1;
                if (diploidNo < 0) {
                    System.err.println("Invalid Diploid Number.");
                    return;
                }
                String filename = args[idx++];
                filename = filename.equals("\"\"") ? "" : filename;
                diploidSNPFiles.set(diploidNo, filename);
            } else if (parameter.equals("-bc")) { // polyploid Base Coverage File
                polyploidBaseCoverageFile = args[idx++];
            } else if (parameter.startsWith("-bc")) { // Diploid Base Coverage Files
                int diploidNo = Integer.parseInt(parameter.substring(3)) - 1;
                if (diploidNo < 0) {
                    System.err.println("Invalid Diploid Number.");
                    return;
                }
                String filename = args[idx++];
                filename = filename.equals("\"\"") ? "" : filename;
                diploidBaseCoverageFiles.set(diploidNo, filename);
            } else if (parameter.startsWith("-out")) { // Diploid Output Files
                int diploidNo = Integer.parseInt(parameter.substring(4)) - 1;
                if (diploidNo < 0) {
                    System.err.println("Invalid Diploid Number.");
                    return;
                }
                String filename = args[idx++];
                filename = filename.equals("\"\"") ? "" : filename;
                outFiles.set(diploidNo, filename);
            } else if (parameter.equals("-p")) { // Position file (optional) to provide positions to check in the pre-processing step
                positionFile = args[idx++];
            } else if (parameter.equals("-sp")) { // Variant Pair Proportion threshold
                SNP_PAIR_PROPORTION_THRESHOLD = Double.valueOf(args[idx++]);
            } else if (parameter.equals("-pm")) { // Variant Pattern Matching threshold
                SNP_PATTERN_MATCHING_THRESHOLD = Double.valueOf(args[idx++]);
            } else if (parameter.equals("-r")) { // rectify using reference genome
                char rectify = args[idx++].charAt(0);
                RECTIFY_USING_REFERENCE = (Character.toUpperCase(rectify) == 'T');
            } else if (parameter.equals("-d")) { //  distant genome number
                DISTANT_GENOME = Integer.valueOf(args[idx++]);
            } else if (parameter.equals("-m")) { //  merge base patterns
                merge = args[idx++].charAt(0);
                MERGE_BASE_PATTERNS = (Character.toUpperCase(merge) == 'T');
            } else if (parameter.equals("-vcf")) { //  Output in VCF format
                vcf = args[idx++].charAt(0);
                OUTPUT_VCF = (Character.toUpperCase(vcf) == 'T');
            } else if (parameter.equals("-pa")) { // Variant Pattern Assignment Mode (Keep Maximum Proportion or Add All Proportions)
                BASE_PATTERN_PROPORTION_MODE = Character.toUpperCase(args[idx++].charAt(0));
                if (BASE_PATTERN_PROPORTION_MODE != BASE_PATTERN_PROPORTION_MODE_MAXIMUM && BASE_PATTERN_PROPORTION_MODE != BASE_PATTERN_PROPORTION_MODE_ADDITIVE) {
                    System.err.println("Invalid SNP Pattern Assignment Mode.");
                    return;
                }
            } else if (parameter.equals("-u")) { // perform base set closure -- assign the left over base to unassigned subgenome
                char closure = args[idx++].charAt(0);
                TRY_BASE_SET_CLOSURE = (Character.toUpperCase(closure) == 'T');
            } else {
                System.err.println("Invalid option: " + args[idx - 1]);
                printHelpAssign();
                return;

            }
        }

        // validate the input 
        if (alignmentFile.isEmpty()) {
            System.err.println("SAM/BAM file must be specified.");
            return;
        } else if (GFFFile.isEmpty()) {
            System.err.println("GFF file must be specified.");
            return;
        } else if (polyploidHSPFile.isEmpty()) {
            System.err.println("Polyploid HSP File must be specified.");
            return;
        }

        // validate diploid Variant files
        for (int d = MAX_DIPLOIDS - 1; d >= 0; d--) {
            if (diploidSNPFiles.get(d) != null) {
                DIPLOID_COUNT = d + 1;
                break;
            }
        }
        if (DIPLOID_COUNT == 0) {
            System.err.println("At least one diploid SNP file must be specified.");
            return;
        }

        int missingDiploidCount = 0;
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            if (diploidSNPFiles.get(d) == null) {
                System.err.println("Diploid " + Integer.toString(d + 1) + " SNP file missing. Use \"\" as the file name if you would like to skip the file.");
                return;
            } else if (diploidSNPFiles.get(d).isEmpty()) {
                missingDiploidCount++;
                MISSING_DIPLOID = d + 1;
            }
        }
        if (missingDiploidCount > 1) {
            System.err.println("Only one diploid SNP file can be skipped.");
            return;
        }
        // shrink the array list
        diploidSNPFiles.subList(DIPLOID_COUNT, diploidSNPFiles.size()).clear();

        if (DISTANT_GENOME > DIPLOID_COUNT) {
            System.err.println("Invalid Distant Genome.");
            return;
        }
        if (DISTANT_GENOME > 0 && MISSING_DIPLOID > 0) {
            System.err.println("Distant Genome cannot be used when a genome is missing.");
            return;
        } else if (MISSING_DIPLOID > 0) {
            DISTANT_GENOME = MISSING_DIPLOID;
        }

        // validate the base distribution and output files       
        for (int d = MAX_DIPLOIDS - 1; d >= 0; d--) {
            if (diploidBaseCoverageFiles.get(d) != null) {
                if (d > DIPLOID_COUNT - 1) {
                    System.err.println("Extra diploid base coverage file(s) supplied.");
                    return;
                } else if (!diploidBaseCoverageFiles.get(d).isEmpty() && d == MISSING_DIPLOID - 1) {
                    System.err.println("Diploid base coverage file for skipped SNP file supplied.");
                    return;
                }
            }
            if (outFiles.get(d) != null && d > DIPLOID_COUNT - 1) {
                System.err.println("Extra output file(s) supplied.");
                return;
            }
        }
        // shrink the array lists
        diploidBaseCoverageFiles.subList(DIPLOID_COUNT, diploidBaseCoverageFiles.size()).clear();
        outFiles.subList(DIPLOID_COUNT, outFiles.size()).clear();

        // replace null with empty string for base coverage and a suitable name for output file
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            if (diploidBaseCoverageFiles.get(d) == null) {
                diploidBaseCoverageFiles.set(d, "");
            }
            if (outFiles.get(d) == null || outFiles.get(d).isEmpty()) {
                if (d + 1 == MISSING_DIPLOID) {
                    outFiles.set(d, "diploid" + Integer.toString(d + 1) + ".out");
                } else {
                    outFiles.set(d, diploidSNPFiles.get(d) + ".out");
                }
            }
        }
        int diploidBaseDistributionFilesCount = 0;
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            if (!diploidBaseCoverageFiles.get(d).isEmpty()) {
                diploidBaseDistributionFilesCount++;
            }
        }

        if (!positionFile.isEmpty() && polyploidBaseCoverageFile.isEmpty() && diploidBaseDistributionFilesCount == 0) {
            // skip pre-processing 
            System.out.println("Position file spcified but base coverage files not specified. Ignoring position file ...");
        }

        System.out.println("Number of diploids: " + Integer.toString(DIPLOID_COUNT));
        if (MISSING_DIPLOID > 0) {
            System.out.println("Treating diploid " + Integer.toString(MISSING_DIPLOID) + " as missing diploid.");
            if (merge == ' ') { // if no input is specified then set 'Merge Base Patterns = True' by default
                MERGE_BASE_PATTERNS = true;
            }
        }

        // set the flag to skip preprocessing if no base coverage files are provided
        boolean skipPreprocessing = (polyploidBaseCoverageFile.isEmpty() && diploidBaseDistributionFilesCount == 0);

        // Read the SAM/BAM File
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader samReader = samReaderFactory.open(new File(alignmentFile));
        SAMFileHeader samFileHeader = samReader.getFileHeader();
        SAMSequenceDictionary sequenceDictionary = samFileHeader.getSequenceDictionary();

        // read the GFF file
        GFF theGFF = new GFF(GFFFile);

        // read Variant lists
        System.out.print("Reading HSP list ... ");
        HashMap<String, LinkedHashMap<Integer, Variant>> referenceBaseList = (skipPreprocessing) ? null : Variant.initialiseVariantList(sequenceDictionary, true);
        HashMap<String, HashMap<Integer, Variant>> polyploidHSPList = Variant.initialiseVariantList(sequenceDictionary, false);
        //int polyploidHSPListSize = Variant.readVariantList(polyploidHSPFile, theSAM.getSequenceLengths(), polyploidHSPList);
        int polyploidHSPListSize = Variant.readVariantList(polyploidHSPFile, polyploidHSPList, referenceBaseList);
        if (polyploidHSPListSize < 0) {
            System.out.println();
            return;
        }
        System.out.println(Integer.toString(polyploidHSPListSize) + " positions.");
        if (DIPLOID_COUNT - missingDiploidCount > 1) {
            System.out.println("Reading Diploid SNP lists ...");
        } else {
            System.out.println("Reading Diploid SNP list ...");
        }
        ArrayList<HashMap<String, HashMap<Integer, Variant>>> diploidSNPLists = new ArrayList<HashMap<String, HashMap<Integer, Variant>>>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            HashMap<String, HashMap<Integer, Variant>> diploidSNPList = Variant.initialiseVariantList(sequenceDictionary, false);
            //int diploidSNPListSize = Variant.readVariantList(diploidSNPFiles.get(d), theSAM.getSequenceLengths(), diploidSNPList);            
            int diploidSNPListSize = Variant.readVariantList(diploidSNPFiles.get(d), diploidSNPList, referenceBaseList);
            if (diploidSNPListSize < 0) {
                return;
            }
            if (d != (MISSING_DIPLOID - 1)) {
                System.out.println("\tDiploid " + Integer.toString(d + 1) + ": " + Integer.toString(diploidSNPListSize) + " positions.");
            }
            diploidSNPLists.add(diploidSNPList);
        }

        if (polyploidBaseCoverageFile.isEmpty() && diploidBaseDistributionFilesCount == 0) {
            // skip pre-processing 
            System.out.println("Base distibution files not specified. Skipping data pre-processing ...");
        } else {
            System.out.println("Data pre-processing ...");
            // perform the data pre-processing
            dataPrepocessing(positionFile, polyploidBaseCoverageFile, diploidBaseCoverageFiles, referenceBaseList, polyploidHSPList, diploidSNPLists, sequenceDictionary);
        }

        //now reference base list is redundant. free the memory
        referenceBaseList = null;

        ArrayList<BufferedWriter> outList = new ArrayList<BufferedWriter>(DIPLOID_COUNT);
        try {
            for (int d = 0; d < DIPLOID_COUNT; d++) {
                //Construct the BufferedWriter object
                outList.add(new BufferedWriter(new FileWriter(outFiles.get(d))));
            }
        } catch (IOException ex) {
            System.err.println("Unable to open output file(s)");
            return;
        }

        if (alignmentFile.endsWith(".bam")) {
            System.out.println("Processing BAM file ...");
        } else if (alignmentFile.endsWith(".sam")) {
            System.out.println("Processing SAM file ...");
        } else {
            System.out.println("Processing alignment file ...");
        }

        // process SAM/BAM file
        processAlignmentFile(samReader, theGFF, polyploidHSPList, diploidSNPLists, outList, OUTPUT_VCF, args);

        try {
            for (int d = 0; d < DIPLOID_COUNT; d++) {
                //close the BufferedWriter object
                outList.get(d).close();
            }
        } catch (IOException ex) {
            System.err.println("Unable to close output file(s)");
        }

        try {
            samReader.close();
        } catch (IOException ex) {
            Logger.getLogger(HANDS2.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private static void dataPrepocessing(String positionFile, String polyploidBaseCoverageFile, ArrayList<String> diploidBaseCoverageFiles,
            HashMap<String, LinkedHashMap<Integer, Variant>> positionList, HashMap<String, HashMap<Integer, Variant>> polyploidHSPList, ArrayList<HashMap<String, HashMap<Integer, Variant>>> diploidSNPLists,
            SAMSequenceDictionary sequenceDictionary) {

        if (!positionFile.isEmpty()) { // position file provided with extra positions to be checked in pre-processing
            // add them to reference base list
            System.out.println("Reading Position file ...");
            int positionListSize = Variant.readPositionList(positionFile, positionList);
            System.out.println("\t Number of positions: " + Integer.toString(positionListSize));
        }

        // Polyploid
        if (!polyploidBaseCoverageFile.isEmpty()) {
            if (new File(polyploidBaseCoverageFile).exists()) {
                System.out.println("Checking HSPs ...");
                // read base distribution files
                HashMap<String, BaseCoverage[]> polyploidBaseDistributionList = BaseCoverage.read(polyploidBaseCoverageFile, sequenceDictionary);
                // check the polyploid Variant list to see if consensus bases have been called correctly
                int polyploidHSPListSize = Variant.checkPolyploidHSPList(polyploidHSPList, polyploidBaseDistributionList, positionList);
                System.out.println("\t Number of positions: " + Integer.toString(polyploidHSPListSize));
            } else {
                System.err.println("Error: Cannot locate the file " + polyploidBaseCoverageFile);
                return;
            }
        }
        // Diploids
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            String diploidBaseDistributionFile = diploidBaseCoverageFiles.get(d);
            if (new File(diploidBaseDistributionFile).exists()) {
                if (!diploidBaseDistributionFile.isEmpty()) {
                    System.out.println("Checking Diploid " + Integer.toString(d + 1) + " SNPs ...");
                    // read base distribution files
                    HashMap<String, BaseCoverage[]> diploidBaseDistributionList = BaseCoverage.read(diploidBaseDistributionFile, sequenceDictionary);
                    // check the diploid Variant list to see if any Variant is missing 
                    Variant.checkDiploidVariantList(diploidSNPLists.get(d), diploidBaseDistributionList, polyploidHSPList, positionList);
                } else {
                    System.err.println("Error: Cannot locate the file " + diploidBaseDistributionFile);
                    return;
                }
            }
        }
    }

    private static boolean processAlignmentFile(SamReader samReader, GFF theGFF, HashMap<String, HashMap<Integer, Variant>> polyploidHSPList,
            ArrayList<HashMap<String, HashMap<Integer, Variant>>> diploidSNPLists, ArrayList<BufferedWriter> outList, boolean outputVCF, String[] args) {

        ArrayList<HashMap<Integer, Variant>> seqDiploidSNPLists = new ArrayList<HashMap<Integer, Variant>>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            seqDiploidSNPLists.add(null);
        }
        ArrayList<LinkedHashMap> geneDiploidSNPLists = new ArrayList<LinkedHashMap>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            geneDiploidSNPLists.add(null);
        }

        // initialise output
        if (outputVCF) {
            int d = 0;
            for (BufferedWriter out : outList) {
                try {
                    writeVCFHeader(out, args, ++d, samReader.getFileHeader().getSequenceDictionary());
                } catch (IOException ex) {
                    System.err.println("Error writing output.");
                    return false;
                }
            }
        }

        GFFEntry theGFFEntry = null;
        ArrayList theSequenceGFFList;
        Iterator<GFFEntry> itGFFEntry = null;
        HashMap<Integer, Variant> seqPolyploidHSPList = null;
//AMT-CHECK        boolean withinValidGene = false;
        LinkedHashMap<Integer, Variant> geneHSPList = null;
        HashMap<String, TreeSet> pendingBasePatterns = new HashMap();
        LinkedHashMap<TreeSet<Base>, Integer> basePatterns = new LinkedHashMap();
        Iterator<SAMRecord> samRecordIterator = samReader.iterator();
        while (samRecordIterator.hasNext()) {
            // get the next alignment
            SAMRecord samRecord = samRecordIterator.next();

            //System.out.println(samRecord.getReadName() + " " + samRecord.getAlignmentStart() + " " + samRecord.getAlignmentEnd());
            // get the GFF list for the reference against which this read is aligned
            if (theGFFEntry == null // new gene data to  be fetched
                    || !samRecord.getReferenceName().equals(theGFFEntry.getSequenceName()) // new reference sequence. new gene data to  be fetched
                    || samRecord.getAlignmentStart() > theGFFEntry.getEnd()) { // all the reads for the current gene are finsihed
//AMT-CHECK                    || samRecord.getAlignmentEnd() < theGFFEntry.getStart()) {

                // all the reads for the current gene are finsihed. process the data for current gene
                if (theGFFEntry != null) {
//AMT-CHECK     if (theGFFEntry != null && withinValidGene) {
                    // process the previous set of basePatterns
                    processBasePatterns(basePatterns, geneHSPList, geneDiploidSNPLists, theGFFEntry.getSequenceName(), outList, outputVCF);

                    // reset the variables
                    basePatterns.clear();
                    pendingBasePatterns.clear();
//AMT-CHECK                    withinValidGene = false;
                }

                // new gene data to  be fetched (from the current or the new reference seqeunce)
                if (theGFFEntry == null || !samRecord.getReferenceName().equals(theGFFEntry.getSequenceName())) {

                    seqPolyploidHSPList = (HashMap<Integer, Variant>) polyploidHSPList.get(samRecord.getReferenceName());
                    // and the SNP lists for diploids
                    for (int d = 0; d < DIPLOID_COUNT; d++) {
                        HashMap<String, HashMap<Integer, Variant>> diploidSNPList = diploidSNPLists.get(d);
                        if (diploidSNPList != null) {
                            seqDiploidSNPLists.set(d, (HashMap<Integer, Variant>) diploidSNPList.get(samRecord.getReferenceName()));
                        } else {
                            seqDiploidSNPLists.set(d, null);
                        }
                    }

                    // get the GFF list for this sequence
                    theSequenceGFFList = (ArrayList) theGFF.getGFFList().get(samRecord.getReferenceName());

                    if (theSequenceGFFList == null) {
                        System.err.println("Sequence " + samRecord.getReferenceName() + " missing from the GFF file. Unable to continue.");
                        return false;
                    }

                    // and the iterator to the GFF entries
                    itGFFEntry = theSequenceGFFList.iterator();
                    if (itGFFEntry.hasNext()) {
                        theGFFEntry = itGFFEntry.next();

                        System.out.print(theGFFEntry.toString());

                        // Now theGFFEntry corresponds to the gene in which the read lies.
                        // Get the SNPs/HSPs within this gene
                        geneHSPList = getGeneSNPs(seqPolyploidHSPList, theGFFEntry.getStart() - GFF.GENE_FLANKING_REGION_LENGTH, theGFFEntry.getEnd() + GFF.GENE_FLANKING_REGION_LENGTH, true);

                        // also output the number of HSPs in this region
                        System.out.println("\t" + Integer.toString(geneHSPList.size()));

                        for (int d = 0; d < DIPLOID_COUNT; d++) {
                            HashMap<Integer, Variant> seqDiploidSNPList = seqDiploidSNPLists.get(d);
                            if (seqDiploidSNPList != null) {
                                geneDiploidSNPLists.set(d, getGeneSNPs(seqDiploidSNPList, theGFFEntry.getStart() - GFF.GENE_FLANKING_REGION_LENGTH, theGFFEntry.getEnd() + GFF.GENE_FLANKING_REGION_LENGTH, false));
                            } else {
                                geneDiploidSNPLists.set(d, null);
                            }
                        }

                    } else {
                        theGFFEntry = null;
                        continue;
                    }
                }

                // move to the gene in which this read lies
                while (theGFFEntry.getEnd() <= samRecord.getAlignmentStart()) {

                    if (itGFFEntry.hasNext()) {
                        theGFFEntry = itGFFEntry.next();

                        System.out.print(theGFFEntry.toString());

                        // Now theGFFEntry corresponds to the gene in which the read lies.
                        // Get the SNPs/HSPs within this gene
                        geneHSPList = getGeneSNPs(seqPolyploidHSPList, theGFFEntry.getStart() - GFF.GENE_FLANKING_REGION_LENGTH, theGFFEntry.getEnd() + GFF.GENE_FLANKING_REGION_LENGTH, true);

                        // also output the number of HSPs in this region
                        System.out.println("\t" + Integer.toString(geneHSPList.size()));

                        for (int d = 0; d < DIPLOID_COUNT; d++) {
                            HashMap<Integer, Variant> seqDiploidSNPList = seqDiploidSNPLists.get(d);
                            if (seqDiploidSNPList != null) {
                                geneDiploidSNPLists.set(d, getGeneSNPs(seqDiploidSNPList, theGFFEntry.getStart() - GFF.GENE_FLANKING_REGION_LENGTH, theGFFEntry.getEnd() + GFF.GENE_FLANKING_REGION_LENGTH, false));
                            } else {
                                geneDiploidSNPLists.set(d, null);
                            }
                        }

                    } else {
                        theGFFEntry = null;
                        break;
                    }
                }

                // the alignment is after the gene or the alignment starts before the end of this gene. for latter, check if it ends before this gene start (outside a gene) or starts within the current gene but ends outside the gene. 
                if (theGFFEntry == null) {
                    // this read pair lies after the last gene of this sequence. Ignore it.
                    continue;
                } else if (samRecord.getAlignmentEnd() < theGFFEntry.getStart()) {
                    // this read starts before a valid gene. Ignore it.
//AMT-CHECK                    withinValidGene = false;
                    continue;
                } else if (samRecord.getAlignmentEnd() > theGFFEntry.getEnd() + GFF.GENE_FLANKING_REGION_LENGTH) { // allow some extra bases after the gene
                    // this read ends outside a valid gene. Ignore it.
//AMT-CHECK                    withinValidGene = false;
                    continue;
                }

            } // end if theGFFEntry == null || ...

            // the alignment starts before the end of this gene. Now check if it ends before this gene start (outside a gene) or starts within the current gene but ends outside the gene. 
            if (samRecord.getAlignmentEnd() < theGFFEntry.getStart()) {
                // this read starts before a valid gene. Ignore it.
//AMT-CHECK                withinValidGene = false;
                continue;
            } else if (samRecord.getAlignmentEnd() > theGFFEntry.getEnd() + GFF.GENE_FLANKING_REGION_LENGTH) { // allow some extra bases after the gene
                // this read ends outside a valid gene. Ignore it.
//AMT-CHECK                withinValidGene = false;
                continue;
            }

            // reset the flag
//AMT-CHECK            withinValidGene = true;
            if (geneHSPList == null || geneHSPList.isEmpty()) { // no HSPs found in the gene. No need to process the alignment
                continue;
            }

            // get the HSPs that lie within the region of this alignment
            LinkedHashSet<Variant> alignmentHSPList = getSNPs(geneHSPList, samRecord.getAlignmentStart(), samRecord.getAlignmentEnd());
            // check if any HSPs are present within the alignment boundaries 
            TreeSet<Base> alignmentBasePattern;
            if (alignmentHSPList.isEmpty()) { // no -- create an empty pattern
                alignmentBasePattern = new TreeSet();
            } else { // yes -- get the SNPs in the alignment read
                alignmentBasePattern = getMismatchingBases(alignmentHSPList, samRecord);
                if (alignmentBasePattern == null) { // poor data present at one of the Variant positions .. ignore the alignment
                    //alignmentBasePattern = new TreeSet();
                    continue;
                }
            }

            // check if the read's mate is mapped or not
            TreeSet<Base> mateBasePattern;
            if (!samRecord.getReadPairedFlag() || samRecord.getMateUnmappedFlag()) { // unpaired data OR paired data with unmapped mate. add to the list of base patterns to be processed

                if (alignmentBasePattern.size() > 0) { // check if we have the list of valid SNPs present in this read                    
                    // Check if this snp pattern has been seen before
                    Integer frequency = basePatterns.get(alignmentBasePattern);
                    if (frequency == null) { // if not, add it to the list of snp patterns
                        basePatterns.put(alignmentBasePattern, 1);
                    } else { // otherwise increment the number of reads containing this Variant pattern
                        basePatterns.put(alignmentBasePattern, frequency + 1);
                    }
                } // if Read Variant list size > 0 
            } else if ((mateBasePattern = pendingBasePatterns.get(samRecord.getReadName())) != null) { // see if the read's mate has been seen before
                // Yes. Add them to the list of base patterns to be processed.

                // remove the base pattern data from the list of pending base patterns
                pendingBasePatterns.remove(samRecord.getReadName());

                // check if there is any descrepancy in the SNPs in the two reads. This will only happen if the reads overlap.
                // the mate Battern goes first because the mate is actually mapped before this alignment                
                if (checkDiscrepancyInReadPair(mateBasePattern, alignmentBasePattern)) { // ignore the base patterns if there is a discrepancy
                    continue;
                }

                // No discrepancy. Add Read 2 SNPs (alignment) to the Read 1 (mate) Variant list to get a combined Variant list
                mateBasePattern.addAll(alignmentBasePattern);

                if (mateBasePattern.size() > 0) {
                    // we have the list of valid SNPs present in this read pair.

                    // Check if this snp pattern has been seen before
                    Integer frequency = basePatterns.get(mateBasePattern);
                    if (frequency == null) { // if not, add it to the list of snp patterns
                        basePatterns.put(mateBasePattern, 1);
                    } else { // otherwise increment the number of reads containing this Variant pattern
                        basePatterns.put(mateBasePattern, frequency + 1);
                    }
                } // if Read Variant list size > 0 

            } else { //no, add it to the list of pending alignments
                pendingBasePatterns.put(samRecord.getReadName(), alignmentBasePattern);
            }
        } // end while

        // process the final set of alignments 
//AMT-CHECK  if (theGFFEntry != null && withinValidGene) {
        if (theGFFEntry != null) {
            processBasePatterns(basePatterns, geneHSPList, geneDiploidSNPLists, theGFFEntry.getSequenceName(), outList, outputVCF);
            // reset the variables
            basePatterns.clear();
            pendingBasePatterns.clear();
        }

        return true;

    }

    private static LinkedHashMap<Integer, Variant> getGeneSNPs(HashMap<Integer, Variant> seqVariantList, int start, int end, boolean onlyHSPs) {

        LinkedHashMap<Integer, Variant> snpList = new LinkedHashMap();

        for (int i = start; i <= end; i++) {
            Variant variant = seqVariantList.get(i);
            if (variant != null) {
                if (onlyHSPs && !DNA.isExtendedNucleotide(variant.getNucleotide())) { // if onlyHSPs are to be returned and this is an unambiguous position then ignore it
                    continue;
                }
                snpList.put(i, variant); //add to the actual position (1-based)
            }
        }

        return snpList;
    }

    public static TreeSet<Base> getMismatchingBases(LinkedHashSet<Variant> readRegionSNPList, SAMRecord samRecord) { //, int readStart, String MD, String readSequence, String baseQualities, String cigar) {

        String readSequence = samRecord.getReadString();
        String baseQualities = samRecord.getBaseQualityString();
        // returns null if bad quality snps are present. also ignores SNPs found in deletion
        TreeSet<Base> mismatchingBases = new TreeSet();
        Iterator<Variant> itReadRegionSNPList = readRegionSNPList.iterator();
        while (itReadRegionSNPList.hasNext()) {
            // get the next Variant
            Variant variant = itReadRegionSNPList.next();

            // Get the relative position of this SNP on the read
            int variantRelativePosition = samRecord.getReadPositionAtReferencePosition(variant.getPosition()) - 1; // We need 0-based indexing            

            // ignore this base if there's is a deletion in this read at SNP position (SNPRelativePosition == -1)
            if (variantRelativePosition == -1) {
                continue;
            }

            // return null this base if the sequencing quality is low
            int baseQuality = (int) (baseQualities.charAt(variantRelativePosition)) - 33; // base qualities are ascii - 33
            if (baseQuality < BASE_QUALITY_THRESHOLD) {
                return null;
            }

            // Get the base on the read
            char readBase = readSequence.charAt(variantRelativePosition);

            TreeSet bases = (TreeSet) DNA.getNucleotideMap().get(variant.getNucleotide());
            // ignore this base if it's not part of the consensus
            if (bases == null || !bases.contains(readBase)) {
                continue;
            }

            // Add the base to this read's list of mismatching bases
            mismatchingBases.add(new Base(variant.getPosition(), readBase));
        }

        return mismatchingBases;
    }

    private static LinkedHashSet<Variant> getSNPs(LinkedHashMap<Integer, Variant> theSNPList, int start, int end) {

        LinkedHashSet<Variant> theSNPSubList = new LinkedHashSet();
        Iterator<Map.Entry<Integer, Variant>> itSNP = theSNPList.entrySet().iterator();
        while (itSNP.hasNext()) {
            Map.Entry<Integer, Variant> entry = itSNP.next();
            int position = entry.getKey();

            if (position < start) {
                // do nothing
            } else if (position > end) {
                break;
            } else {
                theSNPSubList.add(entry.getValue());
            }
        }

        return theSNPSubList;
    }

    private static boolean checkDiscrepancyInReadPair(TreeSet<Base> alignment1BasePattern, TreeSet<Base> alignment2BasePattern) {

        Iterator<Base> itRead1SNPList = alignment1BasePattern.iterator();
        Iterator<Base> itRead2SNPList = alignment2BasePattern.iterator();
        Base base1;
        Base base2;
        if (itRead1SNPList.hasNext()) {
            base1 = itRead1SNPList.next();
        } else {
            return false;
        }
        if (itRead2SNPList.hasNext()) {
            base2 = itRead2SNPList.next();
        } else {
            return false;
        }
        while (true) {

            if (base1.getPosition() < base2.getPosition()) {
                // get next Variant
                if (itRead1SNPList.hasNext()) {
                    base1 = itRead1SNPList.next();
                } else {
                    break;
                }
            } else if (base1.getPosition() == base2.getPosition()) { // same position
                if (base1.getNucleotide() == base2.getNucleotide()) { //same Variant
                    if (itRead1SNPList.hasNext() && itRead2SNPList.hasNext()) {
                        base1 = itRead1SNPList.next();
                        base2 = itRead2SNPList.next();
                    } else {
                        break;
                    }
                } else {
                    return true;
                }
            } else if (base1.getPosition() > base2.getPosition()) {
                if (itRead2SNPList.hasNext()) {
                    base2 = itRead2SNPList.next();
                } else {
                    break;
                }
            }
        } // end while

        return false;
    }

    private static void processBasePatterns(LinkedHashMap<TreeSet<Base>, Integer> basePatterns, LinkedHashMap<Integer, Variant> geneHSPList,
            ArrayList<LinkedHashMap> geneDiploidSNPLists, String sequenceName, ArrayList<BufferedWriter> outList, boolean outputVCF) {

        if (basePatterns.isEmpty()) {
            return;
        }

//         printBasePatterns(basePatterns);
        // filter the base patterns to remove unwanted patterns
        filterBasePatterns(basePatterns, geneHSPList);

        //printBasePatterns(basePatterns);
        // sort the base patterns by size (longest patterns first)
        LinkedHashSet<TreeSet<Base>> sortedBasePatterns = sortBasePatternsBySize(basePatterns);

        // remove base patterns embedded in another pattern
        removeEmbeddedBasePatterns(sortedBasePatterns);

        //       printBasePatterns(sortedBasePatterns);
        if (MERGE_BASE_PATTERNS) {
//       System.out.println ("Merging base patterns");
            sortedBasePatterns = mergeBasePatterns(sortedBasePatterns);
//        printBasePatterns(sortedBasePatterns);
        }

        // assign patterns to 1st, 2nd and/or 3rd genomes
        ArrayList<LinkedHashMap<TreeSet<Base>, Double>> basePatternsList = new ArrayList<LinkedHashMap<TreeSet<Base>, Double>>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            basePatternsList.add(new LinkedHashMap());
        }
        // first false: ignoreUnassingedOrAmbiguousPositions, second false: assignToMissingDiploid (forcefully)
        assignBasePatternsToGenomes(sortedBasePatterns, geneDiploidSNPLists, basePatternsList, false, false);

        // assign genome bases
        ArrayList<LinkedHashMap> baseLists = new ArrayList<LinkedHashMap>(DIPLOID_COUNT);
        for (LinkedHashMap<TreeSet<Base>, Double> genomeBasePatterns : basePatternsList) {
            baseLists.add(assignBasesToGenome(genomeBasePatterns));
        }

        // get the list of base patterns not assigned to any genome
        LinkedHashSet<TreeSet<Base>> unassignedBasePatterns = new LinkedHashSet(sortedBasePatterns);
        for (LinkedHashMap<TreeSet<Base>, Double> genomeBasePatterns : basePatternsList) {
            unassignedBasePatterns.removeAll(genomeBasePatterns.keySet());
        }

        ArrayList<Boolean> genomeFinalised = new ArrayList<Boolean>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            genomeFinalised.add(false);
        }

        if (RECTIFY_USING_REFERENCE) { // correction using reference genome
            // check which genome(s) contributed to the reference sequence
            // we use a set as multiple genome might have the maximum identity with the reference
            HashSet referenceGenome = checkReferenceGenome(geneHSPList, baseLists);

            for (int d = 0; d < DIPLOID_COUNT; d++) {
                if (referenceGenome.contains(d)) { // reference comes from this genome 
                    // rectify the Variant list for this genome using reference base list
                    rectifyBaseListUsingReferenceGenome(geneHSPList, geneDiploidSNPLists.get(d), baseLists.get(d));
                    // remove the patterns which no longer contradict the base list (i.e. no longer unassigned)
                    checkBasePatternsForConsistency(unassignedBasePatterns, baseLists.get(d), true);
                    genomeFinalised.set(d, true);
                }
            }
        }

        // assign unassigned Variant patterns to genomes
        ArrayList<LinkedHashMap<TreeSet<Base>, Double>> unassignedBasePatternsList = new ArrayList<LinkedHashMap<TreeSet<Base>, Double>>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            unassignedBasePatternsList.add(new LinkedHashMap());
        }
        if (MISSING_DIPLOID > 0) { // if a genome is marked as missing, assign all unassigned base patterns to that genome
            // first true: ignoreUnassingedOrAmbiguousPositions (irrelevant here), second true: assignToMissingDiploid (forcefully)
            assignBasePatternsToGenomes(unassignedBasePatterns, baseLists, unassignedBasePatternsList, true, true);
            if (TRY_BASE_SET_CLOSURE) { // perform base set closure, if required
                tryBaseSetClosure(geneHSPList, baseLists, baseLists.get(MISSING_DIPLOID - 1));
            }
        } else { // otherwise, assign them to appropriate genomes based on similarity
            // first true: ignoreUnassingedOrAmbiguousPositions, second false: assignToMissingDiploid (forcefully)
            assignBasePatternsToGenomes(unassignedBasePatterns, baseLists, unassignedBasePatternsList, true, false);
            if (TRY_BASE_SET_CLOSURE) { // perform base set closure, if required
                tryBaseSetClosure(geneHSPList, baseLists);
            }
        }

        // finalise the base lists, if need and then write the base list
        Iterator<LinkedHashMap> itBaseLists = baseLists.iterator();
        Iterator<LinkedHashMap<TreeSet<Base>, Double>> itBasePatternsList = basePatternsList.iterator();
        Iterator<LinkedHashMap<TreeSet<Base>, Double>> itUnassignedBasePatternsList = unassignedBasePatternsList.iterator();
        Iterator<Boolean> itGenomeFinalised = genomeFinalised.iterator();
        Iterator<LinkedHashMap> itGeneDiploidSNPLists = geneDiploidSNPLists.iterator();
        Iterator<BufferedWriter> itOutList = outList.iterator();
        while (itBaseLists.hasNext()) {
            LinkedHashMap<Integer, Base> baseList = itBaseLists.next();
            if (!itGenomeFinalised.next()) {
                finaliseBaseList(itUnassignedBasePatternsList.next(), baseList, itBasePatternsList.next());
            }
            // commented on 2-jan-2015 1330 - writing is now done once all the processing is complete. This will allow futher processing/validation, if required.
            // write the output
            //writeBases(baseList, geneHSPList, geneReferenceBaseList, itGeneDiploidSNPLists.next(), sequenceName, itOutList.next());
        }

        LinkedHashMap<Integer, Character> hspBasesCovered = checkIfAllHSPBasesAreCovered(geneHSPList, baseLists);
        // added on 2-jan-2015 1330 - writing is now done once all the processing is complete. This will allow futher processing/validation, if required.
        // write the output for each base list
        itBaseLists = baseLists.iterator();
        int d = 0;
        boolean isMissingDiploid;
        while (itBaseLists.hasNext()) {
            // write the output for the next base list
            if (MISSING_DIPLOID > 0 && d == (MISSING_DIPLOID - 1)) {
                isMissingDiploid = true;
            } else {
                isMissingDiploid = false;
            }
            writeBases(itBaseLists.next(), geneHSPList, itGeneDiploidSNPLists.next(), sequenceName, itOutList.next(), hspBasesCovered, isMissingDiploid, outputVCF);
            d++;
        }
    }

    private static void filterBasePatterns(LinkedHashMap<TreeSet<Base>, Integer> basePatterns, LinkedHashMap<Integer, Variant> geneHSPList) {

        // get Variant pairs along with the number of reads supporting each pair
        HashMap<BasePair, Integer> basePairList = getNucleotidePairs(basePatterns);
        // get Variant position pairs along with the number of reads supporting each position pair
        HashMap<PositionPair, Integer> positionPairList = getPositionPairs(basePatterns);

        HashSet basePatternsToBeDeleted = new HashSet();
        Iterator<Map.Entry<TreeSet<Base>, Integer>> itBasePatterns = basePatterns.entrySet().iterator();
        while (itBasePatterns.hasNext()) {
            // get the next entry
            Map.Entry<TreeSet<Base>, Integer> entry = itBasePatterns.next();

            // Get the next Variant pattern
            TreeSet basePattern = entry.getKey();

            // check if is a single position pattern
            if (basePattern.size() == 1) {
                // Yes! check the threshold
                if (entry.getValue() < SINGLE_POSITION_SNP_PATTERN_COUNT_THRESHOLD) {
                    basePatternsToBeDeleted.add(basePattern);
                }
                continue;
            }
            // now check the multi position base patterns

            // get the iterator over this Base pattern
            Iterator<Base> itBase = basePattern.iterator();
            // Get the first Variant in the Base Pattern
            Base previousBase = itBase.next();
            // process all remaining SNPs
            while (itBase.hasNext()) {
                Base base = itBase.next();

                // create a new Variant position pair using the previous and current SNPs
                PositionPair positionPair = new PositionPair(previousBase.getPosition(), base.getPosition());

                if (consecutivePositions(positionPair, geneHSPList)) {

                    // create a new Variant pair using the previous and current SNPs
                    BasePair basePair = new BasePair(previousBase, base);

                    // get the readcount for this Variant pair
                    int pairReadCount = basePairList.get(basePair);

                    // create a new Variant position pair using the previous and current SNPs and get the read count for this snp position pair
                    int positionReadCount = positionPairList.get(positionPair);

                    if ((double) pairReadCount / (double) positionReadCount < SNP_PAIR_PROPORTION_THRESHOLD) {
                        basePatternsToBeDeleted.add(basePattern);
                        break;
                    }
                }

                // save this Variant as previous Variant
                previousBase = base;
            } // end while itSNPPattern.hasNext()

        } // it.hasNext()

        // remove the groups which have been marked for deletion
        basePatterns.keySet().removeAll(basePatternsToBeDeleted);

    }

    private static HashMap getNucleotidePairs(LinkedHashMap<TreeSet<Base>, Integer> basePatterns) {

        HashMap<BasePair, Integer> snpPairList = new HashMap();
        Iterator<Map.Entry<TreeSet<Base>, Integer>> itBasePattern = basePatterns.entrySet().iterator();
        while (itBasePattern.hasNext()) {
            // Get the current Variant pattern and corresponding read counts
            Map.Entry<TreeSet<Base>, Integer> entry = itBasePattern.next();

            TreeSet<Base> basePattern = entry.getKey();
            int readCount = entry.getValue();

            Base previousBase = null;
            // get the iterator over this Variant pattern
            Iterator<Base> itBase = basePattern.iterator();
            // Get the first Variant in the Variant Pattern
            if (itBase.hasNext()) {
                previousBase = itBase.next();
            }
            // process all remaining SNPs
            while (itBase.hasNext()) {
                Base base = itBase.next();

                // create a new Variant pair using the previous and current SNPs
                BasePair basePair = new BasePair(previousBase, base);

                // get the entry from the list of snp pairs
                Integer currentReadCount = snpPairList.get(basePair);
                if (currentReadCount == null) {
                    // add it to the list of Variant pairs
                    snpPairList.put(basePair, readCount);
                } else {
                    // add it to the list of Variant pairs
                    snpPairList.put(basePair, currentReadCount + readCount);
                }
                // save this Variant as previous Variant
                previousBase = base;
            }
        } // end while (snpPatterns)

        return snpPairList;
    }

    private static HashMap getPositionPairs(LinkedHashMap<TreeSet<Base>, Integer> basePatterns) {

        HashMap<PositionPair, Integer> positionPairList = new HashMap();
        Iterator<Map.Entry<TreeSet<Base>, Integer>> itBasePattern = basePatterns.entrySet().iterator();
        while (itBasePattern.hasNext()) {
            // Get the current Variant pattern and corresponding read counts
            Map.Entry<TreeSet<Base>, Integer> entry = itBasePattern.next();

            TreeSet<Base> basePattern = entry.getKey();
            int readCount = entry.getValue();

            // get the iterator over this Variant pattern
            Iterator<Base> itBase = basePattern.iterator();
            // Get the first Variant in the Variant Pattern
            Base previousBase = itBase.next();
            // process all remaining SNPs
            while (itBase.hasNext()) {
                Base base = itBase.next();

                // create a new Variant pair using the previous and current SNPs
                PositionPair positionPair = new PositionPair(previousBase.getPosition(), base.getPosition());

                // get the entry from the list of snp pairs
                Integer currentReadCount = positionPairList.get(positionPair);
                if (currentReadCount == null) {
                    // add it to the list of Variant pairs
                    positionPairList.put(positionPair, readCount);
                } else {
                    // add it to the list of Variant pairs
                    positionPairList.put(positionPair, currentReadCount + readCount);
                }

                // save this Variant as previous Variant
                previousBase = base;
            }
        } // end while (snpPatterns)

        return positionPairList;
    }

    private static boolean consecutivePositions(PositionPair positionPair, LinkedHashMap<Integer, Variant> geneHSPList) {
        Iterator<Integer> itPosition = geneHSPList.keySet().iterator();
        while (itPosition.hasNext()) {
            int position = itPosition.next();
            if (positionPair.getPosition1() == position) {
                if (itPosition.hasNext()) {
                    position = itPosition.next();
                    return positionPair.getPosition2() == position;
                } else {
                    return false;
                }
            } else if (positionPair.getPosition1() < position) {
                // ideally this should not happen
                return false;
            }
        }
        return false;
    }

    private static LinkedHashSet<TreeSet<Base>> sortBasePatternsBySize(LinkedHashMap<TreeSet<Base>, Integer> basePatterns) {
        /* Sort the base patterns such that largest pattern is at the begining */

        TreeMap<Integer, LinkedHashSet<TreeSet<Base>>> basePatternsBySize = new TreeMap();
        Iterator<TreeSet<Base>> itBasePatterns = basePatterns.keySet().iterator();
        while (itBasePatterns.hasNext()) {
            // get the next base pattern
            TreeSet basePattern = itBasePatterns.next();

            // get the list of patterns of this size
            LinkedHashSet<TreeSet<Base>> sizeBasePatterns = basePatternsBySize.get(basePattern.size());
            if (sizeBasePatterns == null) { // if no pattern found, create a new list
                sizeBasePatterns = new LinkedHashSet();
            }

            // add this pattern to the list of patterns
            sizeBasePatterns.add(basePattern);
            basePatternsBySize.put(basePattern.size(), sizeBasePatterns);
        } // end while()

        // now order them descendingly according to their size
        LinkedHashSet<TreeSet<Base>> sortedBasePatterns = new LinkedHashSet();
        Iterator<Integer> itSize = basePatternsBySize.descendingKeySet().iterator();
        while (itSize.hasNext()) {
            sortedBasePatterns.addAll(basePatternsBySize.get(itSize.next()));
        }

        return sortedBasePatterns;
    }

    private static LinkedHashSet<TreeSet<Base>> sortBasePatternsBySize(LinkedHashSet<TreeSet<Base>> basePatterns) {
        /* Sort the base patterns such that largest pattern is at the begining */

        TreeMap<Integer, LinkedHashSet<TreeSet<Base>>> basePatternsBySize = new TreeMap();
        Iterator<TreeSet<Base>> itBasePatterns = basePatterns.iterator();
        while (itBasePatterns.hasNext()) {
            // get the next base pattern
            TreeSet basePattern = itBasePatterns.next();

            // get the list of patterns of this size
            LinkedHashSet<TreeSet<Base>> sizeBasePatterns = basePatternsBySize.get(basePattern.size());
            if (sizeBasePatterns == null) { // if no pattern found, create a new list
                sizeBasePatterns = new LinkedHashSet();
            }

            // add this pattern to the list of patterns
            sizeBasePatterns.add(basePattern);
            basePatternsBySize.put(basePattern.size(), sizeBasePatterns);
        } // end while()

        // now order them descendingly according to their size
        LinkedHashSet<TreeSet<Base>> sortedBasePatterns = new LinkedHashSet();
        Iterator<Integer> itSize = basePatternsBySize.descendingKeySet().iterator();
        while (itSize.hasNext()) {
            sortedBasePatterns.addAll(basePatternsBySize.get(itSize.next()));
        }

        return sortedBasePatterns;
    }

    private static void removeEmbeddedBasePatterns(LinkedHashSet<TreeSet<Base>> basePatterns) {

        HashSet basePatternsToBeDeleted = new HashSet();
        Iterator<TreeSet<Base>> itBasePatterns1 = basePatterns.iterator();
        while (itBasePatterns1.hasNext()) {
            // get the next Variant group
            TreeSet<Base> basePattern1 = itBasePatterns1.next();

            // move to the next pattern if this group has been marked as inactive
            if (basePatternsToBeDeleted.contains(basePattern1)) {
                continue;
            }

            // get the first and last Variant position in this group
            int firstBasePositionGroup1 = basePattern1.first().getPosition();
            int lastBasePositionGroup1 = basePattern1.last().getPosition();

            Iterator<TreeSet<Base>> itBasePatterns2 = basePatterns.iterator();
            while (itBasePatterns2.hasNext()) {
                TreeSet<Base> basePattern2 = itBasePatterns2.next();
                if (basePattern1 == basePattern2) {
                    break;
                }
            }
            while (itBasePatterns2.hasNext()) {
                // get the Variant pattern
                TreeSet<Base> basePattern2 = itBasePatterns2.next();

                if (lastBasePositionGroup1 < basePattern2.first().getPosition()) {// if the group 1 ends before group 2 starts, no overlap possible
                    continue;
                } else if (basePattern2.last().getPosition() < firstBasePositionGroup1) {// if the group 2 ends before group 1 starts, no overlap possible
                    continue;
                }

                int overlappingSNPs = 0;
                Iterator<Base> itBasePattern1 = basePattern1.iterator();
                Iterator<Base> itBasePattern2 = basePattern2.iterator();
                Base base1 = itBasePattern1.next();
                Base base2 = itBasePattern2.next();
                while (true) {
                    if (base1.getPosition() == base2.getPosition()) { // same position
                        if (base1.getNucleotide() == base2.getNucleotide()) { //same Variant
                            overlappingSNPs++;
                            // get next SNPs
                            if (itBasePattern1.hasNext() && itBasePattern2.hasNext()) {
                                base1 = itBasePattern1.next();
                                base2 = itBasePattern2.next();
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    } else if (base1.getPosition() < base2.getPosition()) {
                        if (itBasePattern1.hasNext()) {
                            base1 = itBasePattern1.next();
                        } else {
                            break;
                        }
                    } else if (base1.getPosition() > base2.getPosition()) {
                        if (itBasePattern2.hasNext()) {
                            base2 = itBasePattern2.next();
                        } else {
                            break;
                        }
                    }
                } // end while

                if (overlappingSNPs == basePattern2.size()) { // all SNPs are embedded in group 1
                    // add it to the groups to be deleted
                    basePatternsToBeDeleted.add(basePattern2);
                }
            } // end while itSNPGroups2.hasNext()
        }// end while itSNPGroups1.hasNext()

        // remove the groups which were found to be embedded in other groups
        basePatterns.removeAll(basePatternsToBeDeleted);
    }

    private static LinkedHashSet<TreeSet<Base>> mergeBasePatterns(LinkedHashSet<TreeSet<Base>> basePatterns) {

        TreeSet<Base> selectedBasePattern1;
        LinkedHashSet<TreeSet<Base>> selectedBasePatterns2 = new LinkedHashSet();
        HashSet basePatternsToBeDeleted = new HashSet();

        while (true) {
            int maxOverlapSize = 0;
            basePatternsToBeDeleted.clear();
            selectedBasePattern1 = null;
            selectedBasePatterns2.clear();

            //printBasePatterns(basePatterns);
            Iterator<TreeSet<Base>> itBasePatterns1 = basePatterns.iterator();
            while (itBasePatterns1.hasNext()) {
                // get the next Variant group
                TreeSet<Base> basePattern1 = itBasePatterns1.next();

                // move to the next pattern if this group has been marked as inactive
                if (basePatternsToBeDeleted.contains(basePattern1)) {
                    continue;
                }

                // get the first and last Variant position in this group
                int firstSNPPositionGroup1 = basePattern1.first().getPosition();
                int lastBasePositionGroup1 = basePattern1.last().getPosition();

                Iterator<TreeSet<Base>> itBasePatterns2 = basePatterns.iterator();
                while (itBasePatterns2.hasNext()) {
                    TreeSet basePattern2 = itBasePatterns2.next();
                    if (basePattern1 == basePattern2) {
                        break;
                    }
                }
                while (itBasePatterns2.hasNext()) {
                    // get the Variant pattern
                    TreeSet<Base> basePattern2 = itBasePatterns2.next();

                    // get the first Variant position in this group
                    int firstBasePositionGroup2 = basePattern2.first().getPosition();
                    if (lastBasePositionGroup1 < firstBasePositionGroup2) {// if the group 1 ends before group 2 starts, no overlap possible
                        continue;
                    } else if (basePattern2.last().getPosition() < firstSNPPositionGroup1) {// if the group 2 ends before group 1 starts, no overlap possible
                        continue;
                    }

                    int overlapSize = 0;
                    Iterator<Base> itBasePattern1 = basePattern1.iterator();
                    Iterator<Base> itBasePattern2 = basePattern2.iterator();
                    Base base1 = itBasePattern1.next();
                    Base base2 = itBasePattern2.next();
                    while (true) {
                        if (base1.getPosition() == base2.getPosition()) { // same position
                            if (base1.getNucleotide() == base2.getNucleotide()) { //same Variant
                                overlapSize++;
                                // get next SNPs
                                if (itBasePattern1.hasNext() && itBasePattern2.hasNext()) {
                                    base1 = itBasePattern1.next();
                                    base2 = itBasePattern2.next();
                                } else {
                                    break;
                                }
                            } else {
                                overlapSize = 0;
                                break;
                            }
                        } else if (base1.getPosition() < base2.getPosition()) {
                            if (itBasePattern1.hasNext()) {
                                base1 = itBasePattern1.next();
                            } else {
                                break;
                            }
                        } else if (base1.getPosition() > base2.getPosition()) {
                            if (itBasePattern2.hasNext()) {
                                base2 = itBasePattern2.next();
                            } else {
                                break;
                            }
                        }
                    } // end while

                    if (overlapSize > 0) {
                        if (overlapSize == basePattern2.size()) { // all SNPs are embedded in group 1 -- remove this base pattern
                            // add it to the groups to be deleted
                            basePatternsToBeDeleted.add(basePattern2);
                        } else if (overlapSize > maxOverlapSize) { // a better overlap found
                            selectedBasePattern1 = basePattern1;
                            selectedBasePatterns2.clear();
                            selectedBasePatterns2.add(basePattern2);
                            maxOverlapSize = overlapSize;
                        } else if (basePattern1 == selectedBasePattern1 && overlapSize == maxOverlapSize) { // multiple best overlaps
                            selectedBasePatterns2.add(basePattern2);
                        }
                    }
                } // end while itSNPGroups2.hasNext()
            }// end while itSNPGroups1.hasNext()

            if (maxOverlapSize > 0) {
                // merge all group 2's with group 1

                // iterate over all selected base patterns 2
                Iterator<TreeSet<Base>> itBasePatterns2 = selectedBasePatterns2.iterator();
                while (itBasePatterns2.hasNext()) {
                    TreeSet basePattern2 = itBasePatterns2.next();

                    // create a new base pattern using base pattern 1
                    TreeSet newBasePattern = (TreeSet) selectedBasePattern1.clone();
                    // add all SNPs from base pattern 2 to newly created base pattern
                    newBasePattern.addAll(basePattern2);
                    // add the newly created base pattern to the list of base patterns
                    basePatterns.add(newBasePattern);

                    // add base pattern 2 to the groups to be deleted
                    basePatternsToBeDeleted.add(basePattern2);

                }
                // add base pattern 2 to the groups to be deleted
                basePatternsToBeDeleted.add(selectedBasePattern1);

                // remove the groups which were found to be embedded in other groups or merged
                basePatterns.removeAll(basePatternsToBeDeleted);

                // re-sort the base patterns
                basePatterns = sortBasePatternsBySize(basePatterns);

            } else {
                // remove the groups which were found to be embedded in other groups or merged
                basePatterns.removeAll(basePatternsToBeDeleted);
                break;
            }
        } // end while (true)

        return basePatterns;
    }

    private static void assignBasePatternsToGenomes(LinkedHashSet<TreeSet<Base>> basePatterns,
            ArrayList<LinkedHashMap> geneDiploidSNPLists,
            ArrayList<LinkedHashMap<TreeSet<Base>, Double>> basePatternsList,
            boolean ignoreUnassingedOrAmbiguousPositions, boolean assignToMissingDiploid) {

        // iterate through all Variant Patterns
        Iterator<TreeSet<Base>> itBasePattern = basePatterns.iterator();
        while (itBasePattern.hasNext()) {
            int genomeAssignedCount = 0;
            // Get the next Variant pattern
            TreeSet basePattern = itBasePattern.next();

            if (assignToMissingDiploid) {
                if (MISSING_DIPLOID > 0) {
                    basePatternsList.get(MISSING_DIPLOID - 1).put(basePattern, 2.0);
                } else {
                    System.err.println("Error in assingment to missing diploid.");
                }
            } else {
                ArrayList<Double> proportionList = new ArrayList<Double>(DIPLOID_COUNT);
                ArrayList<Boolean> isGenomeList = new ArrayList<Boolean>(DIPLOID_COUNT);

                // check the pattern for each genome
                for (int d = 0; d < DIPLOID_COUNT; d++) {
                    LinkedHashMap<Integer, Base> geneDiploidSNPList = geneDiploidSNPLists.get(d);
                    if (geneDiploidSNPList != null) {
                        double proportion = calculateBasePatternIdentity(basePattern, geneDiploidSNPList, ignoreUnassingedOrAmbiguousPositions);
                        proportionList.add(proportion);
                        // check if the matching proportion meets the threshold
                        boolean isGenome = (proportion >= SNP_PATTERN_MATCHING_THRESHOLD);
                        isGenomeList.add(isGenome);
                        if (isGenome && d != (DISTANT_GENOME - 1)) {
                            genomeAssignedCount++;
                        }
                    } else {
                        proportionList.add(0.0);
                        isGenomeList.add(false);
                    }
                }

                if (DISTANT_GENOME > 0 && genomeAssignedCount == 0 && basePattern.size() > 1) {
                    basePatternsList.get(DISTANT_GENOME - 1).put(basePattern, 2.0);
                } else {
                    for (int d = 0; d < DIPLOID_COUNT; d++) {
                        if (isGenomeList.get(d)) {
                            basePatternsList.get(d).put(basePattern, proportionList.get(d));
                        }
                    }
                }
            }
        } // end while (snpPatterns)

    }

    private static double calculateBasePatternIdentity(TreeSet<Base> basePattern, LinkedHashMap<Integer, Base> snpList, boolean ignoreUnassingedOrAmbigousPositions) {
        int matchingSNPs = 0;
        int positionsToBeIgnored = 0;
        Iterator<Base> itBasePattern = basePattern.iterator();
        while (itBasePattern.hasNext()) {
            // get the next Variant in the pattern
            Base base = itBasePattern.next();
            // get the corresponding Variant in the Variant list
            Base theListSNP = snpList.get(base.getPosition());

            if (theListSNP == null || !DNA.isNucleotide(theListSNP.getNucleotide())) {
                positionsToBeIgnored++;
                continue;
            }

            // see if the snp base in the pattern matches the snp in the list
            if (base.getNucleotide() == theListSNP.getNucleotide()) {
                // if yes, increment the match count
                matchingSNPs++;
            }
        } // end while it.hasNext()

        if (ignoreUnassingedOrAmbigousPositions) {
            return (double) matchingSNPs / (double) basePattern.size();
        } else {
            // return the match proportion
            return (double) matchingSNPs / (double) (basePattern.size() - positionsToBeIgnored);
        }

    }

    private static LinkedHashMap<Integer, Base> assignBasesToGenome(LinkedHashMap<TreeSet<Base>, Double> potentialBasePatterns) {

        LinkedHashMap<Integer, Base> baseList = new LinkedHashMap();

        if (potentialBasePatterns.isEmpty()) {
            return baseList;
        }

        boolean baseListChanged;
        do {

            TreeMap<Integer, HashMap<Base, Double>> potentialBases = new TreeMap();
            // iterate through all Variant Patterns
            Iterator<Map.Entry<TreeSet<Base>, Double>> itSNPPattern = potentialBasePatterns.entrySet().iterator();
            while (itSNPPattern.hasNext()) {
                // get the next entry
                Map.Entry<TreeSet<Base>, Double> entry = itSNPPattern.next();
                // Get the next Variant pattern
                TreeSet snpPattern = entry.getKey();
                // update potential bases
                updatePotentialBases(snpPattern, entry.getValue(), potentialBases, baseList);

            } // end while (snpPatterns)

            // assign bases
            baseListChanged = assignBases(potentialBases, baseList);

            // check the Variant patterns for inconsistency
            checkBasePatternsForConsistency(potentialBasePatterns.keySet(), baseList, false);
        } while (baseListChanged);

        return baseList;
    }

    /* AMT 30-Jul-2015: Use the best proportion or add all the proportions depending on user choice*/
    private static void updatePotentialBases(TreeSet<Base> basePattern, double proportion,
            TreeMap<Integer, HashMap<Base, Double>> potentialBases, HashMap<Integer, Base> baseList) {
        Iterator<Base> itBase = basePattern.iterator();
        while (itBase.hasNext()) {
            Base base = itBase.next();

            if (baseList.containsKey(base.getPosition())) { // position already assigned
                continue;
            }

            // get base proportions for this snp, if present
            HashMap<Base, Double> baseProportions = potentialBases.get(base.getPosition());
            if (baseProportions == null) { // this position has not been seen before
                // create a new map and add the proportion against this base
                baseProportions = new HashMap();
                baseProportions.put(base, proportion);
                potentialBases.put(base.getPosition(), baseProportions);

            } else {
                Double baseProportion = baseProportions.get(base);
                if (BASE_PATTERN_PROPORTION_MODE == BASE_PATTERN_PROPORTION_MODE_ADDITIVE) {
                    if (baseProportion == null) {// this base has not been seen at this position 
                        baseProportion = 0.0;
                    }
                    baseProportions.put(base, proportion + baseProportion);
                    potentialBases.put(base.getPosition(), baseProportions);
                } else if (BASE_PATTERN_PROPORTION_MODE == BASE_PATTERN_PROPORTION_MODE_MAXIMUM) {
                    if (baseProportion == null || proportion > baseProportion) {// this base has not been seen at this position or another proportion found for this base

                        baseProportions.put(base, proportion);
                        potentialBases.put(base.getPosition(), baseProportions);
                    }
                }
            }
        } // end while it.hasNext()

    }

    private static boolean assignBases(TreeMap<Integer, HashMap<Base, Double>> potentialBases, LinkedHashMap<Integer, Base> baseList) {

        boolean baseListChanged = false;
        Iterator<Map.Entry<Integer, HashMap<Base, Double>>> itPotentialBases = potentialBases.entrySet().iterator();
        while (itPotentialBases.hasNext()) {
            // get the next entry
            Map.Entry<Integer, HashMap<Base, Double>> entry = itPotentialBases.next();

            // get the position
            int position = entry.getKey();

            if (baseList.containsKey(position)) { // position already assinged
                continue;
            }

            // and the base proportions
            HashMap baseProportions = entry.getValue();
            Iterator<Map.Entry<Base, Double>> itBaseProportions = baseProportions.entrySet().iterator();

            if (baseProportions.size() == 1) { // only one base found at this position
                // add it to the list
                baseList.put(position, itBaseProportions.next().getKey());
                // mark the base list as changed
                baseListChanged = true;
            } else {
                HashSet<Base> selectedBases = new HashSet();
                double maximumProportion = 0.0;
                // assign the base with maximum  proprotion
                while (itBaseProportions.hasNext()) {
                    Map.Entry<Base, Double> entryBP = itBaseProportions.next();

                    Base base = entryBP.getKey();
                    double proportion = entryBP.getValue();

                    if (proportion > maximumProportion) {
                        // save this base as the selected base;
                        selectedBases.clear();
                        selectedBases.add(base);
                        // also update the maximum proportion
                        maximumProportion = proportion;
                    } else if (proportion == maximumProportion) {
                        selectedBases.add(base);
                    }
                } // end while (itBaseProportions.hasNext())
                // add it to the list

                if (selectedBases.size() == 1) {
                    baseList.put(position, selectedBases.iterator().next());
                    // mark the base list as changed
                    baseListChanged = true;
                }
            }
        } // end while itPotentialBases.hasNext()

        return baseListChanged;
    }

    private static void checkBasePatternsForConsistency(Set<TreeSet<Base>> basePatterns, HashMap<Integer, Base> baseList, boolean checkReverse) {

        if (baseList.isEmpty()) {
            return;
        }

        HashSet basePatternsToBeDeleted = new HashSet();
        // check all snps in this pattern
        Iterator<TreeSet<Base>> itBasePattern = basePatterns.iterator();
        while (itBasePattern.hasNext()) {
            // get the next Variant pattern
            TreeSet basePattern = itBasePattern.next();

            // check if the base pattern is consistent with the already assigned bases
            boolean isConsistent = checkBasePattern(basePattern, baseList);
            if (checkReverse) {
                if (isConsistent) {
                    // Yes but we are checking for reverse. mark it for deletion
                    basePatternsToBeDeleted.add(basePattern);
                }
            } else if (!isConsistent) {
                // if not, mark it for deletion
                basePatternsToBeDeleted.add(basePattern);
            }
        } // end while it.hasNext()

        // remove the Variant patterns marked for deletion
        basePatterns.removeAll(basePatternsToBeDeleted);
    }

    private static boolean checkBasePattern(TreeSet<Base> basePattern, HashMap<Integer, Base> baseList) {

        if (baseList.isEmpty()) {
            return true;
        }
        // check all snps in this pattern
        Iterator<Base> itSNP = basePattern.iterator();
        while (itSNP.hasNext()) {
            // get the next Variant
            Base theBase = itSNP.next();

            // get the Variant assigned at this position, if present
            Base assignedBase = baseList.get(theBase.getPosition());
            if (assignedBase != null) { // only check assigned positions
                if (theBase.getNucleotide() != assignedBase.getNucleotide()) { // the base does not match the assigned base
                    return false;
                }
            }
        } // end while it.hasNext()

        return true;
    }

    private static HashSet<Integer> checkReferenceGenome(LinkedHashMap<Integer, Variant> geneHSPList, ArrayList<LinkedHashMap> baseLists) {

        HashSet<Integer> referenceGenome = new HashSet<Integer>(DIPLOID_COUNT); // we use hashset as multiple genome might have the maximum identity with the reference

        // calculate identity of each genome with the reference
        double maxIdentity = 0.0;
        ArrayList<Double> identityList = new ArrayList<Double>(DIPLOID_COUNT);
        for (int d = 0; d < DIPLOID_COUNT; d++) {
            double identity = calculateReferenceIdentity(geneHSPList, baseLists.get(d));
            identityList.set(d, identity);
            // also check for max identity
            maxIdentity = Math.max(identity, maxIdentity);
        }

        for (int d = 0; d < DIPLOID_COUNT; d++) {
            if (identityList.get(d) == maxIdentity) {
                referenceGenome.add(d);
            }
        }

        return referenceGenome;

    }

    private static double calculateReferenceIdentity(LinkedHashMap<Integer, Variant> geneHSPList, LinkedHashMap<Integer, Base> baseList) {

        int identity = 0;
        // go through each HSP position for this gene
        Iterator<Integer> itHSPPosition = geneHSPList.keySet().iterator();
        while (itHSPPosition.hasNext()) {
            // get the next position
            int position = itHSPPosition.next();
            // get the base from the base list (if any) at the current position
            Base theBase = baseList.get(position);
            if (theBase == null) { // snp was not assigned for this position, move to the next one
                continue;
            }
            // check if the base is same
            if (geneHSPList.get(position).getReferenceBase() == theBase.getNucleotide()) {
                identity++;
            }

        } // while (it.hasNext())

        return (double) identity / (double) geneHSPList.size();

    }

    private static void rectifyBaseListUsingReferenceGenome(LinkedHashMap<Integer, Variant> geneHSPList, LinkedHashMap<Integer, Variant> geneDiploidSNPList,
            LinkedHashMap<Integer, Base> baseList) {

        // go through each HSP position for this gene
        Iterator<Integer> itHSPPosition = geneHSPList.keySet().iterator();
        while (itHSPPosition.hasNext()) {
            // get the next position
            int position = itHSPPosition.next();

            // get the base from the diploid base list at the current position
            Variant theDiploidBase = geneDiploidSNPList.get(position);
            if (theDiploidBase == null || !DNA.isNucleotide(theDiploidBase.getNucleotide())) { // dont correct for positions where there is no coverage
                continue;
            }

            // get the base from the reference base list at the current position
            if (theDiploidBase.getReferenceBase() == DNA.BASE_N) { // ignore positions with N in the reference
                continue;
            }

            Base theBase = baseList.get(position);
            if (theBase == null) {
                // snp was not assigned for this position or incorrectly assigned.
                // Replace it with reference base
                baseList.get(position).setNucleotide(theDiploidBase.getReferenceBase());
            }
        } // while (it.hasNext())

    }

    private static void finaliseBaseList(LinkedHashMap<TreeSet<Base>, Double> basePatterns, HashMap<Integer, Base> baseList, LinkedHashMap<TreeSet<Base>, Double> assignedSNPPatterns) {

        // assign genome SNPs
        LinkedHashMap<Integer, Base> newBaseList = assignBasesToGenome(basePatterns);
        // update the previous base list using the new Base list

        for (Map.Entry entry : newBaseList.entrySet()) {
            // get the position
            int position = (Integer) entry.getKey();
            // and the Variant
            Base newBase = (Base) entry.getValue();
            // get the Variant in the original Variant List
            Base base = (Base) baseList.get(position);
            if (base == null) { // snp not previously assigned
                baseList.put(position, newBase);
            } else if (!base.equals(newBase)) { // snp doesnt match the previously assigned Variant
                double maxIdentityAssigned = getMaximumIdentifyForBasePatterns(assignedSNPPatterns, base);
                double maxIdentityNew = getMaximumIdentifyForBasePatterns(basePatterns, newBase);
                if (maxIdentityNew > maxIdentityAssigned) {
                    baseList.put(position, newBase);
                }

            }
        }
    }

    private static double getMaximumIdentifyForBasePatterns(LinkedHashMap<TreeSet<Base>, Double> assignedPatterns, Base base) {

        double maximumIdentity = 0.0;
        // check all snps in this pattern
        Iterator<Map.Entry<TreeSet<Base>, Double>> itSNPPatterns = assignedPatterns.entrySet().iterator();
        while (itSNPPatterns.hasNext()) {
            Map.Entry<TreeSet<Base>, Double> entry = itSNPPatterns.next();
            // get the next Variant pattern
            TreeSet<Base> basePattern = entry.getKey();

            // check all snps in this pattern
            Iterator<Base> itBasePattern = basePattern.iterator();
            while (itBasePattern.hasNext()) {
                // get the next Variant
                Base currentBase = itBasePattern.next();

                if (base.equals(currentBase)) {
                    maximumIdentity = Math.max(maximumIdentity, entry.getValue());
                    break;
                }
            } // end while itSNPPattern.hasNext()
        } // end while itSNPPatterns.hasNext()

        return maximumIdentity;

    }

    private static LinkedHashMap<Integer, Character> checkIfAllHSPBasesAreCovered(LinkedHashMap<Integer, Variant> geneHSPList, ArrayList<LinkedHashMap> baseLists) {

        // initialise - set the count as zero of the number of sub-genomes assigned at each HSP position
        LinkedHashMap<Integer, Integer> assignmentCountList = new LinkedHashMap<Integer, Integer>();
        Iterator<Integer> itHSPPosition = geneHSPList.keySet().iterator();
        while (itHSPPosition.hasNext()) {
            int position = itHSPPosition.next();
            assignmentCountList.put(position, 0);
        } // while (it.hasNext())

        LinkedHashMap<Integer, Character> hspBasesCovered = new LinkedHashMap<Integer, Character>();

        // get all the bases assigned to sub-genomes at each position
        LinkedHashMap<Integer, TreeSet<Character>> consensusBaseList = new LinkedHashMap<Integer, TreeSet<Character>>();
        Iterator<LinkedHashMap> itBaseLists = baseLists.iterator();
        while (itBaseLists.hasNext()) {

            LinkedHashMap<Integer, Base> baseList = itBaseLists.next();
            Iterator<Integer> itPosition = baseList.keySet().iterator();
            while (itPosition.hasNext()) {
                int position = itPosition.next();

                TreeSet bases = consensusBaseList.get(position);
                if (bases == null) {
                    bases = new TreeSet();
                    consensusBaseList.put(position, bases);
                }
                bases.add(baseList.get(position).getNucleotide());

                // Also, increment the count of the sub-genome assigned at this position
                assignmentCountList.put(position, assignmentCountList.get(position) + 1);

            } // while (it.hasNext())
        }// while (itBaseLists.hasNext())

        // check if the consensus base from all assigned bases is same as HSP base
        Iterator<Integer> itPosition = geneHSPList.keySet().iterator();
        while (itPosition.hasNext()) {
            int position = itPosition.next();

            // get the HSP at this position
            Variant theHSP = (Variant) geneHSPList.get(position);

            // get the consensus base at this position
            TreeSet bases = consensusBaseList.get(position);
            if (bases == null) { //no base assigned at this position
                hspBasesCovered.put(position, 'N');
            } else if (assignmentCountList.get(position) != DIPLOID_COUNT) { // not all sub-genomes have a base assigned
                hspBasesCovered.put(position, 'N');
            } else {
                char consensusBase = (Character) DNA.getReverseNucleotideMap().get(bases.toString());
                if (theHSP.getNucleotide() != consensusBase) {
                    hspBasesCovered.put(position, 'N');
                } else {
                    hspBasesCovered.put(position, 'Y');
                }
            }

        } // while (it.hasNext())
        return hspBasesCovered;
    }

    private static void tryBaseSetClosure(LinkedHashMap<Integer, Variant> geneHSPList, ArrayList<LinkedHashMap> baseLists, LinkedHashMap<Integer, Base> baseListToUpdate) {

        // get all the bases assigned to sub-genomes at each position
        LinkedHashMap<Integer, TreeSet<Character>> consensusBaseList = new LinkedHashMap<Integer, TreeSet<Character>>();
        Iterator<LinkedHashMap> itBaseLists = baseLists.iterator();
        while (itBaseLists.hasNext()) {

            LinkedHashMap<Integer, Base> baseList = itBaseLists.next();
            Iterator<Integer> itPosition = baseList.keySet().iterator();
            while (itPosition.hasNext()) {
                int position = itPosition.next();

                TreeSet bases = consensusBaseList.get(position);
                if (bases == null) {
                    bases = new TreeSet();
                    consensusBaseList.put(position, bases);
                }
                bases.add(baseList.get(position).getNucleotide());
            } // while (it.hasNext())
        }// while (itBaseLists.hasNext())

        // get all bases from the HSP consensus base 
        Iterator<Integer> itPosition = geneHSPList.keySet().iterator();
        while (itPosition.hasNext()) {
            int position = itPosition.next();

            // get the HSP at this position
            Variant theHSP = (Variant) geneHSPList.get(position);
            // get all bases for this HSP 
            TreeSet bases = new TreeSet((TreeSet) DNA.getNucleotideMap().get(theHSP.getNucleotide()));
            // remove the bases assigned to the sub-genomes    
            TreeSet assignedBases = consensusBaseList.get(theHSP.getPosition());
            if (assignedBases != null) {
                bases.removeAll(assignedBases);
            }

            // Check if one base is remaining
            if (bases.size() == 1) {
                Base base = (Base) baseListToUpdate.get(position);
                if (base == null) {
                    base = new Base(position, (Character) bases.first());
                    baseListToUpdate.put(position, base);
                }
            }

        } // while (it.hasNext())

    }

    private static void tryBaseSetClosure(LinkedHashMap<Integer, Variant> geneHSPList, ArrayList<LinkedHashMap> baseLists) {

        // initialise - set the count as zero of the number of sub-genomes assigned at each HSP position
        LinkedHashMap<Integer, Integer> assignmentCountList = new LinkedHashMap<Integer, Integer>();
        Iterator<Integer> itHSPPosition = geneHSPList.keySet().iterator();
        while (itHSPPosition.hasNext()) {
            int position = itHSPPosition.next();
            assignmentCountList.put(position, 0);
        } // while (it.hasNext())

        // get all the bases assigned to sub-genomes at each position
        LinkedHashMap<Integer, TreeSet<Character>> consensusBaseList = new LinkedHashMap<Integer, TreeSet<Character>>();
        Iterator<LinkedHashMap> itBaseLists = baseLists.iterator();
        while (itBaseLists.hasNext()) {

            LinkedHashMap<Integer, Base> baseList = itBaseLists.next();
            Iterator<Integer> itPosition = baseList.keySet().iterator();
            while (itPosition.hasNext()) {
                int position = itPosition.next();

                // increment the count of the sub-genome assigned at this position
                assignmentCountList.put(position, assignmentCountList.get(position) + 1);

                TreeSet bases = consensusBaseList.get(position);
                if (bases == null) {
                    bases = new TreeSet();
                    consensusBaseList.put(position, bases);
                }
                bases.add(baseList.get(position).getNucleotide());
            } // while (it.hasNext())
        }// while (itBaseLists.hasNext())

        // get all bases from the HSP consensus base 
        Iterator<Integer> itPosition = geneHSPList.keySet().iterator();
        while (itPosition.hasNext()) {
            int position = itPosition.next();

            // get the HSP at this position
            Variant theHSP = (Variant) geneHSPList.get(position);
            // get all bases for this HSP 
            TreeSet bases = new TreeSet((TreeSet) DNA.getNucleotideMap().get(theHSP.getNucleotide()));
            // remove the bases assigned to the sub-genomes    
            TreeSet assignedBases = consensusBaseList.get(theHSP.getPosition());
            if (assignedBases != null) {
                bases.removeAll(assignedBases);
            }
            // Check if one base is remaining and only one sub-genome is left unassigned
            if (bases.size() == 1 && DIPLOID_COUNT - assignmentCountList.get(position) == 1) {
                //iterate on all base lists
                itBaseLists = baseLists.iterator();
                while (itBaseLists.hasNext()) {
                    // get the next base list
                    LinkedHashMap<Integer, Base> baseList = itBaseLists.next();
                    // get the snp assigned at this position
                    Base base = baseList.get(position);

                    if (base == null) { // no assingment has been made, this must be the missing subgenome
                        // assign the remaining base to this position on this genome
                        base = new Base(position, (Character) bases.first());
                        baseList.put(position, base);
                        break;
                    }
                }
            }
        } // while (it.hasNext())

    }

    private static void writeBases(LinkedHashMap<Integer, Base> baseList, LinkedHashMap<Integer, Variant> geneHSPList, LinkedHashMap<Integer, Variant> geneDiploidSNPList,
            String sequenceName, BufferedWriter out, LinkedHashMap<Integer, Character> hspBasesCovered, boolean isMissingDiploid, boolean outputVCF) {

        try {
            Iterator<Integer> itPosition = geneHSPList.keySet().iterator();
            while (itPosition.hasNext()) {
                int position = itPosition.next();
                Variant theHSP = geneHSPList.get(position);

                Base base = (Base) baseList.get(position);
                if (base == null) { // snp was not assigned for this genome, create a dummy snp
                    base = new Base(position, DNA.BASE_N);
                }
                Variant theDiploidSNP;
                if (geneDiploidSNPList != null) {
                    theDiploidSNP = geneDiploidSNPList.get(position);

                    if (theDiploidSNP == null) {
//                        if (isMissingDiploid) {
                        theDiploidSNP = new Variant(position, DNA.BASE_MISSING, theHSP.getReferenceBase());
//                        } else {
//                            theDiploidSNP = new Variant(position, theHSP.getReferenceBase(), theHSP.getReferenceBase());
//                        }
                    } else if (theDiploidSNP.getNucleotide() == DNA.BASE_BLANK) {
                        theDiploidSNP.setNucleotide(theHSP.getReferenceBase());
                    }
                } else {
                    theDiploidSNP = new Variant(position, DNA.BASE_AMBIGUOUS, theHSP.getReferenceBase());
                }

                Character posHspBasesCovered = hspBasesCovered.get(position);
                if (outputVCF) {
                    out.write(sequenceName + "\t" + position + "\t.\t" + theHSP.getReferenceBase() + "\t" + (theHSP.getReferenceBase() == base.getNucleotide() ? "." : base.getNucleotide()) + "\t" + "100" + "\t.\t" + "HSP=" + geneHSPList.get(position).getNucleotide() + ";FC=" + posHspBasesCovered.toString() + ";DB=" + theDiploidSNP.getNucleotide());
                } else {
                    out.write(sequenceName + "\t" + position + "\t" + theHSP.getReferenceBase() + "\t" + theDiploidSNP.getNucleotide() + "\t" + base.getNucleotide() + "\t" + hspBasesCovered.get(position));
                }
                out.newLine();

            } // while (it.hasNext())
        } catch (IOException ex) {
            // do noting;
        }
    }

    // UNUSED
    private static double calculateBaseScore(Variant referenceBase, Variant diploidSNP, Character hspBasesCovered) {
        double score = 0.0;

        if (DNA.isNucleotide(referenceBase.getNucleotide())) {
            score += 1.0;
        }
        if (DNA.isNucleotide(diploidSNP.getNucleotide())) {
            score += 1.0;
        }
        if (hspBasesCovered == 'Y') {
            score += 1.0;
        }

        return (double) Math.round(-10 * Math.log10(1 - score / 3.0));
    }

    private static void writeVCFHeader(BufferedWriter out, String[] args, int diploidNo, SAMSequenceDictionary sequenceDictionary) throws IOException {
        out.write("##fileformat=VCFv4.1");
        out.newLine();
        out.write("##source=HANDS2" + HANDS2_VERSION);
        out.newLine();
        out.write("##Genome=" + Integer.toString(diploidNo));
        out.newLine();
        StringBuilder commandLine = new StringBuilder();
        for (String s : args) {
            commandLine.append(s);
            commandLine.append(" ");
        }
        out.write("##HANDS2CommandLine<ID=assign,CommandLine=\"" + commandLine.substring(0, commandLine.length() - 1) + "\">");
        out.newLine();
        out.write("##ALT=<ID=X,Description=\"Homoeallelic base identity assigned at this position. (\".\" for same as the reference, N for unassigned)\">");
        out.newLine();
        out.write("##QUAL=<ID=X,Description=\"Quality of the assignment. Currently unmeaningful. Always 100\">");
        out.newLine();
        out.write("##INFO=<ID=HSP,Number=1,Type=String,Description=\"HSP Base present at this position in IUPAC code\">");
        out.newLine();
        out.write("##INFO=<ID=DB,Number=1,Type=String,Description=\"Diploid Base (< for low coverage, 0 for no coverage, * for heterozygous base, ? for unknown base)\">");
        out.newLine();
        out.write("##INFO=<ID=FC,Number=1,Type=String,Description=\"Fully characterized position. Y if all HSP bases were assigned to subgenomes, N otherwise\">");
        out.newLine();
        // write sequence information
        Iterator<SAMSequenceRecord> sequenceIterator = sequenceDictionary.getSequences().iterator();
        while (sequenceIterator.hasNext()) {
            // get the next sequence
            SAMSequenceRecord sequenceEntry = sequenceIterator.next();
            // Initialise a new variant list
            out.write("##contig=<ID=" + sequenceEntry.getSequenceName() + ",length=" + Integer.toString(sequenceEntry.getSequenceLength()) + ">");
            out.newLine();
        } // end while        
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        out.newLine();

    }

    private static void printBasePatterns(LinkedHashSet<TreeSet<Base>> basePatterns) {
        Iterator<TreeSet<Base>> itBasePatterns = basePatterns.iterator();
        while (itBasePatterns.hasNext()) {

            // get the next Variant pattern
            TreeSet basePattern = itBasePatterns.next();
            System.out.println(basePattern);
        }
    }

    private static void printBasePatterns(LinkedHashMap<TreeSet<Base>, Integer> basePatterns) {
        Iterator<Map.Entry<TreeSet<Base>, Integer>> itBasePatterns = basePatterns.entrySet().iterator();
        while (itBasePatterns.hasNext()) {

            // get the next Variant pattern
            Map.Entry entry = itBasePatterns.next();
            System.out.println(entry.getKey() + ": " + entry.getValue());
        }
    }

    private static void multiFastaToFasta(String[] args) {

        String multiFastaFile = "";
        String outFastaFile = "";
        String fastaHeader = "";
        int gapSize = 200;

        // read input arguments
        int idx = 1;
        int nArgs = args.length;
        while (idx < nArgs) {
            String parameter = args[idx++];
            if (parameter.equals("-h") || parameter.equals("-help")) { // Help
                if (nArgs == 2) {
                    printHelpSeq2Ref();
                } else {
                    System.err.println("Invalid option used with Help!");
                }
                return;
            } else if (parameter.equals("-i")) { // multifasta file
                multiFastaFile = args[idx++];
            } else if (parameter.equals("-o")) { // output file
                outFastaFile = args[idx++];
            } else if (parameter.equals("-n")) { // sequence name
                fastaHeader = args[idx++];
            } else if (parameter.equals("-g")) { // Gap Size
                gapSize = Integer.valueOf(args[idx++]);

            } else {
                System.err.println("Invalid option: " + args[idx - 1]);
                printHelpSeq2Ref();
                return;
            }
        }
        // validate the input 
        if (multiFastaFile.isEmpty()) {
            System.err.println("Multifasta file must be specified.");
            return;
        }
        if (fastaHeader.isEmpty()) {
            fastaHeader = "sequence";
        }
        if (outFastaFile.isEmpty()) {
            outFastaFile = multiFastaFile + ".fa";
        }

        DNA.multiFastaToFasta(new File(multiFastaFile), new File(outFastaFile), fastaHeader, gapSize);
    }

    private static void calculateBaseCoverages(String[] args) {

        String samFile = "";
        String outCoverageFile = "";
        int baseQualityThreshold = 20;

        // read input arguments
        int idx = 1;
        int nArgs = args.length;
        while (idx < nArgs) {
            String parameter = args[idx++];
            if (parameter.equals("-h") || parameter.equals("-help")) { // Help
                if (nArgs == 2) {
                    printHelpCoverage();
                } else {
                    System.err.println("Invalid option used with Help!");
                }
                return;
            } else if (parameter.equals("-i")) { // sam file
                samFile = args[idx++];
            } else if (parameter.equals("-o")) { // output file
                outCoverageFile = args[idx++];
            } else if (parameter.equals("-q")) { // Base quality threshold
                baseQualityThreshold = Integer.valueOf(args[idx++]);

            } else {
                System.err.println("Invalid option: " + args[idx - 1]);
                printHelpCoverage();
                return;

            }
        }

        // validate the input 
        if (samFile.isEmpty()) {
            System.err.println("SAM file must be specified.");
            return;
        }
        if (outCoverageFile.isEmpty()) {
            outCoverageFile = samFile + ".bc";
        }

        BaseCoverage.calculateBaseCoverages(samFile, outCoverageFile, baseQualityThreshold);
    }

    private static void printHelp() {

        System.out.println("HANDS2 " + HANDS2_VERSION);
        System.out.println("Assign homoeallelic base identities in allopolyploids using diploid similarity.");
        System.out.println();
        System.out.println("Usage: java -jar hands2.jar <command> <input parameters>");
        System.out.println();
        System.out.println("Commands");
        System.out.println("\thelp      :   Display this help");
        System.out.println("\tassign    :   Assign homoeallelic base identities");
        System.out.println("\tcoverage  :   Calculate the number of reads supporting a particular base at each position");
        System.out.println("\tseq2ref   :   Create an in silico reference using a set of unigenes, contigs or other sequences");
    }

    private static void printHelpAssign() {

        System.out.println("HANDS2 " + HANDS2_VERSION);
        System.out.println();
        System.out.println("Command: assign - Assign homoeallelic base identities");
        System.out.println("Usage: java -jar hands2.jar assign <input parameters>");
        System.out.println();
        System.out.println("Input Parameters");
        System.out.println("\t-h or -help    :   Display this help");
        System.out.println("\t-i <str>       :   Polyploid SAM/BAM file");
        System.out.println("\t-g <str>       :   GFF3 file containing gene start/end coordinates");
        System.out.println("\t-hsp <str>     :   Polyploid HSP file in VCF Format.");
        System.out.println("\t-snp<n> <str>  :   Diploid # n SNP file in VCF Format.");
        System.out.println("\t-bc <str>      :   Polyploid Base coverage file (optional). See coverage command.");
        System.out.println("\t-bc<n> <str>   :   Diploid # n Base coverage file, e.g. bc1 (optional). See coverage command.");
        System.out.println("\t-out<n> <str>  :   Sub-Genome # n output file, e.g. out1");
        System.out.println("\t-vcf <boolean> :   Generate VCF output (Default: TRUE). When FALSE, tab-delimited output is generated.");
        System.out.println("\t-sp <double>   :   SNP pair proportion threshold (Default: 0.05)");
        System.out.println("\t-pm <double>   :   Base pattern matching threshold (Default: 0.5)");
        System.out.println("\t-pa <char>     :   Base pattern assignment mode (M: Keep maximum proportion for a base or A: Add all proportions; Default: M)");
        System.out.println("\t-r <boolean>   :   Rectify Assignment using reference genome (Default: FALSE)");
        System.out.println("\t-m <boolean>   :   Merge Base Patterns before assignment (Default: FALSE)");
        System.out.println("\t-u <boolean>   :   Assign the unassigned base, if any, to the subgenome to which no base is assigned (Default: TRUE)");
        System.out.println("\t-d <int>       :   Use genome <int> as distant genome (Default: <null>)");
        System.out.println("Note: At most one diploid SNP file can be missing. Use \"\" for the missing file.");
        System.out.println("      HANDS2 supports up to 10 genomes.");
    }

    private static void printHelpSeq2Ref() {

        System.out.println("HANDS2 " + HANDS2_VERSION);
        System.out.println();
        System.out.println("Command: seq2ref - Create an in silico reference from given sequences/contigs");
        System.out.println("Usage: java -jar hands2.jar seq2ref <input parameters>");
        System.out.println();
        System.out.println("Input Parameters");
        System.out.println("\t-h or -help  :   Display this help");
        System.out.println("\t-i <str>     :   Input sequence file (multifasta format)");
        System.out.println("\t-o <str>     :   Output file");
        System.out.println("\t-n <str>     :   Header for the in silico reference");
        System.out.println("\t-g <int>     :   Gap size between two sequences (Default: 200)");
    }

    private static void printHelpCoverage() {

        System.out.println("HANDS2 " + HANDS2_VERSION);
        System.out.println();
        System.out.println("Command: coverage - Calculate base coverage for each position from a SAM file");
        System.out.println("Usage: java -jar hands2.jar coverage <input parameters>");
        System.out.println();
        System.out.println("Input Parameters");
        System.out.println("\t-h or -help  :   Display this help");
        System.out.println("\t-i <str>     :   SAM/BAM file");
        System.out.println("\t-o <str>     :   Output file");
        System.out.println("\t-q <int>     :   Base quality threshold (Default: 20)");
    }
}
