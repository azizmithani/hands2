# hands2
HANDS2: accurate assignment of homoeallelic base-identity in allopolyploids despite missing data

Current version: HANDS2 v1.1.1
Assign homoeallelic base identities in allopolyploids using diploid similarity.

Usage: java -jar hands2.jar <command> <input parameters>

Commands
	help      :   Display this help
	assign    :   Assign homoeallelic base identities
	coverage  :   Calculate the number of reads supporting a particular base at each position
	seq2ref   :   Create an in silico reference using a set of unigenes, contigs or other sequences


Command: assign - Assign homoeallelic base identities
Usage: java -jar hands2.jar assign <input parameters>

Input Parameters
	-h or -help    :   Display this help
	-i <str>       :   Polyploid SAM/BAM file
	-g <str>       :   GFF3 file containing gene start/end coordinates
	-hsp <str>     :   Polyploid HSP file in VCF Format.
	-snp<n> <str>  :   Diploid # n SNP file in VCF Format.
	-bc <str>      :   Polyploid Base coverage file (optional). See coverage command.
	-bc<n> <str>   :   Diploid # n Base coverage file, e.g. bc1 (optional). See coverage command.
	-out<n> <str>  :   Sub-Genome # n output file, e.g. out1
	-vcf <boolean> :   Generate VCF output (Default: TRUE). When FALSE, tab-delimited output is generated.
	-sp <double>   :   SNP pair proportion threshold (Default: 0.05)
	-pm <double>   :   Base pattern matching threshold (Default: 0.5)
	-pa <char>     :   Base pattern assignment mode (M: Keep maximum proportion for a base or A: Add all proportions; Default: M)
	-r <boolean>   :   Rectify Assignment using reference genome (Default: FALSE)
	-m <boolean>   :   Merge Base Patterns before assignment (Default: FALSE)
	-u <boolean>   :   Assign the unassigned base, if any, to the subgenome to which no base is assigned (Default: TRUE)
	-d <int>       :   Use genome <int> as distant genome (Default: <null>)
Note: At most one diploid SNP file can be missing. Use "" for the missing file.
      HANDS2 supports up to 10 genomes.

  
  Command: coverage - Calculate base coverage for each position from a SAM file
Usage: java -jar hands2.jar coverage <input parameters>

Input Parameters
	-h or -help  :   Display this help
	-i <str>     :   SAM/BAM file
	-o <str>     :   Output file
	-q <int>     :   Base quality threshold (Default: 20)

  
  Command: seq2ref - Create an in silico reference from given sequences/contigs
Usage: java -jar hands2.jar seq2ref <input parameters>

Input Parameters
	-h or -help  :   Display this help
	-i <str>     :   Input sequence file (multifasta format)
	-o <str>     :   Output file
	-n <str>     :   Header for the in silico reference
	-g <int>     :   Gap size between two sequences (Default: 200)
