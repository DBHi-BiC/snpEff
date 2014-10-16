snpEff with HGVS
=================

This codebase is branched off of official SnpEff version 3.6

### Setup ###

The snpEff make script assumes a specific Eclipsey directory structure. Let's not muck with that.

```
mkdir ~/workspace ~/snpEff
cd ~/workspace
git clone git@github.com:CBMi-BiG/snpEff.git SnpEff
cd SnpEff
```

Maven/Ant/Ivy have a hard time getting Picard and Samtools jars, but they are in the Picard distribution. Alter the pom.xml versions as necessary.

```
wget http://downloads.sourceforge.net/project/picard/picard-tools/1.84/picard-tools-1.84.zip
unzip picard-tools-1.84.zip
mvn install:install-file -DgroupId=net.sf.samtools -DartifactId=Sam -Dversion=1.84 \
-Dpackaging=jar -Dfile=$PWD/picard-tools-1.84/sam-1.84.jar
mvn install:install-file -DgroupId=net.sf.picard -DartifactId=Picard -Dversion=1.84 \
-Dpackaging=jar -Dfile=$PWD/picard-tools-1.84/picard-1.84.jar
```

CHANGELOG
=========
* 3.1h2 Introduces proper HGVS for insertions and deletions, including those which require "walking and rolling" to identify the correct indel frame.
* 3.6 Merges official SnpEff 3.6

To build snpEff-3.6-jar-with-dependencies.jar in the ~/workspace/SnpEff/target directory:


```
./scripts/make.sh
```

This is copied to ~/workplace/SnpEff/snpEff.jar per the config file

To get GRCh37 annotations:
```
mkdir -p data/GRCh37.64
java -jar snpEff.jar download GRCh37.64
```

To test (assuming you have downloaded all the test cases):

```
java -cp snpEff.jar \
ca.mcgill.mcb.pcingola.snpEffect.testCases.TestSuiteAll
```


```
java -Xmx4g -jar snpEff.jar hg19 -i vcf -o vcf tests/hgvs_test_in.vcf
```
### Output ###
This fork of SnpEff introduces the following annotations to VCF files
*   `Segment` (previously Exon) - a verbose description of  exons (e.g. `NM_152486.2.ex.3`) or introns (`NR_024540.1_intron_7`) for all applicable transcript hits
*   `HGVS_DNA_nomenclature` HGVS coding DNA nomenclature for exonic SNPs e.g. `NM_001005484.1:c.655C>T`. See http://www.hgvs.org/mutnomen/recs-DNA.html for more info.

#### VCF output example ####
```
##SnpEffVersion="SnpEff_cbmi 3.6h"
##SnpEffCmd="SnpEff  hg19 -i vcf -o vcf tests/hgvs_test_in.vcf "
##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Rank | Segment | HGVS_DNA_nomenclature | HGVS_protein_nomenclature [ | ERRORS | WARNINGS ])' ">
#CHROME POS     ID      REF     ALT     QUAL    FILTER  INFO
1	16856	.	A	G	100.0	PASS	DP=100;EFF=DOWNSTREAM(MODIFIER|||||DDX11L1||NON_CODING|NR_046018.2||||),INTRON(MODIFIER|||||WASH7P||NON_CODING|NR_024540.1|7|NR_024540.1_intron_7||),SPLICE_SITE_DONOR(HIGH|||||WASH7P||NON_CODING|NR_024540.1|7|NR_024540.1.ex.5||)
1	69745	.	C	T	100.0	PASS	DP=100;EFF=STOP_GAINED(HIGH|NONSENSE|Cag/Tag|Q219*|305|OR4F5||CODING|NM_001005484.1|1|NM_001005484.1.ex.1|NM_001005484.1:c.655C>T|p.Q219*)
1	865697	.	A	G	100.0	PASS	DP=100;EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Atc/Gtc|I79V|681|SAMD11||CODING|NM_152486.2|3|NM_152486.2.ex.3|NM_152486.2:c.235A>G|p.I79V)
1	1114650	.	C	T	100.0	PASS	DP=100;EFF=STOP_GAINED(HIGH|NONSENSE|Cga/Tga|R19*|673|TTLL10||CODING|NM_001130045.1|4|NM_001130045.1.ex.4|NM_001130045.1:c.55C>T|p.R19*),UPSTREAM(MODIFIER||||404|TTLL10||CODING|NM_153254.2||||)
```
