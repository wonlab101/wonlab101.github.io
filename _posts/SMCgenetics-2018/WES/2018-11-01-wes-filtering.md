---
layout: page
permalink: /SMCgenetics-2018/WES/wes-filtering
---

## **2. Filtering**

This course assume that the input data is "Raw SNPs + Indels VCF" (blue box), the outcome of joint genotyping of GATK ([Genome Anallysis Tool Kit](https://software.broadinstitute.org/gatk/){: target="blank" }), and the final output data would be "Analysis-ready VCF" (red box) as shown in [GATK Best Practice](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145){: target="blank" }.  
After Genotyping, we need to 'filtering process' to make clean variant call set.  

We are going to practice filtering process using ```GATK-VariantRecalibrator```, ```GATK-ApplyRecalibration```, ```GATK-SelectVariants```, ```GATK-VariantFiltration```, ```GATK-VariantAnnotator``` and ```VCFtools``` in certain situation.  

The **input data** is a toy VCF file for _five different samples_ with _biallelic_ variants of only _chromosome 22_.

![gatk_bp](https://user-images.githubusercontent.com/26876362/48392916-bb441000-e750-11e8-8861-5139ba595892.png){: width="800" height="400"}  

credit: [GATK Best Practice](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145){: target="blank" }


#### Link for each tool/module:  
- [GATK-VariantRecalirator](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php){: target="blank" }  
- [GATK-ApplyRecalibration](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php){: target="blank" }   
- [GATK-SelectVariants](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php){: target="blank" }  
- [GATK-VariantFiltration](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php){: target="blank" }  
- [GATK-VariantAnnotator](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator.php){: target="blank" }  
- [VCFtools](https://vcftools.github.io/man_latest.html){: target="blank" }  

- GATK version: 3.7  
- VCFtools version: 0.1.15

---------------------------
<br>

### **<u>1) GATK-VQSR (Variant Quality Score Recalibration)</u>**

As described earlier, our input file is the outcome of joint genotyping. We start the filtering process with VQSR at first.   
VQSR makes a new variant quality score for each variant based on the information of known SNPs and indels such as 1000 genome project. Detail description is in [VQSR](https://gatkforums.broadinstitute.org/gatk/discussion/2805/(https://www.broadinstitute.org/gatk/guide/article?id=1259)){: target="blank" }


#### **1-1) GATK-VariantRecalibrator**

Create a Gaussian mixture model by looking at the annotations values over a high quality subset of the input call set and then evaluate all input variants. This step produces a recalibration file. [read](https://software.broadinstitute.org/gatk/documentation/article.php?id=39){: target="blank" }

- Input: BGDA.jg.chr22.bi.vcf  
- Output:   
1) SNP: BGDA.snp.recal, BGDA.snp.tranches  
2) Indel: BGDA.indel.split.recal, BGDA.indel.split.tranches 
- Script: 01-1_VQSR.sh / **15min (homework)**   

(a) SNPs

```bash
$ java -jar GenomeAnalysisTK.jar  -T VariantRecalibrator  -R ucsc.hg19.fasta  \
 -input BGDA.jg.chr22.vcf  -recalFile BGDA.snp.recal \
 -tranchesFile BGDA.snp.tranches -allPoly \
 -tranche 100.0  -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ  -an InbreedingCoeff \
 -resource:hapmap,known=false,training=true,truth=true,prior=15 $HapMap \
 -resource:omni,known=false,training=true,truth=true,prior=12 $KG_Omni  \
 -resource:1000G,known=false,training=true,truth=false,prior=10 $KG_phase1 \
 -resource:dbsnp138,known=false,training=false,truth=false,prior=7 $dbSNP138 \
 -resource:dbsnp129,known=true,training=false,truth=false,prior=3  $dbSNP138excluding \
 --maxGaussians 4 -mode SNP 
```

(b) Indels

```bash
java -jar GenomeAnalysisTK.jar -T VariantRecalibrator -R ucsc.hg19.fasta \
 -input BGDA.jg.chr22.vcf  -recalFile BGDA.indel.split.recal \
 -tranchesFile  BGDA.indel.split.tranches  -allPoly \
 -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
 -an FS -an ReadPosRankSum -an InbreedingCoeff -an MQRankSum -an QD     \
 -resource:mills,known=false,training=true,truth=true,prior=12 $Mills \
 -resource:axiomPoly,known=false,training=true,truth=false,prior=10 $Axiom_poly \
 -resource:dbsnp138,known=true,training=false,truth=false,prior=2 $dbSNP138  \
 --maxGaussians 4 -mode INDEL 
```
---
<br>

#### **1-2) GATK-ApplyRecalibration**

Apply the model parameters to each variant in input VCF files producing a recalibrated VCF file in which each variant is annotated with its VQSLOD value. In addition, this step will filter the calls based on this new lod score by adding lines to the FILTER column for variants that don't meet the specified lod threshold.[read](https://software.broadinstitute.org/gatk/documentation/article.php?id=39){: target="blank" }


**<u>#NOTE</u>** For convenience, we use **only SNP file** in the following steps. 

  - Input:     
    - BGDA.jg.chr22.bi.vcf  
    - BGDA.snp.recal  
    - BGDA.snp.tranches
  - Output:  
    - BGDA.jg.chr22.bi.VQSR.snp.vcf  
  - Script: 01-2_ApplyVQSR_snp.sh /**1 min**

```bash
java -jar GenomeAnalysisTK.jar -T ApplyRecalibration -R /home/user/Reference/ucsc.hg19.fasta \
 -o BGDA.jg.chr22.bi.VQSR.snp.vcf -input BGDA.jg.chr22.bi.vcf  \
 -recalFile BGDA.snp.recal  -tranchesFile BGDA.snp.tranches \
 -ts_filter_level 99.8  -mode SNP
```

- Let's open the input vcf and output vcf file, then see FILTER column.  

```bash
vi BGDA.jg.chr22.bi.vcf # input data
```  
![before VQSR](https://user-images.githubusercontent.com/26876362/48540741-48789780-e8fe-11e8-8868-e3a424efa910.png)

```bash
vi BGDA.jg.chr22.bi.VQSR.snp.vcf # output data
```  

![after VQSR](https://user-images.githubusercontent.com/26876362/48548560-9e563b00-e910-11e8-8e05-57ac76916dd1.png)
- chr22:16157701:T:C --> VQSRTrancheSNP99.90to100.00



**<u>#NOTE</u> VQSR is not suitable for <30 samples.** [read](https://software.broadinstitute.org/gatk/documentation/article.php?id=3225) 

[VQSR tutorial for step-by-step instructions](https://software.broadinstitute.org/gatk/documentation/article.php?id=2805){: target="blank" }

---
<br>

### **<u>2) GATK-SelectVariants</u>**

Often, a VCF containing many samples and/or variants will need to be subset in order to facilitate certain analyses (e.g. comparing and contrasting cases vs. controls; extracting variant or non-variant loci that meet certain requirements, displaying just a few samples in a browser like IGV, etc.). SelectVariants can be used for this purpose.[read](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php){: target="blank" }

**<u>#NOTE</u>** GATK-SelectVariants applies to **INFO** field of VCF. 

- Basic synopsis of GATK-SelectVariants:  
```bash
java -jar GenomeAnalysisTK.jar \
 -T SelectVariants \
 -R reference.fasta \
 -V input.vcf \
 -o output.vcf \
 options
```  
- Input: VCF
- Output: A new VCF file (filtered VCF)

<br>

#### **#Example 02-1) Select variants of VQSR PASS** 

- Input: BGDA.jg.chr22.VQSR.snp.vcf
- Output: BGDA.jg.chr22.VPASS.snp.vcf
- Script: 02-1_VQSRPASS.sh 

```bash
java -jar GenomeAnalysisTK.jar  -T SelectVariants  -R /home/user/Reference/ucsc.hg19.fasta \
  -V BGDA.jg.chr22.bi.VQSR.snp.vcf  -o BGDA.VPASS.snp.vcf \
  -ef 
```  
  - ```-ef```: same as ```--excluedFiltered``` Don't include filtered sites  


- Let's open the input vcf and output vcf file and confirm the result.

```bash
vi BGDA.VPASS.snp.vcf # output data
```  

![vcf with VPASS](https://user-images.githubusercontent.com/26876362/48548953-c003f200-e911-11e8-9a33-5698e038d3eb.png)  

- chr22:16157701:T:C --> VQSRTrancheSNP99.90to100.00 --> **removed**  

_Q. How do I know the number of filtered variants?_  
```bash
wc -l BGDA.jg.chr22.bi.VQSR.snp.vcf BGDA.VPASS.snp.vcf 
```

![count VPASS filtered vcf](https://user-images.githubusercontent.com/26876362/48549038-06595100-e912-11e8-9f87-456e04acfaad.png)


<br>

#### **#Example 02-2) Select variants of Quality by Depth >= 3.0**

- Input: BGDA.VPASS.snp.vcf
- Output: BGDA.VPASS.snp.QD.vcf
- Script: 02-2_QD.sh

```bash
java -jar GenomeAnalysisTK.jar  -T SelectVariants  -R /home/user/Reference/ucsc.hg19.fasta \
  -V BGDA.VPASS.snp.vcf  -o BGDA.VPASS.snp.QD.vcf \
  -select "QD >= 3.0"
```

- Let's open the output file and confirm the result.

```bash
vi BGDA.VPASS.snp.QD.vcf # output data
```

![02-2 output](https://user-images.githubusercontent.com/26876362/48541887-0bfa6b00-e901-11e8-9da2-67b36563d68c.png)

- Count lines.

```bash
wc -l BGDA.VPASS.snp.vcf  BGDA.VPASS.snp.QD.vcf
```

![count QD filtered vcf](https://user-images.githubusercontent.com/26876362/48549376-f68e3c80-e912-11e8-8b00-21e81fce4e0b.png)


---
<br>

### **<u>3) GATK-VariantFiltration</u>**

This tool is designed for **hard-filtering variant calls** based on certain criteria. Records are hard-filtered by changing the value **in the FILTER field** to something other than PASS. Filtered records will be preserved in the output unless their removal is requested in the command line.[read](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php){: target="blank" }


- Basic synopsis of GATK-VariantFiltration:  

```bash
 java -jar GenomeAnalysisTK.jar \
   -T VariantFiltration \
   -R reference.fasta \
   -o output.vcf \
   --variant input.vcf \
   --filterExpression "AB < 0.2 || MQ0 > 50" \
   --filterName "SomeFilterName" 
```
<br>

#### **#Example 03) Make filter: AC = 0 name as "ACzero"**

- Input: BGDA.VPASS.snp.QD.vcf
- Output: BGDA.VPASS.snp.QD.AC.vcf
- Script: 03_VariantFiltration_ACzero.sh

```bash
java -jar GenomeAnalysisTK.jar  -T VariantFiltration -R /home/user/Reference/ucsc.hg19.fasta \
 -V BGDA.VPASS.snp.QD.vcf  -o BGDA.VPASS.snp.QD.AC.vcf \
 --filterExpression "AC==0" --filterName "ACzero"
```

- Let's open the output file and confirm the result.  

```bash
vi BGDA.VPASS.snp.QD.AC.vcf # output data
```

![after ACzero](https://user-images.githubusercontent.com/26876362/48549543-6c92a380-e913-11e8-9514-8b41276b5ace.png)

_Q. How do I make filtered variant by "ACzero"?_  
A. Same as '#Example 02-1'. Use '```SelectVariants```' and '```-ef```'.
- Output: BGDA.VPASS.snp.QD.ACfilterd.vcf

- Count lines.

```bash
wc -l BGDA.VPASS.snp.QD.AC.vcf  BGDA.VPASS.snp.QD.ACfiltered.vcf
```

![count lines](https://user-images.githubusercontent.com/26876362/48556933-ad48e780-e928-11e8-93da-a7d45247dd1d.png)


---
<br>

### **<u>4) VCFtools</u>**

vcftools is a suite of functions for use on genetic variation data in the form of VCF and BCF files. The tools provided will be used mainly to summarize data, run calculations on data, filter out data, and convert data into other useful file formats. [vcftools](https://vcftools.github.io/man_latest.html#DESCRIPTION){: target="blank" }

- Basic synopsis:  

```perl
vcftools [ --vcf FILE | --gzvcf FILE | --bcf FILE] \
 [ --out OUTPUT PREFIX ] \ 
 [ FILTERING OPTIONS] \ 
 [ OUTPUT OPTIONS ]
```

#### **#Example 04-1) Calculate Hardy-Weinberg p-value for every site in the vcf file**

- Input: BGDA.VPASS.snp.QD.ACfiltered.vcf
- Output: BGDA.VPASS.snp.QD.ACfiltered.hwe
- Script: 04-1_vcftools_HWE_p.sh  
```bash
vcftools --vcf BGDA.VPASS.snp.QD.ACfiltered.vcf \
 --out BGDA.VPASS.snp.QD.ACfiltered  --hardy    #  --out prefix of output
```
- Description:  
```--hardy``` reports a p-value for each site from a HWE test using an exact test, as defined by Wigginton, Cutler and Abecasis (2005). Sites with a p-value below the threshold defined by this option are taken to be out of HWE, and therefore excluded.    
```suffix.hwe``` output file

- Let's open the output file and confirm the result.

```bash
vi BGDA.VPASS.snp.QD.ACfiltered.hwe
```

![BGDA.VPASS.snp.QD.ACfiltered.hwe](https://user-images.githubusercontent.com/26876362/48557269-a40c4a80-e929-11e8-8128-aa506365c957.png)


<br>

#### **#Example 04-2) Filter variants HWE < 1e-06**

- Input: BGDA.VPASS.snp.QD.ACfiltered.vcf
- Output: BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf
- Script: 04-2_vcftools_applyHWE.sh  
```bash
vcftools --vcf BGDA.VPASS.snp.QD.ACfiltered.vcf  \
 --out BGDA.VPASS.snp.QD.AC.hwefiltered \
 --hwe 0.000001 \
 --recode --recode-INFO-all
```

- Count lines.

```bash
wc -l BGDA.VPASS.snp.QD.ACfiltered.vcf  BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf
```

![image](https://user-images.githubusercontent.com/26876362/48557565-7b388500-e92a-11e8-920f-032b78c7004a.png)
- No variant was filtered out by hwe test.

---
<br>

Related articles:  
[Apply hard filters to a call set](https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set/p1){: target="blank" } 

---
<br>

Go to other sections,  
[1. Interpretation of VCF]({% post_url 2018-11-01-wes-vcf%})  
[3. Annotation]({% post_url 2018-11-01-wes-annotation%})  
[4. Homework]({% post_url 2018-11-01-wes-practice%})  

