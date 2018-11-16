---
layout: page
permalink: /WES/wes-practice
---

## **4. Homework: ExAC Project**

Finally, with [ExAC project/2016 Nature](https://www.nature.com/articles/nature19057){: target="blank" }, we would like to practice filtering/annotation steps by using same threshold described in the paper. Detail process of data generation is described in supplementary information and a piece of contents is in below. We plan to follow the highlighted sections for practice: **1.5 Assessment of variant filtering** and **1.10 Functional annotation**.

![Nature 2016, ExAC](https://user-images.githubusercontent.com/26876362/48522638-bb611e80-e8bc-11e8-8210-8c1371aa0b9e.png){: width="750" height="400"}

<br>

![Table of Contents in Supplementary information](https://user-images.githubusercontent.com/26876362/48399323-ac685800-e766-11e8-8484-77e4f431f24d.png){: width="500" height="350"}

_**#Comment:** There may be many ways to set the condition of the paper described. So, my solution is just one of the ways!_

---
<br>

### **<u>1) Assessment of variant filtering</u>**

#### **#Paper description**  
- **Filter out:**   
	- Indels and multiallelic SNPs (biallelic SNP only) # Analysis of multi allele was performed in ExAC project, but we do not deal with that in this section.    
	- InbreedingCoeff < -0.2  
	- AC = 0  
	- DP < 10  
	- GQ < 20  

#### **#Strategy**  

1. Make vcf composed of only biallelic SNPs (using ```bcftools```)  
2. Annotate InbreedingCoeff to vcf (using ```GATK-VariantAnnotator```)  
3. Make genotype filter of condition to 'GQ<20 \|\| DP<10', and set genotype as missing (./.) to the filtered samples (using ```GATK-VariantFiltration```)  
4. Filter out variants with AC = 0 (```GATK-SelectVariants```)  

#### **#Script**

```bash
## Make biallelic SNPs only
bcftools view -m2 -M2 -v snps example.vcf > example.bi.vcf

## Annotate InbreedingCoeff
java -jar GenomeAnalysisTK.jar  -T VariantAnnotator   -R ucsc.hg19.fasta  \
 -V example.bi.vcf   -o example.bi.Inbreed.vcf  -A InbreedingCoeff

## VariantFiltration
java -jar GenomeAnalysisTK.jar  -T VariantFiltration -R ucsc.hg19.fasta \
 -V example.bi.Inbreed.vcf   -o example.bi.Inbreed.GT.vcf \
 -G_filter "GQ<20||DP<10"   -G_filterName 'GQ20DP10' \
 --setFilteredGtToNocall    --missingValuesInExpressionsShouldEvaluateAsFailing

## SelectVariants
java -jar GenomeAnalysisTK.jar  -T SelectVariants  -R ucsc.hg19.fasta \
  -V example.bi.Inbreed.GT.vcf \
  -o example.bi.Inbreed.GTfiltered.vcf \
  -select "AC > 0 " -select "InbreedingCoeff >= -0.2"
```
 - ```--setFilteredGtToNocall``` Set filtered genotypes to no-call  
 - ```--missingValuesInExpressionsShouldEvaluateAsFailing``` When evaluating the JEXL expressions, missing values should be considered failing the expression


---
<br>

### **<u>2) Functional annotation</u>**

#### **#Paper description**  
VEP v81 was used with Gencode v19 on GRCh37  
Protein-truncating variant (PTV) annotation using LOFTEE  

- All stop-gained, splice-disrupting, and frameshift variants
- Ensembl Gene ID
- Gene symbol
- Ensembl Transcript ID
- PolyPhen2
- SIFT
- CADD

#### **#Strategy**

1. Set up dbNSFP and LOFTEE plugin. [Read me](https://github.com/Ensembl/VEP_plugins/blob/release/94/dbNSFP.pm)  
2. Run VEP with Ensembl Gene ID, Transcript ID, Gene symnol and predictors (PPh2, SIFT, and CADD).

### **#Script**

```bash
vep --cache --dir_cache /data/software/VEP/ensembl-vep-release-90/cache \
 --assembly GRCh37 --port 3337 --symbol --pick --vcf \
 -i example.vcf  -o  example.LOFTEE.vcf \
 --plugin LoF,human_ancestor_fa:/data/software/VEP/ensembl-vep-release-90/cache/Plugins/LoF/human_ancestor.fa \
 --dir_plugins /home/jsh/Plugins

vep --cache --dir_cache /data/software/VEP/ensembl-vep-release-90/cache \
 --assembly GRCh37 --port 3337 --symbol --pick \
 -i example.LOFTEE.vcf -o example.LOFTEE.dbNSFP.vcf \
 --plugin dbNSFP,/path/to/dbNSFP.gz,col1,col2 
```

---
<br>

Go to other sections,  
[1. Interpretation of VCF]({% post_url 2018-11-01-wes-vcf%})  
[2. Filtering]({% post_url 2018-11-01-wes-filtering%})  
[3. Annotation]({% post_url 2018-11-01-wes-annotation%})  

or [main page]({{ site.baseurl }})