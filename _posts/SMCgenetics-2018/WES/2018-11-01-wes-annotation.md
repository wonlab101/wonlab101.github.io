---
layout: page
permalink: /SMCgenetics-2018/WES/wes-annotation
---

## **3. Annotation**

In this session, we are going to annotate variant information to each variant usnig VEP (Variant Effect Predictor) version 90.6. With VEP, we are able to add gene symbol, population frequency, prediction score, etc. Details of VEP in [here](http://asia.ensembl.org/info/docs/tools/vep/index.html){: target="blank" }    
(**Note** The meagning of 'annotation' in this section is close to **'functional annotation'** that is different from meaning of annotation used in previous section such as QD, MQ of VCF)


### **<u>1) Run VEP</u>**  
  
VEP can use VCF, variant identifiers and HGVS notations in addition to default format as input. We practice annotation only with **VCF** file. Details in [here: run vep](http://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html){: target="blank" }

- Basic synopsis of vep:

```bash
./vep [options] -i input.txt -o output.txt
```

---
<br>

### **<u>2) Practice: Four examples</u>**

#### **#Example 1) Simple annotation**

- Input: BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf
- Output: BGDA.vcf.vep
- Script: 01_vep.sh

```bash
vep --cache --dir_cache /data/software/VEP/ensembl-vep-release-90/cache \
 --port 3337 --assembly GRCh37 \
 -i BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf  \
 -o BGDA.vcf.vep \
 --no_stats
```

- Let's open the output file.

```bash
vi BGDA.vcf.vep
```

![00 output](https://user-images.githubusercontent.com/26876362/48534518-b8315700-e8eb-11e8-9f44-82271ec07265.png)


- Find variants:  
	- synonymous: rs5747988
	- missense: rs114306778
	- missense_variant,splice_region_variant: rs1129172
	- non_coding_transcript_exon_variant: rs6518594

```bash
# in vi editor,
/rs6518594 
```

![rs6518594 no pick](https://user-images.githubusercontent.com/26876362/48600124-fbec9500-e9ad-11e8-966c-995301deebda.png)

- Q1. How do I pick one of the variants? 
- Q2. How to add gene symbol to variants?

---
<br>

#### **#Example 2) Annotate gene symbols and pick one annotation per variant**

- Input: BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf
- Output: BGDA.vcf.gene.vep
- Script: 02_vep_gene_symbol.sh  

```bash
vep  --cache --dir_cache /data/software/VEP/ensembl-vep-release-90/cache \
 --port 3337 --assembly GRCh37 \
 -i BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf \
 -o BGDA.vcf.gene.vep \
 -symbol --pick \
 --no_stats
```
  
  - ```--pick``` Pick once line or block of consequence data per variant, including transcript-specific columns. Consequences are chosen according to the criteria described here, and the order the criteria are applied may be customised with ```--pick_order```. This is the best method to use if you are interested only in one consequence per variant.  
	- pick_order
		1. canonical status of transcript
		2. APPRIS isoform annotation
		3. transcript support level
		4. biotype of transcript ("protein_coding" preferred)
		5. CCDS status of transcript
		6. consequence rank according to this table
		7. translated, transcript or feature length (longer preferred)

  - ```--symbol``` Adds the gene symbol (e.g. HGNC) (where available) to the output.  
  - Additional description is in [here](http://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html#output){: target="blank" } 
<br>

- Let's open the output file and find rs6518594.  

```bash
vi BGDA.vcf.gene.vep

# in vi editor,
/rs6518594 
```

![02 output](https://user-images.githubusercontent.com/26876362/48600220-656ca380-e9ae-11e8-8502-c965a168b82f.png)

- rs6518594 

![rs6518594](https://user-images.githubusercontent.com/26876362/48600174-33f3d800-e9ae-11e8-9792-aa102fdf3e09.png)

---
<br>

#### **#Example 3) Annotate population frequencies**

- Input: BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf
- Output: BGDA.vcf.gene.popfreq.vep
- Script: 03_vep_gene_popfreq.sh  

```bash
vep  --cache --dir_cache /data/software/VEP/ensembl-vep-release-90/cache \
 --port 3337 --assembly GRCh37 \
 -i BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf  -o BGDA.vcf.gene.popfreq.vep \
 --symbol --pick \
 --af_1kg --af_esp --af_gnomad \
 --no_stats
```

- Let's open the output file and confirm population frequencies.   

```bash
vi BGDA.vcf.gene.popfreq.vep
```   
![image](https://user-images.githubusercontent.com/26876362/48600414-0f4c3000-e9af-11e8-8aef-46f8c219c0c6.png)

- **Description (in header)**
	- AFR_AF/EAS_AF/.../gnomAD_SAS_AF: Frequency of existing variant in the population DB
	- CLIN_SIG: ClinVar clinical significance of the dbSNP variant
	- PHENO: Indicates if existing variant(s) is associated with a phenotype, disease or trait; multiple values correspond to multiple variants

---
<br>

#### **#Example 4) Annotate variant predictor**

To know what functional effect in my variant, we can use VEP with plugin modules. Many plugins in VEP and details in [here.](https://github.com/Ensembl/VEP_plugins){: target="blank" }  
This section uses [**'LOFTEE' module.**](https://github.com/konradjk/loftee){: target="blank" }  LOFTEE is a VEP plugin to identify LoF (loss-of-function) variation.


- Input: BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf
- Output: BGDA.vcf.gene.predict.vep
- Script: 04_vep_gene_predict.sh

```bash
vep --cache --dir_cache /data/software/VEP/ensembl-vep-release-90/cache \
 --port 3337 --assembly GRCh37 \
 -i BGDA.VPASS.snp.QD.AC.hwefiltered.recode.vcf \
 -o BGDA.vcf.gene.predict.vep \
 --symbol --pick \
 --af_1kg --af_esp --af_gnomad \
 --plugin LoF,human_ancestor_fa:/data/software/VEP/ensembl-vep-release-90/cache/Plugins/LoF/human_ancestor.fa \
 --dir_plugins /home/jsh/Plugins \
 --no_stats
```

- Let's open the output file and find rs2301558 and rs885985.

```bash
vi BGDA.vcf.gene.predict.vep

# in vi editor,
/rs2301558
/rs885985
```   

![rs2301558](https://user-images.githubusercontent.com/26876362/48600588-d9f41200-e9af-11e8-9982-1ce72adde37b.png)

![rs885985](https://user-images.githubusercontent.com/26876362/48600618-f001d280-e9af-11e8-959b-b6f21cfcf221.png)

- **Description (details in [here](https://github.com/konradjk/loftee))**  
	- LoF : Loss-of-function annotation (HC = High Confidence; LC = Low Confidence)
	- LoF_filter : Reason for LoF not being HC
	- LoF_flags : Possible warning flags for LoF
	- LoF_info : Info used for LoF annotation

---
<br>

Go to other sections,  
[1. Interpretation of VCF]({% post_url 2018-11-01-wes-vcf%})  
[2. Filtering]({% post_url 2018-11-01-wes-filtering%})  
[4. Homework]({% post_url 2018-11-01-wes-practice%})  


