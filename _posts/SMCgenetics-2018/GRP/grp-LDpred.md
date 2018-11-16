---
layout: page
permalink: /GRP/grp-LDpred
---

## **2. LDpred**

Bayesian Shrinkage
- Shrinkage strategy: Prior distribution, e.g. fraction of causal SNPs
- Handling Linkage Disequilibrium: Shrink effect sizes with respect to LD

---
<br>
LDpred is a Python based software package that adjusts GWAS summary statistics for the effects of linkage disequilibrium (LD). The details of the method is described in [Vilhjalmsson et al. (AJHG 2015)](http://www.cell.com/ajhg/abstract/S0002-9297(15)00365-1){: target="blank" } 

The current version is 0.9.9

Github Page: [https://github.com/bvilhjal/ldpred](https://github.com/bvilhjal/ldpred){: target="blank" } 
<br>
Copyright (C) 2015, Bjarni J. Vilhjalmsson (bjarni.vilhjalmsson@gmail.com)
<br>
<br>
#### **Using LDpred**

To run LDpred, at least two steps are required:

1. The first step is a data synchronization step, where two or three data sets, genotypes and summary statistics are synchronized. This generates a HDF5 file which contains the synchronized genotypes. This step is implemented in the **** coord **** script. This step requires at least one genotype file (the LD reference genotypes), where we recommend at least 1000 unrelated individuals with the same ancestry make-up as the individuals for which summary statistics datasets are obtained from. Another genotype file can also be given if the user intends to validate the predictions using a separate set of genotypes.

2. After generating the coordinated data file then the one can apply LDpred and run it on the synchronized dataset. This step is implemented in **** ldpred **** script. This step generates two files, a LD file with LD information for the given LD radius, and the re-weighted effect estimates. The LD file enables the user to not have to generate the LD file again when trying, e.g., different values of p (the fraction of causal variants). However, it is re-generated if a different LD radius is given. The other file that LDpred generates contains the LDpred-adjusted effect estimates.

#### **Generating individual risk scores**
Individual risk scores can be generated using the **** validate **** script. It calculates polygenic risk scores for the individuals in the validation data if given, otherwise it treats the LD reference genotypes as validation genotypes. A phenotype file can be provided, covariate file, as well as plink-formatted principal components file.
<br>
<br>

#### **Installing LDpred**
```bash
$ pip install ldpred
$ pip install plinkio
$ git clone https://github.com/bvilhjal/ldpred.git
```

---
<br>

### **1) coord_genotypes.py**

Coordinate genotypes and summary statistics datasets for calculating polygenic risk scores. Only SNPs that overlap between the two (or three) different datasets are retained.

**Usage:**
```bash
$ coord --gf=PLINK_LD_REF_GENOTYPE_FILE --ssf=SUM_STATS_FILE --N=SS_SAMPLE_SIZE --out=OUT_COORD_FILE 
```

PLINK_LD_REF_GENOTYPE_FILE (and PLINK_VAL_GENOTYPE_FILE) should be a (full path) filename prefix to a standard PLINK bed file (without .bed) Make sure that the fam and bim files with same names are in the same directory. PLINK_LD_REF_GENOTYPE_FILE refers LD reference genotypes, and PLINK_VAL_GENOTYPE_FILE refers to validation genotypes. It is not necessary to have LD validation genotypes at this stage.

SUM_STATS_FILE should be a (full path) filename prefix for a text file with the GWAS summary statistics. Several formats are supported, see below.

SS_SAMPLE_SIZE should be the approximate number of individuals used for calculating the GWAS summary statistics.

OUT_COORD_FILE is the output file. This file will follow a HDF5 format and contain both LD-reference genotypes and summary statistics.

**Try:**
```bash
$ coord --gf=train --ssf=ss.txt --N=8000 --out=coord 
```
---
<br>

### **2) LDpred.py**

Implements LDpred, an approximate Gibbs sampler that calculate posterior means of effects, conditional on LD information. The method requires the user to have generated a coordinated dataset using coord_genotypes.py

**Usage:**
```bash
$ ldpred --coord=COORD_DATA_FILE --ld_radius=LD_RADIUS --local_ld_file_prefix=LD_FILE_NAME \
--PS=FRACTIONS_CAUSAL --N=SAMPLE_SIZE --out=OUTPUT_FILE_PREFIX
```

COORD_DATA_FILE: The HDF5 file obtained by running the coord_genotypes.py

LD_RADIUS: An integer number which denotes the number of SNPs on each side of the focal SNP for which LD should be adjusted. A value corresponding M/3000, where M is the number of SNPs used for the analysis is reasonable for genome length of 3000Mb. This should result in a LD-radius of about 1Mb on average.

LD_FILE_NAME: A path and filename prefix for the LD file. If it doesn't exist, it will be generated. This can take up to several hours, depending on LD radius number of SNPs, etc. If it does exits, that file will be used.

FRACTION_CAUSAL: A list of comma separated (without space) values between 1 and 0, excluding 0. 1 corresponds to the infinitesimal model and will yield results similar to LDpred-inf. Default is --PS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001

N: This is the sample size which LDpred assumes was used to calculate the GWAS summary statistics.

OUTPUT_FILE_PREFIX: The prefix of output file.

**Try:**
```bash
$ ldpred --coord=coord --ld_radius=1000000 --local_ld_file_prefix=coord \
--PS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001 --N=8000 --out=ldpred
```

```
Calculating LD information w. radius 1000000
Working on chrom_1
1000 400
Done calculating the LD table and LD score, writing to file: /home/ishim1/GRP/coord_ldradius1000000.pickled.gz
Genome-wide average LD score was: 9.162306170225143
LD information is now pickled.
Applying LDpred with LD radius: 1000000
Genome-wide lambda inflation: 1.3248282588753428 Genome-wide mean LD score: 9.162306170225143
Estimated genome-wide heritability: 0.004431584320044624
Calculating LDpred-inf weights
Calculating scores for Chromosome 1

Starting LDpred with p=1.0000
Calculating scores for Chromosome 1
The R2 prediction accuracy of PRS using chrom_1 was: 0.0424
There were 1000 (SNP) effects
Prediction accuracy was assessed using 400 individuals.
The  R2 prediction accuracy (observed scale) for the whole genome was: 0.0424 (0.002292)
203 cases, 197 controls
AUC: 0.6196
AUC for the whole genome was: 0.6196
The slope for predictions with P-value derived  effects is: 6.692116360826834

Starting LDpred with p=0.3000
Calculating scores for Chromosome 1
The R2 prediction accuracy of PRS using chrom_1 was: 0.0422
There were 1000 (SNP) effects
Prediction accuracy was assessed using 400 individuals.
The  R2 prediction accuracy (observed scale) for the whole genome was: 0.0425 (0.002292)
203 cases, 197 controls
AUC: 0.6191
AUC for the whole genome was: 0.6191
The slope for predictions with P-value derived  effects is: 3.067052772232231

(...)

```
---
<br>
### **3) validate.py**

Takes LDpred.py (or LD_pruning_thres.py) effect estimates , and (validation) genotypes in PLINK bed format as input. The script then works out overlap and outputs predictions or risk scores as well as some prediction accuracy statistics.

Note that for maximal accuracy all SNPs with LDpred weights should be included in the validation dataset. If they are a subset of the validation dataset, then we suggest recalculate LDpred for the overlapping SNPs.

**Usage:**
```bash
$ validate --vgf=PLINK_VAL_GENOTYPE_FILE --rf=RESULT_FILE_PREFIX --out=OUTPUT_FILE_PREFIX 
```

PLINK_VAL_GENOTYPE_FILE: PLINK formatted genotypes for which we want to calculate risk scores.

RESULT_FILE_PREFIX: SNP weights file, e.g. LDpred SNP weights.

OUTPUT_FILE_PREFIX: The prefix of output file.


**Try:**
```bash
$ validate --vgf=test --rf=ldpred --out=val
```

```
Calculating LDpred-inf risk scores
Parsing PLINK bed file: /home/ishim1/GRP/test
100 individuals have phenotype and genotype information.
Iterating over BED file to calculate risk scores.
DONE!
Number of non-matching NTs: 0
Number of flipped NTs: 0
Raw effects PRS correlation: 0.0930
Raw effects PRS r2: 0.0086
Weigted effects PRS correlation: -0.0766
Weigted effects PRS r2: 0.0059
Final raw effects PRS correlation: 0.0930
Final raw effects PRS r2: 0.0086
Final weighted effects PRS correlation: -0.0766
Final weighted effects PRS r2: 0.0059
The slope for predictions with raw effects is: 0.0998698504764381
The slope for predictions with weighted effects is: -2.348045404802415

Calculating LDpred risk scores using p=1.000e+00
Parsing PLINK bed file: /home/ishim1/GRP/test
100 individuals have phenotype and genotype information.
Iterating over BED file to calculate risk scores.
DONE!
Number of non-matching NTs: 0
Number of flipped NTs: 0
Raw effects PRS correlation: 0.0930
Raw effects PRS r2: 0.0086
Weigted effects PRS correlation: -0.0775
Weigted effects PRS r2: 0.0060
Final raw effects PRS correlation: 0.0930
Final raw effects PRS r2: 0.0086
Final weighted effects PRS correlation: -0.0775
Final weighted effects PRS r2: 0.0060
The slope for predictions with raw effects is: 0.0998698504764381
The slope for predictions with weighted effects is: -2.364484736020801

Calculating LDpred risk scores using p=3.000e-01
Parsing PLINK bed file: /home/ishim1/GRP/test
100 individuals have phenotype and genotype information.
Iterating over BED file to calculate risk scores.
DONE!
Number of non-matching NTs: 0
Number of flipped NTs: 0
Raw effects PRS correlation: 0.0930
Raw effects PRS r2: 0.0086
Weigted effects PRS correlation: -0.0583
Weigted effects PRS r2: 0.0034
Final raw effects PRS correlation: 0.0930
Final raw effects PRS r2: 0.0086
Final weighted effects PRS correlation: -0.0583
Final weighted effects PRS r2: 0.0034
The slope for predictions with raw effects is: 0.0998698504764381
The slope for predictions with weighted effects is: -1.4324891224974239

(...)
```

**Output**
```bash
$ vi val_LDpred-inf.txt
```
```
IID, true_phens, raw_effects_prs, pval_derived_effects_prs, sex
i0, 1.000000e+00, -1.227203e+00, 3.239673e-02, 2,
i1, 2.000000e+00, 9.557595e-01, -3.363165e-02, 2,
i2, 2.000000e+00, -1.745532e+00, 5.763071e-02, 2,
i3, 2.000000e+00, 1.647665e+00, -5.947706e-02, 2,
i4, 2.000000e+00, -1.196074e+00, 3.953134e-02, 2,
i5, 1.000000e+00, -4.813022e-02, 1.818085e-03, 2,
i6, 1.000000e+00, -9.330380e-02, 2.491208e-03, 2,
i7, 2.000000e+00, 8.156840e-01, -3.784357e-02, 2,
i8, 2.000000e+00, -3.120926e-01, 6.808986e-03, 2,
i9, 2.000000e+00, 5.569845e-01, -2.826669e-02, 2,
i10, 1.000000e+00, 1.127705e-01, -2.645583e-03, 2,
```

---
<br>
### **Comparison of Methods**

![image1](/assets/images/GRP/Figure1.png){: width="60%" height="60%"}

Figure 2. Comparison of Four Prediction Methods Applied to Simulated Traits 
[Vilhjalmsson et al. (AJHG 2015)](http://www.cell.com/ajhg/abstract/S0002-9297(15)00365-1){: target="blank" } 

![image2](/assets/images/GRP/Figure2.png){: width="60%" height="60%"}

Figure 4. Comparison of Methods Training on Large GWAS Summary Statistics for Five Different Diseases 
[Vilhjalmsson et al. (AJHG 2015)](http://www.cell.com/ajhg/abstract/S0002-9297(15)00365-1){: target="blank" } 


---
<br>
### **Genetic Risk Prediction**

![image3](/assets/images/GRP/Figure3.png)

Fig. 2 | Risk for CAD according to GPS 
[Khera et al. (NatGenet 2018)](https://www.nature.com/articles/s41588-018-0183-z){: target="blank" } 


---
<br>

Go to previous step,  
[1. LD clumping + thresholding]({% post_url 2018-11-01-grp-PT %})

or [main page](../)