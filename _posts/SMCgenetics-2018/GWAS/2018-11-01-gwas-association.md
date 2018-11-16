---
layout: page
permalink: /SMCgenetics-2018/GWAS/gwas-association
---

## **3. Association Analysis**

---
<br>

### **Hypothesis**

<br>

#### **> Hypothesis for GWAS**
H<sub>0</sub> : None of the SNP loci in these data are associated with the disease of interest.  
H<sub>1</sub> : At least one of the SNPs is associated with the trait in these data.


#### **> Hypothesis for each test**
H<sub>0</sub> : The SNP is not linearly associated with the trait in these data (effect size (β) = 0).  
H<sub>1</sub> : The SNP is linearly associated with the trait in these data (effect size (β) ≠ 0).

- Because we perform as many association analyses as the number of SNPs, multiple comparison correction is required.  
(e.g. Bonferroni correction)

---
<br>
<br>

### **1) Association Analysis**

---
<br>

#### a. Model

- view phenotype file

```bash
$ head phenotype.gwas
```

- set test model  
**LDL = β<sub>0</sub> + β<sub>1</sub>SNP + β<sub>2</sub>AGE + β<sub>3</sub>SEX + ε**  
: We are interested in the coefficient of the SNP, β<sub>1</sub>.

---
<br>

#### b. Linear regression

- linear regression : ```--linear```. It can be replaced to ```--logistic```, when you perform logistic regression.
- hide result of covariates : ```--linear hide-covar```
- view 95% confidence interval of beta : ```--ci 0.95```
- phenotype file : ```---pheno```, column name of phenotype : ```--pheno-name {column name}```
- covariate file : ```---covar```, column names of covariates : ```--covar-name {column1,column2,column3, ...}```

```bash
$ plink --bfile gwas.qc --linear hide-covar --ci 0.95 \
--pheno phenotype.gwas --pheno-name ldl \
--covar phenotype.gwas --covar-name age,sex \
--allow-no-sex --out gwas
```
```bash
$ less gwas.assoc.linear
```


---
<br>

#### c. Sorting the result file by P-value

- As mentioned above, we have to perform multiple comparison correction like [Bonferroni correction](https://en.wikipedia.org/wiki/Bonferroni_correction){: target="blank" }.
- *Bonferroni correction*  
If we want to test under the significance level=0.05, we have to correct the threshold of *P-value* to 0.05/8,998=5.56x10<sup>-6</sup>.  
(the number of test SNPs in our QCed data = 8,998)
- In GWAS study, 5x10<sup>-8</sup> is usually recommended for the threshold of *P-value* and it is called the genome-wide significance level.

```bash
cat gwas.assoc.linear | sort -k 12 -g > gwas.assoc.linear.sort
```
```bash
$ less gwas.assoc.linear.sort
```

---
<br>
<br>

### **2) Visualization**

We will draw a Quantile-Quantile (Q-Q) plot and Manhattan plot using "[**qqman**](https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html){: target="blank" }" package in R.

---
<br>

#### 0. Install and load packages

We don't need to *install packages* when we practice in the BGDA class.

```bash
$ R
```
```R
# install package
install.packages("qqman")
# load package
library(qqman)
```

---
<br>

#### a. Quantile-Quantile (Q-Q) plot

```R
# load a result file of the PLINK association test
df <- read.table("gwas.assoc.linear", sep="", header=T, na.strings="NA")

# draw and extract Q-Q plot
pdf("gwas_qq.pdf")
qq(df$P)
dev.off()
```

![qqplot](/assets/images/GWAS/gwas_qq.png)

---
<br>

#### b. Manhattan plot

```R
# draw and extract Manhattan plot
pdf("gwas_man.pdf")
manhattan(df, cex=0.6, cex.axis=0.8, ylim=c(0,10))
dev.off()

quit()
```

![manplot](/assets/images/GWAS/gwas_man.png)

---
<br>
<br>

### **3) Additional summary statistics**

We can get more summary statistics using PLINK. ([http://zzz.bwh.harvard.edu/plink/summary.shtml](http://zzz.bwh.harvard.edu/plink/summary.shtml){: target="blank" })

---
<br>

- minor allele frequency (\*.frq) : ```--freq```
- result of Hardy-Weinberg equilibrium test (\*.hwe) : ```--hardy```
- missingness by individual (\*.imiss) and by SNP (\*.lmiss) : ```--missing```

```bash
$ mkdir summary
```
```bash
$ plink --bfile gwas.qc --freq --hardy --missing --out summary/gwas.qc
```

\> Let's view summary files
```bash
$ cd summary
```
```bash
$ less gwas.qc.frq
$ less gwas.qc.hwe
$ less gwas.imiss
$ less gwas.lmiss
```

<br>

**※ *Tips for using PLINK***

- Summary statistics v.s. inclusion criteria

| Feature | As summary statistics | As inclusion criteria |
| ------- | ---------------------:| ---------------------:|
| Missingness per individual | `--missing` | `--mind N` |
| Missingness per marker     | `--missing` | `--geno N` |
| Allele frequency           | `--freq`    | `--maf N`  |
| Hardy-Weinberg equilibrium | `--hardy`   | `--hwe N`  |

---
<br>

Let's go to next step [**4. Browser Searching**]({% post_url 2018-11-01-gwas-browser %}),
or go back to [**GWAS contents list**]({% post_url 2018-11-01-gwas %}).
