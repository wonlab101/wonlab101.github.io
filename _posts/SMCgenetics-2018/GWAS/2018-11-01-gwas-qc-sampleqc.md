---
layout: page
permalink: /SMCgenetics-2018/GWAS/gwas-qc-sampleqc
---

go back to previous page, [**1. Data Quality Control**]({% post_url 2018-11-01-gwas-qc %}).
### **2) Sample Quality Control**

---
<br>

#### a. Filtering by threshold
- remove missingness rate ≥0.05 individuals :  ```--mind 0.05```
- generate a lists of missingness rate statistics :  ```--missing```
- generate a lists of heterozygosity rate statistics :  ```--het```

```bash
$ plink --bfile QC/gwas.1 --mind 0.05 --missing --het --make-bed --out QC/gwas.2
```

```bash
$ head QC/gwas.2.imiss
$ head QC/gwas.2.lmiss
$ head QC/gwas.2.het
```

---
<br>

#### b. Heterozygosity rate & Missingness rate of individuals

- remove outliers of heterozygosity rate or missingness rate (e.g. >\|5σ\|)

![QCflow](/assets/images/GWAS/gwas.2.het.png){: width="60%" height="60%"}

<br>

\> draw plot using R
```bash
$ R
```
```R
#!/bin/R
# draw Heterozygosity X Missingness plot
# and extract outliers (>|5sd|)

library(ggplot2)

file.imiss <- "QC/gwas.2.imiss"
file.het <- "QC/gwas.2.het"
out.prefix <- "QC/gwas.2"

## load data
df.imiss <- read.table(file.imiss, header=T)
df.het <- read.table(file.het, header=T)

## calculate heterozygosity
df.het <- cbind(df.het,
                het_rate=(df.het$N.NM. - df.het$O.HOM.)/df.het$N.NM.)

df <- merge(df.imiss, df.het, by=c("FID", "IID"))


## prepare parameters
df$F_MISS <- as.numeric(df$F_MISS)
df$het_rate <- as.numeric(df$het_rate)

ymax = max(df$F_MISS)*1.01
xmin = min(df$het_rate)*0.99
xmax = max(df$het_rate)*1.01

m.miss = mean(df$F_MISS)
sd.miss = sd(df$F_MISS)
m.het = mean(df$het_rate)
sd.het = sd(df$het_rate)

df$het_DST <- (df$het_rate - m.het)/sd.het


## plot
p <- ggplot()
p <- p + geom_point(data=df, aes(x=het_rate, y=F_MISS), size=0.8, color="gray20")
p <- p + xlab("Heterozygosity") + ylab("Missingness")
p <- p + theme_classic()
# threshold::missingness
p <- p + geom_hline(yintercept=m.miss+3*sd.miss, color="blue", lty=2, size=0.2)
p <- p + geom_hline(yintercept=m.miss+5*sd.miss, color="red", lty=2, size=0.2)
# threshold::heterozygosity
p <- p + geom_vline(xintercept=c(m.het-3*sd.het, m.het+3*sd.het), color="blue", lty=2, size=0.2)
p <- p + geom_vline(xintercept=c(m.het-5*sd.het, m.het+5*sd.het), color="red", lty=2, size=0.2)

pdf(paste0(out.prefix, ".het.pdf"), width=4, height=4)
p
dev.off()


## extract outliers
# missingness
miss_fail <- subset(df, F_MISS > m.miss+5*sd.miss)
write.table(miss_fail, paste0(out.prefix, ".miss_5sdout"), row.names=FALSE, quote=FALSE, sep="\t")

# heterozygosity
het_fail <- subset(df, (het_rate < m.het-5*sd.het) | (df$het_rate > m.het+5*sd.het))
write.table(het_fail, paste0(out.prefix, ".het_5sdout"), row.names=FALSE, quote=FALSE, sep="\t")


quit()
```

\> remove outliers
```bash
$ cat QC/gwas.2.miss_5sdout QC/gwas.2.het_5sdout | cut -f1,2 > QC/rmlist.het
$ plink --bfile QC/gwas.2 --remove QC/rmlist.het --make-bed --out QC/gwas.2.rmhet
```


---
<br>

#### c. Check sex

- Compare reported sex in the **\*.fam** file with estimated sex (given genotype data).
- "STATUS" column in output file (\*.sexcheck) displays "PROBLEM" or "OK" for each individual.
- male: F>0.8, Female: F<0.2

```bash
$ plink --bfile QC/gwas.2.rmhet --check-sex --out QC/gwas.2.rmhet
```

```bash
$ head QC/gwas.2.rmhet.sexcheck
```

```bash
$ awk '{if($5=="PROBLEM") print $1,$2}' QC/gwas.2.rmhet.sexcheck > QC/rmlist.sex
$ plink --bfile QC/gwas.2.rmhet --remove QC/rmlist.sex --make-bed --out QC/gwas.2.rmhet.rmsex
```


---
<br>

#### d. Relatedness

- Using [KING tool](http://people.virginia.edu/~wc9c/KING/), we can estimate the relationships between samples.
- remove one of related pairs of individuals with the 2nd or closer relationship.

\> make directory and go into the directory
```bash
$ mkdir KING
$ cd KING
```

\> make bfiles for running KING tool (leave only common variants (MAF ≥0.05) in autosomal chromosomes)
```bash
$ plink --bfile ../QC/gwas.2.rmhet.rmsex --chr 1-22 --maf 0.05 --geno 0.01 --make-bed --out gwas.king
```

\> run KING tool
```bash
$ king -b gwas.king.bed --unrelated --degree 2 --prefix gwas.king > king.log
$ head gwas.kingunrelated_toberemoved.txt    # one of KING output files: list of samples to be removed
```

\> go to the parent directory and remove relatedness samples
```bash
$ cd ..
$ plink --bfile QC/gwas.2.rmhet.rmsex --remove KING/gwas.kingunrelated_toberemoved.txt --make-bed --out QC/gwas.2.rmhet.rmsex.unrelated
```


---
<br>

#### e. PCA

- Using **smartpca** tool of [EIGENSOFT software](https://www.hsph.harvard.edu/alkes-price/software/){: target="blank" }, we can perform PCA.
- remove outliers (e.g. >\|5σ\|) or check homogeneity of genetic variation (population structure).

![pcaplot](/assets/images/GWAS/gwas.pca.png){: width="60%" height="60%"}

<br>

\> make directory and go into the directory
```bash
$ mkdir PCA
$ cd PCA
```

\> make bfiles for running smartpca (leave only common variants (MAF ≥0.05) in autosomal chromosomes)
```bash
$ plink --bfile ../QC/gwas.2.rmhet.rmsex.unrelated --chr 1-22 --maf 0.05 --geno 0.01 --make-bed --out gwas.pca
```

\> perfom [LD pruning](http://zzz.bwh.harvard.edu/plink/summary.shtml#prune){: target="blank" } to leave only independent variants (smartpca input file)
```bash
# make independent variants list file (*.prune.in)
$ plink --bfile gwas.pca --indep-pairwise 50 5 0.5 --out gwas.pca

# extract independent variants
$ plink --bfile gwas.pca --extract gwas.pca.prune.in --make-bed --out gwas.pca.pruned
```

\> make smartpca parameter file  
&nbsp;&nbsp; ※ *Tips for using vi text editor*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; i : insert mode on (you can type text on the vi text editor)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; esc : insert mode off (replace mode on)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :qw : save file and quit the vi text editor (at replace mode)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :q! : **do not** save file and quit the vi text editor (at replace mode)
```bash
$ vi smartpca.par
```
```bash
# type below text on the vi text editor
genotypename: gwas.pca.pruned.bed
snpname: gwas.pca.pruned.bim
indivname: gwas.pca.pruned.fam
evecoutname: ./gwas.pca.pruned.evec
evaloutname: ./gwas.pca.pruned.eval
altnormstyle: NO
numoutevec: 10
numoutlieriter: 0
numoutlierevec: 10
outliersigmathresh: 6
```
```bash
:qw
```

\> run smartpca
```bash
$ smartpca -p smartpca.par > smartpca.log
```

\> draw plot using R
```bash
$ R
```
```R
#!/bin/R
# draw biplot by smartpca output (*.evec file)
# and extract outliers (>|5sd|)

library(ggplot2)

file = "gwas"
evecfile = paste0(file, ".pca.pruned.evec")
outfile = paste0(file, ".pca")

## load smartpca output (*.evec file)
evec <- read.table(evecfile, sep="", header=FALSE)
names(evec) <- c("ID", paste0("PC", 1:(ncol(evec)-2)), "pheno")


## calculate mean and standard deviation of PC1 and PC2
m1 = mean(evec[,2])
m2 = mean(evec[,3])
s1 = sd(evec[,2])
s2 = sd(evec[,3])


## plotting
plot <- ggplot(data=evec, aes(x=PC1, y=PC2))
plot <- plot + geom_vline(xintercept=c(m1-3*s1, m1+3*s1), col='blue', lty=2)
plot <- plot + geom_hline(yintercept=c(m2-3*s2, m2+3*s2), col='blue', lty=2)
plot <- plot + geom_vline(xintercept=c(m1-5*s1, m1+5*s1), col='red')
plot <- plot + geom_hline(yintercept=c(m2-5*s2, m2+5*s2), col='red')
plot <- plot + geom_point(size=1.5, alpha=0.7)
plot <- plot + theme_bw()


## make outlier list
outlierlist <- function(df, threshold=5){
  # threshold : n sigma (5 sigma by default)
  m1 = mean(df$PC1)
  m2 = mean(df$PC2)
  s1 = sd(df$PC1)
  s2 = sd(df$PC2)

  rslt = subset(df, (df$PC1 < m1-threshold*s1) | (df$PC1 > m1+threshold*s1) | (df$PC2 < m2-threshold*s2) | (df$PC2 > m2+threshold*s2))

  return(rslt)
}

outlier_5sd <- outlierlist(evec, threshold=5)

outlier_5sd <- subset(evec, (evec$PC1 < m1-5*s1) | (evec$PC1 > m1+5*s1) | (evec$PC2 < m2-5*s2) | (evec$PC2 > m2+5*s2))

## writing
write.table(outlier_5sd, file = paste0(outfile, ".out5sd"), sep="\t", quote=FALSE, row.names=FALSE, fileEncoding="utf-8")
ggsave(plot=plot, filename=paste0(outfile, ".pdf"))


quit()
```

\> go to the parent directory and remove outliers
```bash
$ cd ..
$ plink --bfile QC/gwas.2.rmhet.rmsex.unrelated --remove PCA/gwas.pca.out5sd --make-bed --out QC/gwas.2.rmhet.rmsex.unrelated.rmpca
```

<br>

\> Quality controlled file
```bash
$ plink --bfile QC/gwas.2.rmhet.rmsex.unrelated.rmpca --chr 1-22 --make-bed --out gwas.qc
```

---
<br>

Let's go to nex step [**2. Imputation**]({% post_url 2018-11-01-gwas-imputation %}),
go back to previous page [**1. Data Quality Control**]({% post_url 2018-11-01-gwas-qc %}),
or go back to [**GWAS contents list**]({% post_url 2018-11-01-gwas %}).
