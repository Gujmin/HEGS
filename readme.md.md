# Guideline for HEGS
## Run Encrypt and Decrypt Locally
Use `Encrypt.r` to convert VCF→PLINK, align phenotypes to `.fam`, run HIBLUP to compute `ADphe`, and overwrite the target phenotype column with an orthogonal-matrix–encrypted value (outputs: `HE_phenotype.txt`, `encrypter.rds`). Use `Decrypt.r` with the encrypted model output (e.g., HIBLUP `SOL`) and the saved `encrypter.rds` to recover values on the original scale for a given sample size. Run with `Rscript`; quote paths that contain spaces.

```bash
# Encrypt (Encrypt.r)
Rscript /var/www/html/HEGS/extend/Encrypt.r   --workPath=/your/path/   --phenotypePath=/your/path/data/phe.tsv   --genotypePath=/your/path/data/geno.vcf.gz   --pedigreePath=""   --selectedTraits=2   --selectedFixedFactors=""   --selectedCovariates=""   --selectedRandomOtherEffects=""   --selectedRandomGeneticEffects=1   --plink=/your/path/to/plink   --gcta=/your/path/to/gcta64   --hiblup=/your/path/to/hiblup
```

```bash
# Decrypt (Decrypt.r)
Rscript /var/www/html/HEGS/extend/Decrypt.r   --workPath=/your/path/   --resultPath=/your/path/to/SOL   --matrixPath=/your/path/to/encrypter.rds   --SampleSize=xxxx
```


## HELP for HEGS web
### 1. Encryption Genomic Selection

The "Encryption Genomic Selection" page is designed for **independent** encrypted genomic selection. Users need to upload **encrypted** phenotype and genotype files. The phenotype file should include **header information **to identify the factors in the model, while the genotype file must be in **a compressed genotype inverse matrix format**. We provide R code to help users process data locally, and there is also an online function available on the webpage. After uploading the phenotype file with the header, users can select the corresponding model factors, including:**Traits to be selected, Fixed effects(factor), Fixed effects(covarite), Random effect(genetic), Random effect(other)**. Finally, after entering the username and email address, users can click the submit button. Once the analysis is complete, the results will be sent to the users **email**.

![Figure 1: schematic diagram of encrypted genome selection process](<./picture/encryptGS.png>)

### 2. Encryption Joint Genomic Selection

The "Encryption Joint Genomic Selection" page is designed for conducting **joint** encrypted genomic selection. On this page, we have provided three functions for users: **creating a task, joining a public task, and joining a task provided by our laboratory.** The information included in a task consists of: data provider, breed information, sample size, traits, downloadable header files, and SNP Marker files. Users can also create public tasks, and in addition to the aforementioned information, they will need to specify the factors in the model. When performing joint genomic selection, users need to format the data according to the task requirements, **replacing missing information with NA**. If the **headers and SNP numbers do not match during the subsequent merging process**, the program will report an error. Similarly, once the analysis is complete, the results will be sent to the users **email**.

![Figure 2: schematic diagram of Joint Genome selection collaboration](<./picture/encryptJGS.png>)

### 3. Homomorphic Encryption

The "Homomorphic Encryption" page is designed for conducting homomorphic encryption and decryption. Uploaded data is converted into a homomorphic encryption format. Once encrypted, **a download link** for the data is sent to the user via email. The decryption process follows a similar procedure. At the same time, the platform provides **R code** for encryption and decryption, so that users can run the R code locally to encrypt phenotypic and genotypic data. During genomic selection analysis, the platform performs computations directly on the encrypted data, ensuring that data remains encrypted throughout the entire process. The use of homomorphic encryption facilitates this direct computation without requiring decryption. After completing the analysis, the platform deletes the user-uploaded data. This comprehensive approach ensures the security of user data during transmission, storage, and processing, effectively preventing data leakage and unauthorized access.

![Figure 3: schematic diagram of homomorphic encryption technology](<./picture/encrypt.png>)

### 4. Regular Genomic Selection

In this platform, we also provide genome selection functions that do not require encryption. This section also supports operations **BLUP, GBLUP and variance component estimation**. For different functional modules, users need to upload the corresponding files. For example, when performing BLUP analysis, users need to upload the phenotype file and pedigree file, and specify the factors in the model. Please note that encrypted data should not be uploaded in this module to avoid incorrect results. Once the analysis is complete, the uploaded files will be **automatically deleted from the server**.

![Figure 4: routine genome selection analysis process](<./picture/regular.png>)


### 5. Calculate Breeding Value

The "Calculate Breeding Value" page is designed for enabling users to upload genotype data and generate genomic breeding values (GEBV) for individuals. Users upload their genotypic data files via the platform interface, ensuring compatibility with Plink software requirements.Once the data is uploaded, the platform performs automatic quality control by removing SNPs with high deletion rates and conducting Hardy-Weinberg equilibrium tests. Following data preprocessing, the platform calculates breeding values by evaluating the effects of loci on various traits and predicting each individuals GEBV based on these effects. After the calculation is complete, results are presented in an easily interpretable format, for data security reasons, we will provide **a ranking of breeding values** for each trait in the results. This automated process enhances data processing efficiency and ensures accuracy and reliability. The HEGS platform thus provides users with straightforward access to breeding value information, facilitating better breeding decisions.

![Figure 5: Calculation and analysis process of breeding value](<./picture/EBV.png>)
