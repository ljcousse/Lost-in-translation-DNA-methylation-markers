# Lost in translation: DNA methylation markers

This repository contains the key methodological scripts written for our manuscript "Lost in translation: How CpG site selection and assay design determine the clinical value of cancer DNA methylation markers"

Developed and tested on a local Unix-based server, using bash scripting
for workflow automation and conda for systematic package
and dependency administration. Besides stand-alone tools all analyses were performed in R (version 4.3.2).  

### Obtain the code for this project

```bash
git lfs install
git clone https://github.com/ljcousse/Lost-in-translation-DNA-methylation-markers.git
cd Lost-in-translation-DNA-methylation-markers

conda env create -f environment.yml
conda activate lost_in_translation
```

### Data
The targeted bisulfite sequencing data generated in this study, comprising both raw and processed files, have been deposited in the NCBI Gene Expression Omnibus (GEO) under accession number GSE338748. The data are currently private and accessible to reviewers via a secure reviewer token; they will be made publicly available upon publication.

### License
© 2026 Louis Coussement. All rights reserved. This code is provided for review and reproducibility purposes only. No permission is granted to use, copy, modify, or distribute it without the author's written consent.

### Citation
Publication will hopefully be coming soon:

Lost in translation: How CpG site selection and assay design determine the clinical value of cancer DNA methylation markers. Kim Lommen*, Louis Coussement*, Johan Vandersmissen, James G. Herman, Wim Van Criekinge, Kim M. Smits, Tim De Meyer**, Manon van Engeland**

* co-first authors
** co-last authors


