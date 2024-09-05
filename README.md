# Annotation
I recently started working on annotating the genome of an Opossum sample we have. I am still learning which is the correct step by step to use for this, and I have also been using this pipeline (https://github.com/KrabbenhoftLab/genome_annotation_pipeline). However, in parallel, I am trying to run some steps that I will upload here.

Before doing the annotation using RNA-seq data, one needs to trim the RNA, and control for quality. Therefore, this tool called https://github.com/FelixKrueger/TrimGalore uses cutadapt and fastqc to do this, and can easily be installed using conda environments. Therefore, I am currently working on that and will upload the scripts for this here. 

1) cutadapt and fastqc are trimmers for adaptors and quality check, respectively.
2) hisat2 is for mapping your RNA reads into the genome. You can use multiple RNA-seq files or you can do it one by one. The manual recommends doing it per tissue. 
3) RepeatModeler does not use the RNA-seq files, but is important for identifying (and further down masking) the repeats on the genome. 
