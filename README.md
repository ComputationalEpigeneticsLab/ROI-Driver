# ROI-Driver

disorder or domain regions of interest which enriching cancer mutation

## Pan-cancer assessment of mutational landscape in intrinsically disordered hotspots reveals potential driver genes

Haozhe Zou, Tao Pan , Yueying Gao, Renwei Chen, Si Li, Jing Guo, Zhanyu Tian, Gang Xu, Juan Xu, Yanlin Ma, Yongsheng Li.

### Abstract
Large-scale cancer genome sequencing has enabled the catalogs of somatic mutations; however, the mutational impact on intrinsically disordered protein regions (IDRs) has not been systematically investigated to date. Here, we comprehensively characterized the mutational landscapes of IDRs and found that IDRs have higher mutation frequencies across diverse cancers. We thus developed a computational method, ROI-Driver, to identify putative driver genes enriching IDR and domain hotspots in cancer. Numerous well-known cancer-related oncogenes or tumor suppressors that play important roles in cancer signaling regulation, development and immune response were identified at a higher resolution. In particular, the incorporation of IDR structures helps in the identification of novel potential driver genes that play central roles in human protein-protein interaction networks. Interestingly, we found that the putative driver genes with IDR hotspots were significantly enriched with predicted phase separation propensities, suggesting that IDR mutations disrupt phase separation in key cellular pathways. We also identified an appreciable number of clinically relevant genes enriching IDR mutational hotspots that exhibited differential expression patterns and are associated with cancer patient survival. In summary, combinations of mutational effects on IDRs significantly increase the sensitivity of driver detection and are likely to open new therapeutic avenues for various cancers.

![image](https://user-images.githubusercontent.com/91582097/157180700-78111f05-3fcb-49b7-8ee0-1c379642af01.png)


## Running ROI-Driver

#### Running ROI-Driver starting with maf files (typical):

```R
Rscript ROI-Driver.R -i ./input_data/Input_mutation.txt -t domain -o ./output_data/output.txt
```

-i:input file

-t:region type(domain or disorder)

-o:output file

### Cite
Haozhe Zou, Tao Pan, Yueying Gao, Renwei Chen, Si Li, Jing Guo, Zhanyu Tian, Gang Xu, Juan Xu, Yanlin Ma, Yongsheng Li, Pan-cancer assessment of mutational landscape in intrinsically disordered hotspots reveals potential driver genes, Nucleic Acids Research, 2022;, gkac028, https://doi.org/10.1093/nar/gkac028
