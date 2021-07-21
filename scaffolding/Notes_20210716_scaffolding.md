# Brahman and Angus bionano genome maps
---
*16/07/2021*

## Table of Contents
* [Bionano map](#bn)

<a name="bn"></a>
## Brahman and Angus bionano data

I will be generating an upgraded Angus x Brahman assemblies using a lot more
PacBio and Oxford Nanopore reads. The first version of the assembly is
[published](https://www.nature.com/articles/s41467-020-15848-y).
After contig assembly, I'll need to redo the scaffolding with HiC data and Bionano.

Checking that I have the raw data from Bionano to be used for (i) scaffolding (ii)
checking copy number of certain genes.

How do I use the Bionano genome maps and align with the new assemblies? Do I have
the raw data already below?

```console
[a1223107@l01 bionano]$ ls
Angus_selected  Brahman_selected  compressed  Offspring  readme
[a1223107@l01 bionano]$ du -h .
23M	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/ref
44K	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/alignmolvref/copynumber
2.0G	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/alignmolvref/merge
2.0G	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/alignmolvref
80M	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/exp_refineFinal1/alignref_final
120M	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/exp_refineFinal1
1.1G	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise
90M	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/exp_refineFinal1_sv/merged_smaps
258M	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/exp_refineFinal1_sv
3.5G	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs
4.2G	./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results
4.0K	./Angus_selected/cur_results_2750695/fa2cmap
9.9G	./Angus_selected/cur_results_2750695/hybrid_scaffolds
11G	./Angus_selected/cur_results_2750695
15G	./Angus_selected
4.0K	./Offspring/cur_results_2750748/fa2cmap
7.9G	./Offspring/cur_results_2750748/hybrid_scaffolds
7.9G	./Offspring/cur_results_2750748
4.0K	./Offspring/cur_results_2750747/fa2cmap
7.8G	./Offspring/cur_results_2750747/hybrid_scaffolds
7.8G	./Offspring/cur_results_2750747
23M	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/ref
44K	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/alignmolvref/copynumber
3.9G	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/alignmolvref/merge
3.9G	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/alignmolvref
87M	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/exp_refineFinal1/alignref_final
136M	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/exp_refineFinal1
2.2G	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/auto_noise
104M	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/exp_refineFinal1_sv/merged_smaps
341M	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/exp_refineFinal1_sv
6.5G	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs
8.5G	./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results
24G	./Offspring
17G	./compressed
4.0K	./Brahman_selected/cur_results_2750697/fa2cmap
11G	./Brahman_selected/cur_results_2750697/hybrid_scaffolds
12G	./Brahman_selected/cur_results_2750697
24M	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/ref
44K	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/alignmolvref/copynumber
2.1G	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/alignmolvref/merge
2.1G	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/alignmolvref
83M	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/exp_refineFinal1/alignref_final
124M	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/exp_refineFinal1
1.2G	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise
88M	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/exp_refineFinal1_sv/merged_smaps
223M	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/exp_refineFinal1_sv
3.6G	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs
4.3G	./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results
16G	./Brahman_selected
72G	.
```

*16/07/2021*
Joyce: The molecules (.bnx) should be in the auto_noise folder, which is one of the folders you listed, so I think you have the raw data. The name of the molecule file is autoNoise1_rescaled.bnx.

```console
[a1223107@l01 bionano]$ pwd
/uofaresstor/avsci2/Lloyd_a1223107/big_migration_end2020/Angus_x_Brahman_data/bionano

[a1223107@l01 bionano]$ ls -lah ./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx
-rw-r--r-- 1 a1223107 phoenix-hpc-avsci 853M Nov 26  2020 ./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx

[a1223107@l01 bionano]$ ls -lah ./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx
-rw-r--r-- 1 a1223107 phoenix-hpc-avsci 875M Nov 26  2020 ./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx

[a1223107@l01 bionano]$ ls -lah ./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx
-rw-r--r-- 1 a1223107 phoenix-hpc-avsci 1.7G Nov 26  2020 ./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx


[a1223107@l01 bionano]$ md5sum ./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx
fe52ce3ec4e6c29304cc36e2f658ed23  ./Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx

[a1223107@l01 bionano]$ md5sum ./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx
2d3b44c74fb9f4ea035a9824962bdab4  ./Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx

[a1223107@l01 bionano]$ md5sum ./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx
0969987093a7767290177b3512f1cb68  ./Offspring/USDA_bull_fetus_nonhap_ES_sdb_08012018_pipeline_results/contigs/auto_noise/autoNoise1_rescaled.bnx

```
