# Workflow to explore motif SNP in fine mapping result

`motif_analysis/` performs motif scanning using `fimo` on a list of gene candidates (defined in `motif_analysis/scripts/task-name.R` which can be specified in config file) and generate a list of candidate motif binding sites.

`fine_mapping` performs fine mapping on given tissues (specify fine mapping details in `fine_mapping/scripts/`) and output fine mapping result as in `rds` format. Then, `fine_mapping/spotlight.snakemake` can be used to preliminarily visualize the fine mapping result along with motif information.

# Example

```
$ cd motif_analysis
$ snakemake --configfile config.eqtl-skin.yaml all_spotlight
$ cd ../fine_mapping
$ snakemake --configfile config.skin.yaml all_map
$ snakemake -s spotlight.snakemake --configfile config.spotlight-motif_snp.yaml all_spotlight 
```

