# PBI

Phages Bacteria Interaction exploration and overview

## Stage 1 - PhageScope

Chosen for its large number of source database, [PhageScope](https://phagescope.deepomics.org/) contains phage and protein sequences alongside multiple metadata such as virulence factor, protein annotation, anti-crispr, transmembrane protein, tRNA - mRNA or transcription terminator. Phages metadata contain the **Host species** of the phage, necessary information for further work on prediction of these interactions. Moreover, sequences are clustered 

All the data from PhageScope come from a set of 13 databases: 
- RefSeq
- Genbank
- EMBL
- DDBJ
- PhagesDB
- GVD
- GPD
- MGV
- TemPhD
- CHVD
- IGVD
- IMGVR
- GOV2
- STV

![Image from phagescope.deepomics.org](https://phagescope.deepomics.org/png/databasevis.png)
![Image from phagescope.deepomics.org](https://phagescope.deepomics.org/png/analysisvis.png)
![Image from phagescope.deepomics.org](https://phagescope.deepomics.org/png/visualization.png)

### Download  Snakemake pipeline with Pixi

For details check Snakemake 9.8.0 documentation and Pixi documentation. 

For details on multi-environments with Pixi check [here](https://pixi.sh/latest/tutorials/multi_environment/#lets-get-started). 

Snakemake is run from pixi, with pixi run



TODO: show .toml file here

### To Be Moved : 

- Don't forget to run  `pixi install --environment <envname>`
- Add Pixi and Snakemake files to .gitignore


### How to 

Snakemake is launched from Pixi with `pixi run`. 

The Pixi environment needs to be exported first with `pixi workspace export conda-environment -e base envs/pixi_base_enf.yaml`

The Conda environment to use is specified within each rule (if needed) with 

```
    conda:
        "envs/pixi_base_env.yaml"
```

Cache: to use caching, it is first needed to export snakemake cache with `export SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache/` (create the destination directory first). 

- pixi run snakemake all --cache --use-conda --conda-frontend mamba --printshellcmds --cores 2
    - optionally to generate the DAT you can add: --dag | dot -Tsvg > dag.svg


- TODO: command to remove temporary objects








