# PBI
Phages Bacteria Interaction exploration and overview

## First Snakemake pipeline with Pixi

For details check Snakemake 9.8.0 documentation and Pixi documentation. 

For details on multi-environments with Pixi check [here](https://pixi.sh/latest/tutorials/multi_environment/#lets-get-started). 

TODO: show .toml file here

Notes à supprimer: 
- pas oublier de faire pixi install --environment <envname>
- ajouter un gitignore avec les éléments pour pixi et snakemake (done)


## Effectué: 
- snakemake --cores 4 data/merged/merged_annotated_proteins_metadata.csv
    - en option: --dag | dot -Tsvg > dag.svg
- snakemake --cores 4 data/merged/merged_phage_metadata.csv
    - en option: --dag | dot -Tsvg > dag.svg








