FEATURES = [
    "phage_metadata", 
    "annotated_proteins_metadata", 
    "phage_trna_tmrna_metadata", 
    "phage_anti_crispr_metadata", 
    "phage_virulent_factor_metadata", 
    "phage_transmembrane_protein_metadata"
]



configfile: "config.yaml"

import os

# expand fonctionne avec des templates de chemins, pas directement avec des clés formatées à évaluer dynamiquement dans le dictionnaire config.

# ----------------------------------------
# ✅ RULE ALL (tous les rapports)
# ----------------------------------------

rule all:
    input:
        expand(
            "data/merged/merged_{feature}.csv",
            feature=FEATURES
        ),
        expand(
            "reports/{feature}_report.html",
            feature=FEATURES
        )

# ----------------------------------------
# ✅ RULE DOWNLOAD (UNIQUE)
#     download_all_tsvs -> lists explicit filenames
#     download_tsv -> downloads each file
# ----------------------------------------
rule download_all_tsvs:
    input:
        sum(
            [
                [
                    f"data/intermediate_csv/{feature}/{source}.tsv"
                    for source in config[f"{feature}_urls"].keys()
                ]
                for feature in FEATURES
            ],
            []
        )


rule download_tsv:
    output:
        "data/intermediate_csv/{feature}/{source}.tsv"
    params: 
        url = lambda wildcards: config[f"{wildcards.feature}_urls"][wildcards.source]
    shell:
        """
        mkdir -p data/intermediate_csv/{wildcards.feature}
        wget -O {output} {params.url}
        """

# ----------------------------------------
# ✅ RULE MERGE PHAGE METADATA
# ----------------------------------------
rule merge_phage_metadata_tsvs:
    input:
        expand(
            "data/intermediate_csv/phage_metadata/{source}.tsv",
            source=list(config["phage_metadata_urls"].keys())
        )
    output:
        config["phage_metadata_merged_output"]
    script:
        "scripts/merge_phage_metadata.py"

# ----------------------------------------
# ✅ RULE MERGE ANNOTATED PROTEINS METADATA
# ----------------------------------------
rule merge_annotated_proteins_metadata_tsvs:
    input:
        expand(
            "data/intermediate_csv/annotated_proteins_metadata/{source}.tsv",
            source=list(config["annotated_proteins_metadata_urls"].keys())
        )
    output:
        config["annotated_proteins_metadata_merged_output"]
    script:
        "scripts/merge_annotated_proteins_metadata.py"

# ----------------------------------------
# ✅ RULE MERGE PHAGE tRNA/tmRNA METADATA
# ----------------------------------------
rule merge_phage_trna_tmrna_metadata_tsvs:
    input:
        expand(
            "data/intermediate_csv/phage_trna_tmrna_metadata/{source}.tsv",
            source=list(config["phage_trna_tmrna_metadata_urls"].keys())
        )
    output:
        config["phage_trna_tmrna_metadata_merged_output"]
    script:
        "scripts/merge_phage_trna_tmrna_metadata.py"

# ----------------------------------------
# ✅ RULE MERGE PHAGE ANTI-CRISPR METADATA
# ----------------------------------------
rule merge_phage_anti_crispr_metadata_tsvs:
    input:
        expand(
            "data/intermediate_csv/phage_anti_crispr_metadata/{source}.tsv",
            source=list(config["phage_anti_crispr_metadata_urls"].keys())
        )
    output:
        config["phage_anti_crispr_metadata_merged_output"]
    script:
        "scripts/merge_phage_anti_crispr_metadata.py"

# ----------------------------------------
# ✅ RULE MERGE PHAGE VIRULENT FACTOR METADATA
# ----------------------------------------
rule merge_phage_virulent_factor_metadata_tsvs:
    input:
        expand(
            "data/intermediate_csv/phage_virulent_factor_metadata/{source}.tsv",
            source=list(config["phage_virulent_factor_metadata_urls"].keys())
        )
    output:
        config["phage_virulent_factor_metadata_merged_output"]
    script:
        "scripts/merge_phage_virulent_factor_metadata.py"

# ----------------------------------------
# ✅ RULE MERGE PHAGE TRANSMEMBRANE PROTEIN METADATA
# ----------------------------------------
rule merge_phage_transmembrane_protein_metadata_tsvs:
    input:
        expand(
            "data/intermediate_csv/phage_transmembrane_protein_metadata/{source}.tsv",
            source=list(config["phage_transmembrane_protein_metadata_urls"].keys())
        )
    output:
        config["phage_transmembrane_protein_metadata_merged_output"]
    script:
        "scripts/merge_phage_transmembrane_protein_metadata.py"

