
configfile: "config.yaml"

FEATURES = [
    "phage_metadata", 
    "annotated_proteins_metadata", 
    "transcription_terminator_metadata",
    "phage_trna_tmrna_metadata", 
    "phage_anti_crispr_metadata", 
    "phage_virulent_factor_metadata", 
    "phage_transmembrane_protein_metadata"
]

# Récupération des paramètres de config
protein_fasta_urls = config["protein_fasta_urls"]
compressed_dir = config["protein_fasta_compressed_output"]
output_protein_fasta_dir = config["protein_fasta_output"]


# Génération de la liste des dossiers extraits attendus pour les fichiers fasta (protéines)
extracted_dirs_prot = [
    os.path.join(output_protein_fasta_dir, name)
    for name in protein_fasta_urls.keys()
]

# Récupération des paramètres de config pour les fichiers fasta de phages
phage_fasta_urls = config["phage_fasta_urls"]
compressed_phage_dir = config["phage_fasta_compressed_output"]
output_phage_fasta_dir = config["phage_fasta_output"]

# Génération de la liste des dossiers extraits attendus pour les fichiers fasta de phages
extracted_dirs_phages = [
    os.path.join(output_phage_fasta_dir, name)
    for name in phage_fasta_urls.keys()
]

#print(f"Extracted directories: {extracted_dirs}") # Remove for DAG generation

# ['data/protein_fasta/Genbank', 'data/protein_fasta/RefSeq', 'data/protein_fasta/DDBJ', 'data/protein_fasta/EMBL', 
# 'data/protein_fasta/PhagesDB', 'data/protein_fasta/GPD', 'data/protein_fasta/GVD', 'data/protein_fasta/MGV', 'data/protein_fasta/TemPhD',
# 'data/protein_fasta/CHVD', 'data/protein_fasta/IGVD', 'data/protein_fasta/GOV2', 'data/protein_fasta/STV']



import os

# expand fonctionne avec des templates de chemins, pas directement avec des clés formatées à évaluer dynamiquement dans le dictionnaire config.

# ----------------------------------------
# ✅ RULE ALL (tous les rapports)
# Demander de générer les csv merged est redondant
# ----------------------------------------

rule all:
    input:
        expand(
            "reports/{feature}_report.html",
            feature=FEATURES
        ),
        extracted_dirs_prot,
        extracted_dirs_phages

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
        wget -O {output} {params.url} || echo "Failed download for {wildcards.feature}/{wildcards.source}"
        """

# ----------------------------------------
# ✅ RULE MERGE TRANSCRIPTION TERMINATOR METADATA
# ----------------------------------------
rule merge_transcription_terminator_metadata_tsvs:
    input:
        expand(
            "data/intermediate_csv/transcription_terminator_metadata/{source}.tsv",
            source=list(config["transcription_terminator_metadata_urls"].keys())
        )
    output:
        config["transcription_terminator_metadata_merged_output"]
    script:
        "scripts/merge_transcription_terminator_metadata.py"

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
            source=list(config["phage_transmembrane_protein_metadata_urls"].keys()) # e.g. STV_Phage_Metadata_URL
        )
    output:
        config["phage_transmembrane_protein_metadata_merged_output"]
    script:
        "scripts/merge_phage_transmembrane_protein_metadata.py"

rule generate_report:
    input:
        "data/merged/merged_{feature}.csv"
    output:
        "reports/{feature}_report.html"
    shell:
        """
        pixi run -e reporting python scripts/generate_reports.py {input} {output}
        """

# Protein fasta files

rule download_protein_fasta:
    """
    Télécharge un fichier .tar.gz à partir de son URL.
    L'input est dynamique et dépend du nom.
    """
    output:
        os.path.join(compressed_dir, "{dataset}.tar.gz")
    params:
        url = lambda wildcards: protein_fasta_urls[wildcards.dataset]
    shell:
        """
        wget -O {output} {params.url}
        """

rule extract_protein_fasta:
    """
    Extrait le contenu d'une archive .tar.gz dans un dossier dédié par dataset.
    Dépend du téléchargement de l'archive correspondante.
    """
    input:
        os.path.join(compressed_dir, "{dataset}.tar.gz")
    output:
        directory(os.path.join(output_protein_fasta_dir, "{dataset}"))
    shell:
        """
        mkdir -p {output}
        tar -xzf {input} -C {output}
        """

# Phage fasta files

rule download_phage_fasta:
    """
    Télécharge un fichier .tar.gz à partir de son URL.
    L'input est dynamique et dépend du nom.
    """
    output:
        os.path.join(compressed_phage_dir, "{dataset}.tar.gz")
    params:
        url = lambda wildcards: phage_fasta_urls[wildcards.dataset]
    shell:
        """
        wget -O {output} {params.url}
        """

rule extract_phage_fasta:
    """
    Extrait le contenu d'une archive .tar.gz dans un dossier dédié par dataset.
    Dépend du téléchargement de l'archive correspondante.
    """
    input:
        os.path.join(compressed_phage_dir, "{dataset}.tar.gz")
    output:
        directory(os.path.join(output_phage_fasta_dir, "{dataset}"))
    shell:
        """
        mkdir -p {output}
        tar -xzf {input} -C {output}
        """

