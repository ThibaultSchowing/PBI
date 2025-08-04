

# ----------------------------------------
# RULE DOWNLOAD (UNIQUE)
#     download_all_tsvs -> lists explicit filenames
#     download_tsv -> downloads each file
# ----------------------------------------
rule download_all_tsvs:
    input:
        sum(
            [
                [
                    f"../data/intermediate_csv/{feature}/{source}.tsv"
                    for source in config[f"{feature}_urls"].keys()
                ]
                for feature in FEATURES
            ],
            []
        )


rule download_tsv:
    output:
        temp("../data/intermediate_csv/{feature}/{source}.tsv")
    params: 
        url = lambda wildcards: config[f"{wildcards.feature}_urls"][wildcards.source]
    threads: 8
    shell:
        """
        mkdir -p data/intermediate_csv/{wildcards.feature}
        wget -O {output} {params.url} || echo "Failed download for {wildcards.feature}/{wildcards.source}"
        """

# ----------------------------------------
# RULE MERGE TRANSCRIPTION TERMINATOR METADATA
# ----------------------------------------
rule merge_transcription_terminator_metadata_tsvs:
    input:
        expand(
            "../data/intermediate_csv/transcription_terminator_metadata/{source}.tsv",
            source=list(config["transcription_terminator_metadata_urls"].keys())
        )
    output:
        config["transcription_terminator_metadata_merged_output"]
    conda:
        "../envs/pixi_base_env.yaml"
    script:
        "../scripts/merge_transcription_terminator_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE METADATA
# ----------------------------------------
rule merge_phage_metadata_tsvs:
    input:
        expand(
            "../data/intermediate_csv/phage_metadata/{source}.tsv",
            source=list(config["phage_metadata_urls"].keys())
        )
    output:
        config["phage_metadata_merged_output"]
    conda:
        "../envs/pixi_base_env.yaml"
    script:
        "../scripts/merge_phage_metadata.py"

# ----------------------------------------
# RULE MERGE ANNOTATED PROTEINS METADATA
# ----------------------------------------
rule merge_annotated_proteins_metadata_tsvs:
    input:
        expand(
            "../data/intermediate_csv/annotated_proteins_metadata/{source}.tsv",
            source=list(config["annotated_proteins_metadata_urls"].keys())
        )
    output:
        config["annotated_proteins_metadata_merged_output"]
    conda:
        "../envs/pixi_base_env.yaml"
    script:
        "../scripts/merge_annotated_proteins_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE tRNA/tmRNA METADATA
# ----------------------------------------
rule merge_phage_trna_tmrna_metadata_tsvs:
    input:
        expand(
            "../data/intermediate_csv/phage_trna_tmrna_metadata/{source}.tsv",
            source=list(config["phage_trna_tmrna_metadata_urls"].keys())
        )
    output:
        config["phage_trna_tmrna_metadata_merged_output"]
    conda:
        "../envs/pixi_base_env.yaml"
    script:
        "../scripts/merge_phage_trna_tmrna_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE ANTI-CRISPR METADATA
# ----------------------------------------
rule merge_phage_anti_crispr_metadata_tsvs:
    input:
        expand(
            "../data/intermediate_csv/phage_anti_crispr_metadata/{source}.tsv",
            source=list(config["phage_anti_crispr_metadata_urls"].keys())
        )
    output:
        config["phage_anti_crispr_metadata_merged_output"]
    conda:
        "../envs/pixi_base_env.yaml"
    script:
        "../scripts/merge_phage_anti_crispr_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE VIRULENT FACTOR METADATA
# ----------------------------------------
rule merge_phage_virulent_factor_metadata_tsvs:
    input:
        expand(
            "../data/intermediate_csv/phage_virulent_factor_metadata/{source}.tsv",
            source=list(config["phage_virulent_factor_metadata_urls"].keys())
        )
    output:
        config["phage_virulent_factor_metadata_merged_output"]
    conda:
        "../envs/pixi_base_env.yaml"
    
    script:
        "../scripts/merge_phage_virulent_factor_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE TRANSMEMBRANE PROTEIN METADATA
# ----------------------------------------
rule merge_phage_transmembrane_protein_metadata_tsvs:
    input:
        expand(
            "../data/intermediate_csv/phage_transmembrane_protein_metadata/{source}.tsv",
            source=list(config["phage_transmembrane_protein_metadata_urls"].keys()) # e.g. STV_Phage_Metadata_URL
        )
    output:
        config["phage_transmembrane_protein_metadata_merged_output"]
    conda:
        "../envs/pixi_base_env.yaml"
    script:
        "../scripts/merge_phage_transmembrane_protein_metadata.py"


rule generate_report:
    input:
        "../data/merged/merged_{feature}.csv"
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
        temp(os.path.join(compressed_dir, "{dataset}.tar.gz"))
    params:
        url = lambda wildcards: protein_fasta_urls[wildcards.dataset]
    cache: True 
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
        #done_flag = temp(os.path.join(output_protein_fasta_dir, "{dataset}", ".extraction_done"))
        extracted_dir = temp(directory(os.path.join(output_protein_fasta_dir, "{dataset}")))
    shell:
        """
        mkdir -p {output.extracted_dir}
        tar -xzf {input} -C {output.extracted_dir}
        
        """

        #mkdir -p $(dirname {output.done_flag})
        #tar -xzf {input} -C $(dirname {output.done_flag})
        #touch {output.done_flag}
        

# Phage fasta files

rule download_phage_fasta:
    """
    Télécharge un fichier .tar.gz à partir de son URL.
    L'input est dynamique et dépend du nom.
    """
    output:
        temp(os.path.join(compressed_phage_dir, "{dataset}.tar.gz"))
    params:
        url = lambda wildcards: phage_fasta_urls[wildcards.dataset]
    cache: True 
    threads: 8
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
        #done_flag = temp(os.path.join(output_phage_fasta_dir, "{dataset}", ".extraction_done"))
        extracted_dir = temp(directory(os.path.join(output_phage_fasta_dir, "{dataset}")))
    shell:
        """
        mkdir -p {output.extracted_dir}
        tar -xzf {input} -C {output.extracted_dir}
        
        """
        #mkdir -p $(dirname {output.done_flag})
        #tar -xzf {input} -C $(dirname {output.done_flag})
        #touch {output.done_flag}

rule merge_protein_fasta_by_source:
    input:
        #done_flag = os.path.join(output_protein_fasta_dir, "{dataset}", ".extraction_done"),
        source_dir = os.path.join(output_protein_fasta_dir, "{dataset}")
    output:
        merged_fasta = os.path.join("../data/protein_fasta_merged", "{dataset}.fasta")
    params:
        source_dir = lambda wildcards: os.path.join(output_protein_fasta_dir, wildcards.dataset),
        #dataset = lambda wildcards: wildcards.dataset
    shell:
        # If only one fasta is present, just copy and rename. Otherwise, run the Python merge script.
        # This ensures we don’t waste time unnecessarily merging a single file.
        r'''
        mkdir -p data/protein_fasta_merged
        fasta_files=("$(find {params.source_dir} -type f \( -name "*.fasta" -o -name "*.fa" \))")
        if [ $(echo "$fasta_files" | wc -l) -eq 1 ]; then
            cp "$fasta_files" {output.merged_fasta}
        else
            pixi run -e base python scripts/merge_protein_fasta.py "{params.source_dir}" "{output.merged_fasta}"
        fi
        '''
    
rule merge_phage_fasta_by_source:
    input:
        #done_flag = os.path.join(output_phage_fasta_dir, "{dataset}", ".extraction_done"),
        source_dir = os.path.join(output_phage_fasta_dir, "{dataset}")
    output:
        merged_fasta = os.path.join("../data/phage_fasta_merged", "{dataset}.fasta")
    params:
        source_dir = lambda wildcards: os.path.join(output_phage_fasta_dir, wildcards.dataset),
        #dataset = lambda wildcards: wildcards.dataset
    shell:
        # If only one fasta is present, just copy and rename. Otherwise, run the Python merge script.
        # This ensures we don’t waste time unnecessarily merging a single file.
        r'''
        mkdir -p data/phage_fasta_merged
        fasta_files=("$(find {params.source_dir} -type f \( -name "*.fasta" -o -name "*.fa" \))")
        if [ $(echo "$fasta_files" | wc -l) -eq 1 ]; then
            cp "$fasta_files" {output.merged_fasta}
        else
            pixi run -e base python scripts/merge_phage_fasta.py "{params.source_dir}" "{output.merged_fasta}"
        fi
        '''

rule cleanup_extracted_phage_fasta:
    input:
        flag = os.path.join(output_phage_fasta_dir, "{dataset}", ".extraction_done")
    shell:
        """
        rm -rf $(dirname {input.flag})
        """

rule cleanup_extracted_protein_fasta:
    input:
        flag = os.path.join(output_protein_fasta_dir, "{dataset}", ".extraction_done")
    shell:
        """
        rm -rf $(dirname {input.flag})
        """