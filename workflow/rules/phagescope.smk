

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
                    f"{config['intermediate_csv_output']}/{feature}/{source}.tsv"
                    for source in config[f"{feature}_urls"].keys()
                ]
                for feature in FEATURES
            ],
            []
        )


rule download_tsv:
    output:
        tsv=config["intermediate_csv_output"] + "/{feature}/{source}.tsv",
        provenance_sidecar=config["intermediate_csv_output"] + "/{feature}/{source}.provenance.json"
    params: 
        url = lambda wildcards: config[f"{wildcards.feature}_urls"][wildcards.source],
        intermediate_csv_output = config["intermediate_csv_output"]
    threads: 8
    script:
        "../scripts/preprocessing/download_public_file.py"


rule build_public_data_provenance_manifest:
    input:
        sum(
            [
                [
                    f"{config['intermediate_csv_output']}/{feature}/{source}.provenance.json"
                    for source in config[f"{feature}_urls"].keys()
                ]
                for feature in FEATURES
            ],
            []
        )
    output:
        manifest_json=config["public_data_provenance"]["manifest_json_output"],
        manifest_csv=config["public_data_provenance"]["manifest_csv_output"],
        pipeline_run_json=config["public_data_provenance"]["pipeline_run_provenance_json_output"],
        pipeline_run_csv=config["public_data_provenance"]["pipeline_run_provenance_csv_output"]
    script:
        "../scripts/preprocessing/build_public_data_provenance_manifest.py"

# ----------------------------------------
# RULE MERGE TRANSCRIPTION TERMINATOR METADATA
# ----------------------------------------
rule merge_transcription_terminator_metadata_tsvs:
    input:
        expand(
            config["transcription_terminator_metadata_intermediate_output"] + "/{source}.tsv",
            source=list(config["transcription_terminator_metadata_urls"].keys())
        )
    output:
        config["transcription_terminator_metadata_merged_output"]
    params:
        schema_path=str(Path(__file__).resolve().parent.parent / "schemas" / "transcription_terminator_metadata_merged.yaml")
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/preprocessing/mergers/merge_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE METADATA
# ----------------------------------------
rule merge_phage_metadata_tsvs:
    input:
        expand(
            config["phage_metadata_intermediate_output"] + "/{source}.tsv",
            source=list(config["phage_metadata_urls"].keys())
        )
    output:
        config["phage_metadata_merged_output"]
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/preprocessing/mergers/merge_phage_metadata.py"

# ----------------------------------------
# RULE MERGE ANNOTATED PROTEINS METADATA
# ----------------------------------------
rule merge_annotated_proteins_metadata_tsvs:
    input:
        expand(
            config["annotated_proteins_metadata_intermediate_output"] + "/{source}.tsv",
            source=list(config["annotated_proteins_metadata_urls"].keys())
        )
    output:
        config["annotated_proteins_metadata_merged_output"]
    params:
        schema_path=str(Path(__file__).resolve().parent.parent / "schemas" / "annotated_proteins_metadata_merged.yaml")
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/preprocessing/mergers/merge_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE tRNA/tmRNA METADATA
# ----------------------------------------
rule merge_phage_trna_tmrna_metadata_tsvs:
    input:
        expand(
            config["phage_trna_tmrna_metadata_intermediate_output"] + "/{source}.tsv",
            source=list(config["phage_trna_tmrna_metadata_urls"].keys())
        )
    output:
        config["phage_trna_tmrna_metadata_merged_output"]
    params:
        schema_path=str(Path(__file__).resolve().parent.parent / "schemas" / "trna_tmrna_metadata_merged.yaml")
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/preprocessing/mergers/merge_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE ANTI-CRISPR METADATA
# ----------------------------------------
rule merge_phage_anti_crispr_metadata_tsvs:
    input:
        expand(
            config["phage_anti_crispr_metadata_intermediate_output"] + "/{source}.tsv",
            source=list(config["phage_anti_crispr_metadata_urls"].keys())
        )
    output:
        config["phage_anti_crispr_metadata_merged_output"]
    params:
        schema_path=str(Path(__file__).resolve().parent.parent / "schemas" / "anti_crispr_metadata_merged.yaml")
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/preprocessing/mergers/merge_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE VIRULENT FACTOR METADATA
# ----------------------------------------
rule merge_phage_virulent_factor_metadata_tsvs:
    input:
        expand(
            config["phage_virulent_factor_metadata_intermediate_output"] + "/{source}.tsv",
            source=list(config["phage_virulent_factor_metadata_urls"].keys())
        )
    output:
        config["phage_virulent_factor_metadata_merged_output"]
    params:
        schema_path=str(Path(__file__).resolve().parent.parent / "schemas" / "virulent_factor_metadata_merged.yaml")
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/preprocessing/mergers/merge_metadata.py"

# ----------------------------------------
# RULE MERGE PHAGE TRANSMEMBRANE PROTEIN METADATA
# ----------------------------------------
rule merge_phage_transmembrane_protein_metadata_tsvs:
    input:
        expand(
            config["phage_transmembrane_protein_metadata_intermediate_output"] + "/{source}.tsv",
            source=list(config["phage_transmembrane_protein_metadata_urls"].keys())
        )
    output:
        config["phage_transmembrane_protein_metadata_merged_output"]
    params:
        schema_path=str(Path(__file__).resolve().parent.parent / "schemas" / "transmembrane_protein_metadata_merged.yaml")
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/preprocessing/mergers/merge_metadata.py"

# ----------------------------------------
# RULE MERGE ANTIMICROBIAL RESISTANCE GENE METADATA
# ----------------------------------------
rule merge_antimicrobial_resistance_gene_metadata_tsvs:
    input:
        expand(
            config["antimicrobial_resistance_gene_metadata_intermediate_output"] + "/{source}.tsv",
            source=list(config["antimicrobial_resistance_gene_metadata_urls"].keys())
        )
    output:
        config["antimicrobial_resistance_gene_metadata_merged_output"]
    params:
        schema_path=str(Path(__file__).resolve().parent.parent / "schemas" / "antimicrobial_resistance_gene_metadata_merged.yaml")
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/preprocessing/mergers/merge_metadata.py"

# ----------------------------------------
# RULE MERGE CRISPR ARRAY METADATA
# ----------------------------------------
rule merge_crispr_array_metadata_tsvs:
    input:
        expand(
            config["crispr_array_metadata_intermediate_output"] + "/{source}.tsv",
            source=list(config["crispr_array_metadata_urls"].keys())
        )
    output:
        config["crispr_array_metadata_merged_output"]
    params:
        schema_path=str(Path(__file__).resolve().parent.parent / "schemas" / "crispr_array_metadata_merged.yaml")
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/preprocessing/mergers/merge_metadata.py"

rule generate_report:
    input:
        #"../data/merged/merged_{feature}.csv"
        lambda wildcards: config[f"{wildcards.feature}_merged_output"]
    output:
        config["reports_output"] + "{feature}_report.html"
    conda:
        "../envs/reporting.yaml"
    script:
        "../scripts/utils/generate_reports.py"

# Protein fasta files

rule download_protein_fasta:
    """
    Download a .tar.gz archive of protein FASTA files from PhageScope API.
    """
    output:
        os.path.join(config["protein_fasta_compressed_output"], "{dataset}.tar.gz")
    params:
        url=lambda wildcards: config["protein_fasta_urls"][wildcards.dataset]
    cache: True
    shell:
        """
        wget -c -O {output}.tmp {params.url} && mv {output}.tmp {output} || (rm -f {output}.tmp; exit 1)
        """

rule extract_protein_fasta:
    """
    Extract a .tar.gz archive into a dedicated directory per dataset.
    """
    input:
        os.path.join(config["protein_fasta_compressed_output"], "{dataset}.tar.gz")
    output:
        extracted_dir = temp(directory(os.path.join(config["protein_fasta_extracted_output"], "{dataset}")))
    shell:
        """
        mkdir -p {output.extracted_dir}
        tar -xzf {input} -C {output.extracted_dir}
        """

# Phage fasta files

rule download_phage_fasta:
    """
    Download a .tar.gz archive of phage genome FASTA files from PhageScope API.
    """
    output:
        os.path.join(config["phage_fasta_compressed_output"], "{dataset}.tar.gz")
    params:
        url=lambda wildcards: config["phage_fasta_urls"][wildcards.dataset]
    cache: True
    threads: 8
    shell:
        """
        wget -c -O {output}.tmp {params.url} && mv {output}.tmp {output} || (rm -f {output}.tmp; exit 1)
        """

rule extract_phage_fasta:
    """
    Extract a .tar.gz archive into a dedicated directory per dataset.
    """
    input:
        os.path.join(config["phage_fasta_compressed_output"], "{dataset}.tar.gz")
    output:
        extracted_dir = temp(directory(os.path.join(config["phage_fasta_extracted_output"], "{dataset}")))
    shell:
        """
        mkdir -p {output.extracted_dir}
        tar -xzf {input} -C {output.extracted_dir}
        """

rule merge_protein_fasta_by_source:
    input:
        source_dir = os.path.join(config["protein_fasta_extracted_output"], "{dataset}")
    output:
        merged_fasta = os.path.join(config["protein_fasta_merged_output"], "{dataset}.fasta")
    params:
        source_dir = lambda wildcards: os.path.join(config["protein_fasta_extracted_output"], wildcards.dataset),
        merged_fasta_dir = config["protein_fasta_merged_output"]
    shell:
        # If only one fasta is present, just copy and rename. Otherwise, run the Python merge script.
        # This ensures we don’t waste time unnecessarily merging a single file.
        r'''
        mkdir -p {params.merged_fasta_dir}
        
        # Find all fasta files and store in an array
        mapfile -t fasta_files < <(find {params.source_dir} -type f \( -name "*.fasta" -o -name "*.fa" \))
        
        # Check if any files were found
        if [ "${{#fasta_files[@]}}" -eq 0 ]; then
            echo "⚠️ WARNING: No FASTA files found in {params.source_dir} - creating empty file" >&2
            touch {output.merged_fasta}
        elif [ "${{#fasta_files[@]}}" -eq 1 ]; then
            # Only one file, just copy it
            cp "${{fasta_files[0]}}" {output.merged_fasta}
        else
            # Multiple files, merge them
            python scripts/preprocessing/mergers/merge_protein_fasta.py "{params.source_dir}" "{output.merged_fasta}"
        fi
        '''
    
rule merge_phage_fasta_by_source:
    input:
        source_dir = os.path.join(config["phage_fasta_extracted_output"], "{dataset}")
    output:
        merged_fasta = os.path.join(config["phage_fasta_merged_output"], "{dataset}.fasta")
    params:
        source_dir = lambda wildcards: os.path.join(config["phage_fasta_extracted_output"], wildcards.dataset),
        merged_fasta_dir = config["phage_fasta_merged_output"]
    shell:
        # If only one fasta is present, just copy and rename. Otherwise, run the Python merge script.
        # This ensures we don’t waste time unnecessarily merging a single file.
        r'''
        mkdir -p {params.merged_fasta_dir}
        
        # Find all fasta files and store in an array
        mapfile -t fasta_files < <(find {params.source_dir} -type f \( -name "*.fasta" -o -name "*.fa" \))
        
        # Check if any files were found
        if [ "${{#fasta_files[@]}}" -eq 0 ]; then
            echo "⚠️ WARNING: No FASTA files found in {params.source_dir} - creating empty file" >&2
            touch {output.merged_fasta}
        elif [ "${{#fasta_files[@]}}" -eq 1 ]; then
            # Only one file, just copy it
            cp "${{fasta_files[0]}}" {output.merged_fasta}
        else
            # Multiple files, merge them
            python scripts/preprocessing/mergers/merge_phage_fasta.py "{params.source_dir}" "{output.merged_fasta}"
        fi
        '''

rule cleanup_extracted_phage_fasta:
    input:
        flag = os.path.join(config["phage_fasta_extracted_output"], "{dataset}", ".extraction_done")
    shell:
        """
        rm -rf $(dirname {input.flag})
        """

rule cleanup_extracted_protein_fasta:
    input:
        flag = os.path.join(config["protein_fasta_extracted_output"], "{dataset}", ".extraction_done")
    shell:
        """
        rm -rf $(dirname {input.flag})
        """
