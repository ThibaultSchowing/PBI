# Rules are executed from the snakefile in the workflow directory and the paths are relative to this directory
# so, ../data/ points to the data directory at the root of the project

from pathlib import Path


def private_root_dependency(wildcards):
    root = config.get("private_data_root", "")
    if not root:
        return []
    root_path = Path(root)
    return [str(root_path)] if root_path.exists() else []


def private_blast_db_dependency(wildcards):
    """Return the private BLAST DB done marker if it exists, else empty list.

    This makes the dependency optional so the pipeline works even when
    no private data or BLAST DB is present (e.g. in CI).
    """
    blast_db_dir = config.get("blast_db_dir", "/data/processed/blast_db")
    done_marker = Path(blast_db_dir) / "private" / "makeblastdb_private.done"
    return [str(done_marker)] if done_marker.exists() else []


rule create_duckdb:
    input:
        phage_data=config["phage_metadata_merged_output"], 
        protein_data=config["annotated_proteins_metadata_merged_output"], 
        terminator_data=config["transcription_terminator_metadata_merged_output"],
        anti_crispr_data=config["phage_anti_crispr_metadata_merged_output"], 
        virulent_factor_data=config["phage_virulent_factor_metadata_merged_output"],
        transmembrane_data=config["phage_transmembrane_protein_metadata_merged_output"],
        trna_tmrna_data=config["phage_trna_tmrna_metadata_merged_output"],
        antimicrobial_resistance_data=config["antimicrobial_resistance_gene_metadata_merged_output"],
        crispr_array_data=config["crispr_array_metadata_merged_output"],
        # Host metadata files (ensure database is created after host downloads)
        host_metadata=config["host_metadata_output"],
        assembly_metadata=config["assembly_metadata_output"],
        phage_host_links=config["phage_host_links_output"],
        private_manifest=config["private_manifest_output"],
        public_data_manifest=config["public_data_provenance"]["manifest_csv_output"],
        pipeline_run_provenance=config["public_data_provenance"]["pipeline_run_provenance_csv_output"],
        # Private BLAST DB for duplicate detection (optional, non-blocking)
        private_blast_db=private_blast_db_dependency,
    output:
        db=config["duckdb_output"]
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/database/create_duckdb.py"   

rule prepare_private_sources:
    input:
        # list of directories in private_data_root, [] if empty
        private_root=private_root_dependency 
    output:
        manifest=config["private_manifest_output"]
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/database/prepare_private_sources.py"

rule optimize_database:
    input:
        db=config["duckdb_output"]
    output:
        optimized_db=config["optimized_duckdb_output"]
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/database/optimize_db.py"

rule validate_database:
    input:
        db=config["optimized_duckdb_output"],
        per_source_protein_counts=config["per_source_protein_counts"]
    output:
        report=config["database_validation_report_output"]  
    conda:
        "../envs/base_env.yaml"
    script:
        "../scripts/database/validate_db.py"
