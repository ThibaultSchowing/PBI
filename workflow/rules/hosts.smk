"""
Host genome download and processing rules

This Snakefile handles:
1. Downloading host bacterial genomes from NCBI
2. Merging individual genomes into single FASTA
3. Indexing the merged FASTA with pyfaidx
"""

rule download_host_genomes:
    """
    Download host bacterial genomes from NCBI RefSeq
    
    This rule extracts unique host species from the phage metadata CSV and downloads
    reference genomes for each using assembly_resolver.py + download_host_genomes_robust.py.
    
    The robust downloader:
    - Uses NCBI Assembly database as authoritative source
    - Normalizes all identifiers to assembly accessions (GCF_ preferred)
    - Handles ambiguity in species names explicitly
    - Supports metadata-only mode (no downloads)
    - Preserves existing downloads with validation
    - Creates phage-host assembly links
    
    Note: Reads from CSV instead of database to avoid circular dependency
    (database is created after host downloads).
    """
    input:
        phage_csv = config["phage_metadata_merged_output"]
    output:
        metadata = config["host_metadata_output"],
        assembly_metadata = config.get("assembly_metadata_output", 
                                       config["host_metadata_output"].replace('.csv', '_assemblies.csv')),
        phage_host_links = config.get("phage_host_links_output",
                                      config["host_metadata_output"].replace('.csv', '_phage_host_links.csv')),
        phage_host_candidates = config.get("phage_host_candidates_output",
                                           config["host_metadata_output"].replace('.csv', '_host_candidates.csv')),
        phage_host_assemblies = config.get("phage_host_assemblies_output",
                                           config["host_metadata_output"].replace('.csv', '_host_assemblies.csv')),
        host_resolution_cache = config.get("host_resolution_cache_output",
                                           config["host_metadata_output"].replace('.csv', '_token_resolution_cache.json'))
    params:
        output_dir = config["host_genomes_intermediate"],
        limit = None,  # Set to integer for testing with subset of hosts
        metadata_only = config.get("metadata_only_mode", False),
        skip_existing = config.get("skip_existing_downloads", True),
        validate_checksums = config.get("validate_file_checksums", True),
        reuse_resolution_cache = config.get("reuse_host_resolution_cache", True)
    log:
        config["host_download_log"]
    conda:
        "../envs/sequences.yaml"
    script:
        "../scripts/sequences/download_host_genomes_robust.py"


rule create_host_mapping:
    """
    Create mapping file from Host_ID to individual FASTA file paths
    
    Instead of merging all host genomes, this creates a JSON mapping that allows
    on-demand loading of individual host genome files for efficient memory usage.
    """
    input:
        metadata = config["host_metadata_output"],
        private_mapping = config["private_host_mapping"]
    output:
        mapping = config["host_fasta_mapping"]
    params:
        input_dir = config["host_genomes_intermediate"]
    log:
        config["create_host_mapping_log"]
    conda:
        "../envs/sequences.yaml"
    script:
        "../scripts/sequences/create_host_mapping.py"


rule index_individual_host_sequences:
    """
    Create pyfaidx indexes for individual host genome FASTA files

    FASTA quality-control procedure (non-destructive):
    1. Header audit: duplicate sequence identifiers → file rejected, event logged.
    2. Sequence content audit: identical sequences → warning logged, file still indexed.

    Generates .fai index files for all files that pass the header-uniqueness check.
    Produces a DataFrame-loadable CSV QC log for downstream analysis.
    """
    input:
        mapping = config["host_fasta_mapping"]
    output:
        flag    = touch(config["host_index_complete_flag"]),
        qc_log  = config["host_fasta_qc_log"]
    log:
        config["index_individual_host_sequences_log"]
    conda:
        "../envs/sequences.yaml"
    script:
        "../scripts/sequences/index_individual_hosts.py"


rule create_host_status_report:
    """
    Combined per-phage host status report (DataFrame-loadable CSV)

    Joins phage_host_candidates + phage_host_assemblies + assembly_metadata
    + FASTA QC log into one table with one row per (Phage_ID, Host_Token).

    Enables queries like:
    - For all phages, how many have ≥1 resolved host assembly?
    - Of resolved hosts, how many were downloaded / indexed / rejected?
    """
    input:
        candidates        = config.get("phage_host_candidates_output",
                                       config["host_metadata_output"].replace('.csv', '_host_candidates.csv')),
        assemblies        = config.get("phage_host_assemblies_output",
                                       config["host_metadata_output"].replace('.csv', '_host_assemblies.csv')),
        assembly_metadata = config.get("assembly_metadata_output",
                                       config["host_metadata_output"].replace('.csv', '_assemblies.csv')),
        qc_log            = config["host_fasta_qc_log"]
    output:
        status_report = config["host_status_report"]
    log:
        config["create_host_status_report_log"]
    conda:
        "../envs/sequences.yaml"
    script:
        "../scripts/sequences/create_host_status_report.py"


rule all_hosts:
    """
    Target rule for all host genome processing
    """
    input:
        config["host_index_complete_flag"],
        config["host_fasta_qc_log"],
        config["host_status_report"],
        config["host_fasta_mapping"],
        config["host_metadata_output"],
        config.get("phage_host_candidates_output",
                   config["host_metadata_output"].replace('.csv', '_host_candidates.csv')),
        config.get("phage_host_assemblies_output",
                   config["host_metadata_output"].replace('.csv', '_host_assemblies.csv'))
