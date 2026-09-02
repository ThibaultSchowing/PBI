"""
BLAST database build rules

Creates BLAST databases from the merged FASTA files for sequence similarity searches.
Databases are stored under /data/processed/blast_db/ (configurable via blast_db_dir).

Note on large files:
  PhageScope contains millions of phage sequences. makeblastdb handles this by
  splitting into volumes. -max_file_size controls the maximum volume size (default
  1 GiB). Increase via config['blast_max_file_size'] if needed.
"""

BLAST_DB_DIR = config.get("blast_db_dir", "/data/processed/blast_db")
BLAST_MAX_FILE_SIZE = config.get("blast_max_file_size", "1G")


rule makeblastdb_phages:
    """
    Build a nucleotide BLAST database from merged phage genomes.

    PhageScope headers are not NCBI-formatted, so -parse_seqids is omitted.
    """
    input:
        fasta = config["all_phages_fasta"]
    output:
        db = touch(os.path.join(BLAST_DB_DIR, "phages", "makeblastdb_phages.done"))
    params:
        db_dir = os.path.join(BLAST_DB_DIR, "phages"),
        db_name = "all_phages",
        max_file_size = BLAST_MAX_FILE_SIZE
    log:
        config.get("makeblastdb_phages_log", "/pipeline-logs/logs/makeblastdb_phages.log")
    conda:
        "../envs/blast.yaml"
    shell:
        """
        mkdir -p {params.db_dir} $(dirname {log})
        makeblastdb \
            -in {input.fasta} \
            -dbtype nucl \
            -out {params.db_dir}/{params.db_name} \
            -max_file_sz {params.max_file_size} \
            > {log} 2>&1
        """


rule makeblastdb_proteins:
    """
    Build a protein BLAST database from merged phage proteins.

    PhageScope headers are not NCBI-formatted, so -parse_seqids is omitted.
    """
    input:
        fasta = config["all_proteins_fasta"]
    output:
        db = touch(os.path.join(BLAST_DB_DIR, "proteins", "makeblastdb_proteins.done"))
    params:
        db_dir = os.path.join(BLAST_DB_DIR, "proteins"),
        db_name = "all_proteins",
        max_file_size = BLAST_MAX_FILE_SIZE
    log:
        config.get("makeblastdb_proteins_log", "/pipeline-logs/logs/makeblastdb_proteins.log")
    conda:
        "../envs/blast.yaml"
    shell:
        """
        mkdir -p {params.db_dir} $(dirname {log})
        makeblastdb \
            -in {input.fasta} \
            -dbtype prot \
            -out {params.db_dir}/{params.db_name} \
            -max_file_sz {params.max_file_size} \
            > {log} 2>&1
        """


rule makeblastdb_hosts:
    """
    Build a nucleotide BLAST database from concatenated host genomes.
    
    Host genomes are stored as individual files. This rule reads the host
    FASTA mapping JSON to locate and concatenate all host genome FASTA files
    into a single file, then builds the BLAST database from it.
    
    Hosts are from NCBI RefSeq, so -parse_seqids is safe.
    """
    input:
        mapping = config["host_fasta_mapping"]
    output:
        db = touch(os.path.join(BLAST_DB_DIR, "hosts", "makeblastdb_hosts.done"))
    params:
        db_dir = os.path.join(BLAST_DB_DIR, "hosts"),
        db_name = "all_hosts",
        combined_fasta = os.path.join(BLAST_DB_DIR, "hosts", "all_hosts_combined.fasta"),
        max_file_size = BLAST_MAX_FILE_SIZE
    log:
        config.get("makeblastdb_hosts_log", "/pipeline-logs/logs/makeblastdb_hosts.log")
    conda:
        "../envs/blast.yaml"
    script:
        "../scripts/blast/build_host_blast_db.py"


rule makeblastdb_private:
    """
    Build a nucleotide BLAST database from private phage sequences.

    Collects all per-source phage.fasta files produced by prepare_private_sequences
    and builds a BLAST DB from them. This DB is used for duplicate detection
    against the public data and for user queries via the API/CLI.
    """
    input:
        private_phage_mapping = config["private_phage_mapping"]
    output:
        db = touch(os.path.join(BLAST_DB_DIR, "private", "makeblastdb_private.done"))
    params:
        private_phage_dir = config["private_phage_genomes_intermediate"],
        db_dir = os.path.join(BLAST_DB_DIR, "private"),
        db_name = "all_private",
        max_file_size = BLAST_MAX_FILE_SIZE
    log:
        config.get("makeblastdb_private_log", "/pipeline-logs/logs/makeblastdb_private.log")
    conda:
        "../envs/blast.yaml"
    script:
        "../scripts/blast/build_private_blast_db.py"


rule makeblastdb_combined:
    """
    Build a combined nucleotide BLAST database from public + private phage sequences.

    Merges all_phages.fasta (public) with all per-source private phage.fasta files
    into a single BLAST DB. This enables users to search against everything at once
    via the API/CLI.
    """
    input:
        public_fasta = config["all_phages_fasta"],
        private_phage_mapping = config["private_phage_mapping"]
    output:
        db = touch(os.path.join(BLAST_DB_DIR, "combined", "makeblastdb_combined.done"))
    params:
        private_phage_dir = config["private_phage_genomes_intermediate"],
        db_dir = os.path.join(BLAST_DB_DIR, "combined"),
        db_name = "all_combined",
        max_file_size = BLAST_MAX_FILE_SIZE
    log:
        config.get("makeblastdb_combined_log", "/pipeline-logs/logs/makeblastdb_combined.log")
    conda:
        "../envs/blast.yaml"
    script:
        "../scripts/blast/build_combined_blast_db.py"


rule check_private_duplicates:
    """
    Search private phage sequences against the public BLAST DB to detect duplicates.

    Identifies private phages with >99% identity to public data, which may
    indicate that the private data came from an old public database release
    with a modified ID.
    """
    input:
        private_phage_mapping = config["private_phage_mapping"],
        blast_db_prefix = os.path.join(BLAST_DB_DIR, "phages", "all_phages")
    output:
        report = config.get("private_duplicate_report", "/private-data/private_duplicate_report.json")
    params:
        pident_threshold = float(config.get("duplicate_pident_threshold", 99.0)),
        qcovs_threshold = float(config.get("duplicate_qcovs_threshold", 90.0))
    log:
        config.get("check_private_duplicates_log", "/pipeline-logs/logs/check_private_duplicates.log")
    conda:
        "../envs/blast.yaml"
    script:
        "../scripts/blast/check_duplicates.py"


rule all_blast_dbs:
    """
    Target rule for building all BLAST databases.
    """
    input:
        os.path.join(BLAST_DB_DIR, "phages", "makeblastdb_phages.done"),
        os.path.join(BLAST_DB_DIR, "proteins", "makeblastdb_proteins.done"),
        os.path.join(BLAST_DB_DIR, "hosts", "makeblastdb_hosts.done"),
        os.path.join(BLAST_DB_DIR, "private", "makeblastdb_private.done"),
        os.path.join(BLAST_DB_DIR, "combined", "makeblastdb_combined.done")
