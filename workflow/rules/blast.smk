"""
BLAST database build rules

Creates BLAST databases from the merged FASTA files for sequence similarity searches.
Databases are stored under /data/processed/blast_db/ (configurable via blast_db_dir).
"""

BLAST_DB_DIR = config.get("blast_db_dir", "/data/processed/blast_db")


rule makeblastdb_phages:
    """
    Build a nucleotide BLAST database from merged phage genomes.
    
    Uses makeblastdb with -parse_seqids to preserve sequence identifiers
    for result interpretation.
    """
    input:
        fasta = config["all_phages_fasta"]
    output:
        db = touch(os.path.join(BLAST_DB_DIR, "phages", "makeblastdb_phages.done"))
    params:
        db_dir = os.path.join(BLAST_DB_DIR, "phages"),
        db_name = "all_phages"
    log:
        config.get("makeblastdb_phages_log", "/pipeline-logs/logs/makeblastdb_phages.log")
    conda:
        "../envs/blast.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        makeblastdb \
            -in {input.fasta} \
            -dbtype nucl \
            -out {params.db_dir}/{params.db_name} \
            -parse_seqids \
            -title "PBI Phage Genomes" \
            2> {log}
        """


rule makeblastdb_proteins:
    """
    Build a protein BLAST database from merged phage proteins.
    """
    input:
        fasta = config["all_proteins_fasta"]
    output:
        db = touch(os.path.join(BLAST_DB_DIR, "proteins", "makeblastdb_proteins.done"))
    params:
        db_dir = os.path.join(BLAST_DB_DIR, "proteins"),
        db_name = "all_proteins"
    log:
        config.get("makeblastdb_proteins_log", "/pipeline-logs/logs/makeblastdb_proteins.log")
    conda:
        "../envs/blast.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        makeblastdb \
            -in {input.fasta} \
            -dbtype prot \
            -out {params.db_dir}/{params.db_name} \
            -parse_seqids \
            -title "PBI Phage Proteins" \
            2> {log}
        """


rule makeblastdb_hosts:
    """
    Build a nucleotide BLAST database from concatenated host genomes.
    
    Host genomes are stored as individual files. This rule reads the host
    FASTA mapping JSON to locate and concatenate all host genome FASTA files
    into a single file, then builds the BLAST database from it.
    """
    input:
        mapping = config["host_fasta_mapping"]
    output:
        db = touch(os.path.join(BLAST_DB_DIR, "hosts", "makeblastdb_hosts.done"))
    params:
        db_dir = os.path.join(BLAST_DB_DIR, "hosts"),
        db_name = "all_hosts",
        combined_fasta = os.path.join(BLAST_DB_DIR, "hosts", "all_hosts_combined.fasta")
    log:
        config.get("makeblastdb_hosts_log", "/pipeline-logs/logs/makeblastdb_hosts.log")
    conda:
        "../envs/blast.yaml"
    script:
        "../scripts/blast/build_host_blast_db.py"


rule all_blast_dbs:
    """
    Target rule for building all BLAST databases.
    """
    input:
        os.path.join(BLAST_DB_DIR, "phages", "makeblastdb_phages.done"),
        os.path.join(BLAST_DB_DIR, "proteins", "makeblastdb_proteins.done"),
        os.path.join(BLAST_DB_DIR, "hosts", "makeblastdb_hosts.done")
