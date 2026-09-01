#!/usr/bin/env python

import duckdb
import pandas as pd
import os
import logging
import sys
import csv
import json
from pathlib import Path

logging.basicConfig(level=logging.INFO)

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_PATH = REPO_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from pbi.private_data import ingest_private_sources_into_db


def _table_exists(conn, table_name: str) -> bool:
    return bool(
        conn.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = ? LIMIT 1",
            [table_name],
        ).fetchone()
    )


def _csv_header_columns(path: str) -> set[str]:
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return set()
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, [])
        return {str(name).strip() for name in header if str(name).strip()}
    except Exception:
        return set()


def _create_dataset_provenance_table(conn, manifest_csv_path: str):
    if manifest_csv_path and os.path.exists(manifest_csv_path) and os.path.getsize(manifest_csv_path) > 0:
        conn.execute("""
        CREATE TABLE dataset_provenance AS
        SELECT *
        FROM read_csv(?,
                      header=true,
                      all_varchar=true,
                      ignore_errors=true,
                      null_padding=true)
        """, [manifest_csv_path])
    else:
        conn.execute("""
        CREATE TABLE dataset_provenance (
            provider_name VARCHAR,
            provider_release VARCHAR,
            provider_snapshot_date VARCHAR,
            provider_schema_profile VARCHAR,
            feature VARCHAR,
            source_key VARCHAR,
            normalized_source_db VARCHAR,
            source_url VARCHAR,
            local_path VARCHAR,
            retrieved_at VARCHAR,
            file_size VARCHAR,
            sha256 VARCHAR,
            etag VARCHAR,
            last_modified VARCHAR,
            detected_tabular_columns VARCHAR,
            schema_fingerprint VARCHAR,
            status VARCHAR,
            error_message VARCHAR
        )
        """)


def _create_pipeline_run_provenance_table(conn, run_csv_path: str):
    if run_csv_path and os.path.exists(run_csv_path) and os.path.getsize(run_csv_path) > 0:
        conn.execute("""
        CREATE TABLE pipeline_run_provenance AS
        SELECT *
        FROM read_csv(?,
                      header=true,
                      all_varchar=true,
                      ignore_errors=true,
                      null_padding=true)
        """, [run_csv_path])
    else:
        conn.execute("""
        CREATE TABLE pipeline_run_provenance (
            pipeline_run_timestamp VARCHAR,
            provider_name VARCHAR,
            provider_release VARCHAR,
            provider_snapshot_date VARCHAR,
            provider_schema_profile VARCHAR,
            provider_api_base_url VARCHAR,
            provider_provenance_mode VARCHAR,
            pbi_version VARCHAR,
            git_commit VARCHAR,
            download_records_count VARCHAR
        )
        """)

def create_star_schema_duckdb():
    """Create DuckDB with star schema from PhageScope data"""
    
    # Get inputs and outputs from Snakemake
    phage_data = snakemake.input.phage_data
    protein_data = snakemake.input.protein_data
    terminator_data = snakemake.input.terminator_data
    anti_crispr_data = snakemake.input.anti_crispr_data
    virulent_factor_data = snakemake.input.virulent_factor_data
    transmembrane_data = snakemake.input.transmembrane_data
    trna_tmrna_data = snakemake.input.trna_tmrna_data
    antimicrobial_resistance_data = snakemake.input.antimicrobial_resistance_data
    crispr_array_data = snakemake.input.crispr_array_data
    db_path = snakemake.output.db
    
    # Host metadata files (now provided as inputs from Snakemake rule)
    host_metadata_path = snakemake.input.get('host_metadata', None)
    assembly_metadata_path = snakemake.input.get('assembly_metadata', None)
    phage_host_links_path = snakemake.input.get('phage_host_links', None)
    private_manifest_path = snakemake.input.get('private_manifest', None)
    provenance_cfg = snakemake.config.get("public_data_provenance", {})
    dataset_provenance_manifest_csv = snakemake.input.get(
        "public_data_manifest",
        provenance_cfg.get("manifest_csv_output", ""),
    )
    pipeline_run_provenance_csv = snakemake.input.get(
        "pipeline_run_provenance",
        provenance_cfg.get("pipeline_run_provenance_csv_output", ""),
    )
    host_count = 0
    dataset_provenance_count = 0
    pipeline_run_provenance_count = 0
    
    # Create output directory
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    
    logging.info(f"Creating DuckDB database at {db_path}")
    
    conn = duckdb.connect(db_path)
    
    # 1. CREATE FACT TABLE - PHAGES
    logging.info("Creating fact_phages table")
    phage_columns = _csv_header_columns(phage_data)
    optional_fact_columns = [
        "Provider_Name",
        "Provider_Release",
        "Provider_Snapshot_Date",
        "Provider_Schema_Profile",
        "Input_Source_Key",
        "Input_File",
        "Input_Retrieved_At",
    ]

    optional_select_parts = []
    for column_name in optional_fact_columns:
        if column_name in phage_columns:
            optional_select_parts.append(column_name)
        else:
            optional_select_parts.append(f"NULL::VARCHAR AS {column_name}")

    optional_select_sql = ",\n        ".join(optional_select_parts)
    conn.execute(f"""
    CREATE TABLE fact_phages AS
    SELECT
        Phage_ID,
        Source_DB,
        TRY_CAST(NULLIF(Length, '-') AS INTEGER) as Length,
        TRY_CAST(NULLIF(GC_content, '-') AS DOUBLE) as GC_content,
        Taxonomy,
        Completeness,
        Host,
        Lifestyle,
        Cluster,
        Subcluster,
        {optional_select_sql},
        'public' as source_type,
        FALSE as is_duplicate_of_public,
        NULL::VARCHAR as duplicate_public_id,
        NULL::DOUBLE as duplicate_pident,
        NULL::DOUBLE as duplicate_qcovs
    FROM read_csv(?,
                  header=true,
                  all_varchar=true,
                  ignore_errors=true,
                  null_padding=true)
    WHERE Phage_ID IS NOT NULL
    """, [phage_data])
    
    phage_count = conn.execute("SELECT COUNT(*) FROM fact_phages").fetchone()[0]
    logging.info(f"✅ Created fact_phages: {phage_count:,} rows")
    
    # 2. CREATE DIM_PROTEINS TABLE
    logging.info("Creating dim_proteins table")
    conn.execute(f"""
    CREATE TABLE dim_proteins AS
    SELECT 
        Phage_ID,
        Protein_ID,
        Protein_source, 
        Function_prediction_source, 
        TRY_CAST(NULLIF(Start, '-') AS INTEGER) as Start,
        TRY_CAST(NULLIF(Stop, '-') AS INTEGER) as Stop,
        Strand, 
        Product, 
        Protein_classification,
        TRY_CAST(NULLIF(Molecular_weight, '-') AS DOUBLE) as Molecular_weight,
        TRY_CAST(NULLIF(Aromaticity, '-') AS DOUBLE) as Aromaticity,
        TRY_CAST(NULLIF(Instability_index, '-') AS DOUBLE) as Instability_index,
        TRY_CAST(NULLIF(Isoelectric_point, '-') AS DOUBLE) as Isoelectric_point,
        TRY_CAST(NULLIF(Helix_fraction, '-') AS DOUBLE) as Helix_fraction,
        TRY_CAST(NULLIF(Turn_fraction, '-') AS DOUBLE) as Turn_fraction,
        TRY_CAST(NULLIF(Sheet_fraction, '-') AS DOUBLE) as Sheet_fraction,
        TRY_CAST(NULLIF(Reduced_coefficient, '-') AS DOUBLE) as Reduced_coefficient,
        TRY_CAST(NULLIF(Oxidized_coefficient, '-') AS DOUBLE) as Oxidized_coefficient,
        Source_DB,
        'phagescope' as source_type
    FROM read_csv('{protein_data}', 
                  header=true, 
                  all_varchar=true, 
                  ignore_errors=true,
                  null_padding=true)
    WHERE Protein_ID IS NOT NULL AND Phage_ID IS NOT NULL
    """)
    
    protein_count = conn.execute("SELECT COUNT(*) FROM dim_proteins").fetchone()[0]
    logging.info(f"✅ Created dim_proteins: {protein_count:,} rows")
    
    # 3. CREATE DIM_TERMINATORS TABLE
    logging.info("Creating dim_terminators table") 
    conn.execute(f"""
    CREATE TABLE dim_terminators AS
    SELECT 
        Phage_ID,
        Terminator as terminator_type,
        TRY_CAST(NULLIF(Start, '-') AS INTEGER) as terminator_start,
        TRY_CAST(NULLIF(Stop, '-') AS INTEGER) as terminator_end,
        TRY_CAST(NULLIF(Confidence, '-') AS DOUBLE) as confidence_score,
        Sense,
        Loc,
        Source_DB,
        'phagescope' as source_type
    FROM read_csv('{terminator_data}', 
                  header=true, 
                  all_varchar=true, 
                  ignore_errors=true,
                  null_padding=true)
    WHERE Phage_ID IS NOT NULL
    """)
    
    terminator_count = conn.execute("SELECT COUNT(*) FROM dim_terminators").fetchone()[0]
    logging.info(f"✅ Created dim_terminators: {terminator_count:,} rows")

    # 4. CREATE DIM_ANTI_CRISPR TABLE
    logging.info("Creating dim_anti_crispr table")
    conn.execute(f"""
    CREATE TABLE dim_anti_crispr AS 
    SELECT 
        Phage_ID,
        Protein_ID,
        Source,
        Source_DB,
        'phagescope' as source_type
    FROM read_csv_auto('{anti_crispr_data}')
    WHERE Phage_ID IS NOT NULL
    """)
    
    anti_crispr_count = conn.execute("SELECT COUNT(*) FROM dim_anti_crispr").fetchone()[0]
    logging.info(f"✅ Created dim_anti_crispr: {anti_crispr_count:,} rows")
    
    # 5. CREATE DIM_VIRULENT_FACTORS TABLE
    # Columns: Protein_ID, Aligned_Protein_in_VFDB, Phage_ID, Source_DB
    logging.info("Creating dim_virulent_factors table")
    conn.execute(f"""
    CREATE TABLE dim_virulent_factors AS 
    SELECT 
        Phage_ID,
        Protein_ID,
        Aligned_Protein_in_VFDB as aligned_protein_vfdb,
        Source_DB,
        'phagescope' as source_type
    FROM read_csv_auto('{virulent_factor_data}')
    WHERE Phage_ID IS NOT NULL
    """)
    
    virulent_count = conn.execute("SELECT COUNT(*) FROM dim_virulent_factors").fetchone()[0]
    logging.info(f"✅ Created dim_virulent_factors: {virulent_count:,} rows")
    
    # 6. CREATE DIM_TRANSMEMBRANE_PROTEINS TABLE
    # Columns: Phage_ID, Protein_ID, Length, PredictedTMHsNumber, ExpnumberofAAsinTMHs, 
    #          Expnumberfirst60AAs, TotalprobofNin, POSSIBLENterm, Insidesource, Insidestart, 
    #          Insideend, TMhelixsource, TMhelixstart, TMhelixend, Outsidesource, 
    #          Outsidestart, Outsideend, Source_DB
    logging.info("Creating dim_transmembrane_proteins table")
    conn.execute(f"""
    CREATE TABLE dim_transmembrane_proteins AS 
    SELECT 
        Phage_ID,
        Protein_ID,
        TRY_CAST(NULLIF(Length, '-') AS INTEGER) as protein_length,
        TRY_CAST(NULLIF(PredictedTMHsNumber, '-') AS INTEGER) as predicted_tmhs_number,
        TRY_CAST(NULLIF(ExpnumberofAAsinTMHs, '-') AS DOUBLE) as exp_aa_in_tmhs,
        TRY_CAST(NULLIF(Expnumberfirst60AAs, '-') AS DOUBLE) as exp_first_60_aa,
        TRY_CAST(NULLIF(TotalprobofNin, '-') AS DOUBLE) as total_prob_n_in,
        POSSIBLENterm as possible_n_term,
        Insidesource as inside_source,
        TRY_CAST(NULLIF(Insidestart, '-') AS INTEGER) as inside_start,
        TRY_CAST(NULLIF(Insideend, '-') AS INTEGER) as inside_end,
        TMhelixsource as tm_helix_source,
        TRY_CAST(NULLIF(TMhelixstart, '-') AS INTEGER) as tm_helix_start,
        TRY_CAST(NULLIF(TMhelixend, '-') AS INTEGER) as tm_helix_end,
        Outsidesource as outside_source,
        TRY_CAST(NULLIF(Outsidestart, '-') AS INTEGER) as outside_start,
        TRY_CAST(NULLIF(Outsideend, '-') AS INTEGER) as outside_end,
        Source_DB,
        'phagescope' as source_type
    FROM read_csv('{transmembrane_data}',
                  header=true,
                  all_varchar=true,
                  ignore_errors=true,
                  null_padding=true)
    WHERE Phage_ID IS NOT NULL
    """)
    
    transmembrane_count = conn.execute("SELECT COUNT(*) FROM dim_transmembrane_proteins").fetchone()[0]
    logging.info(f"✅ Created dim_transmembrane_proteins: {transmembrane_count:,} rows")
    
    # 7. CREATE DIM_TRNA_TMRNA TABLE
    # Columns: t(m)RNA_ID, Source, t(m)RNA, Start, Stop, Strand, Length, Permuted, 
    #          Sequence, Phage_ID, Source_DB
    logging.info("Creating dim_trna_tmrna table")
    conn.execute(f"""
    CREATE TABLE dim_trna_tmrna AS 
    SELECT 
        Phage_ID,
        "t(m)RNA_ID" as trna_tmrna_id,
        Source as source,
        "t(m)RNA" as trna_type,
        TRY_CAST(NULLIF(Start, '-') AS INTEGER) as start_pos,
        TRY_CAST(NULLIF(Stop, '-') AS INTEGER) as stop_pos,
        Strand,
        TRY_CAST(NULLIF(Length, '-') AS INTEGER) as length,
        Permuted as permuted,
        Sequence as sequence,
        Source_DB,
        'phagescope' as source_type
    FROM read_csv('{trna_tmrna_data}',
                  header=true,
                  all_varchar=true,
                  ignore_errors=true,
                  null_padding=true)
    WHERE Phage_ID IS NOT NULL
    """)
    
    trna_count = conn.execute("SELECT COUNT(*) FROM dim_trna_tmrna").fetchone()[0]
    logging.info(f"✅ Created dim_trna_tmrna: {trna_count:,} rows")

    # 8. CREATE DIM_ANTIMICROBIAL_RESISTANCE_GENES TABLE
    logging.info("Creating dim_antimicrobial_resistance_genes table")
    
    conn.execute(f"""
    CREATE TABLE dim_antimicrobial_resistance_genes AS 
    SELECT 
        Phage_ID,
        Protein_ID,
        Aligned_Protein_in_CARD as aligned_protein_card,
        Source_DB,
        'phagescope' as source_type
    FROM read_csv_auto('{antimicrobial_resistance_data}')
    WHERE Phage_ID IS NOT NULL
    """)
    amr_count = conn.execute("SELECT COUNT(*) FROM dim_antimicrobial_resistance_genes").fetchone()[0]
    logging.info(f"✅ Created dim_antimicrobial_resistance_genes: {amr_count:,} rows")

    # 9. CREATE DIM_CRISPR_ARRAYS TABLE
    logging.info("Creating dim_crispr_arrays table")
    conn.execute(f"""
    CREATE TABLE dim_crispr_arrays AS 
    SELECT 
        Phage_ID,
        CRISPR_ID as crispr_id,
        TRY_CAST(NULLIF(Duplicated_Spacers, '-') AS INTEGER) as duplicated_spacers,
        TRY_CAST(NULLIF(CRISPR_Start, '-') AS INTEGER) as crispr_start,
        TRY_CAST(NULLIF(CRISPR_End, '-') AS INTEGER) as crispr_end,
        TRY_CAST(NULLIF(CRISPR_Length, '-') AS INTEGER) as crispr_length,
        "Potential_Orientation (AT%)" as potential_orientation_at_percent,
        CRISPRDirection as crispr_direction,
        Consensus_Repeat as consensus_repeat,
        "Repeat_ID (CRISPRdb)" as repeat_id_crisprdb,
        TRY_CAST(NULLIF("Nb_CRISPRs_with_same_Repeat (CRISPRdb)", '-') AS INTEGER) as nb_crisprs_same_repeat,
        TRY_CAST(NULLIF(Repeat_Length, '-') AS INTEGER) as repeat_length,
        TRY_CAST(NULLIF(Spacers_Nb, '-') AS INTEGER) as spacers_count,
        TRY_CAST(NULLIF(Mean_size_Spacers, '-') AS DOUBLE) as mean_spacer_size,
        TRY_CAST(NULLIF(Standard_Deviation_Spacers, '-') AS DOUBLE) as std_dev_spacers,
        TRY_CAST(NULLIF(Nb_Repeats_matching_Consensus, '-') AS INTEGER) as repeats_matching_consensus,
        TRY_CAST(NULLIF("Ratio_Repeats_match/TotalRepeat", '-') AS DOUBLE) as ratio_repeats_match_total,
        TRY_CAST(NULLIF("Conservation_Repeats (% identity)", '-') AS DOUBLE) as conservation_repeats_pct,
        TRY_CAST(NULLIF(EBcons_Repeats, '-') AS DOUBLE) as ebcons_repeats,
        TRY_CAST(NULLIF("Conservation_Spacers (% identity)", '-') AS DOUBLE) as conservation_spacers_pct,
        TRY_CAST(NULLIF(EBcons_Spacers, '-') AS DOUBLE) as ebcons_spacers,
        TRY_CAST(NULLIF(Repeat_Length_plus_mean_size_Spacers, '-') AS DOUBLE) as repeat_plus_spacer_length,
        TRY_CAST(NULLIF("Ratio_Repeat/mean_Spacers_Length", '-') AS DOUBLE) as ratio_repeat_spacer_length,
        "CRISPR_found_in_DB (if sequence IDs are similar)" as crispr_found_in_db,
        TRY_CAST(NULLIF(Evidence_Level, '-') AS INTEGER) as evidence_level,
        Source_DB,
        'phagescope' as source_type
    FROM read_csv('{crispr_array_data}',
                  header=true,
                  all_varchar=true,
                  ignore_errors=true,
                  null_padding=true)
    WHERE Phage_ID IS NOT NULL
    """)
    
    crispr_count = conn.execute("SELECT COUNT(*) FROM dim_crispr_arrays").fetchone()[0]
    logging.info(f"✅ Created dim_crispr_arrays: {crispr_count:,} rows")

    # 10. CREATE DIM_HOSTS TABLE (OPTIONAL - may not exist on first run)
    # This table contains host bacterial genome metadata
    # It's created from host_metadata.csv which is generated by download_host_genomes_robust.py
    # If the file doesn't exist yet, we skip this table to avoid circular dependency
    if host_metadata_path and str(host_metadata_path).strip() and os.path.exists(host_metadata_path):
        logging.info("Creating dim_hosts table")
        try:
            conn.execute(f"""
            CREATE TABLE dim_hosts AS 
            SELECT 
                Host_ID,
                Species_Name,
                Strain_Name,
                Assembly_Accession,
                Assembly_Name,
                Assembly_Level,
                TRY_CAST(NULLIF(Genome_Length, '-') AS BIGINT) as Genome_Length,
                TRY_CAST(NULLIF(GC_Content, '-') AS DOUBLE) as GC_Content,
                RefSeq_Category,
                Download_Date,
                Source,
                'phagescope' as source_type
            FROM read_csv('{host_metadata_path}',
                          header=true,
                          all_varchar=true,
                          ignore_errors=true,
                          null_padding=true)
            WHERE Host_ID IS NOT NULL
            """)
            
            host_count = conn.execute("SELECT COUNT(*) FROM dim_hosts").fetchone()[0]
            logging.info(f"✅ Created dim_hosts: {host_count:,} rows")
        except Exception as e:
            logging.warning(f"⚠️  Could not create dim_hosts table: {e}")
            logging.warning("   This is expected if host genomes haven't been downloaded yet")
    else:
        logging.info("⚠️  Skipping dim_hosts table (host metadata not available yet)")
        logging.info("   Host genomes can be downloaded after database creation")
    
    # 10a. CREATE DIM_ASSEMBLY_METADATA TABLE (NEW - comprehensive assembly info)
    # This table contains detailed assembly metadata from the robust downloader
    assembly_metadata_count = 0
    if assembly_metadata_path and os.path.exists(assembly_metadata_path):
        logging.info("Creating dim_assembly_metadata table")
        try:
            conn.execute(f"""
            CREATE TABLE dim_assembly_metadata AS
            SELECT
                Assembly_Accession,
                Assembly_Name,
                Organism_Name,
                TRY_CAST(NULLIF(Species_TaxID, '-') AS INTEGER) as Species_TaxID,
                Strain,
                Assembly_Level,
                RefSeq_Category,
                BioSample,
                BioProject,
                FTP_Path,
                Submission_Date,
                TRY_CAST(Is_Latest AS BOOLEAN) as Is_Latest,
                TRY_CAST(Quality_Score AS INTEGER) as Quality_Score,
                TRY_CAST(Is_RefSeq AS BOOLEAN) as Is_RefSeq,
                Download_Status,
                Download_Date,
                TRY_CAST(Metadata_Only AS BOOLEAN) as Metadata_Only
            FROM read_csv('{assembly_metadata_path}',
                          header=true,
                          all_varchar=true,
                          ignore_errors=true,
                          null_padding=true)
            WHERE Assembly_Accession IS NOT NULL
            """)
            
            assembly_metadata_count = conn.execute("SELECT COUNT(*) FROM dim_assembly_metadata").fetchone()[0]
            logging.info(f"✅ Created dim_assembly_metadata: {assembly_metadata_count:,} rows")
        except Exception as e:
            logging.warning(f"⚠️  Could not create dim_assembly_metadata table: {e}")
    
    # 10b. CREATE DIM_PHAGE_HOST_LINKS TABLE (NEW - phage to host assembly links)
    # This table links phages to their host assembly accessions
    phage_host_links_count = 0
    if phage_host_links_path and os.path.exists(phage_host_links_path):
        logging.info("Creating dim_phage_host_links table")
        try:
            conn.execute(f"""
            CREATE TABLE dim_phage_host_links AS
            SELECT
                Phage_ID,
                Host_Species,
                Host_Full_Name,
                Assembly_Accession,
                Assembly_Level,
                RefSeq_Category,
                Link_Quality
            FROM read_csv('{phage_host_links_path}',
                          header=true,
                          all_varchar=true,
                          ignore_errors=true,
                          null_padding=true)
            WHERE Phage_ID IS NOT NULL AND Assembly_Accession IS NOT NULL
            """)
            
            phage_host_links_count = conn.execute("SELECT COUNT(*) FROM dim_phage_host_links").fetchone()[0]
            logging.info(f"✅ Created dim_phage_host_links: {phage_host_links_count:,} rows")
        except Exception as e:
            logging.warning(f"⚠️  Could not create dim_phage_host_links table: {e}")

    # 11. INGEST PRIVATE SOURCES (NON-BLOCKING)
    # Driven entirely by the manifest produced by prepare_private_sources.
    # If the manifest is empty (no valid sources) nothing is ingested and
    # the PhageScope-only database is produced unchanged.
    private_ingestion_summary = {"ingested": [], "skipped": []}
    if private_manifest_path and os.path.exists(private_manifest_path):
        try:
            with open(private_manifest_path, "r", encoding="utf-8") as handle:
                manifest = json.load(handle) or {}
            # Only ingest validated sources from the manifest; invalid entries
            # are kept in the manifest for diagnostics but skipped here.
            valid_source_dirs = [
                source.get("source_dir")
                for source in manifest.get("sources", [])
                if source.get("is_valid", False) and source.get("source_dir")
            ]
            if valid_source_dirs:
                # Resolve BLAST DB dir and private phage mapping for duplicate detection
                blast_db_dir = snakemake.config.get("blast_db_dir", "")
                if blast_db_dir and not os.path.isabs(blast_db_dir):
                    blast_db_dir = str(REPO_ROOT / blast_db_dir.lstrip("/"))
                private_phage_mapping = snakemake.config.get("private_phage_mapping", "")
                pident_threshold = float(snakemake.config.get("duplicate_pident_threshold", 99.0))
                qcovs_threshold = float(snakemake.config.get("duplicate_qcovs_threshold", 90.0))

                private_ingestion_summary = ingest_private_sources_into_db(
                    conn,
                    valid_source_dirs,
                    blast_db_dir=blast_db_dir if blast_db_dir and os.path.exists(blast_db_dir) else None,
                    private_phage_mapping_path=private_phage_mapping if private_phage_mapping and os.path.exists(private_phage_mapping) else None,
                    pident_threshold=pident_threshold,
                    qcovs_threshold=qcovs_threshold,
                )
                dup_stats = private_ingestion_summary.get("duplicate_check", {})
                if dup_stats.get("duplicates_found", 0) > 0:
                    logging.info(
                        "⚠️  Duplicate check: %d/%d private phages match public data",
                        dup_stats["duplicates_found"], dup_stats["checked"],
                    )
                for skipped_source in private_ingestion_summary["skipped"]:
                    logging.warning(
                        "⚠️  Skipped private source '%s': %s",
                        skipped_source.get("source_db", "unknown"),
                        "; ".join(skipped_source.get("errors", [])),
                    )
                if _table_exists(conn, "dim_hosts"):
                    host_count = conn.execute("SELECT COUNT(*) FROM dim_hosts").fetchone()[0]
        except Exception as e:
            logging.warning(f"⚠️  Private ingestion failed, continuing with PhageScope-only database: {e}")

    # 12. PROVENANCE TABLES (NON-BLOCKING)
    try:
        _create_dataset_provenance_table(conn, str(dataset_provenance_manifest_csv or ""))
        dataset_provenance_count = conn.execute("SELECT COUNT(*) FROM dataset_provenance").fetchone()[0]
        logging.info(f"✅ Created dataset_provenance: {dataset_provenance_count:,} rows")
    except Exception as e:
        logging.warning(f"⚠️  Could not create dataset_provenance table: {e}")

    try:
        _create_pipeline_run_provenance_table(conn, str(pipeline_run_provenance_csv or ""))
        pipeline_run_provenance_count = conn.execute("SELECT COUNT(*) FROM pipeline_run_provenance").fetchone()[0]
        logging.info(f"✅ Created pipeline_run_provenance: {pipeline_run_provenance_count:,} rows")
    except Exception as e:
        logging.warning(f"⚠️  Could not create pipeline_run_provenance table: {e}")

    # Normalize provenance labels to avoid NULL/blank source_type values.
    # Public rows default to 'public'; rows linked to private interactions are
    # enforced as 'private'.
    # Note: some CSV ingestion paths may materialize missing values as the
    # literal string 'nan', so that sentinel is normalized as well.
    has_private_interactions = _table_exists(conn, "private_interactions")
    if has_private_interactions:
        conn.execute(
            """
            WITH private_keys AS (
                SELECT DISTINCT Phage_ID, Source_DB
                FROM private_interactions
            )
            UPDATE fact_phages AS fp
            SET source_type = CASE
                WHEN (fp.Phage_ID, fp.Source_DB) IN (
                    SELECT Phage_ID, Source_DB FROM private_keys
                ) THEN 'private'
                ELSE 'public'
            END
            WHERE source_type IS NULL
               OR trim(source_type) = ''
               OR lower(trim(source_type)) = 'nan'
            """
        )
    else:
        conn.execute(
            """
            UPDATE fact_phages
            SET source_type = 'public'
            WHERE source_type IS NULL
               OR trim(source_type) = ''
               OR lower(trim(source_type)) = 'nan'
            """
        )

    # CREATE PERFORMANCE INDEXES
    logging.info("Creating indexes")

    # Indexes for fact_phages table
    conn.execute("CREATE INDEX idx_phages_id ON fact_phages(Phage_ID)")
    conn.execute("CREATE INDEX idx_phages_source ON fact_phages(Source_DB)")
    conn.execute("CREATE INDEX idx_phages_source_type ON fact_phages(source_type)")
    conn.execute("CREATE INDEX idx_phages_dup ON fact_phages(is_duplicate_of_public)")
    
    # Indexes for dim_proteins table
    conn.execute("CREATE INDEX idx_proteins_phage ON dim_proteins(Phage_ID)")
    conn.execute("CREATE INDEX idx_proteins_source ON dim_proteins(Source_DB)")
    conn.execute("CREATE INDEX idx_proteins_source_type ON dim_proteins(source_type)")
    
    # Indexes for terminator table
    conn.execute("CREATE INDEX idx_terminators_phage ON dim_terminators(Phage_ID)")
    conn.execute("CREATE INDEX idx_terminators_source ON dim_terminators(Source_DB)")
    conn.execute("CREATE INDEX idx_terminators_source_type ON dim_terminators(source_type)")
    
    # Indexes for anti_crispr table
    conn.execute("CREATE INDEX idx_anti_crispr_phage ON dim_anti_crispr(Phage_ID)")
    conn.execute("CREATE INDEX idx_anti_crispr_source ON dim_anti_crispr(Source_DB)")
    conn.execute("CREATE INDEX idx_anti_crispr_source_type ON dim_anti_crispr(source_type)")
    
    # Indexes for virulent_factors table
    conn.execute("CREATE INDEX idx_virulent_phage ON dim_virulent_factors(Phage_ID)")
    conn.execute("CREATE INDEX idx_virulent_source ON dim_virulent_factors(Source_DB)")
    conn.execute("CREATE INDEX idx_virulent_source_type ON dim_virulent_factors(source_type)")
    
    # Indexes for transmembrane_proteins table
    conn.execute("CREATE INDEX idx_transmembrane_phage ON dim_transmembrane_proteins(Phage_ID)")
    conn.execute("CREATE INDEX idx_transmembrane_source ON dim_transmembrane_proteins(Source_DB)")
    conn.execute("CREATE INDEX idx_transmembrane_source_type ON dim_transmembrane_proteins(source_type)")
    
    # Indexes for trna_tmrna table
    conn.execute("CREATE INDEX idx_trna_phage ON dim_trna_tmrna(Phage_ID)")
    conn.execute("CREATE INDEX idx_trna_source ON dim_trna_tmrna(Source_DB)")
    conn.execute("CREATE INDEX idx_trna_source_type ON dim_trna_tmrna(source_type)")

    # Indexes for AMR table
    conn.execute("CREATE INDEX idx_amr_phage ON dim_antimicrobial_resistance_genes(Phage_ID)")
    conn.execute("CREATE INDEX idx_amr_protein ON dim_antimicrobial_resistance_genes(Protein_ID)")
    conn.execute("CREATE INDEX idx_amr_source ON dim_antimicrobial_resistance_genes(Source_DB)")
    conn.execute("CREATE INDEX idx_amr_source_type ON dim_antimicrobial_resistance_genes(source_type)")
    
    # Indexes for CRISPR table
    conn.execute("CREATE INDEX idx_crispr_phage ON dim_crispr_arrays(Phage_ID)")
    conn.execute("CREATE INDEX idx_crispr_id ON dim_crispr_arrays(crispr_id)")
    conn.execute("CREATE INDEX idx_crispr_source ON dim_crispr_arrays(Source_DB)")
    conn.execute("CREATE INDEX idx_crispr_evidence ON dim_crispr_arrays(evidence_level)")
    conn.execute("CREATE INDEX idx_crispr_source_type ON dim_crispr_arrays(source_type)")
    
    # Indexes for hosts table (if it exists)
    if host_count > 0:
        conn.execute("CREATE INDEX idx_hosts_id ON dim_hosts(Host_ID)")
        conn.execute("CREATE INDEX idx_hosts_species ON dim_hosts(Species_Name)")
        conn.execute("CREATE INDEX idx_hosts_accession ON dim_hosts(Assembly_Accession)")
        conn.execute("CREATE INDEX idx_hosts_source_type ON dim_hosts(source_type)")
    
    # Indexes for assembly metadata table (if it exists)
    if assembly_metadata_count > 0:
        conn.execute("CREATE INDEX idx_assembly_accession ON dim_assembly_metadata(Assembly_Accession)")
        conn.execute("CREATE INDEX idx_assembly_organism ON dim_assembly_metadata(Organism_Name)")
        conn.execute("CREATE INDEX idx_assembly_taxid ON dim_assembly_metadata(Species_TaxID)")
    
    # Indexes for phage-host links table (if it exists)
    if phage_host_links_count > 0:
        conn.execute("CREATE INDEX idx_phage_host_phage ON dim_phage_host_links(Phage_ID)")
        conn.execute("CREATE INDEX idx_phage_host_assembly ON dim_phage_host_links(Assembly_Accession)")

    if _table_exists(conn, "dataset_provenance") and dataset_provenance_count > 0:
        conn.execute("CREATE INDEX idx_dataset_provenance_feature ON dataset_provenance(feature)")
        conn.execute("CREATE INDEX idx_dataset_provenance_source_key ON dataset_provenance(source_key)")
    
    # CREATE ANALYTICAL VIEWS
    logging.info("Creating analytical views")
    conn.execute("""
    CREATE VIEW phage_summary AS
    SELECT 
        Source_DB,
        COUNT(*) as total_phages,
        AVG(Length) as avg_length,
        AVG(GC_content) as avg_gc_content,
        MIN(Length) as min_length,
        MAX(Length) as max_length
    FROM fact_phages
    WHERE Length IS NOT NULL AND GC_content IS NOT NULL
    GROUP BY Source_DB
    ORDER BY total_phages DESC
    """)
    
    conn.execute("""
    CREATE VIEW phage_size_distribution AS
    SELECT 
        Source_DB,
        CASE 
            WHEN Length < 10000 THEN 'Small (<10kb)'
            WHEN Length < 100000 THEN 'Medium (10-100kb)'
            ELSE 'Large (>100kb)'
        END as size_category,
        COUNT(*) as count,
        AVG(Length) as avg_length
    FROM fact_phages
    WHERE Length IS NOT NULL
    GROUP BY Source_DB, size_category
    """)
    
    # CREATE COMPREHENSIVE PHAGE VIEW (updated with AMR and CRISPR)
    conn.execute("""
    CREATE VIEW phage_complete_profile AS
    SELECT 
        f.Phage_ID,
        f.Source_DB,
        f.Length,
        f.GC_content,
        f.Host,
        f.Lifestyle,
        COUNT(DISTINCT p.Protein_ID) as protein_count,
        COUNT(DISTINCT t.terminator_type) as terminator_types,
        COUNT(DISTINCT a.Protein_ID) as anti_crispr_count,
        COUNT(DISTINCT v.Protein_ID) as virulent_factor_count,
        COUNT(DISTINCT tm.Protein_ID) as transmembrane_protein_count,
        COUNT(DISTINCT tr.trna_tmrna_id) as trna_tmrna_count,
        COUNT(DISTINCT amr.Protein_ID) as amr_gene_count,
        COUNT(DISTINCT cr.crispr_id) as crispr_array_count,
        AVG(cr.spacers_count) as avg_spacers_per_array,
        MAX(cr.evidence_level) as max_crispr_evidence
    FROM fact_phages f
    LEFT JOIN dim_proteins p ON f.Phage_ID = p.Phage_ID
    LEFT JOIN dim_terminators t ON f.Phage_ID = t.Phage_ID
    LEFT JOIN dim_anti_crispr a ON f.Phage_ID = a.Phage_ID
    LEFT JOIN dim_virulent_factors v ON f.Phage_ID = v.Phage_ID
    LEFT JOIN dim_transmembrane_proteins tm ON f.Phage_ID = tm.Phage_ID
    LEFT JOIN dim_trna_tmrna tr ON f.Phage_ID = tr.Phage_ID
    LEFT JOIN dim_antimicrobial_resistance_genes amr ON f.Phage_ID = amr.Phage_ID
    LEFT JOIN dim_crispr_arrays cr ON f.Phage_ID = cr.Phage_ID
    GROUP BY f.Phage_ID, f.Source_DB, f.Length, f.GC_content, f.Host, f.Lifestyle
    """)
    
    # AMR GENE SUMMARY VIEW
    conn.execute("""
    CREATE VIEW amr_gene_summary AS
    SELECT 
        Source_DB,
        COUNT(DISTINCT Phage_ID) as phages_with_amr,
        COUNT(DISTINCT Protein_ID) as total_amr_proteins,
        COUNT(DISTINCT aligned_protein_card) as unique_card_matches,
        ROUND(COUNT(DISTINCT Phage_ID) * 100.0 / 
              (SELECT COUNT(DISTINCT Phage_ID) FROM fact_phages WHERE Source_DB = amr.Source_DB), 2) 
              as percentage_phages_with_amr
    FROM dim_antimicrobial_resistance_genes amr
    GROUP BY Source_DB
    ORDER BY phages_with_amr DESC
    """)
    
    # CRISPR ARRAY SUMMARY VIEW
    conn.execute("""
    CREATE VIEW crispr_array_summary AS
    SELECT 
        Source_DB,
        COUNT(DISTINCT Phage_ID) as phages_with_crispr,
        COUNT(DISTINCT crispr_id) as total_crispr_arrays,
        AVG(crispr_length) as avg_array_length,
        AVG(spacers_count) as avg_spacers_per_array,
        AVG(repeat_length) as avg_repeat_length,
        AVG(mean_spacer_size) as avg_spacer_size,
        AVG(conservation_repeats_pct) as avg_repeat_conservation,
        ROUND(COUNT(DISTINCT Phage_ID) * 100.0 / 
              (SELECT COUNT(DISTINCT Phage_ID) FROM fact_phages WHERE Source_DB = cr.Source_DB), 2) 
              as percentage_phages_with_crispr
    FROM dim_crispr_arrays cr
    WHERE crispr_length IS NOT NULL
    GROUP BY Source_DB
    ORDER BY phages_with_crispr DESC
    """)
    
    # CRISPR EVIDENCE LEVEL DISTRIBUTION
    conn.execute("""
    CREATE VIEW crispr_evidence_distribution AS
    SELECT 
        evidence_level,
        CASE 
            WHEN evidence_level = 4 THEN 'Very High'
            WHEN evidence_level = 3 THEN 'High'
            WHEN evidence_level = 2 THEN 'Medium'
            WHEN evidence_level = 1 THEN 'Low'
            ELSE 'Unknown'
        END as evidence_category,
        COUNT(*) as array_count,
        AVG(spacers_count) as avg_spacers,
        AVG(crispr_length) as avg_length,
        AVG(conservation_repeats_pct) as avg_conservation
    FROM dim_crispr_arrays
    WHERE evidence_level IS NOT NULL
    GROUP BY evidence_level
    ORDER BY evidence_level DESC
    """)
    
    # PHAGE DEFENSE SYSTEMS PROFILE
    conn.execute("""
    CREATE VIEW phage_defense_profile AS
    SELECT 
        f.Phage_ID,
        f.Source_DB,
        f.Host,
        f.Lifestyle,
        COUNT(DISTINCT cr.crispr_id) as crispr_arrays,
        COUNT(DISTINCT ac.Protein_ID) as anti_crispr_proteins,
        COUNT(DISTINCT amr.Protein_ID) as amr_genes,
        CASE 
            WHEN COUNT(DISTINCT cr.crispr_id) > 0 AND COUNT(DISTINCT ac.Protein_ID) > 0 
                THEN 'CRISPR + Anti-CRISPR'
            WHEN COUNT(DISTINCT cr.crispr_id) > 0 THEN 'CRISPR only'
            WHEN COUNT(DISTINCT ac.Protein_ID) > 0 THEN 'Anti-CRISPR only'
            WHEN COUNT(DISTINCT amr.Protein_ID) > 0 THEN 'AMR only'
            ELSE 'None detected'
        END as defense_category,
        MAX(cr.evidence_level) as max_crispr_evidence
    FROM fact_phages f
    LEFT JOIN dim_crispr_arrays cr ON f.Phage_ID = cr.Phage_ID
    LEFT JOIN dim_anti_crispr ac ON f.Phage_ID = ac.Phage_ID
    LEFT JOIN dim_antimicrobial_resistance_genes amr ON f.Phage_ID = amr.Phage_ID
    GROUP BY f.Phage_ID, f.Source_DB, f.Host, f.Lifestyle
    """)
    
    # HOST-SPECIFIC AMR PROFILE
    conn.execute("""
    CREATE VIEW host_amr_profile AS
    SELECT 
        f.Host,
        f.Source_DB,
        COUNT(DISTINCT f.Phage_ID) as total_phages,
        COUNT(DISTINCT amr.Phage_ID) as phages_with_amr,
        COUNT(DISTINCT amr.Protein_ID) as total_amr_genes,
        ROUND(COUNT(DISTINCT amr.Phage_ID) * 100.0 / COUNT(DISTINCT f.Phage_ID), 2) as amr_prevalence_pct
    FROM fact_phages f
    LEFT JOIN dim_antimicrobial_resistance_genes amr ON f.Phage_ID = amr.Phage_ID
    WHERE f.Host IS NOT NULL AND f.Host != '-'
    GROUP BY f.Host, f.Source_DB
    HAVING COUNT(DISTINCT f.Phage_ID) > 10
    ORDER BY amr_prevalence_pct DESC
    """)
    
    # CRISPR SPACER STATISTICS BY HOST
    conn.execute("""
    CREATE VIEW host_crispr_profile AS
    SELECT 
        f.Host,
        f.Source_DB,
        COUNT(DISTINCT f.Phage_ID) as total_phages,
        COUNT(DISTINCT cr.Phage_ID) as phages_with_crispr,
        AVG(cr.spacers_count) as avg_spacers,
        MAX(cr.spacers_count) as max_spacers,
        AVG(cr.conservation_repeats_pct) as avg_repeat_conservation,
        ROUND(COUNT(DISTINCT cr.Phage_ID) * 100.0 / COUNT(DISTINCT f.Phage_ID), 2) as crispr_prevalence_pct
    FROM fact_phages f
    LEFT JOIN dim_crispr_arrays cr ON f.Phage_ID = cr.Phage_ID
    WHERE f.Host IS NOT NULL AND f.Host != '-'
    GROUP BY f.Host, f.Source_DB
    HAVING COUNT(DISTINCT f.Phage_ID) > 10
    ORDER BY crispr_prevalence_pct DESC
    """)
    
    # HOST-RELATED VIEWS (only if host table exists)
    if host_count > 0:
        logging.info("Creating host-related analytical views")
        
        # View linking phages to their downloaded host genomes
        conn.execute("""
        CREATE VIEW phage_host_genomes AS
        SELECT 
            f.Phage_ID,
            f.Host as Host_Name,
            h.Host_ID,
            h.Species_Name,
            h.Strain_Name,
            h.Assembly_Accession,
            h.Genome_Length as Host_Genome_Length,
            h.GC_Content as Host_GC_Content,
            f.Length as Phage_Length,
            f.GC_content as Phage_GC_Content,
            f.Source_DB
        FROM fact_phages f
        INNER JOIN dim_hosts h ON f.Host LIKE h.Species_Name || '%'
        WHERE f.Host IS NOT NULL AND f.Host != '-'
        """)
        
        # Summary of host genome statistics
        conn.execute("""
        CREATE VIEW host_genome_summary AS
        SELECT 
            Species_Name,
            COUNT(*) as genome_count,
            AVG(Genome_Length) as avg_genome_length,
            AVG(GC_Content) as avg_gc_content,
            Assembly_Level,
            RefSeq_Category
        FROM dim_hosts
        GROUP BY Species_Name, Assembly_Level, RefSeq_Category
        ORDER BY Species_Name
        """)

    conn.execute("DROP VIEW IF EXISTS phage_host_associations")
    if _table_exists(conn, "dim_hosts"):
        has_public_links = _table_exists(conn, "dim_phage_host_links")
        has_private_links = _table_exists(conn, "private_phage_host_associations")
        if has_public_links and has_private_links:
            conn.execute("""
            CREATE VIEW phage_host_associations AS
            SELECT DISTINCT
                phl.Phage_ID,
                h.Host_ID
            FROM dim_phage_host_links phl
            JOIN dim_hosts h ON phl.Assembly_Accession = h.Assembly_Accession
            WHERE phl.Phage_ID IS NOT NULL AND h.Host_ID IS NOT NULL
            UNION
            SELECT DISTINCT
                Phage_ID,
                Host_ID
            FROM private_phage_host_associations
            """)
        elif has_public_links:
            conn.execute("""
            CREATE VIEW phage_host_associations AS
            SELECT DISTINCT
                phl.Phage_ID,
                h.Host_ID
            FROM dim_phage_host_links phl
            JOIN dim_hosts h ON phl.Assembly_Accession = h.Assembly_Accession
            WHERE phl.Phage_ID IS NOT NULL AND h.Host_ID IS NOT NULL
            """)
        elif has_private_links:
            conn.execute("""
            CREATE VIEW phage_host_associations AS
            SELECT DISTINCT
                Phage_ID,
                Host_ID
            FROM private_phage_host_associations
            """)
        else:
            conn.execute("""
            CREATE VIEW phage_host_associations AS
            SELECT NULL::VARCHAR AS Phage_ID, NULL::VARCHAR AS Host_ID
            WHERE FALSE
            """)
            # Placeholder empty view keeps downstream queries stable when no associations exist.
        logging.info("✅ Created phage_host_associations view (public + private)")

    conn.close()
    
    logging.info(f"✅ Database created successfully!")
    logging.info(f"   📊 Summary Statistics:")
    logging.info(f"   • Phages: {phage_count:,}")
    logging.info(f"   • Proteins: {protein_count:,}")
    logging.info(f"   • Terminators: {terminator_count:,}")
    logging.info(f"   • Anti-CRISPR: {anti_crispr_count:,}")
    logging.info(f"   • Virulent Factors: {virulent_count:,}")
    logging.info(f"   • Transmembrane Proteins: {transmembrane_count:,}")
    logging.info(f"   • tRNA/tmRNA: {trna_count:,}")
    logging.info(f"   • AMR Genes: {amr_count:,}")
    logging.info(f"   • CRISPR Arrays: {crispr_count:,}")
    if host_count > 0:
        logging.info(f"   • Host Genomes: {host_count:,}")
    if assembly_metadata_count > 0:
        logging.info(f"   • Assembly Metadata: {assembly_metadata_count:,}")
    if phage_host_links_count > 0:
        logging.info(f"   • Phage-Host Links: {phage_host_links_count:,}")
    logging.info(f"   • Dataset Provenance Rows: {dataset_provenance_count:,}")
    logging.info(f"   • Pipeline Run Provenance Rows: {pipeline_run_provenance_count:,}")
    n_ingested = len(private_ingestion_summary["ingested"])
    n_skipped = len(private_ingestion_summary["skipped"])
    if n_ingested or n_skipped:
        logging.info(
            "   • Private Sources: %d ingested / %d skipped",
            n_ingested,
            n_skipped,
        )

if __name__ == "__main__":
    create_star_schema_duckdb()
