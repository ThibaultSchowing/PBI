#!/usr/bin/env python

import duckdb
import shutil
import logging

logging.basicConfig(level=logging.INFO)

def optimize_database():
    """Optimize database for performance"""
    
    input_db = snakemake.input.db
    output_db = snakemake.output.optimized_db
    
    logging.info(f"Optimizing database from {input_db} to {output_db}")
    
    # Copy the database
    shutil.copy2(input_db, output_db)
    
    # Open and optimize
    conn = duckdb.connect(output_db)
    
    # Performance settings
    conn.execute("PRAGMA memory_limit='4GB'")
    conn.execute("PRAGMA threads=4")
    
    # Analyze tables for query optimization
    logging.info("Analyzing tables")
    conn.execute("ANALYZE fact_phages")
    conn.execute("ANALYZE dim_proteins") 
    conn.execute("ANALYZE dim_terminators")
    conn.execute("ANALYZE dim_anti_crispr")
    conn.execute("ANALYZE dim_virulent_factors")
    conn.execute("ANALYZE dim_transmembrane_proteins")
    conn.execute("ANALYZE dim_trna_tmrna")
    conn.execute("ANALYZE dim_antimicrobial_resistance_genes")
    conn.execute("ANALYZE dim_crispr_arrays")
    
    # Vacuum for cleanup
    logging.info("Vacuuming database")
    conn.execute("VACUUM")
    
    conn.close()
    
    logging.info("✅ Database optimization completed")

if __name__ == "__main__":
    optimize_database()