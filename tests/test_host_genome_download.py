#!/usr/bin/env python3
"""
Test the host genome download CSV extraction logic.
This test verifies that unique hosts are correctly extracted from the phage metadata CSV.
"""

import os
import sys
import tempfile
import pandas as pd
from pathlib import Path


def extract_unique_hosts_from_csv(csv_path):
    """
    Extract unique host species from phage metadata CSV
    This mirrors the logic in download_host_genomes_robust.py
    """
    df = pd.read_csv(csv_path)
    
    if 'Host' not in df.columns:
        raise ValueError("Host column not found in phage metadata CSV")
    
    # Filter for valid hosts using a single boolean mask for efficiency
    valid_mask = (
        df['Host'].notna() &                                                    # Not null
        (df['Host'] != '-') &                                                    # Not dash
        (df['Host'] != '') &                                                     # Not empty
        (~df['Host'].str.contains('unknown', case=False, na=False)) &           # Not unknown
        (~df['Host'].str.contains('unidentified', case=False, na=False))       # Not unidentified
    )
    
    # Get unique hosts
    unique_hosts = df.loc[valid_mask, 'Host'].unique()
    
    # Extract species names (Genus species format)
    species_names = set()
    for host in unique_hosts:
        parts = str(host).strip().split()
        if len(parts) >= 2:
            if parts[0][0].isupper():
                species = f"{parts[0]} {parts[1]}"
                species_names.add(species)
        elif len(parts) == 1:
            if parts[0][0].isupper():
                species_names.add(parts[0])
    
    return sorted(list(species_names))


def test_host_extraction_from_csv():
    """Test that unique hosts are correctly extracted from CSV"""
    
    # Create a temporary CSV with sample phage data
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        csv_path = f.name
        
        # Write sample data
        f.write("Phage_ID,Host,Source_DB\n")
        f.write("phage1,Escherichia coli,RefSeq\n")
        f.write("phage2,Escherichia coli K12,RefSeq\n")
        f.write("phage3,Staphylococcus aureus,Genbank\n")
        f.write("phage4,Pseudomonas aeruginosa,RefSeq\n")
        f.write("phage5,-,RefSeq\n")  # Should be filtered
        f.write("phage6,unknown host,RefSeq\n")  # Should be filtered
        f.write("phage7,Escherichia coli,RefSeq\n")  # Duplicate
    
    try:
        # Test extraction
        species = extract_unique_hosts_from_csv(csv_path)
        
        # Verify results
        print(f"✅ Extracted species: {species}")
        
        # Expected: Escherichia coli, Staphylococcus aureus, Pseudomonas aeruginosa
        assert len(species) == 3, f"Expected 3 species, got {len(species)}: {species}"
        assert "Escherichia coli" in species, "Escherichia coli not found"
        assert "Staphylococcus aureus" in species, "Staphylococcus aureus not found"
        assert "Pseudomonas aeruginosa" in species, "Pseudomonas aeruginosa not found"
        
        print("✅ All assertions passed!")
        return True
        
    finally:
        # Clean up
        if os.path.exists(csv_path):
            os.unlink(csv_path)


def test_csv_missing_host_column():
    """Test that script handles missing Host column gracefully"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        csv_path = f.name
        f.write("Phage_ID,Source_DB\n")
        f.write("phage1,RefSeq\n")
    
    try:
        try:
            species = extract_unique_hosts_from_csv(csv_path)
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "Host column not found" in str(e)
            print("✅ Correctly raised ValueError for missing Host column")
            return True
            
    finally:
        if os.path.exists(csv_path):
            os.unlink(csv_path)


def test_edge_cases():
    """Test edge cases in host extraction"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        csv_path = f.name
        
        # Write edge cases
        f.write("Phage_ID,Host,Source_DB\n")
        f.write("phage1,Bacillus,RefSeq\n")  # Single word (genus only)
        f.write("phage2,escherichia coli,RefSeq\n")  # lowercase - should be filtered
        f.write("phage3,Mycobacterium tuberculosis H37Rv,RefSeq\n")  # Extra strain info
        f.write("phage4,,RefSeq\n")  # Empty host
        f.write("phage5,Salmonella enterica serovar Typhimurium,RefSeq\n")  # Complex name
    
    try:
        species = extract_unique_hosts_from_csv(csv_path)
        
        print(f"✅ Extracted species from edge cases: {species}")
        
        # Should extract: Bacillus, Mycobacterium tuberculosis, Salmonella enterica
        assert "Bacillus" in species, "Bacillus (single word) not found"
        assert "Mycobacterium tuberculosis" in species, "Mycobacterium tuberculosis not found"
        assert "Salmonella enterica" in species, "Salmonella enterica not found"
        
        # Should NOT extract lowercase genus
        assert len([s for s in species if s.lower().startswith('escherichia')]) == 0, \
            "Lowercase genus should not be extracted"
        
        print("✅ Edge cases handled correctly!")
        return True
        
    finally:
        if os.path.exists(csv_path):
            os.unlink(csv_path)


if __name__ == "__main__":
    print("Testing host genome download CSV extraction logic...")
    print("=" * 70)
    
    print("\n=== Test 1: Extract unique hosts from CSV ===")
    test_host_extraction_from_csv()
    
    print("\n=== Test 2: Handle missing Host column ===")
    test_csv_missing_host_column()
    
    print("\n=== Test 3: Edge cases ===")
    test_edge_cases()
    
    print("\n" + "=" * 70)
    print("✅ All tests passed!")

