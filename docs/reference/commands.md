# Command Reference

Quick reference for common PBI-Scope operations and useful commands.

## Docker Commands

### Pipeline

```bash
# Build pipeline image
docker compose build pipeline

# Run pipeline
docker compose run --rm pipeline

# View pipeline logs
docker compose logs pipeline

# Run pipeline with custom cores
docker compose run --rm pipeline snakemake --cores 2 --use-conda

# Force host resolution in Docker
docker compose run --rm pipeline \
  snakemake --cores all --use-conda --printshellcmds \
  --directory /app/workflow --snakefile /app/workflow/Snakefile \
  --forcerun download_host_genomes

# Force host resolution in Docker and disable token-cache reuse
docker compose run --rm pipeline \
  snakemake --cores all --use-conda --printshellcmds \
  --directory /app/workflow --snakefile /app/workflow/Snakefile \
  --forcerun download_host_genomes \
  --config reuse_host_resolution_cache=false
```

### API

```bash
# Build API image
docker compose build api

# Start API (detached)
docker compose up -d api

# Start API (foreground with logs)
docker compose up api

# View API logs
docker compose logs -f api

# Stop API
docker compose down

# Restart API
docker compose restart api
```

### Data Access

```bash
# List database files
docker run --rm -v pbi-data:/data alpine ls -lh /data/processed/databases/

# List reports
docker run --rm -v pbi-data:/data alpine ls -lh /data/processed/reports/

# Copy database to host
docker run --rm -v pbi-data:/data -v $(pwd):/backup alpine \
  cp /data/processed/databases/phage_database_optimized.duckdb /backup/

# Copy reports to host
docker run --rm -v pbi-data:/data -v $(pwd):/backup alpine \
  cp -r /data/processed/reports /backup/

# Copy specific report
docker run --rm -v pbi-data:/data -v $(pwd):/backup alpine \
  cp /data/processed/reports/database_validation.html /backup/
```

### Volume Management

```bash
# List volumes
docker volume ls

# Inspect volume
docker volume inspect pbi-data

# Remove cache volume (keeps data)
docker compose down
docker volume rm pbi-cache

# Remove all volumes (CAREFUL - deletes all data!)
docker compose down -v

# Clean up using provided script
./cleanup_cache.sh
```

### Troubleshooting

```bash
# Check running containers
docker ps

# Check all containers (including stopped)
docker ps -a

# View container resource usage
docker stats

# Clean up Docker system
docker system prune -a --volumes

# Rebuild without cache
docker compose build --no-cache pipeline
docker compose build --no-cache api
```

## Pipeline Execution

The basic pipeline execution is executed from the pipeline docker container ('/workflow/Dockerfile') and (re)execute all required tasks. Should you need to force or re-execute specific tasks, you can run specific pipeline tasks within the pipeline container. 

### Basic Execution

```bash
# Run full pipeline (first run: use 2-4 cores due to I/O)
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda --printshellcmds

# Run with caching (recommended for development)
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cache --cores 4 --use-conda --printshellcmds

# Keep temporary files (useful for debugging)
snakemake --directory workflow --snakefile workflow/Snakefile \
  --notemp --cores 4 --use-conda --printshellcmds

# Dry run (see what would execute)
snakemake --directory workflow --snakefile workflow/Snakefile -n

# Use all available cores
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores all --use-conda
```

### Specific Targets

```bash
# Create database only
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda \
  ../data/databases/phage_database.duckdb

# Create optimized database
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda \
  ../data/databases/phage_database_optimized.duckdb

# Generate validation reports only
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda \
  reports/database_validation.html
```

### Force Re-run (Snakemake)

```bash
# Force host resolution/download rule
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda \
  --forcerun download_host_genomes

# Force host resolution and ignore persisted token cache for this run
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda \
  --forcerun download_host_genomes \
  --config reuse_host_resolution_cache=false

# Force CSV download + merge rule examples
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda \
  --forcerun download_all_tsvs merge_phage_metadata_tsvs
```

### Workflow Visualization

```bash
# Generate DAG (Directed Acyclic Graph) visualization
cd workflow
snakemake --dag | dot -Tsvg > dag/workflow.svg

# Generate rulegraph (simplified)
snakemake --rulegraph | dot -Tsvg > dag/rulegraph.svg

# View in browser
xdg-open dag/workflow.svg  # Linux
open dag/workflow.svg      # macOS
```

### Cleanup

```bash
# Remove temporary files after execution
snakemake --delete-temp-output

# Clean all generated files (CAREFUL!)
snakemake --directory workflow --snakefile workflow/Snakefile --delete-all-output

# Remove conda environments
rm -rf workflow/.snakemake/conda/

# Remove Snakemake metadata
rm -rf workflow/.snakemake/
```

## Database Operations

Similarly as for the Snakemake workflow, the DuckDB is created by the pipeline itself. The Python 'pbi' package or the API can then access it in read-only. 

### DuckDB CLI

```bash
# Connect to database
duckdb data/databases/phage_database_optimized.duckdb

# Inside DuckDB:
# Show all tables
.tables

# Describe table schema
DESCRIBE fact_phages;

# Show table row count
SELECT COUNT(*) FROM fact_phages;

# Exit
.quit
```

### Python/DuckDB

```python
import duckdb

# Connect to database
conn = duckdb.connect('data/databases/phage_database_optimized.duckdb')

# Simple query
result = conn.execute("SELECT COUNT(*) FROM fact_phages").fetchone()
print(f"Total phages: {result[0]:,}")

# Query to DataFrame
df = conn.execute("""
    SELECT Source_DB, COUNT(*) as count 
    FROM fact_phages 
    GROUP BY Source_DB
""").df()
print(df)

# Close connection
conn.close()
```

### Common Queries

```sql
-- Count phages by source
SELECT Source_DB, COUNT(*) as count
FROM fact_phages
GROUP BY Source_DB
ORDER BY count DESC;

-- Get phages by host
SELECT Phage_ID, Host, Length, GC_content
FROM fact_phages
WHERE Host LIKE '%Escherichia%'
LIMIT 10;

-- Find large phages
SELECT Phage_ID, Length, Host, Lifestyle
FROM fact_phages
WHERE Length > 200000
ORDER BY Length DESC
LIMIT 20;

-- Phage with protein count
SELECT 
    f.Phage_ID,
    f.Length,
    f.Host,
    COUNT(p.Protein_ID) as protein_count
FROM fact_phages f
LEFT JOIN dim_proteins p ON f.Phage_ID = p.Phage_ID
GROUP BY f.Phage_ID, f.Length, f.Host
ORDER BY protein_count DESC
LIMIT 10;

-- tRNA type distribution
SELECT 
    trna_type,
    COUNT(*) as count
FROM dim_trna_tmrna
WHERE trna_type IS NOT NULL
GROUP BY trna_type
ORDER BY count DESC;
```

## Environment Setup

### Conda

```bash
# Create environment
conda create -n pbi python=3.10

# Activate environment
conda activate pbi

# Deactivate environment
conda deactivate

# Remove environment
conda env remove -n pbi

# Export environment
conda env export > environment.yml

# Create from file
conda env create -f environment.yml
```

### PBI Package

```bash
# Install in development mode
pip install -e .

# Install with specific extras
pip install -e ".[dev]"

# Verify installation
python -c "import pbi; print(pbi.__version__)"

# Validate private dataset roots before full pipeline run
# (run from the repository root — defaults to ./private_data/)
pbi validate-private

# Or point to an explicit root
pbi validate-private --path /path/to/private_root

# Uninstall
pip uninstall pbi
```

### Snakemake Cache

```bash
# Create cache directory
mkdir -p /mnt/snakemake-cache

# Set cache location (add to ~/.bashrc for persistence)
export SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache/

# Verify cache is set
echo $SNAKEMAKE_OUTPUT_CACHE
```

## API Operations

Here as well, the normal usecase is to start the docker container but you can still start the API locally. 

### Starting API

```bash
# Local (development with auto-reload)
cd api
uvicorn app:app --reload --host 0.0.0.0 --port 8000

# Local (production)
cd api
uvicorn app:app --host 0.0.0.0 --port 8000 --workers 4

# Docker
docker compose up api
```

### Testing API

```bash
# Health check
curl http://localhost:8000/health

# Get statistics
curl http://localhost:8000/stats

# Get phage metadata
curl "http://localhost:8000/phage-metadata?limit=10"

# Filtered phage metadata
curl "http://localhost:8000/phage-metadata?where=Source_DB%20%3D%20%27RefSeq%27&limit=50"

# Get host metadata
curl "http://localhost:8000/host-metadata?limit=10"

# Get protein metadata
curl "http://localhost:8000/protein-metadata?limit=10"

# Single phage sequence
curl http://localhost:8000/phage/NC_001330.1/sequence

# SQL query
curl -X POST http://localhost:8000/query \
  -H "Content-Type: application/json" \
  -d '{"query": "SELECT Source_DB, COUNT(*) FROM fact_phages GROUP BY Source_DB"}'
```

## File Operations

### Check File Sizes

```bash
# Database size
du -h data/databases/phage_database_optimized.duckdb

# All databases
du -sh data/databases/*

# Sequence files
du -sh data/sequences/*

# Total data directory size
du -sh data/
```

### Find Files

```bash
# Find all FASTA files
find data/ -name "*.fasta"

# Find all databases
find data/ -name "*.duckdb"

# Find all reports
find workflow/reports/ -name "*.html"

# Find large files (>1GB)
find data/ -type f -size +1G
```

## Git Operations

```bash
# Check what would be committed
git status

# View changes
git diff

# View file history
git log --oneline -- path/to/file

# Discard changes to specific file
git checkout -- path/to/file

# Pull latest changes
git pull origin main

# Create feature branch
git checkout -b feature/my-feature
```

## Jupyter

```bash
# Start Jupyter Lab
jupyter lab

# Start on specific port
jupyter lab --port 8889

# Start without browser
jupyter lab --no-browser

# List running servers
jupyter lab list

# Stop server
jupyter lab stop 8888
```

## Process Management

```bash
# Find process by port
lsof -i :8000

# Kill process by PID
kill <PID>

# Kill by name (use with caution)
pkill -f "uvicorn"

# Monitor system resources
htop

# Check disk space
df -h

# Check memory usage
free -h
```

## Debugging

### Snakemake Debug Mode

```bash
# Verbose output
snakemake --directory workflow --snakefile workflow/Snakefile \
  -v --printshellcmds --cores 4 --use-conda

# Keep going despite errors
snakemake --directory workflow --snakefile workflow/Snakefile \
  --keep-going --cores 4 --use-conda

# Print execution reason
snakemake --directory workflow --snakefile workflow/Snakefile \
  -p -r --cores 4 --use-conda
```

### Python Debugging

```python
# Enable debugging
import logging
logging.basicConfig(level=logging.DEBUG)

# In code
import pdb; pdb.set_trace()  # Set breakpoint

# Better debugging with ipdb
import ipdb; ipdb.set_trace()
```

### Docker Debugging

```bash
# Interactive shell in container
docker compose run --rm pipeline /bin/bash

# Check container logs
docker compose logs --tail 100 api

# Inspect running container
docker exec -it <container_id> /bin/bash

# View container processes
docker top <container_id>
```

## Performance Monitoring

```bash
# Monitor pipeline execution
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda --printshellcmds \
  --benchmark-repeats 3

# Time command execution
time snakemake --directory workflow --snakefile workflow/Snakefile --cores 4

# Profile Python script
python -m cProfile -o output.prof script.py

# View profile
python -m pstats output.prof
```

## Quick Workflows

### Complete Fresh Installation without Docker (not a good idea unless you know what you are doing)

```bash
# 1. Clone repository
git clone https://github.com/ThibaultSchowing/PBI-Scope.git
cd PBI

# 2. Create conda environment
conda create -n pbi python=3.10
conda activate pbi

# 3. Install PBI package
pip install -e .

# 4. Run pipeline
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda --printshellcmds

# 5. Verify
ls -lh data/databases/
ls -lh workflow/reports/
```

### Docker Fresh Start

```bash
# 1. Build images
docker compose build

# 2. Run pipeline
docker compose run --rm pipeline

# 3. Start API (-d flags runs the container **detached** which will bring up the container in the background)
docker compose up -d api

# 4. Test
curl http://localhost:8000/health
```

### Update Existing Installation

```bash
# 1. Pull latest code
git pull origin main

# 2. Update package
pip install -e . --upgrade

# 3. Re-run pipeline
snakemake --directory workflow --snakefile workflow/Snakefile \
  --cores 4 --use-conda

# 4. Verify
curl http://localhost:8000/stats
```

---

**Tip**: Bookmark this page for quick reference when working with PBI-Scope!
