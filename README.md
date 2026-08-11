![alt](https://github.com/ThibaultSchowing/PBI/blob/main/docs/img/PBI_Schema_Note.png)

# PBI-Scope
## Dockerized Phage Bacteria Interactions toolkit based on PhageScope

> A proof-of-concept dockerized bioinformatics pipeline that makes phage genomic data from [PhageScope](https://phagescope.deepomics.org/database) and their hosts available in an efficient, structured format for training neural networks and AI models for phage-host interaction prediction. 

**Install - Wait - Work** The pipeline takes care of everything within Docker !

[![Documentation](https://img.shields.io/badge/docs-github%20pages-blue)](https://thibaultschowing.github.io/PBI-Scope/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo-blue.svg)](https://doi.org/10.5281/zenodo.21469490)
[![CI Pipeline and DB Tests](https://github.com/ThibaultSchowing/PBI-Scope/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ThibaultSchowing/PBI-Scope/actions/workflows/ci.yml)
[![ExPASy SIB](https://img.shields.io/badge/ExPASy-SIB_Resource-E2001A)](https://www.expasy.org/resources/pbi-scope)

## 🎯 What is PBI-Scope?

PBI-Scope builds a unified data product from:

- Public phage metadata/sequences from PhageScope
- Optional private datasets from `private_data/` (validated and merged by source)
- Host genomes resolved from NCBI RefSeq

> **Note**: PBI-Scope is dependent on PhageScope as its primary data source. Regarding data such as Host range or lifestyle, unavailable data were predicted using various tools (e.g. DeepHost). Refer to [the publication](https://academic.oup.com/nar/article/52/D1/D756/7334092) for more information.

Outputs are stored in a shared Docker volume and exposed through:

- DuckDB metadata database
- Indexed FASTA files (phage/protein/host)
- `pbi` Python package (recommended access path)
- Analysis container (Jupyter Lab + VS Code Dev Containers)

> The REST API is now supported for database exploration ! For sequence-heavy usage, load the database's sequence retriever directly from the analysis container. **Check Notebook examples [in the notebooks folder !](https://github.com/ThibaultSchowing/PBI-Scope/tree/main/notebooks)**

![](https://github.com/ThibaultSchowing/PBI/blob/main/docs/img/PBI_Schema_Note_HighRes.jpg)

## 📚 Documentation

- **Home**: https://thibaultschowing.github.io/PBI/
- **Installation**: https://thibaultschowing.github.io/PBI/guides/installation/
- **Story (one-read walkthrough)**: https://thibaultschowing.github.io/PBI/guides/storytelling/
- **Private data handling**: https://thibaultschowing.github.io/PBI/guides/private-data-ingestion/
- **Analysis container**: https://thibaultschowing.github.io/PBI/guides/analysis-guide/

**Check Notebook examples [in the notebooks folder !](https://github.com/ThibaultSchowing/PBI/tree/main/notebooks)**

## 🚀 Quick Start

```bash
git clone https://github.com/ThibaultSchowing/PBI.git
cd PBI

# Configure credentials and set your host UID/GID so containers
# write files as your user instead of root:
cp .env.example .env
echo "UID=$(id -u)" >> .env
echo "GID=$(id -g)" >> .env
# Then edit .env and fill in NCBI_EMAIL (and NCBI_API_KEY if you have one).

# Set up SSH port forwarding first (on your local machine):
# ssh -L 8888:localhost:8888 username@your-server

tmux new -s pbi

docker compose build pipeline
docker compose run --rm pipeline
# ~4 hours for PhageScope data, ~12-18 hours for host genomes retrieval

# analysis container
docker compose build analysis
docker compose up -d analysis


# analysis container (run in a dedicated terminal e.g. tmux session)
docker compose build api
docker compose up api

```

Open `http://localhost:8888` (with SSH tunnel: `ssh -L 8888:localhost:8888 user@server`).

## 🏗️ Infrastructure overview

```text
                  +--------------------------------------+
                  |         shared bind mounts           |
                  | ./private_data  -> /private-data     |
                  | ./pipeline_logs -> /pipeline-logs    |
                  | (pipeline: rw / analysis: ro)        |
                  +-------------------+------------------+
                                      |
         +----------------------------+----------------------------+
         |                                                         |
+--------v---------+      +-----------v-----------+      +--------v---------+
| pipeline         |----->|      pbi-data         |<-----| analysis         |
| (rw: /data,      |      | named volume (/data)  |      | (ro: /data)      |
|      /cache)     |      +-----------+-----------+      +---------+--------+
+--------+---------+                  |                            |
         |                    +-------v-------+          +---------v---------+
         |                    | api           |          | bind mounts       |
         |                    | (ro on /data) |          | [analysis only]   |
         |                    +---------------+          | ./notebooks       |
         |                                               |  -> /workspace(rw)|
+--------v---------+                                     | ./outputs         |
| named volume     |                                     |  -> /results (rw) |
| pbi-cache        |                                     | ./src             |
| -> /cache (rw)   |                                     |  -> /app/src (ro) |
| [pipeline only]  |                                     +-------------------+
+------------------+
```

Additional mounts/volumes currently used in `docker-compose.yml`:

- **Named volumes**
  - `pbi-data` → `/data` (pipeline: rw, analysis/api: ro)
  - `pbi-cache` → `/cache` (pipeline: rw)
- **Bind mounts**
  - `./private_data` → `/private-data` (pipeline: rw, analysis: ro)
  - `./pipeline_logs` → `/pipeline-logs` (pipeline: rw, analysis: ro)
  - `./notebooks` → `/workspace` (analysis: rw)
  - `./outputs` → `/results` (analysis: rw)
  - `./src` → `/app/src` (analysis: ro)

## License

MIT
