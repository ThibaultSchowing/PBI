# PBI-Scope for PERPHECT

![alt](https://github.com/ThibaultSchowing/PBI/blob/main/docs/img/PBI_Schema_Note.png)

# PBI-Scope
## Dockerized Phage Bacteria Interactions toolkit based on PhageScope

> A proof-of-concept dockerized bioinformatics pipeline that makes phage genomic data from [PhageScope](https://phagescope.deepomics.org/database) and their hosts available in an efficient, structured format for training neural networks and AI models for phage-host interaction prediction. 

**Install - Wait - Work** The pipeline takes care of everything within Docker !

[![Documentation](https://img.shields.io/badge/docs-github%20pages-blue)](https://thibaultschowing.github.io/PBI-Scope/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo-blue.svg)](https://doi.org/10.5281/zenodo.21469490)
[![CI Pipeline and DB Tests](https://github.com/ThibaultSchowing/PBI-Scope/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ThibaultSchowing/PBI-Scope/actions/workflows/ci.yml)
[![ExPASy SIB](https://img.shields.io/badge/ExPASy-SIB_Resource-E2001A)](https://www.expasy.org/resources/pbi-scope)

> Real-life usage example and guide on how to use PBI-Scope for Phage Bacteria Interaction model training. 


1. Preprocessing and integrating you private data
2. Run the pipeline to include these data into the PBI-Scope database
3. Use the _pbi_ Python package to generate datasets or stream data into your model training. 

## Workflow

You can first find the PBI-Scope readme [here](https://github.com/CI4CB-lab/PBI-Scope-PERPHECT/blob/main/README.md) for all information about PBI-Scope itself. 

### Integrate your private data

First we need to include our own data which will be stored within the 'private_data/' directory under a 'source' folder. Follow the given examples format and don't forget to remove the examples before running your analysis !



> Real-life usage example and guide on how to use PBI-Scope for Phage Bacteria Interaction model training. 


1. Preprocessing and integrating you private data
2. Run the pipeline to include these data into the PBI-Scope database
3. Use the _pbi_ Python package to generate datasets or stream data into your model training. 

## Workflow

You can first find the PBI-Scope readme [here](https://github.com/CI4CB-lab/PBI-Scope-PERPHECT/blob/main/README.md) for all information about PBI-Scope itself. 

### Integrate your private data

First we need to include our own data which will be stored within the 'private_data/' directory under a 'source' folder. Follow the given examples format and don't forget to remove the examples before running your analysis !

Here we do not show our data however we can give a small piece of advice: 
    - If your data are fixed in time and not evolving, you can pre-process it and upload them in the _private\_data_ folder.
    - If your data are regularly changing, it would be a good idea to implement a data version control with, for instance, DVC in order to pull your most recent data easilly. 



