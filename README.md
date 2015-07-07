#NetView P
##v.0.7

NetView P is an implementation of the graph-based population structure visualization NetView (Neuditschko et al. 2012) in Python. The pipeline generates networks based on a distance matrix from genome-wide SNPs using a minimum spanning tree and mutual kNN. Networks can visualize large- and fine-scale patterns of population structure, including individual and family-level relationships in natural and captive populations.

NetView P is under development and now uses scientific computing packages (`numpy, scipy, scikit-learn`) to provide additional configurations and operate efficiently on large data sets. Please note that the use of these packages (in particular the choice of nearest neighbour algorithm) can produce slightly different network topologies than previous versions of NetView P. Installation of the appropriate environment and configuration of the pipeline are outlined below, with more detailed descriptions available in the manual and walktrough. The project has moved to GitHub to facilitate access to the source code and encourage improvements from the community.

## Installation

The pipeline was tested in Python 3.4 on Ubuntu 14.04. Please see the manual for installation on Windows.

We highly recommend installing [Anaconda](http://continuum.io/downloads#py34) for Python 3.4. The distribution includes packages used in Netview P and runs under both Windows and Linux. Once you have downloaded Anaconda:
```
bash Anaconda3-2.3.0-Linux-x86_64.sh
```

## Dependencies

**1. PLINK**

Quality Control and ASD are computed in [PLINK v1.07](http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml). Download the program to your home directory (or elsewhere) and update PATH:

```
cd $HOME
wget -O $HOME/plink.zip "http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-x86_64.zip"
unzip plink.zip
echo 'export PATH=$PATH:$HOME/plink-1.07-x86_64'  >> ~/.profile
```

**2. Cytoscape**

Graph visualization and analysis within the pipeline will be included in the next release. Visualizations can be done manually in [Cytoscape v3](http://www.cytoscape.org/download.php).

## Input

**Data:**

* Files for [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml) (.ped/.map) 
* SNP (N x SNP) Matrix
* Distance (N x N) Matrix.

**Node Attributes:**

Data attributes (**.csv**) ordered by N as in Data. One row per sample, column headers, first column contains IDs.
Node attributes can be any of the following:

* Node descriptor (ID, Population, Location...)
* Node [shapes](http://js.cytoscape.org/#style/node-body) for Cytoscape: 
* Node [colour names](http://www.w3schools.com/html/html_colornames.asp) for Cytoscape
              
## Configuration

##### Required
```
-f            Name of input file, format can be:

-p            PLINK (.ped/.map)
-s            SNP matrix (N x SNPs)
-m            Distance Matrix (N x N)

-a            Name of node attribute file (.csv), column headers, same order and 
              first column ID as IDs in PLINK / SNP matrix
```
##### Options
```
--quality     Switch on Quality control in PLINK (default: OFF):

--mind        Filter samples with missing data > (default: 0.1)
--maf         Filter SNPs with minor allele frequency < (default: 0.01)
--geno        Filter SNPs with missing data > (default: 0.1)
--hwe         Filter SNPs failing HWE at P > (default: 0.001)

--distance    Distance calculated from SNP data (default: asd):
              asd, hamming

--algorithm   Choice of algorithm for nearest neighbour search (default: auto):
              auto, kd_tree, ball_tree, brute
              
--start       Start pipeline at k = (default: 10)
--stop        Stop pipeline at k = (default: 40)
--step        Step pipeline at k = (default: 10)

--mst-off     Switch off minimum spanning tree (default: ON)

--project     Name of project, otherwise time stamped directory for outputs
--prefix      Prefix for output files (default: project)
--sep         Delimiter of input files (default: \t)

--missing     Set missing character (default: 0)
--ploidy      Set ploidy of input data (default: diploid):
              haploid, diploid
--visual      Generate node attribute files, skips pipeline

```
## Outputs

##### **`project/networks`:**

`.edges`

Network files formatted as weighted edges.

`.nodes`

Node attribute files for networks.

##### **`project/plink`:**

`_plink_qc.ped/.map`

Data files after QC.

`_plink_qc.log`

Information on excluded samples and SNPs.

`_plink_qc.irem`

Excluded sample IDs.

`_plink.mdist`

ASD (1-IBS) distance matrix (if enabled, after QC).

##### **`project/other`:**

`_qc.csv`

Updated node attributes (if samples were excluded in QC)

## Examples

Run default NetView P with (tab-delimited) input file for PLINK:

```
netview.py -f oyster.ped -a oyster.csv -p
```

Set project and prefix:

```
netview.py -f oyster.ped -a oyster.csv -p --project oyster --prefix test
```

Run with default quality control in new project directory:
```
netview.py -f oyster.ped -a oyster.csv -p --project oyster_qc --prefix test_qc 
           --quality
```

Run with quality control, but different parameters for MIND and MAF:
```
netview.py -f oyster.ped -a oyster.csv -p --project oyster_qc --prefix test_qc 
           --quality --mind 0.2 --maf 0.05
```

Run without quality control, use different algorithm for nearest neighbour search and stop at k=70:
```
netview.py -f oyster.ped -a oyster.csv -p --project oyster --prefix test
           --algorithm kd_tree --stop 70
```

Run without minimum spanning tree and use different distance measure (Hamming):
```
netview.py -f oyster.ped -a oyster.csv -p --project oyster --prefix test
           --algorithm kd_tree --stop 70 --mst-off --distance hamming
```

Generate node attribute files only (`--visual` overrides all other options, requires `-a`):
```
netview.py -f oyster.ped -a oyster.csv -p --project oyster --prefix test
           --algorithm kd_tree --stop 70 --visual
```

Run haploid data (no PLINK), with missing character (NA):

```
netview.py -f staph.snps -a staph.csv -s --project staph --prefix staph
           --ploidy haploid --distance hamming --missing NA
```

Run pre-calculated distance matrix with brute force nearest neighbours search:

```
netview.py -f oyster.dist -a oyster.csv -m --project oyster_matrix --prefix oyster_matrix
           --algorithm brute
```

##Citations:

Steinig et al. (2015) - NetView P: A network visualization tool to unravel complex population structure using genome-wide SNPs, Molecular Ecology Resources

Neuditschko et al. (2012) - NetView: A High-Definition Network-Visualization Approach to Detect Fine-Scale Population Structures from Genome-Wide Patterns of Variation, PLoS One
