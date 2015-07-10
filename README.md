#NetView P
###v.0.7

NetView P is an implementation of the graph-based population structure visualization NetView (Neuditschko et al. 2012) in Python. The pipeline generates networks based on a distance matrix derived from genome-wide SNPs using a minimum spanning tree and mutual k-NN. Networks can visualize large- and fine-scale patterns of population structure, including individual and family-level relationships in natural and captive populations.

NetView P now uses scientific computing packages (`numpy` `scipy` `scikit-learn`) to provide additional configurations and operate efficiently on large data sets. Installation of the appropriate environment and configuration of the pipeline are outlined below. The project has moved to Github to facilitate access to the source code and implement improvements from the community. If you find a bug or have any other questions, feel free to send us a message or use the issues function on Github.

###Table of Contents

- [NetView P](#)
	- [Installation](#)
		- [PLINK](#)
		- [Cytoscape](#)
	- [Input](#)
	- [Configuration](#)
		- [Required](#)
		- [Options](#)
	- [Outputs](#)
	- [Examples](#)
	- [Considerations](#)
	- [Future Implementations](#)
	- [Updates](#)
	- [Citations:](#)
	- [Contact](#)

## Installation  

We recommend installing [Anaconda](http://continuum.io/downloads#py34) for Python 3.4. The distribution includes all packages used in Netview P and runs under both Linux and Windows.

###**PLINK**

Quality Control and ASD are computed in [PLINK v1.07](http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml).

Download the program to your home directory (or elsewhere) and update PATH:

**Linux**

```
cd $HOME
wget -O $HOME/plink.zip "http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-x86_64.zip"
unzip plink.zip
echo 'export PATH=$PATH:$HOME/plink-1.07-x86_64'  >> ~/.profile
```

**Windows**

```
System -> Advanced Sytem Settings -> Advanced -> Environmental Variables

Add your directory to systems variable 'Path' (separated by a semicolon): 

... ;C:\Users\UserName\plink-1.07-x86_64

```
###**Cytoscape**

Visualizations can be constructed in most graph visualization software. We frequently use the open-source platform [Cytoscape v.3.2](http://www.cytoscape.org/download.php). Visualization and community clustering will be available in the next release of NetView P. 

## Input

**Data:**

* Files for [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml) (`.ped` `.map`) 
* SNP (N x SNP) matrix
* Distance (N x N) matrix


**Node Attributes:**

Data attributes (`.csv`) ordered by N as in Data. One row per sample, column headers, first column contains IDs.
Node attributes can be any of the following:

* Node descriptor string (ID, Population, Location...)
* Node [shapes](http://js.cytoscape.org/#style/node-body) for Cytoscape 
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
--quality     Switch on quality control in PLINK (default: OFF):

--mind        Filter samples with missing data > (default: 0.1)
--maf         Filter SNPs with minor allele frequency < (default: 0.01)
--geno        Filter SNPs with missing data > (default: 0.1)
--hwe         Filter SNPs failing HWE at P > (default: 0.001)

--distance    Distance calculated from SNP data (default: asd):
              asd, hamming, correlation...

--algorithm   Choice of algorithm for nearest neighbour search (default: brute):
              auto, kd_tree, ball_tree, brute
              
--start       Start pipeline at k = (default: 10)
--stop        Stop pipeline at k = (default: 40)
--step        Step pipeline at k = (default: 10)

--mst-off     Switch off minimum spanning tree (default: ON)

--project     Name of project (default: timestamp)
--prefix      Prefix for output files (default: project)
--sep         Delimiter of input files (default: \t)

--missing     Set missing character (default: 0)
--ploidy      Set ploidy of input data (default: diploid):
              haploid, diploid
              
--visual      Generate node attribute files (skips NetView and QC)
--off         Switch off Netview and run only QC (requires --quality)
--dist        Calculate distance matrix only (if --quality, after QC)

```
## Outputs

#### **/networks:**

Network files formatted as weighted edges:
`.edges`

#### **/nodes:**

Node attribute files for networks:
`.nat`
#### **/plink:**

Data files after QC: 
`_plink_qc.ped` `_plink_qc.map`

Information on excluded samples and SNPs: 
`_plink_qc.log`

Excluded sample IDs: 
`_plink_qc.irem`

ASD (1-IBS) distance matrix (if enabled, after QC): 
`_plink.mdist`

#### **/other:**

Updated node attributes (if samples were excluded in QC): 
`_qc.csv`

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

##Considerations

**Data**

NetView is designed for large data sets (n > 100-1k, SNPs > 10k-100k, see Neuditschko et al. 2012)  and has been tested on smaller data sets (n > 82, SNPs > 300). Care needs to be taken for heterogeneous population sizes and populations with n < 10 (see Neuditschko et al. 2012). A stepwise reduction in *k* (usually *k* < 10) can sometimes recover resolution for such datasets. Very small data sets should be avoided as the effect of a small number of individuals or populations has not been properly assessed.

**Distance Metrics**

The choice of genetic distance measure is crucial to the analysis, since it determines which nodes are connected by the minimum spanning tree and mk-NN. The default distance is the allele sharing distance (1 - Identity by Similarity) (`--distance-matrix` in PLINK). This does not work for haploid genotypes or other non-genetic data, for instance when constructing networks for ecological or biogeographical data. Hamming distance or any of the [supported metrics](http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.spatial.distance.pdist.html) in `pdist` from `scipy` can be calculated instead. A pre-computed distance matrix can be used as input with `-m`. In this case, the pipeline will use the matrix directly for NetView. Future updates will implement additional genetic distance measures.

**Nearest Neighbour**

The user-defined parameter *k* determines the number of mutual nearest neighbours and is essential to the resolution of the networks. There is currently no optimisation for *k* (Neuditschko et al. 2012). However the effect of the parameter on the network topologies offers a possibilty to investigate the networks at different levels of genetic similarity. A small value of *k* focuses on fine-scale structure (connecting fewer, genetically similar individuals), whereas a large value of *k* focuses on large-scale patterns of similarity (connecting more, genetically distant individuals) (see Neuditscho et al. 2012). 

We therefore highly recommend to explore the data in a range of *k*, although a working empirical value has been suggested at *k* = 10 (see Neuditscho et al. 2012). It should be noted that the application of a simple machine learning algorithm such as the mutual nearest neighbour search is not rooted in established population models.

Netview P implements the nearest neighbour search from `scikit-learn`. The choice of the algorithm can be relevant for large data sets and potentially affect the final network topology. The default mode in the pipeline is the brute force algorithm, for a detailed description and choice of K-D Tree or Ball Tree, an excellent [documentation] (http://scikit-learn.org/stable/modules/neighbors.html) is available for `NearestNeighbors`.

**Minimum Spanning Tree**

The inclusion of a minimum spanning tree is the default mode of NetView P. This is the standard version of NetView implemented by Neuditschko et al. (2012). However, there are some caveats to consider. Older versions (< v.0.7) implemented Prims, wheras later versions now implement `minimum_spanning_tree` from `scipy` which uses Kruskal. Neither algorithm is necessarily deterministic when not all the edge weights in the graph are different and can therefore result in differences in the final network topology.

In a connected network, the edges of the minimum spanning tree that connect communities can drag nodes away from these communities. This is only the case for some visualizations, for instance when using the organic layout in Cytoscape. The edge files contain one column of colour attributes so that edges inluded from the minimum spanning tree can be visualized in the networks and their effect on the network topology considered. 

The minimum spanning tree can be switched off to yield only communities connected at a particular *k* by mk-NN. This can result in the exclusion of nodes that are not considered to be mutual nearest neighbours to any other node in the graph - these are visualized as unconnected singletons.

**Visualization**

The choice of the visualization algorithm can strongly influence the final representation of the network topology. The layout should therefore be consistent between visualizations at different *k* and in comparative network analyses. We generally use the circular or organic layouts provided by yFiles in Cytoscape as recommended by Neuditschko et al. (2012).

##Future Implementations

* Sun Grid Engine support for batch job submission
* Additional genetic distance measures with `scipy.spatial.distance.pdist`
* Summary statistics module
* `BioPython` and [`Bio.PopGen.Genepop`](http://biopython.org/wiki/PopGen_Genepop)
* PCA and [`Admixture`](https://www.genetics.ucla.edu/software/admixture/) for [PCA-InformativeIndividuals](https://livestockgenomics.files.wordpress.com/2014/04/markus-talk.pdf) (PCA-IIs)
* Network visualization and analysis with [`cyREST`](https://github.com/idekerlab/cyREST/wiki) and [`iGraph`](http://igraph.org/python/)
* Interactive data visualization using `MongoDB` and `D3` (Linux)


##Updates

NetView P is under development and is subject to changes.

**v.0.7**

The update changes the architecture of the program considerably. NetView P is now using scientific computing packages and two main components, a data storage class (`Data`) and an analysis class (`Analysis`). These changes have significantly decreased runtime, but the multiprocessing module is still implemented for each iteration of mk-NN.

Data structures processed with:

* `numpy`
* `scipy`
* `scikit-learn`


New options:

```
--algorithm   Choice of algorithm for nearest neighbour search (default: brute):
              auto, kd_tree, ball_tree, brute
--project     Name of project (default: timestamp)
--prefix      Prefix for output files (default: project)
--sep         Delimiter of input files (default: \t)
--missing     Set missing character (default: 0)
--ploidy      Set ploidy of input data (default: diploid):
              haploid, diploid
--off         Switch off Netview and run only QC (requires --quality)
--dist        Calculate distance matrix only (if --quality, after QC)
```

Input:

* Attribute file is now required and must be comma-delimited (`.csv`)
* SNP matrix does no longer require population and individual columns (N x SNP).


Errors:

We found a bug in the calculation of Prims algorithm, which may produce slightly different topologies in older versions (< v.0.7). Since the algorithm is no longer supported and replaced by Kruskal, we recommend using the latest version of NetView P.


##Citations

Steinig et al. (2015) - NetView P: A network visualization tool to unravel complex population structure using genome-wide SNPs, Molecular Ecology Resources

[Neuditschko et al. (2012)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0048375) - NetView: A High-Definition Network-Visualization Approach to Detect Fine-Scale Population Structures from Genome-Wide Patterns of Variation, PLoS One

##Contact

<eikejoachim.steinig@my.jcu.edu.au>
