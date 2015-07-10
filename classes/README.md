#Classes in Netview P

###`Data`


#### Attributes

`prefix`    prefix for output files       `str`

`ploidy`    ploidy of data                `str`

`missing`   missing character             `str`

`n`         number of samples             `int`

`nSNP`      number of snps                `int`         stored as number of snp columns regardless of ploidy

`ids`       list of ids                   `list`

`alleles`   list of unique snp alleles    `list`

`biodata`   list of bipython objects      `list`

`snps`      array of n x snps             `ndarray`

`matrix`    matrix                        `ndarray`     current matrix

`meta_data` dictionary of metadata        `dict`        stores metadata lists of length n

`snp_data`  dictionary of snp data        `dict`        stores snp data lists of length n

`matrices`  dictionary of matrices        `dict`        stores distance matrices by metric

`networks`  dictionary of networks        `dict`        stores networks by k and run (netview_k*k*_*run*)

`netview_runs`            accessory       `int`

`filetype`                accessory       `str`

#### Methods

**`readData(file, f, sep='\t', header=False, add_col=0)`**

`file`      in file name                `str`

`f`         file type                   `str`           plink, nexus, raxml, snp_mat, matrix, attributes

`sep`       delimiter                   `str`

`header`    header in matrix            `bool`

`add_col`   initial columns in matrix   `int`


**`writeData(file, f, sep='\t')`**

`file`      out file name                 `str`

`f`         file type                     `str`         plink, nexus, raxml, matrix, attributes, meta, snp

`sep`       delimiter                     `str`

###`Analysis`

Analysis class for distance calculation, NetView and PLINK.

**Attributes**

`data`      data class object             `Data`

**Methods**

**`getDistance(target='snps', distance='hamming')`**

`target`    targets data.snps or .matrix  `str`         matrix, snps

`distance`  distance to be calculated     `str`         asd (plink), pdist distances

**`runPLINK(qc_parameters={}, commandstring='', asd=False, quality=False)`**

`qc_parameters`  dictionary of qc parameters  `dict`    --mind, --maf, --geno, --hwe

`commandstring`  commands for plink           `str`     overwrites asd and qc

`asd`            calculate 1-IBS              `bool`    

`quality`        run  quality control         `bool`

**`runNetView(tree=True, start=10, stop=40, step=10, algorithm='auto')`**

`tree`           include mst edges            `bool`

`start`          start iterations k           `int`     netview operates on data.matrix

`stop`           stop iterations k            `int`

`step`           iterations at intervals      `int`

`algorithm`      nearest neighbour search     `str`     auto, brute, kd_tree, ball_tree


**`updateNodeAttributes(attribute_file)`**

`attribute_file`  file containing node attributes `str`   updates file and data after qc with plink

###`CommandLine`

Builds the command line parser and runs NetView P.
