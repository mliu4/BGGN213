Class11&12: Structural Bionformatics (Part2) for drug discovery
================

Clean up our protein target structure
-------------------------------------

First we download a target(i.e. receptor) structure from the pdb database. We will pick PDB ID "1hsg"

``` r
#Loading bio3d tools
library(bio3d)
```

``` r
#Getting target structure
pdb.code <- "1hsg"
file.name <- get.pdb(pdb.code)
```

    ## Warning in get.pdb(pdb.code): ./1hsg.pdb exists. Skipping download

``` r
#Reading pdb structure into pdb file.
hiv <- read.pdb(file.name)
```

Extract the protein only segment of this PDB entry and write out a new PDB format file. We will also do the same for the ligand.

``` r
#extracting the protein section only.
protein <- trim(hiv, "protein")
ligand <- trim(hiv, "ligand")

#wrting to new pdb files.
write.pdb(protein, file="1hsg_protein.pdb")
write.pdb(ligand, file="1hsg_ligand.pdb")
```

Convert our docking results for viewing in VMD
----------------------------------------------

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
res
```

    ## 
    ##  Call:  read.pdb(file = "all.pdbqt", multi = TRUE)
    ## 
    ##    Total Models#: 18
    ##      Total Atoms#: 50,  XYZs#: 2700  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 50  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, xyz, calpha, call

``` r
write.pdb(res,file = "results.pdb")
```

Quantitiatively determining the goodness of fit for our ligand results
----------------------------------------------------------------------

``` r
res2 <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("1hsg_ligand.pdbqt")
rmsd(ori, res2)
```

    ##  [1]  0.649  4.206 11.110 10.529  4.840 10.932 10.993  3.655 10.996 11.222
    ## [11] 10.567 10.372 11.019 11.338  8.390  9.063  8.254  8.978
