# pbi
Microbiome analysis for phosphate-defense interaction

The code in this repository was used to analyze the microbial community assembled on roots of
plants growing in a wild soil (census) or on plates with a synthetic community (SynCom).

## Directories

The `census` directory contains code used to analyze the microbial community composition of
plants gorwing in wild Mason Farm soil.

The `syncom` directory contains code used to analyze the microbial community composition,
and the plant phenotypes, of plants growing on agar plates with a well-defined but complex
synthetic bacterial community, under different phosphate conditions.

The `misc` directory contains intermediate files that have been made available after
publication upon request. The `misc/unifrac` subdirectory contains UniFrac distance
matrices as calculated in QIIME/1.5.0 (described in [Castrillo *et al*. 2017](http://www.nature.com/doifinder/10.1038/nature21417)).

## Data

The raw sequencing reads that underlie the data in these analysis are publicly available in
the European Nucleotide Archive's Short Read Archive under study accession: [PRJEB15671](http://www.ebi.ac.uk/ena/data/view/PRJEB15671).


## Functions files

Most scripts load one or two different functions.r files. Those two files are located at the main
directory.

The file that is referenced as "~/rhizogenomics/src/trunk/phosphate_code/functions.r" is
simply named `functions.r`.

The file that is referenced as "~/rhizogenomics/src/trunk/immune_code/PRR_propep/functions.r"
is named `functions_prrpropep.r`.

## Referencing

If you use this code or the data associated with it please cite:

Castrillo G, Teixeira PJPL, Herrera Paredes S, Law TF, de Lorenzo L, Feltcher ME, *et al*.
"Root microbiota drive direct integration of phosphate stress and immunity" (2017). *Nature*
543(7646):513–8. Available from: http://www.nature.com/doifinder/10.1038/nature21417

A BibTex entry for this publication is:

```tex
@article{Castrillo2017,
author = {Castrillo, Gabriel and Teixeira, Paulo Jos{\'{e}} Pereira Lima and {Herrera Paredes}, Sur and Law, Theresa F. and de Lorenzo, Laura and Feltcher, Meghan E. and Finkel, Omri M. and Breakfield, Natalie W. and Mieczkowski, Piotr and Jones, Corbin D. and Paz-Ares, Javier and Dangl, Jeffery L},
doi = {10.1038/nature21417},
issn = {0028-0836},
journal = {Nature},
month = {mar},
title = {{Root microbiota drive direct integration of phosphate stress and immunity}},
url = {http://www.nature.com/doifinder/10.1038/nature21417},
year = {2017}
}
```

## Copyright

(C) Copyright 2016 Sur Herrera Paredes

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
