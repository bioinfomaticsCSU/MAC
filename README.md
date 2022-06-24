# MAC2.0

Merging assemblies by using adjacency algebraic model and classification.

## About Updating

MAC2.0 uses C++ to refactor the main contents of the original MAC, fixing errors in dealing large files.
MAC2.0 removes the step of evalution, which does't need raw reads as input any more.

## Pre-installed:

MUMmer: http://mummer.sourceforge.net/ 

Please make sure the path of MUMmer has been added into system variable.

## Usage: 

`g++ MAC2.0.cpp -o MAC2.0`

`MAC2.0 <contig_query.fa> <contig_reference.fa>`


## Notification:

- The input files need to be placed in the ./input folder.

- The output file will be output to the ./output folder.

- Because of some implementation issues, the MAs may be a little more than the original approach.

- If you have multiple contig files (over 3) to merge, please try to merge any two of them first, then merge the result with the other files iteratively.
