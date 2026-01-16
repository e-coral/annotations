# Genomic features annotations

Package for annotating various BED-like formats of genomic regions of interest, with information including:

* Genes (names, positions, lengths)
* Fragile sites
* Repeats
* Distance between region and telomere

## Disclaimer
This is an in-progress project, developed to accept many different input types, and annotate genomic loci with various features of interest within our lab. It is not currently well standardised, fully comprehensive, or available through pypi/in a packaged format. Questions about using any part of this project are welcome, via the issues/discussion features on this page, as appropriate. However, there are no guarantees this can or will work for all annotation purposes.

## Using the package
### A basic example

```
from annotate import annotate   # point to wherever you have installed your copy of the annotate package within this import statement

infile = your_file.csv  # full path to your input file
outfile = your_annotated_file.csv   # full path for your output file

annotate.annotate_standard_csv_input_file(infile, outfile)  # more options and settings available in example scripts
```

### Other examples
Some example annotation scripts are available in [example_annotation_scripts](example_annotation_scripts/README). Depending on the input format of your data and the required output columns, the example scripts in example_annotation_scripts may be used to annotate input files with different types of information, as required. Otherwise, the example scripts demonstrate how to call various functions and settings of the annotations module in order to achieve various outputs, and so can be used as a basis for generating similar, short annotation scripts specific to other use cases.

## Other information
[ref_files](https://github.com/e-coral/annotations/tree/master/src/annotate/ref_files) contains the reference files required for annotation, except for the repeatmasker reference, which is too large to be hosted here. It can be downloaded from [the T2T GitHub page](https://github.com/marbl/CHM13) (> repeat annotation > repeatmasker), and placed as-is into your copy of the ref_files directory, as 'chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed'.

See the README files in relevant directories for specific information regarding collection and generation of annotation reference resources.

## Licence
MIT (see LICENCE.txt)


[//]: # (<h2>Installation</h2>)

[//]: # (<p>download the wheel from ___, then do:</p>)

[//]: # ()
[//]: # (        pip install <full filepath to downloaded wheel>)