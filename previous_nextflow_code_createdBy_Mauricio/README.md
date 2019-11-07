# VD_pipeline

New version of the Viral Discovery pipeline from B. Andersson's lab @ Karolinska Institutet.  Implemented using [Nextflow](https://www.nextflow.io/)

The pipeline is split in two parts:

1) Preprocessing:  Does QC on the sequences, removes adapters (including SISPA adapters from our viral discovery protocol)  and removes human by mapping to a reference.

2) Discovery:  Search the sequences with as much stuff as we can. Hope we can develop a way of nicely summarizing the results so that the information is easy to use and hopefully discover new stuff !

## Install

In theory, it should be as _simple_ (lol) as installing all the dependencies, and creating a suitable config file with the paths to them.

TODO: In practice, automate this in some way! maybe conda?

## More information

Check the docs/ folder!
