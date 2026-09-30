# Single read recombination detection
> TODO

## Why?
No one seems to care about building decent, isolated CLI tools...
* Porsborg et al.
    * [Code](https://github.com/PeterSoerud/recombination_calling/blob/main/scripts/calling/recombination.py#L5-32) w/27 required arguments and no documentation/defaults

* Schweiger et al.
    * Integrated in Snakemake workflow
        * In run [directives](https://github.com/regevs/recombination/blob/main/snakefiles/read_analysis.snk#L50) so if things crash, you're S.O.L. debugging
        * No dependencies pinned
        * Assumes human and has extraneous alignment [steps](https://github.com/regevs/recombination/blob/main/snakefiles/read_analysis.snk#L648)

I don't want to deal with this so better to start from scratch.

## Sources
* https://www.nature.com/articles/s41467-025-65248-3
* https://www.nature.com/articles/s41586-026-10901-0
