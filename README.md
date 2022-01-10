# UniteM

[![version status](https://img.shields.io/pypi/v/unitem.svg)](https://pypi.python.org/pypi/unitem)
[![downloads](https://img.shields.io/pypi/dm/unitem.svg)](https://pypi.python.org/pypi/unitem)

UniteM is a software toolkit implementing different ensemble binning strategies for producing a non-redundant set of bins from the output of multiple binning methods.

## Major version history

- **1.1.0**: UniteM updated to include HMMs from closely related protein families to ensure accurate annotation of marker genes
- **>=1.0.0**: removed CheckM; genome quality estimated using marker genes identified across [GTDB](https://gtdb.ecogenomic.org/) species representatives using [Prodigal](https://github.com/hyattpd/Prodigal) and [HMMER](http://hmmer.org/); see [UniteM manual](https://github.com/dparks1134/UniteM/blob/master/unitem_manual.pdf) for details
- **<1.0.0**: used CheckM and CheckM marker sets to evaluate the quality of genomes during ensembly binning; used in the [draft UniteM](https://github.com/dparks1134/UniteM/blob/master/unitem_ms.draft.pdf) manuscript

## Installation

```
> pip install unitem
```

## Quick Start

The functionality provided by UniteM can be accessed through the help menu:
```
> unitem -h
```

Usage information about specific functions can also be accessed through the help menu, e.g.:
```
> unitem bin –h
```

Detailed information regarding the use of UniteM can be found in the User's Guide (unitem_manual.pdf).


## Cite

If you find this package useful, please cite this git repository (https://github.com/dparks1134/UniteM)


## Copyright

Copyright © 2017 Donovan Parks. See LICENSE for further details.
