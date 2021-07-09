# Analysis scripts for Cox and Haebig (2021)
**Child-Oriented Word Associations Improve Models of Early Word Learning**
*Christopher R. Cox and Eileen K. Haebig*

## What is this work about?
Over the past decade, researches have begun to model the development of early vocabularies as network that acquire new nodes (words) over time. For these models to accurately predict word learning, the connectivity of the network should reflect the child's environment. In this work, we use a word assocation task with modified instructions to inform network connectivity. This network predicts word learning between 16 and 30 months better than a standard word association task.

## How to use these scripts
The scripts can be run in order to replicate the major analyses of the paper. In addition to `igraph` and packages in the `tidyverse`, some of the scripts rely on new packages written to support this research:

 * [netbuildr](https://github.com/crcox/netbuildr)
 * [netgrowr](https://github.com/crcox/netgrowr)

These can be installed directly from GitHub using the `devtools` package:

```
devtools::install_github("crcox/netbuildr")
devtools::install_github("crcox/netgrowr")
```
