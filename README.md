# ABsuperdimers
This repository contains the code used to generate the results in the paper ["Machine learning assisted discovery of exotic criticality in a planar quasicrystal"](https://arxiv.org/abs/2301.11934). It uses the [`rsmine`](https://github.com/RSMI-NE/RSMI-NE) package, which discovers the emergent degrees of freedom of dimers on AB tilings to be *superdimers* on the deflated tiling by solving an information theoretic variational problem.

## Main dependencies
* `rsmine`
* `networkx`
* `tensorflow`
* `numpy`

## Citation
```bibtex
@article{2023arXiv230111934,
  title = {Machine learning assisted discovery of exotic criticality in a planar quasicrystal},
  author = {G\"okmen, Doruk Efe and 
  			Biswas, Sounak and 
  			Huber, Sebastian D. and 
  			Ringel, Zohar and 
  			Flicker, Felix and 
  			Koch-Janusz, Maciej},
  journal = {arXiv e-prints},
  year = {2023},
  month = {Jan},
  eid = {arXiv:2301.11934},
  archivePrefix = {arXiv},
  eprint = {2301.11934},
  primaryClass = {cond-mat.stat-mech}
}
```