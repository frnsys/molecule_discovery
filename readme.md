View some results at [matter.farm](http://matter.farm/)

# Setup

```
pip install -r requirements.txt
git submodule update --remote
```

# Usage

1. Download necessary data: `cd data; bash download.sh`
2. Generate compound clusters: `python cluster.py`
3. Generate vocabulary for the JTNN VAE model: `python vocab.py`
4. Train the models:
    - JTNN VAE model
        - First, pre-train: `python train.py jtnn_pre`
        - Then train: `python train.py jtnn`
    - ATC code predictor: `python train.py atc`
    - 3N-MCTS retrosynthesis planner (TODO)
5. Sample compounds from the JTNN VAE model: `python sample.py`
    - Samples the JTNN VAE model
    - Generates synthesis plans for the samples
    - Predicts the ATC codes for the samples

# Visualization

To explore the generated clusters, see `viz/`.
To explore the generated compounds, see `site/`.

# Data sources

See `data/readme.md`.

# JTNN VAE model

The JTNN VAE model is adapted from <https://github.com/wengong-jin/icml18-jtnn>, an implementation of [12] licensed under the MIT license - the changes here add [a conditional VAE version](https://github.com/frnsys/icml18-jtnn/tree/cvae).

# References

1. Botev, Viktor, Kaloyan Marinov, and Florian Schäfer. "Word importance-based similarity of documents metric (WISDM): Fast and scalable document similarity metric for analysis of scientific documents." Proceedings of the 6th International Workshop on Mining Scientific Publications. ACM, 2017
2. Jin, Wengong, Regina Barzilay, and Tommi Jaakkola. "Junction Tree Variational Autoencoder for Molecular Graph Generation." arXiv preprint arXiv:1802.04364 (2018).
3. Kusner, Matt J., Brooks Paige, and José Miguel Hernández-Lobato. "Grammar variational autoencoder." arXiv preprint arXiv:1703.01925 (2017).
4. Goh, Garrett B., Nathan O. Hodas, and Abhinav Vishnu. "Deep learning for computational chemistry." Journal of computational chemistry 38.16 (2017): 1291-1307.
5. Yang, Xiufeng, et al. "ChemTS: an efficient python library for de novo molecular generation." Science and technology of advanced materials 18.1 (2017): 972-976.
6. Liu, Yue, et al. "Materials discovery and design using machine learning." Journal of Materiomics 3.3 (2017): 159-177.
7. Segler, Marwin HS, Mike Preuss, and Mark P. Waller. "Planning chemical syntheses with deep neural networks and symbolic AI." Nature 555.7698 (2018): 604.
8. Kim, Edward, et al. "Virtual screening of inorganic materials synthesis parameters with deep learning." npj Computational Materials 3.1 (2017): 53.
9. Josse, Julie, Jérome Pagès, and François Husson. "Testing the significance of the RV coefficient." Computational Statistics & Data Analysis 53.1 (2008): 82-91.
10. Cordasco, Gennaro, and Luisa Gargano. "Community detection via semi-synchronous label propagation algorithms." Business Applications of Social Network Analysis (BASNA), 2010 IEEE International Workshop on. IEEE, 2010.
11. [Community Detection in Python](https://yoyoinwanderland.github.io/2017/08/08/Community-Detection-in-Python/#No-of-Community-Detection-Algorithms)
12. Jin, Wengong, Regina Barzilay, and Tommi Jaakkola. "Junction Tree Variational Autoencoder for Molecular Graph Generation." arXiv preprint arXiv:1802.04364 (2018).
13. [What are the differences between community detection algorithms in igraph?](https://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph/)
14. [Summary of community detection algorithms in igraph 0.6](https://www.r-bloggers.com/summary-of-community-detection-algorithms-in-igraph-0-6/)
15. Wang, Yong-Cui, et al. "Network predicting drug’s anatomical therapeutic chemical code." Bioinformatics 29.10 (2013): 1317-1324.
16. Liu, Zhongyang, et al. "Similarity-based prediction for Anatomical Therapeutic Chemical classification of drugs by integrating multiple data sources." Bioinformatics 31.11 (2015): 1788-1795.
17. Cheng, Xiang, et al. "iATC-mHyb: a hybrid multi-label classifier for predicting the classification of anatomical therapeutic chemicals." Oncotarget 8.35 (2017): 58494.
18. Szklarczyk D, Santos A, von Mering C, Jensen LJ, Bork P, Kuhn M. STITCH 5: augmenting protein-chemical interaction networks with tissue and affinity data. Nucleic Acids Res. 2016 Jan 4;44(D1):D380-4.
19. Wishart DS, Feunang YD, Guo AC, Lo EJ, Marcu A, Grant JR, Sajed T, Johnson D, Li C, Sayeeda Z, Assempour N, Iynkkaran I, Liu Y, Maciejewski A, Gale N, Wilson A, Chin L, Cummings R, Le D, Pon A, Knox C, Wilson M. DrugBank 5.0: a major update to the DrugBank database for 2018. Nucleic Acids Res. 2017 Nov 8. doi: 10.1093/nar/gkx1037.
20. Kim S, Thiessen PA, Bolton EE, Chen J, Fu G, Gindulyte A, Han L, He J, He S, Shoemaker BA, Wang J, Yu B, Zhang J, Bryant SH. PubChem Substance and Compound databases. Nucleic Acids Res. 2016 Jan 4; 44(D1):D1202-13. Epub 2015 Sep 22 [PubMed PMID: 26400175] doi: 10.1093/nar/gkv951.
21. Gilson, Michael K., et al. "BindingDB in 2015: a public database for medicinal chemistry, computational chemistry and systems pharmacology." Nucleic acids research 44.D1 (2015): D1045-D1053.
22. The UniProt Consortium. UniProt: the universal protein knowledgebase. Nucleic Acids Res. 45: D158-D169 (2017)
23. Gaulton A, Hersey A, Nowotka M, Bento AP, Chambers J, Mendez D, Mutowo P, Atkinson F, Bellis LJ, Cibrián-Uhalte E,
Davies M, Dedman N, Karlsson A, Magariños MP, Overington JP, Papadatos G, Smit I, Leach AR. (2017)
'The ChEMBL database in 2017.' Nucleic Acids Res., 45(D1) D945-D954.
24. Papadatos, George, et al. "SureChEMBL: a large-scale, chemically annotated patent document database." Nucleic acids research 44.D1 (2015): D1220-D1228.
25. Lowe, Daniel Mark. Extraction of chemical structures and reactions from the literature. Diss. University of Cambridge, 2012.
26. Lowe, Daniel (2017): [Chemical reactions from US patents (1976-Sep2016)](https://figshare.com/articles/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873). figshare. Fileset. [CC0 License](https://creativecommons.org/publicdomain/zero/1.0/).
27. Liu, Bowen, et al. "Retrosynthetic reaction prediction using neural sequence-to-sequence models." ACS central science 3.10 (2017): 1103-1113.
28. Klucznik, Tomasz, et al. "Efficient syntheses of diverse, medicinally relevant targets planned by computer and executed in the laboratory." Chem 4.3 (2018): 522-532.
29. Segler, Marwin HS, and Mark P. Waller. "Neural‐Symbolic Machine Learning for Retrosynthesis and Reaction Prediction." Chemistry–A European Journal 23.25 (2017): 5966-5971.
30. Law, James, et al. "Route designer: a retrosynthetic analysis tool utilizing automated retrosynthetic rule generation." Journal of chemical information and modeling 49.3 (2009): 593-602.
31. Schwaller, Philippe, et al. "“Found in Translation”: predicting outcomes of complex organic chemistry reactions using neural sequence-to-sequence models." Chemical science 9.28 (2018): 6091-6098.
32. Kipf, Thomas N., and Max Welling. "Semi-supervised classification with graph convolutional networks." arXiv preprint arXiv:1609.02907 (2016).
33. Yang, Zhilin, William W. Cohen, and Ruslan Salakhutdinov. "Revisiting semi-supervised learning with graph embeddings." arXiv preprint arXiv:1603.08861 (2016).
34. Kipf, Thomas, et al. "Neural relational inference for interacting systems." arXiv preprint arXiv:1802.04687 (2018).
35. Wei, Jennifer N., David Duvenaud, and Alán Aspuru-Guzik. "Neural networks for the prediction of organic chemistry reactions." ACS central science 2.10 (2016): 725-732.
36. Coley, Connor W., et al. "Prediction of organic reaction outcomes using machine learning." ACS central science 3.5 (2017): 434-443.
37. Gupta, Anvita. "Predicting Chemical Reaction Type and Reaction Products with Recurrent Neural Networks."
38. Plehiers, Pieter P., et al. "Automated reaction database and reaction network analysis: extraction of reaction templates using cheminformatics." Journal of cheminformatics 10.1 (2018): 11.

---

# RDKit setup

Instructions for a pyenv-virtualenv

Preparation:

```
PYENV_PYTHON=~/.pyenv/versions/3.6.6
PYTHON_EXECUTABLE=$PYENV_PYTHON/bin/python
PYTHON_LIBRARY=$PYENV_PYTHON/lib/libpython3.6m.a
PYTHON_INCLUDE_DIR=$PYENV_PYTHON/include/python3.6m/
pip install numpy # do this with 3.6.6 active
pyenv activate data
```

Get and install Boost 1.58:

```
wget http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.bz2
tar jxvf boost_1_58_0.tar.bz2
cd boost_1_58_0
export CPLUS_INCLUDE_PATH="$CPLUS_INCLUDE_PATH:$PYTHON_INCLUDE_DIR"
./bootstrap.sh --with-libraries=python,serialization --with-python-root=$PYENV_PYTHON --prefix=$PYENV_PYTHON
./b2 link=shared install
cd ..
```

Get and install RDKit:
(note: build this with a `pyenv` `virtualenv` activated to build for that environment)

```
# http://www.rdkit.org/docs/Install.html
RDKIT=Release_2018_03_4
sudo apt install libeigen3-dev libsqlite3-dev libpython3-dev libcairo2-dev
wget https://github.com/rdkit/rdkit/archive/${RDKIT}.tar.gz
tar xzvf ${RDKIT}.tar.gz
cd rdkit-$RDKIT
mkdir build
cd build
cmake -D RDK_BUILD_CAIRO_SUPPORT=ON -D PYTHON_LIBRARY=$PYTHON_LIBRARY -D PYTHON_EXECUTABLE=$PYTHON_EXECUTABLE $PYTHON_INCLUDE_DIR=$PYTHON_INCLUDE_DIR -D BOOST_ROOT=$PYENV_PYTHON -D Boost_NO_SYSTEM_PATHS=ON ..
make -j4
make install
```

Install components:

```
# Hacky, there has to be an easier way to specify where RDKit is built to?
cp -r ../rdkit $PYENV_PYTHON/envs/data/lib/python3.6/site-packages/rdkit

# We could copy to the pyenv lib folder but
# would require settings LD_LIBRARY_PATH for each new terminal session
# cp -r ../lib/libRDKit* $PYENV_PYTHON/lib/
# export LD_LIBRARY_PATH=$PYENV_PYTHON/lib:$LD_LIBRARY_PATH
# Alternatively, easier to just copy to /usr/lib which
# is where Python automatically looks for shared libraries
sudo cp -r ../lib/libRDKit* /usr/local/lib/
sudo cp -r $PYENV_PYTHON/lib/libboost* /usr/local/lib/

for f in /usr/local/lib/libRD*; do
    t=${f##*/}
    sudo ln -s $f /usr/lib/$t
done
for f in /usr/local/lib/libboost*; do
    t=${f##*/}
    sudo ln -s $f /usr/lib/$t
done

```

# Tensorflow setup

```
# With your virtualenv activated
cd /tmp
sudo apt install unzip
wget https://github.com/bazelbuild/bazel/releases/download/0.16.0/bazel-0.16.0-installer-linux-x86_64.sh
chmod +x bazel-0.16.0-installer-linux-x86_64.sh
sudo ./bazel-0.16.0-installer-linux-x86_64.sh
pip install -U keras wheel

git clone https://github.com/tensorflow/tensorflow
cd tensorflow
./configure

# optionally add --config=cuda for CUDA
bazel build -c opt --copt=-mavx512f --copt=-mavx --copt=-mavx2 --copt=-mfma --copt=-mfpmath=both -k //tensorflow/tools/pip_package:build_pip_package
bazel-bin/tensorflow/tools/pip_package/build_pip_package /tmp/tensorflow_pkg
pip install /tmp/tensorflow_pkg/tensorflow-*-linux_x86_64.whl
```
