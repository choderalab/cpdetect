cpDetect
========

Bayesian change point detection

Installation
-------------
Install directly from source directory.
`python setup.py install`

Requirement
-------------
* python 2.7
* NumPy
* Pandas
* SciPy

Usage
-----

To run `cpDetect`, the timeseries data needs to be a list of 1-D numpy arrays. They do not have to be of the same size
First, instantiate the detector. Choose the underlying distribution (normal or log normal) and the log odds threshold
(default is 0). 

```
detector = cpDetector(trajs, distribution='log_normal', log_odds_threshold=0)
detector.detect_cp()
```

`cpDetect` can sometimes miss fast transitions. If you find that this is the case for your data, you can try the refinement step
(see `refinement.ipynb` in `examples/` for an illustration how this step works.

```
detector.refinement(threshold=-2, reject_window=20, split_windor=50)
```

The results of the refinement step are stored in the `.refined_change_point` dictionary. You can also regenerate the
step function:

```
detector.regenerate_step_function()
```

Save the change points, log odds and the start, end for each segment:

```
detector.to_csv('filename.csv')
```

You can save the step function to a pandas data frame:

```
df = pd.DataFrame.from_dict(detector.step_function, orient='index')
df.to_csv('step_function.csv')
```

See `examples/` for confusion matrices and more on the refinement step.

Filtering
---------

`nonlinear_filter.py` implements the non-linear filter from Chung and Kennedy [DOI](https://www.ncbi.nlm.nih.gov/pubmed/1795554)

Reference
---------
* Ensign DL and Pande VS. Bayesian Detection of Intensity Changes in Single Molecule and Molecular Dynamics Trajectories.
J. Phys. Chem B 114:280 (2010) [DOI](http://pubs.acs.org/doi/abs/10.1021/jp906786b)







