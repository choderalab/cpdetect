Examples
========

All examples are on randomly generated synthetic data

Manifest
--------
* **Simple example:**


 
   * `example.ipynb` - simple example on 2 trajectories
   * `ts_log_odss.csv` - saved change points from simple example

`synthetic_data.ipynb` - notebook that generated synthetic trajectories

`synthetic.pdf` - plots of all synthetic trajectories and true step functions (in red).

* **Confusion matrices:**


    * `example_nonfiltered.py` - run cpDetect with several different thresholds
    to calculate confusion matrices for each. Data was not filtered.
    * `example_filtered.py` - run cpDetect on filtered data with several different
    thresholds.
    * `confusion_matrix.ipynb` - ipython notebook that calculates the confusion matrix
    for each threshold. Unfiltered data
    * `confusion_matrix_filtered.ipynb` - ipython notebook that calculates the confusion
    matrix for each threshold. Data is filtered.
    * `confusion_matrix.pdf` - confusion matrices for all 10 runs. Unfiltered data
    * `confusion_matrix.pdf` - confusion matrices for all 10 runs. Filtered data
    
    While cpDetect on filtered data misses less change points, it also founds many
    false positives. Use a higher threshold when data is filtered. 

* **Clean up step (refinement):**


    * `refinement.ipynb` - ipython notebook illustrating refinement step
    * `confusion_matrix.png` - confusion matrix before refinement
    * `refined_confusion_matrix.png` - confusion matrix after refinement
    * `test_refinement.pdf` - results from refinement step on synthetic trajectories. 
    red - true step function; black - predicted step function from initial cpDetect run;
    green lines - change points found with refinement step. 
    
    