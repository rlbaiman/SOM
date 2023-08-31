# Self Organizing Map (SOM)
Self Organizing Map Resources

* example_data files include a small initial data set 'SOM_example_data.csv' and 'SOM_example_data.nc4'. This is identical data but the csv has flattened the 2d data into single rows
* SOM_run.R shows basic runs of SOM using the Kohonen package and the MASS package for sammon mapping. This includes choosing parameters for the SOM, training the model, some basic plots assessing the model, and output options for analysis
* Pearson_Anlysis.ipynb shows pearson correlation coefficient analysis of SOM. This uses the som output from SOM_run.R (example_node_assignments.csv, example_nodes.csv) as well as the original spatial data SOM_example_data.nc4

These resources can be used to perform SOM analysis on any spatial data with multiple timesteps and a list of variable values. If you are using spatial data, you will need to convert a .nc file to a .csv file similar to SOM_example_data.csv where each row is a timestep and the columns are flattened spatial data.

Please email rebecca.baiman@colorado.edu with any questions!
