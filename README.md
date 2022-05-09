# Constant-Beamwidth Kronecker Product Beamforming with Nonuniform Planar Arrays  
Please cite our paper that is available at <https://www.frontiersin.org/article/10.3389/frsip.2022.829463>.  

We included two main scripts to show how to create constant-beamwidth beamformers for symmetric nonuniform linear and planar arrays:  
* plotsForLinearArray.m - for a linear array  
* plotsForPlanarArray.m - for a planar array  

### Resulting plots for the linear array:  
<img src="saved figs/linear_array_sensor_positions.jpg" width="450">  
<img src="saved figs/linear_array_kaiser_window_parameters.jpg" width="450">  

<img src="saved figs/linear_array_weights_and_beampatterns.jpg" width="450">  

<img src="saved figs/linear_array_beamwidth_WNG_DF.jpg" width="450">  

### Resulting plots for the planar array:  
<img src="saved figs/planar_array_sensor_positions.jpg" width="450">  
<img src="saved figs/planar_array_beampatterns.jpg" width="450">  
<img src="saved figs/planar_array_beamwidths.jpg" width="450">  

<img src="saved figs/planar_array_WNG_and_DF.jpg" width="450">  

#### Additional Details:
The symmetric linear array’s weights are computed with Algorithm1_AttainingTheWeights.m.  
A more robust code that can be used: Algorithm1_AttainingTheWeights_robust.m.  
The robust version works correctly for low resolutions of beta.  

The symmetric linear array’s sensor positions can be found with two different algorithms: the iterative and non-iterative algorithms:  
* iterativeAlgorithm_sensorPositions_main.m - the iterative algorithm  
* nonIterativeAlgorithm_sensorPositions.m - the non-iterative algorithm  

The iterative algorithm achieves superior performance but has high computational complexity.  
The non-iterative algorithm has low complexity yet achieves on-par performance.
