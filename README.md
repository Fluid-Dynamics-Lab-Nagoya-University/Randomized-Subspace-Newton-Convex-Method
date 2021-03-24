<!-- Last modified: 2021/01/06 -->
# Randomized Subspace Newton Convex Method  
---
This repository contains Matlab R2020a code to reproduce results for a manuscript entitled __"Randomized Subspace Newton Convex Method Applied to Data-Driven Sensor Selection Problem"__ to be published in IEEE Signal Processing Letters.  
The sparse sensor selection problem is solved by the rabdomized subspace Newton convex method.  
To run the program, excute `P_demo`.  

## Directory  
---
- src: source code is stored  
- work: calculation results are stored (created automatically by running the program)  
- data: __NOAA Optimum Interpolation (OI) Sea Surface Temperature (SST) V2__ data is stored  
  - sst.wkmean.1990-present.nc  
  - lsmask.nc  
NOAA_OI_SST_V2 is provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their Web site at https://www.esrl.noaa.gov/psd/.  
Due to GitHub file size limitations, a dataset is linked online: [NOAA Optimum Interpolation (OI) Sea Surface Temperature (SST) V2](https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html)  

## Code  
---
### Main program  
- P_demo.m  

### Function  
#### Preprocessing  
- F_pre_read_NOAA_SST.m  
- F_pre_SVD_NOAA_SST.m  
- F_pre_truncatedSVD.m  

#### Sensor selection  
- F_sensor_random.m  
- F_sensor_DC.m  
  - F_sensor_DC_sub.m  
    - F_sensor_DC_approxnt_vec.m  
    - F_sensor_DC_approxnt.m  
    - F_sensor_DC_loc.m  
    - F_sensor_DC_locr.m  
- F_sensor_RSNC.m  
  - F_sensor_RSNC_approxnt.m  
- F_sensor_CRSNC.m  
  - F_sensor_CRSNC_approxnt.m  

#### Calculation
- F_calc_det.m  
- F_calc_trace.m  
- F_calc_eigen.m  
- F_calc_sensormatrix.  
- F_calc_error.m  
  - F_calc_reconst.m  
  - F_calc_reconst_error.m  

#### Data organization  
- F_data_ave1.m  
- F_data_ave2.m  
- F_data_arrange1.m  
- F_data_arrange2.m  
- F_data_arrange3.m  
- F_data_normalize.m  

#### Mapping
- F_map_original.m  
	- F_map_videowriter.m  
		- F_map_plot_sensors_forvideo.m  
- F_map_reconst.m  
	- F_map_plot_sensors.m  

### Function  
#### Preprocessing  
- F_pre_read_NOAA_SST.m  
- F_pre_SVD_NOAA_SST.m  
- F_pre_truncatedSVD.m  

## How to cite  
---
If you use `Randomized Subspace Newton Convex Method` code in your work, please cite the software itself and relevent papers.  
### General software reference:  
``` bibtex
@misc{nakai2021github,
      author = {Taku Nonomura},
      title = {Randomized Subspace Newton Convex Method},
      howpublished = {Available online},
      year = {2021},
      url = {https://github.com/KumiNakai/Randomized-Subspace-Newton-Convex-Method}
}
```  

### Randomized subspace Newton algorithm:  
``` bibtex
@misc{nonomura2021randomized,
      title={Randomized Subspace Newton Convex Method Applied to Data-Driven Sensor Selection Problem}, 
      author={Taku Nonomura and Shunsuke Ono and Kumi Nakai and Yuji Saito},
      journal={IEEE Signal Processing Letters},
      <!-- volume={}, -->
      <!-- number={}, -->
      <!-- pages={}, -->
      year={2021},
      publisher={IEEE},
      <!-- doi={}, -->
      <!-- eprint={} -->
}
```

### Greedy algorithm based on D-optimality:  
``` bibtex
@misc{saito2020determinantbased,
      title={Determinant-based Fast Greedy Sensor Selection Algorithm}, 
      author={Yuji Saito and Taku Nonomura and Keigo Yamada and Keisuke Asai and Yasuo Sasaki and Daisuke Tsubakino},
      year={2019},
      eprint={1911.08757},
      archivePrefix={arXiv},
      primaryClass={eess.SP}
}
```

## License  
---
[MIT-License](https://github.com/KumiNakai/Randomized-Subspace-Newton-Convex-Method/blob/master/LICENSE)

## Author
---
Taku Nonomura  
[Experimental Aerodynamics Laboratory](http://www.aero.mech.tohoku.ac.jp/eng/)  
Department of Aerospace Engineering, Tohoku University  
Sendai, JAPAN  
E-mail: nonomura@tohoku.ac.jp
