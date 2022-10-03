README

This set of scripts seeks to determine optimal sampling of particles from microplastics datasets.
They are associated with the manuscript:
"Representative subsampling methods for the chemical identification of microplastic particles in environmental samples"
By Hannah De Frond, Anna M. O'Brien, and Chelsea M. Rochman

When subsampling, the scientist may wish to characterize the proportion of plastic, anthropogenic, and natural particles (broad types).
The scientist may instead wish to characterize the proportions of all material types (narrow types).

The scientist may have further goals when considering narrow types:
	-to know something about the number of material types in the environment.
	-to know how well the sampled particles represent the environmental particles.

Datasets may be one of two kinds:
	-one large sample from the environment (in main folder)
	-many individual samples from the environment (in folder "Subsampling from a group of samples").

The first dataset type can be subsampled by drawing from the single sample.
The second dataset type can be subsampled by:
	-drawing evenly across samples
	-pooling all particles and drawing from the pool
Subsampling may be achieved by the scientist by drawing randomly.
Subsampling may alternately be achieved by drawing particles according to their size, shape, or color.

One script includes the simulation functions for drawing randomly from datasets in the different ways.
sample plastics functions.R

The other scripts run the simulation functions on the datasets and evaluate the results.
	***Each of these other scripts has either "_fulldataset" or "_example" at the end of the filename***
	"_fulldataset" indicates the script that produced the results in the paper, using data provided by a number of authors
	"_example" indicates a nearly identical script that only uses the data provided as example dataset/s
	"_example" scripts will run when downloaded by a user along with all input files.
	"_fulldataset" scripts ***will not run***, as the data provided by outside authors are not included
	"_example" scripts produce all the files detailed below, but with "_example" appended to the name before the extension
	"_example" files are generated from an independent run of the simulations, and therefore may have slight deviations from results reported in the associated manuscript
	For "_fulldataset" scripts, simulations files are not included as these contain the data provided by other authors

***
The following list outputs for the "_fulldataset" scripts, rather than "_example" scripts
Outputs for "_example" scripts are the same, but with "_example" appended to the name
***

The simplest script considers datasets of one single sample, and considers only the proportion of plastic, anthropogenic, and natural particles.
sample plastics polymer diversity btype.R
Outputs:
	Rand_i_b.Rdata - simulations file 
	Col_i_b.Rdata - simulations file 
	Shp_i_b.Rdata - simulations file 
	ShpCol_i_b.Rdata - simulations file 
	minimums_SummedErrorLessThan20_broad_onejar.csv - table for making recommendations
	SummedError3panel_broad_onejar.pdf - summary figure
	SummedError_broad_onejar.pdf - summary figure

The next simplest script again considers datasets of one single sample.
This script now evaluates how well subsampling strategies characterize the proportions of all material types.
This script also considers how well the sample characterizes material types in the environment using metrics of sample coverage and Chao's estimator.
sample plastics polymer diversity.R
Outputs:
	descriptive stats_narrow_onejar.csv - descriptive table of how well the samples are expected to characterize the environment
	Rand_i.Rdata - simulations file 
	Col_i.Rdata - simulations file 
	Shp_i.Rdata - simulations file 
	ShpCol_i.Rdata - simulations file 
	minimums_SummedErrorLessThan20_narrow_onejar.csv - table for making recommendations
	minimums_EstSampleRichnessGreaterThan80_narrow_onejar.csv - table for making recommendations
	minimums_EstEnvRichnessGreaterThan80_narrow_onejar.csv - table for making recommendations
	minimums_EstEnvCoverageGreaterThan80_narrow_onejar.csv - table for making recommendations
	SummedError3Panel_narrow_onejar.pdf - summary figure
	SummedError_narrow_onejar.pdf - summary figure
	EstSampleRichness_narrow_onejar.pdf - summary figure
	EstEnvRichness_narrow_onejar.pdf - summary figure
	EstEnvCoverage_narrow_onejar.pdf - summary figure
	
These next script considers methods for subsampling datasets with multiple samples.
sample plastics structured datasets btype and full.R
Outputs:
	descriptive stats_narrow_manyjars.csv - descriptive table of how well the datasets are expected to characterize the environment
	Rand_str_dp_narrow.Rdata - simulations file 
	Shp_str_dp_narrow.Rdata - simulations file 
	Col_str_dp_narrow.Rdata - simulations file 
	ShpCol_str_dp_narrow.Rdata - simulations file 
	Rand_str_pool_narrow.Rdata - simulations file 
	Shp_str_pool_narrow.Rdata - simulations file 
	Col_str_pool_narrow.Rdata - simulations file 
	ShpCol_str_pool_narrow.Rdata - simulations file 
	Rand_str_dp_b.Rdata - simulations file 
	Shp_str_dp_b.Rdata - simulations file 
	Col_str_dp_b.Rdata - simulations file 
	ShpCol_str_dp_b.Rdata - simulations file 
	Rand_str_pool_b.Rdata - simulations file 
	Shp_str_pool_b.Rdata - simulations file 
	Col_str_pool_b.Rdata - simulations file 
	ShpCol_str_pool_b.Rdata - simulations file 
	minimums_SummedStrErrorLessThan20_narrow_manyjarsEqual.csv - table for making recommendations
	minimums_SummedStrErrorLessThan20_broad_manyjarsEqual.csv - table for making recommendations
	minimums_SummedErrorLessThan20_narrow_manyjarsPool.csv - table for making recommendations
	minimums_SummedErrorLessThan20_broad_manyjarsPool.csv - table for making recommendations
	minimums_EstSammpleRichnessGreaterThan80_narrow_manyjarsEqual.csv - table for making recommendations
	minimums_EstSammpleRichnessGreaterThan80_narrow_manyjarsPool.csv - table for making recommendations
	minimums_EstEnvRichnessgreaterthan80_narrow_manyjarsEqual.csv - table for making recommendations
	minimums_EstEnvRichnessgreaterthan80_narrow_manyjarsPool.csv - table for making recommendations
	minimums_EstEnvCoverageGreaterThan80_narrow_manyjarsEqual.csv - table for making recommendations
	minimums_EstEnvCoverageGreaterThan80_narrow_manyjarsPool.csv - table for making recommendations
	SummedStrError_narrow_manyjarsEqual_3panel.pdf - summary figure
	SummedStrError_broad_manyjarsEqual_3panel.pdf - summary figure

Finally this script uses output from other scripts to generate subfigures for panel figures.
sample plastics grouped figures.R
Outputs: (ALL are  summary figures)
	Broad_all_panel_a_indSumE.pdf
	Broad_all_panel_b_equalSumStrE.pdf
	Broad_all_panel_c_pooledSumE.pdf
	Narrow_all_panel_a_indSumE.pdf
	Narrow_all_panel_b_equalSumStrE.pdf
	Narrow_all_panel_c_poolSumE.pdf
	Narrow_panel_a_indEstSampleRichness.pdf
	Narrow_panel_b_equalEstDatasetRichnes.pdf
	Narrow_panel_c_poolEstDatasetRichness.pdf
	Narrow_panel_a_indEstEnvRichness.pdf
	Narrow_panel_b_equalEstEnvRichness.pdf
	Narrow_panel_c_poolEstEnvRichness.pdf
	Narrow_panel_a_indEstEnvCoverage.pdf
	Narrow_panel_b_equalEstEnvCoverage.pdf
	Narrow_panel_c_poolEstEnvCoverage.pdf
	