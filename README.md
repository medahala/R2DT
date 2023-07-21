# R2DT

Description and purpose of programs that were used to create simulation data contained within R2DT manuscript 
################################################################
#Programs
#################################################################

Scenerio_Data_generation.r
#########################################
creates simulated_data.rds, a large array of random patient outcomes for given scenarios 
also creates a list of fixed scenario probabilities 

F_simulations_cache_sep.r
#########################################
large simulation function read in by output_sims_alpha.r. Premise is to run all decision methods and scenerios at once.

output_sims_alpha.r
########################################
creates sim_allcohorts.rds a list of 20 dataframes, each line of a datataset 
gives a dose recomendation for a given simulation, scenario and decision method following completion of each cohort  
dependences are simulated_data.rds, F_simulations_cache_sep.r, stan model files
and Dosedet functions (This needs a package creating, see below)

inputs for F_simulations_cache_sep.r are also generated and saved:
decision_constants.rds - fixed parameters asociated with decision models
mod_constants.rds - constants and hyper paramters for proability model

output_tables
########################################
creates outputs.html (a searcheable table of summary stats) and paper plots by samplesize


################################################
Dosedet and Stanfiles
################################################
Eff_logisitc.stan and Tox_logisitc.stan
stan logistic model files for efficacy and toxicity

################################################
DoseDet folder
#######################################
The decsision method functions are utilised via a user package to facilitate
parralel programing. PACKAGE NEEDS INSTALLING on local machine via following R code
#install package from local directory where doseDet folder is saved
#install.packages("XXXXXXX/doseDet", repos=NULL, type="source")

brief description of functions used: taking proabilities of efficacy and toxicty from stan and outputing a dose 
recomendation
"doseDet::prospect_pow_addC" - R2DT method with independent addmisiibility rules
"doseDet::prospect_pow_add_u_C" - R2DT method with utility addmisiibility rule applied to each dose at each stage
"doseDet::prospect_pow_u_C" - R2DT method with single trial stopping rule based upon utility contour
"doseDet::efftox_add" - original Efftox method with independent admisibility rules

following function takes efficacy and toxicity probobailities and converts them to a utility value according to
R2DT method
"doseDet::UdetC"







