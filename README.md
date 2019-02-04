# Agent Based Modelling framework for Influenza. 

An age-structured dynamical agent-based model for seasonal influenza. 

  - Age-specific contact patterns in an urban setting
  - Vaccination scenarios
  - Rapid mutation of virus

# API
### Basic Usage
The main framework is built into the module `Influenza`. Basic usage includes:
```
include("Influenza.jl")
using .InfluenzaModel
P = InfluenzaParameters() ## edit entries as required, see below. 
results = main(1, P)
```
This runs one Monte Carlo simulation by running `main`, using the default values of `InfluenzaParameters` (see below) and stores the results in `results`. The signature of the main function is `main(simnum::Int64, P::InfluenzaParameters)` where `simnum` is used to set the seed of the simulation for reproducibility purposes. 

### Influenza Parameters
`Parameters` can be initalized with the following information (default values are shown).  
```
@with_kw struct InfluenzaParameters @deftype Int64
    sim_time = 200
    grid_size_human = 10000
    initial_infected = 1
    ## latent periods
    Latent_period_Max = 2
    Latent_period_Min = 1
    ## probability of going to latent -> asymptomatic/symptomatic
    ProbAsympMin::Float64 = 0.3
    ProbAsympMax::Float64 = 0.7
    reduction_factor::Float64 = 0 ## for now we set this to zero: this means prob of infection is the same whether symp or asymp
    ## if a person is vaccinated, what is their vaccine efficacy 
    ## this number eventually gets multiplied by the (calculated) frailty index.
    vaccine_efficacy::Float64 = 0.0 # change to different values 0.8, 0.6, 0.4
    transmission_beta::Float64 = 0.079 # we need to calibrate this for 1%, 4% and 8% clinical attack rate (symptomatic cases)
    log_normal_mean::Float64 = 1.0 # distribution for infectious period ## scale value
    log_normal_shape::Float64 = 0.4356
end
```
### Parallel runs. 
It is recommended to use multiple cores to run simulations in an embarrasingly parallel methodology. Basic usage as follows.
```
addprocs(64)
@everywhere include("Influenza.jl")
@everywhere using .InfluenzaModel

NUMOFSIMS = 64
@everywhere P = InfluenzaParameters(sim_time = 250, vaccine_efficacy = $ve, transmission_beta=$β)          
results = pmap(x -> main(x, P), 1:NUMOFSIMS)
## saving the results
fn = "foldername/fileappend"
dataprocess(results, P, fileappend=fn)  ## dataprocess uses NUMOFSIMS internally. Make sure its defined. 
```
### Results of each simulation 
The results of a simulation is written to disk via the `dataprocess` function. The file names and their descriptions are as follows:
TODO

### Use Case in the ABM LAB. 
The ABM-LAB comes equipped with a high performance cluster. Large number of simulations (2000, ~7-9 hours) 
The following function sweeps over varying attack rates and vaccine efficacies. 
```
function run_attackrates(ARS, VES)
    ## runs the specified number of simulations for beta (or a range of betas) 
    ## Attack rates are stored in a global vector: ARS
    ## we use regression analysis to calculate a formula for beta. 
    ## the beta values are calculated on the fly for each attack rate
    # common ars: 0.04, 0.08, 0.12, 0.20, 0.30, 0.40
    # common ve: 30%, 40%, 50%, 60%, 70%, 80%
    _clustercheck()
    RF = create_folder()
    #f(y) = round((y + 0.4931677)/24.53868186, digits = 6)  
    #f(t) = round((t + 1.09056771093182)/  58.2096402005676, digits=6)
    #f(t) = -1.09056771093182 + 58.2096402005676t
    
    #betas = [0.0193,  0.0204,  0.0209,  0.0222,  0.0237,  0.0255]
    arb = Dict(0.04=>0.0193, 0.08=>0.0204, 0.12=>0.0209, 0.20=>0.0222, 0.30=>0.0237, 0.40=>0.0255)
    f(t) = arb[t]
    prgess = Progress(length(ARS)*length(VES), 1)   # minimum update interval: 1 second

    for ar in ARS, ve in VES        
        β = f(ar)         
        @everywhere P = InfluenzaParameters(sim_time = 250, vaccine_efficacy = $ve, transmission_beta=$β)          
        results = pmap(x -> main(x, P), 1:NUMOFSIMS)
        dname = "$RF/$(create_fn(ar, ve))"
        dataprocess(results, P, fileappend=dname)
        create_RES([1,2,3,4,5], RF, ar=ar, ve=ve, resultformat="arve");
        create_RES([1], RF, ar=ar, ve=ve, resultformat="arve");
        create_RES([2], RF, ar=ar, ve=ve, resultformat="arve");
        create_RES([3], RF, ar=ar, ve=ve, resultformat="arve");
        create_RES([4], RF, ar=ar, ve=ve, resultformat="arve");
        create_RES([5], RF, ar=ar, ve=ve, resultformat="arve");
        next!(prgess; showvalues = [(:β, β), (:attackrate,ar), (:efficacy,ve)])
    end
    println("\n")   
end
```

### Unit Testing 
TODO. 

### Development
Want to contribute? Great! PRs are welcome. 
#### TO DO LIST: 
- Create a dashboard to visualize and provide information of running simulations. 
 - Write MORE Tests


