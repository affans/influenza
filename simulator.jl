## STDLIB
using Random
using DelimitedFiles
using Distributed
using Base.Filesystem
## ADDED PACKAGES
using ClusterManagers
using DataFrames
using CSV
using JSON
using Statistics
using VegaLite
import REPL
using REPL.TerminalMenus

## call all the using in main package to trigger precompilation 
## the precompilation files will get shared amongst all nodes so won't clash with `@everywhere` triggering precompilation
using Parameters      ## with julia 1.1 this is now built in.
using ProgressMeter   ## can now handle parallel with progress_pmap
#using PmapProgressMeter
#using DataFrames
using Distributions
using StatsBase
using StaticArrays

options = ["yes", "no"]
menu = RadioMenu(options, pagesize=4)
choice = request("Do you want to connect to cluster?", menu)
if options[choice] == "yes"
    addprocs(SlurmManager(544), partition="defq", N=17)
end
@everywhere include("Influenza.jl")
@everywhere using .InfluenzaModel

## simulation global Parameters
const NUMOFSIMS = 3000

function _clustercheck()
    if nprocs() == 1 && NUMOFSIMS > 10
        error("Cluster not initialized, can't run more than 10 simulations. Make sure to `@everywhere include(InfluenzaModel)` before running.")
    end 
end

## given an attack rate and vaccine efficacy, it creates a string out of that information. 
create_fn(ar, ve) = "ar_$(replace(string(ar), "." => "_"))_ve_$(replace(string(ve), "." => "_"))"

function create_folder()
  RF = string("results_", randstring()) ## 
  if !Base.Filesystem.isdir(RF)
      Base.Filesystem.mkpath(RF)
  end
  return RF
end

function run_beta(β_range, ve)
    ## runs the specified number of simulations for beta (or a range of betas) 
    _clustercheck()
    RF = create_folder()  ## to store the results
    prgess = Progress(length(β_range), 1)  ## a progressbar, minimum update interval: 1 second    
    for β in β_range
        dname = "$RF/beta_$(replace(string(β), "." => "_"))"
        @everywhere P = InfluenzaParameters(sim_time = 250, vaccine_efficacy = $ve, transmission_beta=$β)        
        results = pmap(x -> main(x, P), 1:NUMOFSIMS)
        dataprocess(results, P, fileappend=dname)   
        create_RES([1,2,3,4,5], RF, β=β, resultformat="beta");
        create_RES([1], RF, β=β, resultformat="beta");
        create_RES([2], RF, β=β, resultformat="beta");
        create_RES([3], RF, β=β, resultformat="beta");
        create_RES([4], RF, β=β, resultformat="beta");
        create_RES([5], RF, β=β, resultformat="beta");
        next!(prgess; showvalues = [(:β, β)])
    end    

end

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

function dataprocess(results, P::InfluenzaParameters; fileappend="./")     
    ## takes the results of the pmap and stores it to file. 
    
    ## create empty vectors to store the results
    resultsL      = zeros(Int64, P.sim_time, NUMOFSIMS)
    resultsA      = zeros(Int64, P.sim_time, NUMOFSIMS)
    resultsS      = zeros(Int64, P.sim_time, NUMOFSIMS)
    resultsR0     = zeros(Int64, NUMOFSIMS)
    resultsSymp   = zeros(Int64, NUMOFSIMS)
    resultsAsymp  = zeros(Int64, NUMOFSIMS)
    resultsNumAge = zeros(Int64, P.grid_size_human, NUMOFSIMS)
    resultsDemoGroups = zeros(Int64, P.grid_size_human, NUMOFSIMS)
    resultsFailVector = zeros(Int64, P.grid_size_human, NUMOFSIMS)    
    resultsInfOrNot  = zeros(Int64, P.grid_size_human, NUMOFSIMS)   
    resultsSympOrNot = zeros(Int64, P.grid_size_human, NUMOFSIMS)    
    VacStatus = zeros(Int64, P.grid_size_human, NUMOFSIMS)
    Infection_Matrix = zeros(Int64, 15, 15)
    Fail_Matrix = zeros(Int64, 15, 15)
    Infection_Matrix_average = zeros(Float64, 15, 15)
    Contact_Matrix_General = zeros(Float64, 15, 15)
  
    for i=1:NUMOFSIMS
        resultsL[:,i] = results[i][1]
        resultsS[:,i] = results[i][2]
        resultsA[:,i] = results[i][3]
      
        resultsR0[i] = results[i][4]
        resultsSymp[i] = results[i][5]
        resultsAsymp[i] = results[i][6]

        Infection_Matrix = Infection_Matrix + results[i][7]
        Fail_Matrix =  Fail_Matrix + results[i][8]
        Contact_Matrix_General = Contact_Matrix_General + results[i][9]

        resultsNumAge[:,i] = results[i][10]
        resultsFailVector[:,i] = results[i][11]
        resultsInfOrNot[:,i] = results[i][12]
        VacStatus[:,i] = results[i][13]
        resultsSympOrNot[:, i] = results[i][14]
        resultsDemoGroups[:, i] = results[i][15]
    end
    Infection_Matrix = Infection_Matrix/NUMOFSIMS
    Fail_Matrix =  Fail_Matrix/NUMOFSIMS
    Contact_Matrix_General = Contact_Matrix_General/NUMOFSIMS
    
    writedlm(string("$fileappend", "_latent.dat"), resultsL)
    writedlm(string("$fileappend", "_symp.dat"),resultsS)
    writedlm(string("$fileappend", "_asymp.dat"),resultsA)
    writedlm(string("$fileappend", "_R0.dat"),resultsR0)
    writedlm(string("$fileappend", "_SympInf.dat"),resultsSymp)
    writedlm(string("$fileappend", "_AsympInf.dat"),resultsAsymp)
    writedlm(string("$fileappend", "_InfMatrix.dat"),Infection_Matrix)
    writedlm(string("$fileappend", "_FailMatrix.dat"),Fail_Matrix)
    writedlm(string("$fileappend", "_ContactMatrixGeneral.dat"),Contact_Matrix_General)
    writedlm(string("$fileappend", "_NumAgeGroup.dat"),resultsNumAge)
    writedlm(string("$fileappend", "_FailVector.dat"),resultsFailVector)
    writedlm(string("$fileappend", "_InfOrNot.dat"),resultsInfOrNot)
    writedlm(string("$fileappend", "_VacStatus.dat"),VacStatus)
    writedlm(string("$fileappend", "_SympOrNot.dat"),resultsSympOrNot)
    writedlm(string("$fileappend", "_DemoGroups.dat"),resultsDemoGroups)
    JSON.print(open(string("$fileappend", "_parameters.dat"), "w"), P, 4)
end


function create_RES(agegroup, folder; ar=0.0, ve=0.0, β=0.0, resultformat="arve")    
    ## This function takes the raw data by dataprocess() and processes it more (using age groups) for the purposes of MedicagoProject
    ## agegroup = [1, 2, 3, 4], [1], [2], [3], [4] 
    !isdir(folder) && error("input folder is not a valid folder")
    #resultformat in ["arve", "beta"] && error("resultformat it not a valid format: enter either arve or beta")

    ## cfn is the file name structure to read the raw data, resname is the name of the file that will be saved after processing. 
    if resultformat == "arve"
        cfn = "$folder/$(create_fn(ar, ve))"
        if length(agegroup) == 1  ## only dealing with one agegroup
            resname = "$folder/RES_$(create_fn(ar, ve))_ag_$(agegroup[1]).dat"
        else 
            resname = "$folder/RES_$(create_fn(ar, ve))_ag_0.dat"
        end  
    else 
        cfn = "$folder/beta_$(replace(string(β), "." => "_"))"
        if length(agegroup) == 1  ## only dealing with one agegroup
            resname = "$folder/RES_beta_$(replace(string(β), "." => "_"))_ag_$(agegroup[1]).dat"
        else 
            resname = "$folder/RES_beta_$(replace(string(β), "." => "_"))_ag_0.dat"
        end  
    end 
   
    ## create the metadata for the DataFrame we read in. 
    headers=["sim$i" for i = 1:NUMOFSIMS];
    inttypes = [Int64 for i = 1:NUMOFSIMS];
    floattypes = [Float64 for i = 1:NUMOFSIMS];     
    
    ## read the dataframes
    Latent =    CSV.File(string(cfn, "_latent.dat"),      delim='\t', header=headers, types=inttypes) |> DataFrame
    DemoGroup = CSV.File(string(cfn, "_DemoGroups.dat"),  delim='\t', header=headers, types=inttypes) |> DataFrame
    InfOrNot =  CSV.File(string(cfn, "_InfOrNot.dat"),    delim='\t', header=headers, types=inttypes) |> DataFrame
    SympOrNot = CSV.File(string(cfn, "_SympOrNot.dat"),   delim='\t', header=headers, types=inttypes) |> DataFrame
    VacStatus = CSV.File(string(cfn, "_VacStatus.dat"),   delim='\t', header=headers, types=inttypes) |> DataFrame
    
    ## the results are written to this file name. depends on the agegroup
    RES = DataFrame()

    ## number in each agegroup
    _numinagegroup = zeros(Int64, NUMOFSIMS)
    for (i, (colname, coldata)) in zip(1:NUMOFSIMS, eachcol(InfOrNot, true))
        indices = findall(x -> x in agegroup, DemoGroup[colname])
        _numinagegroup[i] = length(indices)
    end      
    RES[:num_in_agegroup] = _numinagegroup
    
    # latentsums = zeros(Int64, numofsimulations)
    # for (i, (colname, coldata)) in zip(1:numofsimulations, eachcol(Latent, true))
    #     indices = findall(x -> x in agegroup, DemoGroup[colname])
    #     latentsums[i] = sum(coldata[indices])
    # end 
    # RES[:latent] = latentsums
    RES[:latent] = sum.(eachcol(Latent, false))

    colsums = zeros(Int64, NUMOFSIMS)
    for (i, (colname, coldata)) in zip(1:NUMOFSIMS, eachcol(InfOrNot, true))
        indices = findall(x -> x in agegroup, DemoGroup[colname])
        colsums[i] = sum(coldata[indices])
    end           
    RES[:total_inf] = colsums
    
    _n_symp  = zeros(Int64, NUMOFSIMS)
    _n_asymp = zeros(Int64, NUMOFSIMS)    

    for (i, (colname, coldata)) in zip(1:NUMOFSIMS, eachcol(SympOrNot, true))
        ## first get the indices of the humans that belong to the correct age group.
        indices = findall(x -> x in agegroup, DemoGroup[colname])
        ## get their sickness data by subsetting the column. 
        _tmparr = coldata[indices]
         ## then from this filtered column data, count how many 3s and 4s there are
        _n_symp[i] = length(findall(x -> x == 3, _tmparr))
        _n_asymp[i] = length(findall(x -> x == 4, _tmparr))    
    end
  
    RES[:total_symp]  = _n_symp
    RES[:total_asymp] = _n_asymp
    
    _totalvac   = zeros(Int64, NUMOFSIMS)
    _sympvax    = zeros(Int64, NUMOFSIMS)
    _asympvax   = zeros(Int64, NUMOFSIMS)
    _sympunvax  = zeros(Int64, NUMOFSIMS)
    _asympunvax = zeros(Int64, NUMOFSIMS)
    for (i, (colname, coldata)) in zip(1:NUMOFSIMS, eachcol(VacStatus, true))
        indices = findall(x -> x in agegroup, DemoGroup[colname])
        ## total vaccination
        _totalvac[i] = sum(coldata[indices])

        ## out of total vaccination how many were sympotmatic   
        map(indices) do x
            if SympOrNot[x, colname] == 3 
                if coldata[x] == 1
                    _sympvax[i] += 1
                else 
                    _sympunvax[i] += 1
                end            
            end   
            if SympOrNot[x, colname] == 4
                if coldata[x] == 1
                    _asympvax[i] += 1
                else
                    _asympunvax[i] += 1
                end
            end  
        end    
    end
    RES[:total_vac] = _totalvac
    RES[:total_vac_symp] = _sympvax    
    RES[:total_vac_asymp] = _asympvax
    RES[:total_vac_inf] = RES[:total_vac_symp] + RES[:total_vac_asymp]
    RES[:total_vac_healthy] = RES[:total_vac] - RES[:total_vac_inf]
    RES[:total_unvac_symp] = _sympunvax
    RES[:total_unvac_asymp] = _asympunvax
    RES[:total_unvac_inf] = RES[:total_unvac_symp] + RES[:total_unvac_asymp]
    RES[:total_unvac_healthy] = (RES[:num_in_agegroup] - RES[:total_vac]) - RES[:total_unvac_inf]
        
    CSV.write(resname, RES)
    return nothing    
end
