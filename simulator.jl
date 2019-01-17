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


## call all the using in main package to trigger precompilation 
## the precompilation files will get shared amongst all nodes so won't clash with `@everywhere` triggering precompilation
using Parameters      ## with julia 1.1 this is now built in.
using ProgressMeter   ## can now handle parallel with progress_pmap
#using PmapProgressMeter
#using DataFrames
using Distributions
using StatsBase
using StaticArrays
#using BenchmarkTools
#include("Influenza.jl")
#using .InfluenzaModel

#addprocs(SlurmManager(544), partition="defq", N=17)
@everywhere include("Influenza.jl")
@everywhere using .InfluenzaModel

function dataprocess(results, P::InfluenzaParameters, numberofsims; fileappend="")    
    ## if the directory argument is passed without a trailing slash, 
    ## it becomes a prepend to the filename.
    
    resultsL      = zeros(Int64, P.sim_time, numberofsims)
    resultsA      = zeros(Int64, P.sim_time, numberofsims)
    resultsS      = zeros(Int64, P.sim_time, numberofsims)
    resultsR0     = zeros(Int64, numberofsims)
    resultsSymp   = zeros(Int64, numberofsims)
    resultsAsymp  = zeros(Int64, numberofsims)
    resultsNumAge = zeros(Int64, P.grid_size_human, numberofsims)
    resultsDemoGroups = zeros(Int64, P.grid_size_human, numberofsims)
    resultsFailVector = zeros(Int64, P.grid_size_human, numberofsims)    
    resultsInfOrNot  = zeros(Int64, P.grid_size_human, numberofsims)   
    resultsSympOrNot = zeros(Int64, P.grid_size_human, numberofsims)    
    VacStatus = zeros(Int64, P.grid_size_human, numberofsims)

    Infection_Matrix = zeros(Int64, 15, 15)
    Fail_Matrix = zeros(Int64, 15, 15)
    Infection_Matrix_average = zeros(Float64, 15, 15)
    Contact_Matrix_General = zeros(Float64, 15, 15)
  
    for i=1:numberofsims
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

    
    Infection_Matrix = Infection_Matrix/numberofsims
    Fail_Matrix =  Fail_Matrix/numberofsims
    Contact_Matrix_General = Contact_Matrix_General/numberofsims
    

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


function run_calibration_R0()
    error("not implemented")
end

function run_calibration_attackrate(beta_range; numofsims = 500)
    ## calibrating to attack rate simply is running the full simulations 
    ## and then seeing the number of symptomatics at the end
    ## the full simulations are run over a range of beta values.     
    #beta_range = 0.01:0.001:0.05
    #beta_range = [0.015 0.02]
    if nprocs() == 1 && numofsims > 10
        error("Cluster not initialized, can't run more than 10 simulations. Make sure to `@everywhere include(InfluenzaModel)` before running.")
    end
    
    for i in beta_range
        println("starting simulation for i: $i")
        ## do not add the trailing slash to have "beta_0_0x" appended to the filename
        dname = "./calibration/beta_$(replace(string(i), "." => "_"))"
        if !Base.Filesystem.isdir("./calibration")
            Base.Filesystem.mkpath("./calibration")
        end
        #println(dname)
        @everywhere P = InfluenzaParameters(sim_time = 250, vaccine_efficacy = 0.0, transmission_beta=$i)        
        results = pmap(x -> main(x, P), 1:numofsims)
        println("... processing results")
        dataprocess(results, P, numofsims, fileappend=dname)   
    end
    println("calibration finished")
end

function run_single_beta(beta; process=true, numofsims=500) 
    ### runs 500 simulations with a particular beta value.    
    if nprocs() == 1 && numofsims > 10
        error("Cluster not initialized, can't run more than 10 simulations. Make sure to `@everywhere include(InfluenzaModel)` before running.")
    end
    @everywhere P = InfluenzaParameters(vaccine_efficacy = 0.0, transmission_beta=$beta)  
    results = pmap(x -> main(x, P), 1:numofsims)
    rns = randstring()
    if !Base.Filesystem.isdir("./results")
        Base.Filesystem.mkpath("./results")
    end
    dname = "./results/$(rns)_beta_$(replace(string(beta), "." => "_"))"
    if process
        dataprocess(results, P, numofsims, fileappend=dname)
    end
    return results
end

function run_attackrates(;process=true, numofsims=500)
    ### runs 500 simulations with a particular attack rate. 
    # we use the regression formula (y + 0.234616)/11.668 (from calibration) to estimate the beta value.

    if nprocs() == 1 && numofsims > 10
        error("Cluster not initialized, can't run more than 10 simulations. Make sure to `@everywhere include(InfluenzaModel)` before running.")
    end

    rns = randstring()
    resultsfolder = string("./results_", rns)
    if !Base.Filesystem.isdir(resultsfolder)
        Base.Filesystem.mkpath(resultsfolder)
    end
    println("results folder: $resultsfolder")
    #ars = [0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26]
    ars = 0.02:0.02:0.40
    ves = [0.4 0.5 0.6 0.7 0.8]
    f(y) = round((y + 0.4931677)/24.53868186, digits = 6)    
    for ar in ars, ve in ves        
        β = f(ar)
        dname =  "$(resultsfolder)/ar_$(replace(string(ar), "." => "_"))_ve_$(replace(string(ve), "." => "_"))"
        println("starting simulations for β=$β, ar=$ar, ve=$ve...")  
        @everywhere P = InfluenzaParameters(sim_time = 250, vaccine_efficacy = $ve, transmission_beta=$β)          
        results = pmap(x -> main(x, P), 1:numofsims)
        println("simulations for β=$β ended...")    
        if process
            println("starting dataprocess β=$β, ar=$ar, ve=$ve ended...")   
            dataprocess(results, P, numofsims, fileappend=dname)
        end
    end      
end

function run_samebeta()
    println("starting same-beta function")
    for _ = 1:30
        run_single_beta(0.027, numofsims=1000)
    end
    println("finished")
end

## Take all of the above and put in a function
function create_RES(agegroup, ar, ve, numofsimulations; abs_path="/results")    
    #ar: attack rate, ve: vaccine efficacy
    cd(abs_path)
    #println("working directory: $(pwd())")
  
    ## create the metadata for the DataFrame we read in. 
    headers=["sim$i" for i = 1:numofsimulations];
    inttypes = [Int64 for i = 1:numofsimulations];
    floattypes = [Float64 for i = 1:numofsimulations];     
    
    ## are we working with one agegroup or multiple agegroups
    if length(agegroup) == 1  ## only dealing with one agegroup
        resname = "RES_ar_$(replace(string(ar), "." => "_"))_ve_$(replace(string(ve), "." => "_"))_ag_$(agegroup[1]).dat"
    else 
        resname = "RES_ar_$(replace(string(ar), "." => "_"))_ve_$(replace(string(ve), "." => "_"))_ag_0.dat"
    end    
    
    ## dynamically create the file name. 
    foldername = "ar_$(replace(string(ar), "." => "_"))_ve_$(replace(string(ve), "." => "_"))"
    #println("file append: $foldername")    

    Latent = CSV.File( string(foldername, "_latent.dat"), delim='\t', header=headers, types=inttypes) |> DataFrame
    DemoGroup = CSV.File(string(foldername, "_DemoGroups.dat"), delim='\t', header=headers, types=inttypes) |> DataFrame
    InfOrNot =  CSV.File(string(foldername, "_InfOrNot.dat"),    delim='\t', header=headers, types=inttypes) |> DataFrame
    SympOrNot = CSV.File(string(foldername, "_SympOrNot.dat"),   delim='\t', header=headers, types=inttypes) |> DataFrame
    VacStatus = CSV.File(string(foldername, "_VacStatus.dat"),   delim='\t', header=headers, types=inttypes) |> DataFrame
    
    RES = DataFrame()
    _numinagegroup = zeros(Int64, numofsimulations)
    for (i, (colname, coldata)) in zip(1:numofsimulations, eachcol(InfOrNot, true))
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

    colsums = zeros(Int64, numofsimulations)
    for (i, (colname, coldata)) in zip(1:numofsimulations, eachcol(InfOrNot, true))
        indices = findall(x -> x in agegroup, DemoGroup[colname])
        colsums[i] = sum(coldata[indices])
    end           
    RES[:total_inf] = colsums
    
    _n_symp  = zeros(Int64, numofsimulations)
    _n_asymp = zeros(Int64, numofsimulations)    

    for (i, (colname, coldata)) in zip(1:numofsimulations, eachcol(SympOrNot, true))
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
    
    _totalvac = zeros(Int64, numofsimulations)
    _sympvax = zeros(Int64, numofsimulations)
    _asympvax = zeros(Int64, numofsimulations)
    _sympunvax = zeros(Int64, numofsimulations)
    _asympunvax = zeros(Int64, numofsimulations)
    for (i, (colname, coldata)) in zip(1:numofsimulations, eachcol(VacStatus, true))
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
        
    #if all the columns are in, add the age group, ve, ar
    #RES[:desc] = resname
         
    CSV.write(resname, RES)
    return nothing    
end

function create_all_RES()
    ars = 0.02:0.02:0.40
    ves = [0.4 0.5 0.6 0.7 0.8]
    p = Progress(length(ars)*length(ves), 1)   # minimum update interval: 1 second
    ag = [1, 2, 3, 4]  
    for ar in ars, ve in ves 
        #println("starting process for ar=$ar, ve=$ve...")
        create_RES(ag, ar, ve, 1000, abs_path="/home/affans/influenza/results_RmEpwH1v");
        create_RES([1], ar, ve, 1000, abs_path="/home/affans/influenza/results_RmEpwH1v");
        create_RES([2], ar, ve, 1000, abs_path="/home/affans/influenza/results_RmEpwH1v");
        create_RES([3], ar, ve, 1000, abs_path="/home/affans/influenza/results_RmEpwH1v");
        create_RES([4], ar, ve, 1000, abs_path="/home/affans/influenza/results_RmEpwH1v");
        next!(p; showvalues = [(:attackrate,ar), (:efficacy,ve)])
    end
    println("results processing finished")
end



