module InfluenzaModel


using Parameters      ## with julia 1.1 this is now built in.
using ProgressMeter   ## can now handle parallel with progress_pmap
#using PmapProgressMeter
#using DataFrames
using Distributions
using StatsBase
using StaticArrays
#using BenchmarkTools

export init, main, Human, InfluenzaParameters

include("parameters.jl")
include("population.jl")
include("functions.jl")

function init()
    [Human(i) for i = 1:10000]
end

function main(simnum::Int64, P::InfluenzaParameters)
    #println("starting simulation number: $simnum")
    humans = init()
    setup_demographic(humans)  ## TO DO, unit tests, plotting  
    apply_vaccination(humans,P)    ## TO DO, unit tests, plotting
    
    initial = setup_rand_initial_latent(humans,P) ## returns the ID of the initial person

    ## data collection variables
    latent_ctr = zeros(Int64, P.sim_time)   #vector for results latent
    symp_ctr =   zeros(Int64, P.sim_time)   #vector for results symp
    asymp_ctr =  zeros(Int64, P.sim_time)   #vector for results asymp
 
    ## contact matrix calculations, create 15x15 matrices
    Fail_Contact_Matrix    = zeros(Int64, 15, 15)
    Contact_Matrix_General = zeros(Int64, 15, 15)
    Number_in_age_group    = zeros(Int64, 15)                      # vector that tells us number of people in each age group.
    Age_group_Matrix       = zeros(Int64, 15, P.grid_size_human)   # don't know what this does. half the columns are empty.
       
    setup_contact_matrix(humans, Age_group_Matrix, Number_in_age_group)

    ## main simulation loop.
    for t=1:P.sim_time
        #@time contact_dynamic2(humans, P, Fail_Contact_Matrix, Age_group_Matrix, Number_in_age_group, Contact_Matrix_General)
        contact_dynamic2(humans, P, Fail_Contact_Matrix, Age_group_Matrix, Number_in_age_group, Contact_Matrix_General)
        for i=1:10000
           increase_timestate(humans[i], P)
        end
      
        latent_ctr[t], symp_ctr[t], asymp_ctr[t] = update_human(humans,P)
    end

    first_inf = findall(x-> x.WhoInf == initial && x.WentTo == SYMP, humans)
    symp_inf = findall(x -> x.WhoInf>0 && humans[x.WhoInf].WentTo == SYMP, humans)
    asymp_inf = findall(x -> x.WhoInf>0 && humans[x.WhoInf].WentTo == ASYMP, humans)

    numb_symp_inf = length(symp_inf)
    numb_asymp_inf = length(asymp_inf)
    numb_first_inf = length(first_inf)

    Number_in_age_group_two = zeros(Int64,P.grid_size_human)
    NumberFailsAge = zeros(Int64,P.grid_size_human)
    InfOrNot = zeros(Int64,P.grid_size_human)
    VacStatus = zeros(Int64,P.grid_size_human)
    Infection_Matrix = zeros(Int64,15,15)

    for i = 1:length(humans)
        Number_in_age_group_two[i] = humans[i].contact_group
        NumberFailsAge[i] = humans[i].NumberFails
        VacStatus[i] = humans[i].vaccineEfficacy > 0 ? 1 : 0 ## TO DO FIX
        if humans[i].WhoInf > 0
            Infection_Matrix[humans[i].contact_group, humans[humans[i].WhoInf].contact_group] += 1 
            if humans[i].health == REC || humans[i].health == SYMP || humans[i].health == ASYMP || humans[i].health == LAT
                InfOrNot[i] = 1
            end
        end
    end

    return latent_ctr, symp_ctr, asymp_ctr, 
    numb_first_inf, numb_symp_inf, numb_asymp_inf, 
    Infection_Matrix, Fail_Contact_Matrix, Contact_Matrix_General, 
    Number_in_age_group_two, NumberFailsAge, InfOrNot, VacStatus
end

end

