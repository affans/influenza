using Parameters
using ProgressMeter
#using PmapProgressMeter
using DataFrames
using QuadGK
using Distributions
using StatsBase
using StaticArrays

include("parameters.jl")
include("population.jl")
include("functions.jl")

function main(cb,simulationNumber::Int64,P::InfluenzaParameters)

    ## set up a grid of humans
    humans =  [Human() for i = 1:10000]

    setup_demographic(humans,P)
    if P.GeneralCoverage == 1
        vaccination(humans,P)
    end
    latent_ctr = zeros(Int64,P.sim_time) #vector of results latent
    symp_ctr = zeros(Int64,P.sim_time) #vector for results symp
    asymp_ctr = zeros(Int64,P.sim_time) #vector for results asymp

    Fail_Contact_Matrix = zeros(Int64,15,15)
    Contact_Matrix_General = zeros(Int64,15,15)



    initial=setup_rand_initial_latent(humans,P) # for now, we are using only 1

    Number_in_age_group = zeros(Int64,15)
    Age_group_Matrix = Matrix{Int64}(15,P.grid_size_human)
    for i = 1:P.grid_size_human
        Age_group_Matrix[humans[i].contact_group,(Number_in_age_group[humans[i].contact_group]+1)] = humans[i].index
        Number_in_age_group[humans[i].contact_group] += 1
    end

    for t=1:P.sim_time

        contact_dynamic2(humans,P,Fail_Contact_Matrix,Age_group_Matrix,Number_in_age_group,Contact_Matrix_General)

        for i=1:P.grid_size_human
            increase_timestate(humans[i],P)
        end
        latent_ctr[t],symp_ctr[t],asymp_ctr[t]=update_human(humans,P)
        cb(1) # increase the progress metre by 1.. callback function
    end
    first_inf = find(x-> x.WhoInf == initial && x.WentTo == SYMP,humans)

    symp_inf = find(x -> x.WhoInf>0 && humans[x.WhoInf].WentTo == SYMP,humans)
    numb_symp_inf = length(symp_inf)

    asymp_inf = find(x -> x.WhoInf>0 && humans[x.WhoInf].WentTo == ASYMP,humans)
    numb_asymp_inf = length(asymp_inf)

    numb_first_inf = length(first_inf)

    Number_in_age_group = zeros(Int64,P.grid_size_human)
    NumberFailsAge = zeros(Int64,P.grid_size_human)
    InfOrNot = zeros(Int64,P.grid_size_human)
    VacStatus = zeros(Int64,P.grid_size_human)
    Infection_Matrix = zeros(Int64,15,15)

    for i = 1:P.grid_size_human

        Number_in_age_group[i] = humans[i].contact_group
        NumberFailsAge[i] = humans[i].NumberFails
        VacStatus[i] = humans[i].vaccinationStatus
        if humans[i].WhoInf > 0
            Infection_Matrix[humans[i].contact_group,humans[humans[i].WhoInf].contact_group]+=1
            if humans[i].health == REC || humans[i].health == SYMP || humans[i].health == ASYMP || humans[i].health == LAT
                InfOrNot[i] = 1
            end
        end
    end


    #return latent_ctr,symp_ctr,asymp_ctr,numb_first_inf,numb_symp_inf,numb_asymp_inf,Infection_Matrix,Fail_Contact_Matrix,Contact_Matrix_General,Number_in_age_group,NumberFailsAge,Risk_Contact2,Risk_Contact3,Proportion
    return latent_ctr,symp_ctr,asymp_ctr,numb_first_inf,numb_symp_inf,numb_asymp_inf,Infection_Matrix,Fail_Contact_Matrix,Contact_Matrix_General,Number_in_age_group,NumberFailsAge,InfOrNot,VacStatus

end
