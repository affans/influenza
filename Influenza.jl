module InfluenzaModel

using Parameters      ## with julia 1.1 this is now built in.
using ProgressMeter   ## can now handle parallel with progress_pmap
#using PmapProgressMeter
#using DataFrames
using Distributions
using StatsBase
using StaticArrays
using Random
#using BenchmarkTools

export init, main, Human, InfluenzaParameters

include("parameters.jl")
include("population.jl")
include("functions.jl")

function init()
    [Human(i) for i = 1:10000]
end

function main(simnum::Int64, P::InfluenzaParameters)
    Random.seed!(simnum)
    #println("starting simulation number: $simnum")
    #println("transmission: $(P.transmission_beta)")
    humans = init()
    setup_demographic(humans)  ## TO DO, unit tests, plotting  
    apply_vaccination(humans,P)    ## TO DO, unit tests, plotting
    
    initial = setup_rand_initial_latent(humans,P) ## returns the ID of the initial person

    ## data collection variables = number of elements is the time units. 
    ## so the vector collects number of latent/symp/asymp at time t. 
    ## it does not collect the initial latent case.
    latent_ctr = zeros(Int64, P.sim_time)   
    symp_ctr =   zeros(Int64, P.sim_time)   
    asymp_ctr =  zeros(Int64, P.sim_time)   
 
    ## contact matrices
    ## these matrices are used to calculate contact patterns
    Fail_Contact_Matrix    = zeros(Int64, 15, 15)    ## how many times did susc/sick contact group i meet contact group j meet but failed to infect.
    Contact_Matrix_General = zeros(Int64, 15, 15)    ## how many times did contact group i meet with contact group j
    Number_in_age_group    = zeros(Int64, 15)                      # vector that tells us number of people in each age group.
    Age_group_Matrix       = zeros(Int64, 15, P.grid_size_human)   # a matrix representation of who is inside that age group (ie. row 1 has all the people that have group 1). 
       
    ## this function just fills in the empty matrices as defined above.
    setup_contact_matrix(humans, Age_group_Matrix, Number_in_age_group)

    NB = N_Binomial()
    CM = ContactMatrixFunc()

    ## main simulation loop.
    for t=1:P.sim_time        
        contact_dynamic2(humans, P, NB, CM, Fail_Contact_Matrix, Age_group_Matrix, Number_in_age_group, Contact_Matrix_General)
        for i=1:P.grid_size_human
           increase_timestate(humans[i], P)
        end      
        latent_ctr[t], symp_ctr[t], asymp_ctr[t] = update_human(humans,P)
    end

    ## find all the humans that went to symptomatic after being infected by the initial latent case.
    first_inf = findall(x-> x.WhoInf == initial && x.WentTo == SYMP, humans)
    first_inf_asymp =  findall(x-> x.WhoInf == initial && x.WentTo == ASYMP, humans)
    
    ## how many total sickness did sympomatics produce?
    ## how many total sickness did asymptomatics produce? 
    ##  -- this will be less 1 from the total sickness.. for some reason
    symp_inf  = findall(x -> x.WhoInf > 0 && humans[x.WhoInf].WentTo == SYMP,  humans)
    asymp_inf = findall(x -> x.WhoInf > 0 && humans[x.WhoInf].WentTo == ASYMP, humans)
    
    numb_symp_inf = length(symp_inf)   ## the total number of people all symptomatics made sick.
    numb_asymp_inf = length(asymp_inf) ## the total number of people all asymptomatics made sick.
    numb_first_inf = length(first_inf) ## the number of people infected by the initial latent case.
    numb_first_inf_asymp = length(first_inf_asymp)



    contact_groups = zeros(Int64, P.grid_size_human)   ## just the contact groups of everyone.
    demographic_group = zeros(Int64, P.grid_size_human)
    number_of_fails = zeros(Int64, P.grid_size_human)  ## this property counts how many times susc i met a sick person and failed to get sick.
    vax_status = zeros(Int64,P.grid_size_human)        ## the vaccination status of individual i. 
    infection_matrix = zeros(Int64, 15, 15)          
    InfOrNot = zeros(Int64, P.grid_size_human)         ## at the end of the simulation, is the person still susceptible?
    SympOrNot = zeros(Int64, P.grid_size_human)        ## at the end of the simulation, if the person was infected, whether they went symptomatic/asymptomatic

    for i = 1:length(humans)        
        contact_groups[i] = humans[i].contact_group           
        demographic_group[i] = humans[i].group
        number_of_fails[i] = humans[i].NumberFails              
        vax_status[i] = humans[i].vaccineEfficacy > 0 ? 1 : 0  
        ## if humans[i] was infected (another way to check is for WentTo == SYMP/ASYMP) -- good way to test the model.
        if !(humans[i].health == SUSC)
            if humans[i].WhoInf > 0 
                infection_matrix[humans[i].contact_group, humans[humans[i].WhoInf].contact_group] += 1 
            end                        
            InfOrNot[i] = 1
            SympOrNot[i] = Int(humans[i].WentTo)            
        end
    end

    return latent_ctr, symp_ctr, asymp_ctr, 
    numb_first_inf, numb_symp_inf, numb_asymp_inf, 
    infection_matrix, Fail_Contact_Matrix, Contact_Matrix_General, 
    contact_groups, number_of_fails, InfOrNot, 
    vax_status, SympOrNot, demographic_group
end

end

