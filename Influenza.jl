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
include("mutation.jl")

function init(P::InfluenzaParameters)
    [Human(i,P) for i = 1:P.grid_size_human]
end

function main(simnum::Int64, P::InfluenzaParameters)
    Random.seed!(simnum)
    rng1 = MersenneTwister(simnum*100)
    #println("starting simulation number: $simnum")
    #println("transmission: $(P.transmission_beta)")
    humans = init(P)
    setup_demographic(humans)  ## TO DO, unit tests, plotting  
    apply_vaccination(humans,P)    ## TO DO, unit tests, plotting
    
    Vaccine_Strain = Vector{Int8}(undef,P.sequence_size)
    Creating_Vaccine_Vector(Vaccine_Strain,P,rng1)

    initial = setup_rand_initial_latent(humans,P,Vaccine_Strain,0,rng1) ## returns the ID of the initial person

    ## data collection variables = number of elements is the time units. 
    ## so the vector collects number of latent/symp/asymp at time t. 
    ## it does not collect the initial latent case.
    latent_ctr = zeros(Int64, P.sim_time)   
    symp_ctr =   zeros(Int64, P.sim_time)   
    asymp_ctr =  zeros(Int64, P.sim_time)   
    latent_vac_ctr = zeros(Int64, P.sim_time) 
    latent_nvac_ctr = zeros(Int64, P.sim_time) 

    exposures_vac_ctr = zeros(Float64, P.sim_time) 
    exposures_nvac_ctr = zeros(Float64, P.sim_time) 
   

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

    vac_prop_ctr = zeros(Float64, P.sim_time) 
    nvac_prop_ctr = zeros(Float64, P.sim_time) 
    total_prop_ctr = zeros(Float64, P.sim_time) 
    ## main simulation loop.
    for t=1:P.sim_time        
        contact_dynamic2(humans, P, NB, CM, Fail_Contact_Matrix, Age_group_Matrix, Number_in_age_group, Contact_Matrix_General,Vaccine_Strain,rng1)
        for i=1:P.grid_size_human
           increase_timestate(humans[i], P)
        end      
        latent_ctr[t], symp_ctr[t], asymp_ctr[t], latent_vac_ctr[t], latent_nvac_ctr[t], exposures_vac_ctr[t], exposures_nvac_ctr[t],vac_prop_ctr[t],
        nvac_prop_ctr[t],total_prop_ctr[t] = update_human(humans,P,rng1,t)
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
    pv = zeros(Int64,P.matrix_strain_lines)
    pnv = zeros(Int64,P.matrix_strain_lines)
    Ef = zeros(Float64,P.matrix_strain_lines)

   
    
   
    for i = 1:length(humans)        
        contact_groups[i] = humans[i].contact_group           
        demographic_group[i] = humans[i].group
        number_of_fails[i] = humans[i].NumberFails              
        vax_status[i] = humans[i].vaccineEfficacy > 0 ? 1 : 0 ##I don't like this system of checking
        #because humans[i].vaccineEfficacy is a float and, sometimes, the compiler might understand it > 0
        auxP1::Int64 = -1
        ## if humans[i] was infected (another way to check is for WentTo == SYMP/ASYMP) -- good way to test the model.
        if !(humans[i].health == SUSC)
            if humans[i].WhoInf > 0 
                
                infection_matrix[humans[i].contact_group, humans[humans[i].WhoInf].contact_group] += 1 
                auxP = Int64(Calculating_Distance_Two_Strains(Vaccine_Strain,humans[i].strains_matrix[1,:]))+1
               # println([i auxP])
               
                if humans[i].vaccinationStatus == 1
                    pv[auxP]+=1
                else
                    pnv[auxP]+=1
                end
                Ef[auxP] += humans[i].EfficacyVS
            end                        
            InfOrNot[i] = 1
            SympOrNot[i] = Int(humans[i].WentTo)            
        end


    end

    DistTime = Vector{Union{Nothing,Int64}}(nothing,P.sim_time)

    for t = 1:P.sim_time
        aux = findall(x-> x.TimeGotInf<t && x.recoveredOn>t , humans)

        if length(aux) > 0
            DistTime[t] = 0
            for i in aux
                for j = 1:humans[i].NumberStrains
                    if (humans[i].Vector_time[j]+humans[i].TimeGotInf) <= t
                        Distance = Calculating_Distance_Two_Strains(Vaccine_Strain,humans[i].strains_matrix[j,:])
                        if DistTime[t] < Distance
                            DistTime[t] = Distance
                        end
                    end
                end

            end
        else
            DistTime[t] = -1
        end

    end

   

    return latent_ctr, symp_ctr, asymp_ctr, 
    numb_first_inf, numb_symp_inf, numb_asymp_inf, 
    infection_matrix, Fail_Contact_Matrix, Contact_Matrix_General, 
    contact_groups, number_of_fails, InfOrNot, 
    vax_status, SympOrNot, demographic_group,
    Ef,pv,pnv,DistTime, latent_vac_ctr, latent_nvac_ctr, 
    exposures_vac_ctr, exposures_nvac_ctr,vac_prop_ctr,nvac_prop_ctr,total_prop_ctr
end

end

