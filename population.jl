mutable struct Human{T <: Number} ## mutable structs are stored on the heap 
    strains_matrix::Array{Int8,2}
    Vector_time::Array{Int64,1}
    NumberStrains::Int64
    EfficacyVS::Float64

    index::T
    health::HEALTH
    swap::HEALTH #do we need  this? We can do a sequential atualization
    timeinstate::T
    statetime::T
    latenttime::Int64
    recoveredOn::Int64

    vaccinationStatus::Int64
    vaccineEfficacy::Float64
    

    WhoInf::T
    WentTo::HEALTH  # we want to know whether human was asymptomatic or symptomatic.. need this information for other calculations
    TimeGotInf::Int64

    age::Union{T, Nothing}    
    group::Union{T, Nothing}  
    contact_group::Union{T, Nothing}      

    daily_contacts::T
    Coverage::Float64
    NumberFails::T
    NumberExposures::Int64
    vac_cont::Int64
    unvac_cont::Int64
    vac_infec::Int64
    unvac_infec::Int64
    #Here I added the strain matrix, a vector for the time the strain showed up, number of strains in the body, and
    #the efficacy of the vaccine against the strain that was transmitted (necessary for the asymp-symp trial)
    Human{Int64}(P::InfluenzaParameters) = new(zeros(Int8,P.matrix_strain_lines,P.sequence_size),zeros(Int64,P.matrix_strain_lines),0,0,-1, SUSC, UNDEF, 0, 999,0,-1,0, 0.0, -1, UNDEF,-1,
                            nothing, nothing, nothing, 0, 0.0, 0, 0, 0, 0, 0, 0)
    function Human(idx,P::InfluenzaParameters)
        h = Human{Int64}(P)
        h.index = idx
        return h
    end
end

function setup_demographic(h)
    dist, AgeMin, AgeMax = distribution_age()
    for i = 1:length(h)
        if h[i].age == nothing
            rn = rand()
            g = findfirst(x -> rn <= x, dist) ## THERE IS A CLOSURE BUG HERE. DO NOT USE INBOUNDS. SEE FTEST() testing. This was fixed. See my issue on github
            if g !== nothing
                h[i].age = rand(AgeMin[g] : AgeMax[g])
                h[i].contact_group = min(g, 15)  ## maximum of 15 contact groups
            end
        end

        ###Group 1 - young child, 2 - school child, 3 - working adult, 4 - elderly
        if h[i].age <= 8
            h[i].group = 1
        elseif h[i].age <= 17
            h[i].group = 2
        elseif h[i].age <= 49
            h[i].group = 3
        elseif h[i].age <= 64
            h[i].group = 4
        else
            h[i].group = 5
        end
    end
end


# # INBOUNDS BUG
# function ftest()
#     dist_local = cumsum([rand() for i = 1:17])
#     dist_local = dist_local/maximum(dist_local)
#     @inbounds for i=1:10000
#         rn_local = rand()
#         newlocalvar = findfirst(x -> rn_local <= x, dist_local)
#     end
#     return nothing
# end
# @btime ftest()

# function ftest_noinbounds()
#     dist_local = cumsum([rand() for i = 1:17])
#     dist_local = dist_local/maximum(dist_local)
#     for i=1:10000
#         rn_local = rand()
#         newlocalvar = findfirst(x -> rn_local <= x, dist_local)
#     end
#     return nothing
# end
# @btime ftest_noinbounds()

# function ftest_suggestion()
#     dist_local = cumsum([rand() for i = 1:17])
#     dist_local = dist_local/maximum(dist_local)
#     @inbounds for i=1:10000
#         rn = rand()
#         let rn = rn
#              newlocalvar = findfirst(x -> rn <= x, dist_local)
#         end
#     end
#     return nothing
# end
# @btime ftest_suggestion()

# function bar()
#     x=[2]
#     @inbounds for i in 1:10_000 #Add @inbounds to this line to trrigger it
#     (m->2i)(x)
#     end
# end
# @btime bar()

function setup_rand_initial_latent(h, P::InfluenzaParameters,Original_Strain::Array{Int8,1},t::Int64,rng1)
    randperson = rand(1:P.grid_size_human)
    make_human_latent(h[randperson], P)
    h[randperson].WhoInf = 0 ## no one infected this person
    for i = 1:P.Number_of_initial_strains
        h[randperson].strains_matrix[i,:] = Original_Strain
        h[randperson].strains_matrix[i,:] = Changing_a_proportion(h[randperson].strains_matrix[i,:],(i-1)*P.AGD_Step,P,rng1)
    end

    h[randperson].NumberStrains = P.Number_of_initial_strains-1
    return randperson
end

@inline function make_human_latent(h::Human, P::InfluenzaParameters)
    ## make the i'th human infected
    h.health = LAT
    h.swap = UNDEF
    h.statetime = rand(P.Latent_period_Min:P.Latent_period_Max)
    h.latenttime = h.statetime
    h.timeinstate = 0    
end

@inline function make_human_asymp(h::Human, P::InfluenzaParameters,rng1)
    ## make the i'th human infected
    h.health = ASYMP    # make the health ->inf
    h.swap = UNDEF
    d = LogNormal(P.log_normal_mean,sqrt(P.log_normal_shape))
    h.statetime = min(P.max_infectious_period,ceil(rand(d)))
    h.timeinstate = 0
    h.WentTo = ASYMP

    h.NumberStrains = h.NumberStrains + 1
    h.strains_matrix,h.Vector_time,h.NumberStrains = mutation(h.strains_matrix,P,h.NumberStrains,h.statetime,h.latenttime,rng1)
end

@inline function make_human_symp(h::Human, P::InfluenzaParameters,rng1)
    ## make the i'th human infected
    h.health = SYMP    # make the health ->inf
    h.swap = UNDEF
    d = LogNormal(P.log_normal_mean,sqrt(P.log_normal_shape))
    h.statetime = min(P.max_infectious_period,ceil(rand(d)))
    h.timeinstate = 0
    h.WentTo = SYMP

    h.NumberStrains = h.NumberStrains + 1
    h.strains_matrix,h.Vector_time,h.NumberStrains = mutation(h.strains_matrix,P,h.NumberStrains,h.statetime,h.latenttime,rng1)
end

@inline function make_human_recovered(h::Human, P::InfluenzaParameters)
    ## make the i'th human recovered
    h.health = REC    # make the health -> latent
    h.swap = UNDEF
    h.statetime = 999
    h.timeinstate = 0
end


function increase_timestate(h::Human,P::InfluenzaParameters)
    h.timeinstate+=1
    if h.timeinstate >= h.statetime
        if h.health == ASYMP || h.health == SYMP
            h.swap = REC
        elseif h.health == LAT
            prob = (P.ProbAsympMax - P.ProbAsympMin)*rand()+P.ProbAsympMin
            if rand() < (1 - prob) * (1 - h.EfficacyVS)
                h.swap = SYMP
            else h.swap = ASYMP
            end
        end
    end
end

function update_human(h, P::InfluenzaParameters,rng1,t::Int64)
    n1::Int64 = 0
    n2::Int64 = 0
    n3::Int64 = 0
    
    n_lat_vac::Int64 = 0
    n_lat_nvac::Int64 = 0

    mean_exp_vac::Float64 = 0.0
    mean_exp_nvac::Float64 = 0.0

    total_contact_vac::Int64 = 0
    total_contact_nvac::Int64 = 0
    total_contact::Int64 = 0
    total_infections_vac::Float64 = 0.0
    total_infections_nvac::Float64 = 0.0
    total_infections::Float64 = 0.0

    for i=1:length(h)
        
        if h[i].swap == LAT
            make_human_latent(h[i],P)
            h[i].TimeGotInf = t
            n1+=1
            if h[i].vaccinationStatus == 1
                n_lat_vac+=1
                mean_exp_vac += h[i].NumberExposures 
            else
                n_lat_nvac+=1
                mean_exp_nvac += h[i].NumberExposures 
            end

        elseif h[i].swap == SYMP
            
            make_human_symp(h[i],P,rng1)
            n2+=1
        elseif h[i].swap == ASYMP
             
            make_human_asymp(h[i],P,rng1)
            n3+=1
        elseif h[i].swap == REC
            make_human_recovered(h[i],P)
            h[i].recoveredOn = t
            total_contact += h[i].vac_cont+h[i].unvac_cont
            total_contact_nvac += h[i].unvac_cont
            total_contact_vac += h[i].vac_cont
            total_infections += h[i].vac_infec+h[i].unvac_infec
            total_infections_nvac += h[i].unvac_infec
            total_infections_vac += h[i].vac_infec

        end
    end

    total_infections = total_infections/total_contact
    total_infections_nvac = total_infections_nvac/total_contact_nvac
    total_infections_vac = total_infections_vac/total_contact_vac


    mean_exp_nvac = mean_exp_nvac/n_lat_nvac
    mean_exp_vac = mean_exp_vac/n_lat_vac

    return n1, n2, n3, n_lat_vac, n_lat_nvac, mean_exp_vac, mean_exp_nvac, total_infections_vac,total_infections_nvac,total_infections #corresponds to latent, symp, asymp
end

function frailty(age)
    if age <= 34
        y = 0.26875-0.00435*age
    elseif age <= 69
        y = 0.01282*age-0.30658
    else
        y = 0.396+0.0039*age
    end
    min1 = max((y-0.05),0.0)

    if y+0.05 > 1.0 ##function min was not working, so I did it manually
        max1 = 1.0
    else
        max1 = y+0.05
    end
    return max1,min1
end

@inline function _apply_vax(bracket, conf, h, P)
    ## this is a helper function. 
    if !(typeof(bracket) == Tuple{Int64, Int64})
        error("age must be a tuple")
    end
    if !(typeof(conf) == Tuple{Int64, Int64})
        error("confidence interval must be a tuple")
    end
    g = findall(x -> x.age > bracket[1] && x.age <= bracket[2], h)    
    howmany = Int(round(rand(conf[1]:conf[2])/100*length(g)))
    whotovax = sample(g, howmany, replace=false) #rand(g, howmany) DO NOT USE RAND.. it samples with replacement
    for i in whotovax
        MaxFra,MinFra = frailty(h[i].age)
        FrIndex = rand()*(MaxFra-MinFra)+MinFra
        h[i].vaccineEfficacy = round(P.vaccine_efficacy*(1.0-FrIndex), digits = 3)
        h[i].vaccinationStatus = 1
    end
end

function apply_vaccination(h, P)
    ## first tuple argument is the age bracket (0, 4): 0 is not included, 4 is included in the if statement
    ## the second tuple is the coverage confidence interval.. select a number between 20-32% of coverage for 0-4 year olds.
    #for (ag,cov) in zip([(0,4), (4,49), (49,64), (64,100)], [(20,32),(18,27),(34,42),(65,73)])
    #    _apply_vax(ag, cov, h, P)
    #end
    _apply_vax((0, 4), (20, 32), h, P)
    _apply_vax((4, 49), (18, 27), h, P)
    _apply_vax((49, 64), (34, 42), h, P)
    _apply_vax((64, 100), (65, 73), h, P)
end

function setup_contact_matrix(h, agm, nagm)
    for i = 1:length(h)
        cg = h[i].contact_group
        nagm[cg] += 1
        agm[cg, nagm[cg]] = h[i].index        
    end
end

