mutable struct Human{T <: Number} ## mutable structs are stored on the heap 
    index::T
    health::HEALTH
    swap::HEALTH #do we need  this? We can do a sequential atualization
    timeinstate::T
    statetime::T

    #vaccinationStatus::T
    vaccineEfficacy::Union{Float64, Nothing}

    WhoInf::T
    WentTo::HEALTH  # we want to know whether human was asymptomatic or symptomatic.. need this information for other calculations

    age::Union{T, Nothing}    ##used
    group::Union{T, Nothing}  ## used
    contact_group::Union{T, Nothing}        ##used

    daily_contacts::T
    Coverage::Float64
    NumberFails::T

    Human{Int64}() = new(-1, SUSC, UNDEF, 0, 999, nothing, -1, UNDEF,
                            nothing, nothing, nothing, 0, 0.0, 0)
    function Human(idx)
        h = Human{Int64}()
        h.index = idx
        return h
    end
end

function setup_demographic(h)
    dist, AgeMin, AgeMax = distribution_age()
    for i = 1:length(h)
        if h[i].age == nothing
            rn = rand()
            g = findfirst(x -> rn <= x, dist) ## THERE IS A CLOSURE BUG HERE. DO NOT USE INBOUNDS. SEE FTEST() testing.
            if g !== nothing
                h[i].age = rand(AgeMin[g] : AgeMax[g])
                h[i].contact_group = min(g, 15)  ## maximum of 15 contact groups
            end
        end

        ###Group 1 - young child, 2 - school child, 3 - working adult, 4 - elderly
        if h[i].age<=4
            h[i].group = 1
        elseif h[i].age<=19
            h[i].group = 2
        elseif h[i].age<= 59
            h[i].group = 3
        else
            h[i].group = 4
        end
        ## TO DO Remove.
        #h[i].Coverage = get_vaccine_coverage(h[i])
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

##
# function get_vaccine_coverage(h::Human)
#     ## DEPRECATED
#     ## vaccine coverage at the age group level.
#     ## may not need this assigned at the individual level.
#     ## TODO: change this and remove the struct field
#     if h.age<=4
#         Coverage = 0.26
#     elseif h.age<= 49
#         Coverage = 0.23
#     elseif h.age<= 64
#         Coverage = 0.38
#     else
#         Coverage = 0.70
#     end
#     return Coverage
# end

function setup_rand_initial_latent(h, P::InfluenzaParameters)
    randperson = rand(1:P.grid_size_human)
    make_human_latent(h[randperson], P)
    return randperson
end

function make_human_latent(h::Human, P::InfluenzaParameters)
    ## make the i'th human infected
    h.health = LAT
    h.swap = UNDEF
    h.statetime = rand(P.Latent_period_Min:P.Latent_period_Max)
    h.timeinstate = 0
end

function make_human_asymp(h::Human, P::InfluenzaParameters)
    ## make the i'th human infected
    h.health = ASYMP    # make the health ->inf
    h.swap = UNDEF
    d = LogNormal(P.log_normal_mean,sqrt(P.log_normal_shape))
    h.statetime = min(15,ceil(rand(d)))
    h.timeinstate = 0
    h.WentTo = ASYMP
end

function make_human_symp(h::Human, P::InfluenzaParameters)
    ## make the i'th human infected
    h.health = SYMP    # make the health ->inf
    h.swap = UNDEF
    d = LogNormal(P.log_normal_mean,sqrt(P.log_normal_shape))
    h.statetime = min(15,ceil(rand(d)))
    h.timeinstate = 0
    h.WentTo = SYMP
end

function make_human_recovered(h::Human, P::InfluenzaParameters)
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
            if rand() < (1 - prob) * (1 - h.vaccineEfficacy)
                h.swap = SYMP
            else h.swap = ASYMP
            end
        end
    end
end

function update_human(h, P::InfluenzaParameters)
    n1::Int64 = 0
    n2::Int64 = 0
    n3::Int64 = 0
    for i=1:length(h)
        if h[i].swap == LAT
            make_human_latent(h[i],P)
            n1+=1
        elseif h[i].swap == SYMP
            make_human_symp(h[i],P)
            n2+=1
        elseif h[i].swap == ASYMP
            make_human_asymp(h[i],P)
            n3+=1
        elseif h[i].swap == REC
            make_human_recovered(h[i],P)
        end
    end
    return n1, n2, n3 #corresponds to latent, symp, asymp
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

function apply_vaccination(h, P)
    numofvax = 0
    for i = 1:length(h)
        rn = rand()
        lage = h[i].age
        vaxon = false
        if lage <= 4 && rn < 0.26
            vaxon = true
        elseif (lage > 4 && lage <= 49) && rn < 0.23
            vaxon = true
        elseif (lage > 49 && lage <= 64) && rn < 0.38
            vaxon = true
        elseif rn < 0.70
            vaxon = true
        end
        MaxFra,MinFra = frailty(lage)
        FrIndex = rand()*(MaxFra-MinFra)+MinFra
        if vaxon
            h[i].vaccineEfficacy = P.vaccine_efficacy*(1.0-FrIndex)
            numofvax = numofvax + 1
        else  ## else statement is neccessary since initially the parameter is initi
            h[i].vaccineEfficacy = 0.0
        end
    end
    return numofvax
end


function setup_contact_matrix(h, agm, nagm)
    for i = 1:length(h)
        cg = h[i].contact_group
        agm[cg, (nagm[cg]+1)] = h[i].index
        nagm[cg] += 1
    end
end