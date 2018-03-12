using Parameters
using Distributions
using StatsBase

#@enum AGENT_TYPE UNDEF=0 PATIENT=1 HCW=2 

abstract type Agent end
mutable struct Patient{T<:Integer} <: Agent
    id::T
    Patient(id = 0) = new(id)
end
mutable struct HCW{T<:Integer} <: Agent
    id::T
    HCW(id = 0) = new(id)
end

## return type Core.Compiler.return_type(read, (String, Type{String}))

## main system parameters
@with_kw mutable struct ModelParameters     
    gridsize::Int64 = 100
end

MP = ModelParameters()
psize = MP.gridsize
## generate patient array 
patients = Array{Patient}(psize)
hcws = Array{HCW}(psize)
[patients[i] = Patient{Int64}(i) for i=1:psize]
[hcws[i] = HCW{Int64}(i) for i=1:psize]

## generate patient contact matrix, memory structure is column wise
pmm = Matrix{Vector{Int64}}(psize, psize)
for i = eachindex(pmm)
    pmm[i] = zeros(Int8, 14)
end
## OR
reshape([zeros(Int8, 14) for i=1:100], 10, 10)

## or use BIT ARRAY
pmm = Matrix{BitArray}(psize,psize)
for i = eachindex(pmm)
    pmm[i] = BitArray(14)
end

pmm = reshape([BitArray(14) for i=1:*(psize, psize)], psize, psize)


## now go through each patient, get the number of distinct contacts, and repeated contacts per patient
for patientid =1:100
    dpc = Int(round(rand(PP_DD))) ## get number of distinct contacts of this patient
    elig = sample(1:psize, dpc, replace=false) ## who will patient meet (IDs)
    deleteat!(elig, elig .== patientid )       ## remove from himself, cant meet himself

    ## for each contact patient will meet, get the repeat contacts and assign to main pmm matrix
    map(x -> begin
        pmm[patientid, x] = repeat_contacts()
        end, elig)
end

function repeat_contacts()
    rcon = Int(round(rand(PP_RD))) ##get the numer of repeated
    m = zeros(Int8, 14)  ## an array of zeros of length 14. counts the number of times contact will meet everyday over 14 days
    randcols = rand(1:14, rcon) ## randomly pick 'rcon' days to meet
    map(j->m[j] += 1, randcols) ## increment the number of meeting times
    return m
end



MP = ModelParameters()
psize = MP.gridsize


## calculate for each distinct contact, the number of repeated contacts over two weeks


dpc = [Int(round(rand(PP_DD))) for i=1:psize];

for i=1:psize
    ppm = zeros(Int64, dpc[i], 14)
    randcols = rand(1:14, dpc[i])
    [ppm[j, randcols[j]] = 1 for j=1:dpc[i]]
end

## contact sampling  - patient/patient
ppc_dist = Poisson(2)
#numcts = Array{Tuple{Int64, Int64}}(psize)  ## for each patient, get the number of distinct contacts and repeated: tuple is (distinct, repeated).
#[numcts[i] = (rand(ppc_dist), rand(ppc_dist)) for i=1:psize]

pp_dist = LogNormal(0.590104, 1.2629)
distinctcts = [rand(ppc_dist) for i = 1:psize]

#distinctmat = Array{Tuple{Int64, Float64}}(psize, psize)
timemat = zeros(Float64, psize, psize)
for i=1:psize
    numofcontacts = rand(ppc_dist)
    who = rand(0:psize, numofcontacts)
    
    
end



tts = Array{Tuple{Int64, Int64}}(psize)

## everyday sample contacts and times.

rand(dist)
a = Array{Float64}(10000)
[a[i] = rand(dist) for i = 1:10000]


function main()    
    
end


a = Patient(2, 4)

b = Array{Agent}(4)
b[1] = Patient(1, 2)
b[2] = Worker(1, 2, 3)