using Parameters
using Distributions

#@enum AGENT_TYPE UNDEF=0 PATIENT=1 HCW=2 

## main system parameters
@with_kw type ModelParameters     
    gridsize::Int64 = 100
end

abstract type Agent end
mutable struct Patient{T<:Integer} <: Agent
    id::T
    Patient(id = 0) = new(id)
end
mutable struct HCW{T<:Integer} <: Agent
    id::T
    HCW(id = 0) = new(id)
end

psize = 90

patients = Array{Patient}(psize)
hcws = Array{HCW}(psize)
[patients[i] = Patient{Int64}(i) for i=1:psize]
[hcws[i] = HCW{Int64}(i) for i=1:psize]

## contact sampling  - patient/patient
ppc_dist = Poisson(2)
cts = Array{Tuple{Int64, Int64}}(psize)
tts = Array{Tuple{Int64, Int64}}(psize)
[cts[i] = (rand(ppc_dist), rand(ppc_dist)) for i=1:psize]

## everyday sample contacts and times.
pp_dist = LogNormal(0.590104, 1.2629)
rand(dist)
a = Array{Float64}(10000)
[a[i] = rand(dist) for i = 1:10000]


function main()    
    
end


a = Patient(2, 4)

b = Array{Agent}(4)
b[1] = Patient(1, 2)
b[2] = Worker(1, 2, 3)