using Parameters
using Distributions
using StatsBase

#@enum AGENT_TYPE UNDEF=0 PATIENT=1 HCW=2 

## create entities
abstract type Agent end
mutable struct Patient{T<:Integer} <: Agent
    id::T
    Patient(id = 0) = new(id)
end
mutable struct HCW{T<:Integer} <: Agent
    id::T
    HCW(id = 0) = new(id)
end

## main system parameters
@with_kw struct ModelParameters     
    numpatients::Int64 = 100   ## number of patients in the system
    numhcws::Int64 = 20        ## number of healthcare workers in the system
end

function main(MP::ModelParameters)
  ## TODO divide n to number of patients, number of hcws
  np = MP.numpatients
  nh = MP.numhcws

  ## generate and allocate patient array 
  patients = [Patient{Int64}(i) for i=1:np]
  hcws = [HCW{Int64}(i) for i=1:nh]

  ## create distinct contacts distributions
  const PP_DD = Normal(4.6, 0.78)
  const HP_DD = Normal(5.3885, 0.73151)
  const HH_DD = Normal(3.28, 0.53)

  ## repeat distributions
  const PP_RD = Normal(5.3885, 0.73151)
  const HP_RD = Normal(9.259, 1.0)
  const HH_RD = Normal(8.321, 1.519)
  
  ## contact matrix allocation
  ## this is a matrix where each row number represents an agent. each column number represents the contact the agent will meet. 
  ## for every agent, the number of distinct contacts he will meet comes from the distributions defined above. 
  ## each element of the matrix is a vector of length 14. Each element in this vector is the number of times the agent will meet contact over 14 days.   
  ppm = reshape([zeros(Int8, 14) for i=1:*(np,np)], np, np)
  hpm = reshape([zeros(Int8, 14) for i=1:*(nh,np)], nh, np)
  hhm = reshape([zeros(Int8, 14) for i=1:*(nh,nh)], nh, nh)

  

end

function generate_contact_matrix(c_array, DD::Distribution, RD::Distribution)
  nrows = size(c_array)[1]
  ncols = size(c_array)[2]
  ## go through each agent
  for id=1:nrows 
    dpc = Int(round(rand(DD)))                    ## get number of distinct contacts of this agent    
    dpc = dpc > ncols ? ncols : dpc               ## check if dpc is bigger then ncols 
    elig = sample(1:ncols, dpc, replace=false)    ## sample the actual array to get IDs of  the 
    deleteat!(elig, elig .== id )                 ## remove from himself, cant meet himself
    ## get the repeated contacts per distinct contact and allocate over 14 days
    map(x -> c_array[id, x] = repeat_contacts(RD), elig)
  end
end

function repeat_contacts(D::Distribution)
  # for every distinct contact for a patient, how many times will this happen over 14 days? 
  rcon = Int(round(rand(D))) 
  # distribute the repeated contact over 2 weeks
  m = zeros(Int8, 14)         ## allocate \empty array
  randcols = rand(1:14, rcon) ## randomly pick 'rcon' days to meet
  map(j->m[j] += 1, randcols) ## increment the number of meeting times
  return m
end
