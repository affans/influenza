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
  np = 100#MP.numpatients
  nh = 10 #MP.numhcws

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

  ppm = reshape([zeros(Int8, 14) for i=1:*(np,np)], np, np);
  hhm = reshape([zeros(Int8, 14) for i=1:*(nh,nh)], nh, nh);
  hpm = reshape([zeros(Int8, 14) for i=1:*(nh,np)], nh, np);
  
  ## the following code will repeat every two weeks
  pp_tens = zeros(Int64, np, np, 14);
  hh_tens = zeros(Int64, nh, nh, 14);
  hp_tens = zeros(Int64, nh, np, 14);

  gen_contact_mat(ppm, pp_tens, PP_DD, PP_RD)
  gen_contact_mat(hhm, hh_tens, HH_DD, HH_RD)
  gen_contact_mat(hpm, hp_tens, HP_RD, HP_RD)

  get_daily_contact(pp_tens, 1)
  #ppm = zeros(Bool, 10, 10)
  #ppm = zeros(Bool, np, np)
  #hpm = zeros(Bool, nh, np)

  #ppm = Array{Tuple{Int64, Int64}}(10,10)
  #cartref = [(i, j)  for i=1:100, j=1:100]
end


function get_daily_contact(c_array, daynumber)
  for i=1:14
    slice = c_array[:, :, i]  ## get the day slice
  end
end

## to do remove the ppm array when testing is complete
function gen_contact_mat(c_array_ppm, c_array_tensor,  DD::Distribution, RD::Distribution)
  ## get the size of the incoming matrix 
  nrows = size(c_array_tensor)[1]
  ncols = size(c_array_tensor)[2]
  ntime = size(c_array_tensor)[3]

  ## are we working with a patient/hcw matrix? If so the logic changes a bit - and matrix does not have to be symmetrical. 
  hcw_patient = nrows == ncols ? false : true
  
  ## go through each agent
  for id=1:nrows
    elig = hcw_patient ? collect(1:ncols) : elig = collect((id+1):ncols)      
    dpc = min(length(elig), Int(round(rand(DD)))) 
    who = sample(elig, dpc, replace=false)    ## sample the actual array to get IDs of  the 
    ## get the repeated contacts per distinct contact and allocate over 14 days
    map(x -> begin a = repeat_contacts(RD); c_array_ppm[id, x] = a; c_array_tensor[id, x, :] = a; end, who)
    ## this ends up creating an upper triangular matrix
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

function contact_statistics_ppm(c_array)
  nrows = size(c_array)[1]
  ncols = size(c_array)[2]

  contact = zeros(Bool, nrows, ncols)
  counts = Array{Int64}(nrows)
  sums = zeros(Int64, nrows, ncols)
  for i = eachindex(c_array)
    if sum(c_array[i]) > 0 
      contact[i] = true
      sums[i] = sum(c_array[i])
    end
  end
  for j = 1:nrows
    counts[j] = length(find(contact[j, :] .== true))
  end 
  ## find(a[1, :] .== true) gets the actual IDs
  
  return contact, counts, sums
  ## LowerTriangular(sums) == 0 for testing purposes.
end

function contact_statistics_tensor(carr)
  # c[1, :, :] ## horizontal slice - who 1 will meet over 14 days
  # c[:, 2, :] ## vertical slice - who is meeting 2
  # c[:, :, 1] ## who is meeting who at time one
  ## we want to see count the number of distinct contacts 
  nrow = size(carr)[1]
  ncol = size(carr)[2]
  ntime = size(carr)[3]
  contacts = zeros(Bool, nrow, ncol)
  for i=1:nrow
    a = find(x -> x > 0, sum(carr[i, :, :], 2))  ## sum() (over the 2nd dimension) looks at each row and then sums up the meeting times. if > 0, then i will meet the person ID atleast once
    contacts[i, a] = true    
  end
  repeatcounts = zeros(Int64, nrow, ncol)
  for i=1:nrow
    a = sum(carr[i, :, :], 2)  ## sum() (over the 2nd dimension) looks at each row and then sums up the meeting times. if > 0, then i will meet the person ID atleast once
    repeatcounts[i, :] = a    
  end
  # lower triangular of this must all be false
  return contacts, repeatcounts
end