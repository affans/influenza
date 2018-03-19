using Parameters
using Distributions
using StatsBase
using BenchmarkTools

@enum HEALTH UNDEF=0 SUSC=1 EXP=2 ASYMP=3 SYMP=4 REC=5 

## create entities
abstract type Agent end
mutable struct Patient{T<:Integer} <: Agent
    id::T
    health::T
    swap::T
    hosp::T
    Patient(id=0, health=SUSC, swap=UNDEF, hosp=0) = new(id, health, swap, hosp)
end
mutable struct HCW{T<:Integer} <: Agent
    id::T
    health::T
    swap::T
    HCW(id=0, health=SUSC, swap=UNDEF) = new(id, health, swap)
end

function gen_contact_mat(c_array_tensor, c_type, DD::Distribution, RD::Distribution)
  ## c_array_tensor: the incoming tensor to be populated
  ## c_type:         do we have Patient/Patient (PP) HCW/HCW (HH) or Patient/HCW (PH)
  ## DD:             number of distinct contacts - distribution
  ## RD:             for each distinct contact, the distribution of how many times this contact will be repeated.

  ## get the size of the incoming matrix dimensions -
  nrows = size(c_array_tensor)[1]
  ncols = size(c_array_tensor)[2]
  ntime = size(c_array_tensor)[3] ## should be 14
  
  ## go through each agent
  for id=1:nrows
    ## are we working with a patient/hcw matrix? If so the logic changes a bit - and matrix does not have to be symmetrical.
    # eligible selects the columns that the current row (id) can possibly meet
    elig = c_type == "PH" ? collect(1:ncols) : elig = collect((id+1):ncols)        
    # sample a random number from the distribution for distinct contacts. ofcourse the minimum number of distinct contacts is the number of elig. 
    dpc = min(length(elig), Int(round(rand(DD)))) 
    # now that we know how many distinct contacts, sample from elig without replacement (otherwise you are repeating)
    who = sample(elig, dpc, replace=false)    ## sample the actual array to get IDs of  the 
    ## get the repeated contacts per distinct contact and allocate over 14 days
    map(x -> begin a = repeat_contacts(RD); c_array_tensor[id, x, :] = a; end, who)    
  end
end

function repeat_contacts(D::Distribution)
  ## this function samples a number from a distribution and distributes this number of an arrary of 14 elements
  rcon = Int(round(rand(D))) 
  m = zeros(Int8, 14)         ## allocate \empty array
  randcols = rand(1:14, rcon) ## randomly pick 'rcon' days to meet
  map(j->m[j] += 1, randcols) ## increment the number of meeting times
  return m
end

function disease_dynamics(c_array, ctype, rowagents, colagents, dn) ##dn = daynumber (1-14)
  ## this function is run everyday.
  ## it takes the contact tensor and calculates disease probability for every contact
  ## it goes down the rows of rowagents and makes them meet colagents
  
  ## error check
  (dn > 14 || dn < 1) && error("dn must be from 1-14")

  ## debug info
  println("Debug Info")
  println("Mem of rowagents: $(pointer_from_objref(rowagents))")
  println("Mem of rowagents: $(pointer_from_objref(colagents))")

  slice = c_array[:, :, dn]  ## get the day slice
  nrows = size(slice)[1]
  ncols = size(slice)[2]
  ## go through the rows/columns. if number > 0, thats how many times they are meeting on this particular day. 
  for j=1:ncols, i=1:nrows    
    agent = rowagents[i]
    whotomeet = colagents[j]
    ## check if the agent captured from the row/col iteration matches what we picked up 
    agent.id == i || error("agent id didn't match")
    whotomeet.id == j || error("whotomeet id didn't match")
  
    ## meeting will happen this many times
    mt = slice[i, j]
    if mt > 0       
      transmission_result = transmission(mt, agent, whotomeet, ctype)     
      
    end    
  end
end


function transmission(mt::Int64, a1::Agent, a2::Agent, ctype)  
  ## calculates the transmission between two agents 
  ## if successfull transmission, the agent is switched to the compartment. 
  ## it considers all possible reductions to transmission values.
  ## mt - # of times contact a1, a2 will meet in a single day. 
  ## ctype = "PP", "HH", "PH"
  
  beta = 0.0
  prob = 0.0
  tflag = false
 
  if ctype == "HP" || ctype == "PH"
    ## we have patient/hcw contact

    ## get the right distribution
    D = LogNormal(0.232197, 0.987036)

    ## h1 determines if one of them is sick, h2 determines if one of them is susceptible 
    h1 = HEALTH(a1.health) == ASYMP || HEALTH(a1.health) == SYMP || HEALTH(a2.health) == ASYMP || HEALTH(a2.health) == SYMP
    h2 = HEALTH(a1.health) == SUSC || HEALTH(a2.health) == SUSC
    
    # note:
    # if the health care worker ends up being symptomatic, nothing happens as they go home
    ## this is accounted for when healthcare worker goes to symptomatic - they are automatically put into recovered case
    ## so h1 can only be true if the sympomatic agent is a patient

    tflag = h1 && h2
    if tflag      
      if HEALTH(a1.health) == ASYMP || HEALTH(a1.health) == SYMP
        sick = a1
        susc = a2
      else
        sick = a2
        susc = a1
      end
    end
  else 
    ## we are in a PP/HHsituation
    
    ## get the right distribution
    D = ctype == "PP" ? LogNormal(0.590104, 1.2629) : LogNormal(0.433809, 1.07664)

    h1 = HEALTH(a1.health) == ASYMP || HEALTH(a2.health) == ASYMP 
    h2 = HEALTH(a1.health) == SUSC || HEALTH(a2.health) == SUSC
    tflag = h1 && h2
    if tflag 
      if HEALTH(a1.health) == ASYMP 
        sick = a1
        susc = a2
      else
        sick = a2
        susc = a1
      end
    end    
  end

  ## debug 
  #println("Contact is between $HP")
  #println("Distribution is $D")
  #tflag && println("sick is $sick, susc is $susc")
  
  if tflag
    unittimes = sum(rand(D, mt))
    prob = 1 - (1 - beta)^unittimes
  end
  return prob
end

function main()
  ## PARAMETERS
  ## TODO divide n to number of patients, number of hcws
  np = 20
  nh = 10 
  
  ## create distinct contacts distributions
  const PP_DD = Normal(4.6, 0.78)
  const HP_DD = Normal(5.3885, 0.73151)
  const HH_DD = Normal(3.28, 0.53)

  ## repeat distributions
  const PP_RD = Normal(5.3885, 0.73151)
  const HP_RD = Normal(9.259, 1.0)
  const HH_RD = Normal(8.321, 1.519)

  ## generate and allocate patient array 
  patients = [Patient{Int64}(i) for i=1:np]
  hcws = [HCW{Int64}(i) for i=1:nh]
  
  systime = 0 ## counter for system time - needed for keeping track of statetimes
  
  ## contact matrix allocation
  ## this is a matrix where each row number represents an agent. each column number represents the contact the agent will meet. 
  ## for every agent, the number of distinct contacts he will meet comes from the distributions defined above. 
  ## each element of the matrix is a vector of length 14. Each element in this vector is the number of times the agent will meet contact over 14 days.  

  newweek = true
  for systime = 1:100 ## 100 days total
    ## the following code will repeat every two weeks
    if newweek
      ## create contact tensors. 
      ## first dimension/second dimension = individuals
      ## third dimension - slices of time, 14 slices for 14 days.       
      pp_tens = zeros(Int64, np, np, 14);
      hh_tens = zeros(Int64, nh, nh, 14);
      hp_tens = zeros(Int64, nh, np, 14);

      gen_contact_mat(pp_tens, "PP", PP_DD, PP_RD)
      gen_contact_mat(hh_tens, "HH", HH_DD, HH_RD)
      gen_contact_mat(hp_tens, "PH", HP_RD, HP_RD)
    end

    ## we have to do this because systime % 14  at the 14th day = 0
    day = systime % 14 > 0 ? systime : 14 

    ## daily disease dynamics between the different groups. 
    disease_dynamics(pp_tens, "PP", patients, patients,  day)
    disease_dynamics(hh_tens, "HH", hcws, hcws, day)
    disease_dynamics(hp_tens, "PH", hcws, patients, day)
    
    ## end of the day, determine if we are at a new week
    newweek = systime % 14 == 0 
  end
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