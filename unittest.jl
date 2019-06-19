using DataArrays, DataFrames
using Match
using ParallelDataTransfer
using QuadGK
using Parameters #module
using Distributions
using StatsBase
using Base.Profile
include("parameters.jl")
include("PopStruct.jl")
include("functions.jl")


ht = init()
P = InfluenzaParameters()
@btime init()
@btime setup_demographic($ht)
@btime apply_vaccination($ht, $P)
function time_setup_contact(h)
    Fail_Contact_Matrix    = zeros(Int64, 15, 15)
    Contact_Matrix_General = zeros(Int64, 15, 15)
    Number_in_age_group    = zeros(Int64, 15)
    Age_group_Matrix       = zeros(Int64, 15, 10000)
    setup_contact_matrix(h, Age_group_Matrix, Number_in_age_group)
end
@btime time_setup_contact($ht)
@btime  increase_timestate(ht[2],P)
@btime main_test()

sigma1 = 1

ef1 = 0.4



humans = Array{Human}(P.grid_size_human)

setup_human(humans)
setup_demographic(humans,P)
vaccination(humans,P)

Number_in_age_group = zeros(Int64,15)
Age_group_Matrix = Matrix{Int64}(15,P.grid_size_human)
for i = 1:P.grid_size_human
    Age_group_Matrix[humans[i].contact_group,(Number_in_age_group[humans[i].contact_group]+1)] = humans[i].index
    Number_in_age_group[humans[i].contact_group] += 1
end

########################
f = open(string("VaccineEffectiveness/VaccineEffec","$ef1",".dat"),"w")
for j = 1:15

    A = find(x -> x.vaccinationStatus == 1 && x.contact_group == j,humans)

    soma = 0.0
    for k = A
        soma += humans[k].vaccineEfficacy

    end
    soma = soma/length(A)
    println(f,"$j $soma")
end
close(f)
end

#################
######erasing the age groups

for i=1:P.grid_size_human
    humans[i].contact_group = 1
end

function test_humanage()
    ag = map(x -> x.age, humans)
    plot(x = ag, Geom.histogram)

    a = find(x -> x.gender == FEMALE, humans)
    ag = zeros(Int64, length(a))
    for i=1:length(a)
        ag[i] = humans[a[i]].agegroup
    end

    ag = map(x -> x.agegroup, humans)

#    ag = map(x -> x.age, humans)
    plot(x = ag, Geom.histogram)

    find(x -> x.age >= 15, humans)
    find(x -> x.age >= 15 && x.gender == MALE, humans)
    find(x -> x.age >= 15 && x.gender == FEMALE, humans)

end

#### Testing LogNormal
d = LogNormal(1,sqrt(0.4356))
d1 = Vector{Float64}(3000)
for i = 1:length(d1)
d1[i] = min(15,rand(d))

end

writedlm("testeLogNormal.dat",d1)

##Contact distribution
X = Matrix{Float64}(15,15)
X[1,:] = ContactMatrix[1,:]
for i = 2:15
    for j = 1:15
        X[i,j] = ContactMatrix[i,j]-ContactMatrix[i-1,j]

    end

end
writedlm("MatrixContact2.dat",X)
########


NB = N_Binomial()
ContactMatrix = ContactMatrixFunc()
ContactMatrix2 = ContactMatrixFunc2()
for i=1:P.grid_size_human
    humans[i].daily_contacts = rand(NB[humans[i].contact_group])
end

 @profile finding_contact(humans,1,ContactMatrix)

## to do
## if there is a protection level, the
function myfunc()
    A = rand(200, 200, 400)
    maximum(A)
end

#############################################################
#### Testing Frailty indexx, vaccine coverage etc
FrIndex = Matrix{Float64}(length(humans),2)
for i = 1:length(humans)
rd = rand()
MaxFra,MinFra = frailty(humans[i])
FrIndex[i,2] = rd*(MaxFra-MinFra)+MinFra
FrIndex[i,1] = humans[i].age
end

writedlm("Frailty.dat",FrIndex)

for i = 1:length(humans)
    FrIndex[i,2] = humans[i].vaccineEfficacy
    FrIndex[i,1] = humans[i].age
end


writedlm("VaccineEffic.dat",FrIndex)

for i = 1:length(humans)
    FrIndex[i,2] = humans[i].Coverage
    FrIndex[i,1] = humans[i].age
end


writedlm("VaccineCov.dat",FrIndex)

find(x -> x.age <= 64 && x.age >= 50, humans)

find(x -> x.age <= 64 && x.age >= 50 && x.vaccinationStatus == 1, humans)

###############################################################33
######################### Testing number of daily contacts

TestMatrix = zeros(Int64,15,15)
NB = N_Binomial()
ContactMatrix = ContactMatrixFunc()
for t = 1:P.sim_time
    for i=1:P.grid_size_human

        humans[i].daily_contacts = rand(NB[humans[i].contact_group])
        for j=1:humans[i].daily_contacts
            r =finding_contact2(humans,i,ContactMatrix,Age_group_Matrix,Number_in_age_group)# rand(1:P.grid_size_human)#
            TestMatrix[humans[i].contact_group,humans[r].contact_group]+=1
        end
    end
end
TestMatrix2 = zeros(Float64,15,15)

for i = 1:15
 TestMatrix2[:,i] = TestMatrix[:,i]/Number_in_age_group[i]
end

TestMatrix2 = TestMatrix2/P.sim_time
writedlm("DailyContactsP10000.dat",TestMatrix2)


## AFFAN'S UNIT TESTS.
using Base.Filesystem

using Revise
include("Influenza.jl")
using .InfluenzaModel
humans = init();
p = InfluenzaParameters(vaccine_efficacy = 0.4);
InfluenzaModel.setup_demographic(humans)
InfluenzaModel.apply_vaccination(humans,p) 


function check_vaccine_coverage_levels(brak, h)
    ## find everyone less than 4 years of age
    lfour = findall(x -> x.age > brak[1] && x.age <= brak[2], h)
    println("lfour: $(length(lfour))")
    lfour_v = findall(x -> x.age > brak[1] && x.age <= brak[2] && x.vaccineEfficacy > 0.0, h)
    println("vax: $(length(lfour_v))")
    cov = length(lfour_v)/length(lfour)
    println("coverage: $cov")
    println("")
end
check_vaccine_coverage_levels((0, 4), humans)
check_vaccine_coverage_levels((4, 49), humans)
check_vaccine_coverage_levels((49, 64), humans)
check_vaccine_coverage_levels((64, 100), humans)

Fail_Contact_Matrix    = zeros(Int64, 15, 15)    ## how many times did susc/sick contact group i meet contact group j meet but failed to infect.
Contact_Matrix_General = zeros(Int64, 15, 15)    ## how many times did contact group i meet with contact group j
Number_in_age_group    = zeros(Int64, 15)                      # vector that tells us number of people in each age group.
Age_group_Matrix       = zeros(Int64, 15, p.grid_size_human)   # a matrix representation of who is inside that age group (ie. row 1 has all the people that have group 1). 

## check if the Age_group_Matrix works 

## this function just fills in the empty matrices as defined above.
InfluenzaModel.setup_contact_matrix(humans, Age_group_Matrix, Number_in_age_group)
ctr = 0
for i = 1:15
    global ctr
    ctr += length(findall(x -> x> 0, Age_group_Matrix[i, :])) 
end
@assert ctr == p.grid_size_human
@assert sum(Number_in_age_group) == p.grid_size_human

for t=1:p.sim_time        
    InfluenzaModel.contact_dynamic2(humans, p, Fail_Contact_Matrix, Age_group_Matrix, Number_in_age_group, Contact_Matrix_General)
    for i=1:p.grid_size_human
       InfluenzaModel.increase_timestate(humans[i], p)
    end          
end
ContactMatrix = InfluenzaModel.ContactMatrixFunc()
r = finding_contact2(humans, 41, ContactMatrix, Age_group_Matrix, Number_in_age_group)


function finding_contact2(h, idx, M, agm, nagm)
    rd = rand()
    g = h[idx].contact_group
    println(g)
    g2 = findfirst(x -> rd <= x, M[:,g])         
    println(g2)
    people_to_search = agm[g2,1:nagm[g2]]
    filter!(x -> x != 18, people_to_search)
    return rand(people_to_search)
end 
function check_group(h)
    ar = zeros(Int64, 4)
    for i = 1:length(h)        
        ar[h[i].group] += 1
    end
    return ar
end
demographic_group = zeros(Int64, 10000)
for i = 1:length(humans)              
    demographic_group[i] = humans[i].group  
end

aux = zeros(Int64,P.grid_size_human)
aux2 = zeros(Float64,P.grid_size_human)
aux3 = zeros(Float64,P.grid_size_human)
matrix_aux = Array{Int8,2}(undef,1,566)
matrix_aux[1,:] = Vaccine_Strain

for i=1:P.grid_size_human
    humans[i].EfficacyVS = Calculating_Efficacy(matrix_aux,1,Vaccine_Strain,humans[i].vaccineEfficacy,P)[1]
    aux[i] = humans[i].age
    aux2[i] = humans[i].vaccineEfficacy
    aux3[i] = humans[i].EfficacyVS
end

A = [aux aux2 aux3]

importants = findall(x -> x>0.0,A[:,2])
A = A[importants,:]

writedlm("testeEff.dat",A)
Random.seed!(1234599)

rng1 = MersenneTwister(1234)

rng2 = MersenneTwister(123024)

rand()
rand(rng1)

rand!(rng2,zeros(1))

rand!(rng1,zeros(1))

rand!(rng2,zeros(1))[1]



Strain = zeros(Int8,2,length(Vaccine_Strain))
Strain[1,:] = Vaccine_Strain
Strain[2,:] = Vaccine_Strain

for i = 1:30
    r = rand(1:length(Vaccine_Strain))
    Strain[2,r] = 0

end

VaccineEfVector = Calculating_Efficacy(Strain,2,Vaccine_Strain,0.6,P)

vacPeople = findall(x->x.vaccinationStatus==1,humans)
Eff = zeros(Float64,length(vacPeople),3)

for i = 1:length(vacPeople)

    Eff[i,1] = humans[vacPeople[i]].vaccineEfficacy
    Eff[i,2] = humans[vacPeople[i]].EfficacyVS
    Eff[i,3] = humans[vacPeople[i]].WhoInf
end
Eff
Dif = Eff[:,1]-Eff[:,2]
findall(x->Eff[x,1]-Eff[x,2]>0 && Eff[x,3]>0,1:length(vacPeople))

for i = 1:10
    r = rand(1:length(Vaccine_Strain))
    Strain[r] = 0

end

Eff = zeros(Float64,length(vacPeople),2)

for i = 1:length(vacPeople)

    Eff[i,1] = humans[vacPeople[i]].vaccineEfficacy
    Eff[i,2] = Calculating_Efficacy(Strain,1,Vaccine_Strain,humans[vacPeople[i]].vaccineEfficacy,P)[1]

end
humans[vacPeople].vaccineEfficacy
NS = humans[WhoGotInf].NumberStrains

auxSum = zeros(Int64,length(WhoGotInf))

for i = 1:length(WhoGotInf)
    auxSum[i] = humans[WhoGotInf[i]].NumberStrains
end
total_number_of_strains = sum(auxSum)

auxSum = zeros(Int64,total_number_of_strains)


count = 0
for i = 1:10
    
    for j = 1:humans[i].NumberStrains
        global count = count+1
        auxSum[count] = Calculating_Distance_Two_Strains(humans[i].strains_matrix[j,:],Vaccine_Strain)
        
        println(count)
    end
end

for i = 1:80

    if i%7 == 0
        println("$i")
    end
end