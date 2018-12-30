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


sigma1 = 0.4

ef1 = 0.4
P=InfluenzaParameters(
    precaution_factorS = sigma1,
    precaution_factorV = 0.0,
    VaccineEfficacy = ef1,
    GeneralCoverage = 1,
    Prob_transmission = 0.079,
    sim_time = 365,
    grid_size_human = 10000
)


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
MaxFra,MinFra = FrailtyIndex(humans[i])
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



