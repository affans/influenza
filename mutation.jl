function mutation(Original_Strain1::Array{Int8,1},P::InfluenzaParameters,n::Int64,statetime::Int64,latenttime::Int64,rng1)
    #this function receives the strain that infected the individual, the parameters, the number of strains (usually 1),
    #and the infected period
    #it will return a matrix with the strains that appeared before 1/2 infected period
    t = statetime+latenttime
    j = n-1 ###strain index
    Matrix_Of_Strains = zeros(Int8,P.matrix_strain_lines,P.sequence_size) ###Returning matrix with strains
    Matrix_Of_Strains[1,:] = Original_Strain1 #setting the first one as the original infection
    Time_Strain = zeros(Int64,P.matrix_strain_lines) ### the time in which the strain was generated
    Time_Strain[1] = 0 ### the original strain is generated at time 0
    

    while j != n ##Do this until there is no more strains to search
        j+=1 ##going to the next strain
        Time_Vector_Size = zeros(Int64,(P.max_infectious_period+P.Latent_period_Max)) ##This vector saves the number of sites changing at that time step
        Time_Matrix = zeros(Int64,100,(P.max_infectious_period+P.Latent_period_Max)) ##This matrix saves which sites changed at a given time step
        if Time_Strain[j] < (latenttime + statetime/2) ###strains after the half of infectious period do not matter
            probability = (1-exp(-1.0*P.mutation_rate*(t-Time_Strain[j])/365.0))
            for i = 1:P.sequence_size ###running the sites
                if rand(rng1) < probability#CumProb[t-Time_Strain[j]] ##checking if it will change
                    time = rand(rng1,(Time_Strain[j]+1):t) ###selecting a random time step to change it
                    Time_Vector_Size[time]+=1 ###increasing the number of sites that changed at time step "time"
                    Time_Matrix[Time_Vector_Size[time],time] = i ##saving the site
                end
            end
        end

        for time = (Time_Strain[j]+1):(Int64(floor(latenttime + statetime/2.0))-1) #running the time (each strain is able to create t-Time_Strain[j] strains)
            ###strains after the half of infectious period do not matter
            if Time_Vector_Size[time] > 0 ##if it has sites that changed
                n+=1
                Matrix_Of_Strains[n,:] = Matrix_Of_Strains[j,:] ###Copying the generating strain
                for time_count = 1:Time_Vector_Size[time] #running the sites that changed

                    ##this loop guarantees that the new state is different from the old one
                    change::Int8 = Matrix_Of_Strains[n,Time_Matrix[time_count,time]] 
                    while change == Matrix_Of_Strains[n,Time_Matrix[time_count,time]]
                        change = rand(rng1,1:P.number_of_states)
                    end
                    ###########################
                    Matrix_Of_Strains[n,Time_Matrix[time_count,time]] = change ##changing the site's state
                end
                Time_Strain[n] = time #saving the time the strain appeared
            end
        end ###Close for time
    
    end

    return Matrix_Of_Strains,Time_Strain,n ##return the strain matrix, when the strains were generated and how many they are

end


function Calculating_Distance(humans::Array{Human},P::InfluenzaParameters)
    #calculating the distance between all strains that appeared in the population (before 1/2 infected period)
    
    n::Int64 = 0
    for i = 1:P.grid_size_human
        n+=humans[i].NumberStrains
    end

    A = zeros(Int64,P.sequence_size,n)
    count::Int64 = 0
    for i = 1:P.grid_size_human
        for j = 1:humans[i].NumberStrains
            count+=1
            A[:,count] = humans[i].strains_matrix[j,:]
        end
    end

    DistMatrix = zeros(Int64,n,n)
    for i = 1:(n-1)
        for j = i+1:n
           DistMatrix[i,j] = Calculating_Distance_Two_Strains(A[:,i],A[:,j])
        end
    end

    return DistMatrix
end

function Calculating_Distance_Two_Strains(A::Array{Int8,1},B::Array{Int8,1})
    #calculating the antigenic distance between strains A and B
    #info: in this model, each strain corresponds to 566 a.a. with 20 possible states each one.
    soma::Int64 = 0
    for  k = 1:length(A)
        if A[k] != B[k]
            soma+=1
        end
    end
    return soma
end

function ProbOfTransmission(ProbTrans::Float64,VaccineEfVector::Array{Float64,1})
    #the prob of transmission is the Basic prob times the average probability of vaccine fails
    Mean1::Float64 = sum(VaccineEfVector)/length(VaccineEfVector)
    prob::Float64 = ProbTrans*(1-Mean1)
    return prob
end

function Calculating_Efficacy(strains_matrix::Array{Int8,2},NumberStrains::Int64,Vaccine_Strain::Array{Int8,1},vaccineEfficacy::Float64,P::InfluenzaParameters)
    ##Vaccine effectiveness is calculated based on the antigenic distance.
    
    VaccineEfVector = zeros(Float64,NumberStrains)
    p::Float64 = 0.0
    for i = 1:NumberStrains
        p = Calculating_Distance_Two_Strains(strains_matrix[i,:],Vaccine_Strain)/P.sequence_size
        VaccineEfVector[i] = vaccineEfficacy-7.37*p
        if VaccineEfVector[i] < 0
            VaccineEfVector[i] = 0.0
        end
    end
    return VaccineEfVector
end

function Which_One_Will_Transmit(VaccineEfVector::Array{Float64,1},Vector_time::Array{Int64,1},timeinstate::Int64,latenttime::Int64,rng1)
## The strain is chosen based on a Fitness function that depends on the time it is in the body and the V.Ef.
    probs = zeros(Float64,length(VaccineEfVector))#ones(Float64,length(VaccineEfVector))# 
    for i = 1:length(VaccineEfVector)
        probs[i] = (1-VaccineEfVector[i])*(timeinstate+latenttime-Vector_time[i])
    end
  
    probs = probs/sum(probs)
    probs = cumsum(probs)
    r = rand(rng1)
    retorno = findfirst(x->x>r,probs)
    return retorno
end

function Creating_Vaccine_Vector(Vaccine_Strain::Array{Int8,1},P::InfluenzaParameters,rng1)
    ##This function just creates a random strain to begin the simulations.
    for i = 1:P.sequence_size
        Vaccine_Strain[i] = Int8(rand(rng1,1:P.number_of_states))
    end
end