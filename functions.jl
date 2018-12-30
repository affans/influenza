function make_human_symp(h::Human, P::InfluenzaParameters)
    ## make the i'th human infected
  h.health = SYMP    # make the health ->inf
  h.swap = UNDEF
  d = LogNormal(P.log_normal_mean,sqrt(P.log_normal_shape))
  h.statetime = min(15,ceil(rand(d)))
  h.timeinstate = 0
  h.WentTo = SYMP
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



function make_human_recovered(h::Human, P::InfluenzaParameters)
    ## make the i'th human recovered
  h.health = REC    # make the health -> latent
  h.swap = UNDEF
  h.statetime = 999
  h.timeinstate = 0
end

function make_human_latent(h::Human, P::InfluenzaParameters)
    ## make the i'th human infected
  h.health = LAT    # make the health ->inf
  h.swap = UNDEF
  h.statetime = rand(P.Latent_period_Min:P.Latent_period_Max)
  h.timeinstate = 0
end

function setup_rand_initial_latent(h::Array{Human}, P::InfluenzaParameters)
   
    randperson = rand(1:P.grid_size_human)
      
    make_human_latent(h[randperson], P)
   
    return randperson
end

function increase_timestate(h::Human,P::InfluenzaParameters)

    h.timeinstate+=1

    if h.timeinstate >= h.statetime
        if h.health == ASYMP || h.health == SYMP
            h.swap = REC

        elseif h.health == LAT
            prob = (P.ProbAsympMax - P.ProbAsympMin)*rand()+P.ProbAsympMin
            if rand() < (1-prob)*(1-h.vaccineEfficacy)
                h.swap = SYMP
            else h.swap = ASYMP
            end
        end
    end

end

function update_human(h::Array{Human},P::InfluenzaParameters)
    n1::Int64 = 0
    n2::Int64 = 0
    n3::Int64 = 0
    for i=1:P.grid_size_human
        if h[i].swap == LAT
            make_human_latent(h[i],P)
            n1+=1
        else
            if h[i].swap == SYMP
            make_human_symp(h[i],P)
                n2+=1
            else
                if h[i].swap == ASYMP
                    make_human_asymp(h[i],P)
                    n3+=1
                elseif h[i].swap == REC
                    make_human_recovered(h[i],P)
                end

            end
        end
    end
    return n1,n2,n3
end

function vaccination(h::Array{Human},P::InfluenzaParameters)
   
   
    for i=1:length(h)
        if rand()< h[i].Coverage
            h[i].vaccinationStatus = 1
            rd = rand()
            MaxFra,MinFra = FrailtyIndex(h[i])
            FrIndex = rd*(MaxFra-MinFra)+MinFra
            h[i].vaccineEfficacy = P.VaccineEfficacy*(1.0-FrIndex)
        end
    end

end

function FrailtyIndex(h::Human)

    if h.age <= 34
        y = 0.26875-0.00435*h.age
    elseif h.age <= 69
        y = 0.01282*h.age-0.30658
    else
        y = 0.396+0.0039*h.age
    end

    min1 = max((y-0.05),0.0)

    if y+0.05 > 1.0 ##function min was not working, so I did it manually
        max1 = 1.0
    else max1 = y+0.05 
    end

    return max1,min1
end



function N_Binomial()

    ##Age group's mean
    AgeMean = Vector{Float64}(15)
    AgeSD = Vector{Float64}(15)
    
    AgeMean[1] = 10.21
    AgeMean[2] = 14.81
    AgeMean[3] = 18.22
    AgeMean[4] = 17.58
    AgeMean[5] = 13.57
    AgeMean[6] = 13.57
    AgeMean[7] = 14.14
    AgeMean[8] = 14.14
    AgeMean[9] = 13.83
    AgeMean[10] = 13.83
    AgeMean[11] = 12.30
    AgeMean[12] = 12.30
    AgeMean[13] = 9.21
    AgeMean[14] = 9.21
    AgeMean[15] = 6.89
    
    AgeSD[1] = 7.65
    AgeSD[2] = 10.09
    AgeSD[3] = 12.27
    AgeSD[4] = 12.03
    AgeSD[5] = 10.60
    AgeSD[6] = 10.60
    AgeSD[7] = 10.15
    AgeSD[8] = 10.15
    AgeSD[9] = 10.86
    AgeSD[10] = 10.86
    AgeSD[11] = 10.23
    AgeSD[12] = 10.23
    AgeSD[13] = 7.96
    AgeSD[14] = 7.96
    AgeSD[15] = 5.83
    
    p = 1 - (AgeSD[1]^2-AgeMean[1])/(AgeSD[1]^2)
    r = AgeMean[1]^2/(AgeSD[1]^2-AgeMean[1])
    
    nc1 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[2]^2-AgeMean[2])/(AgeSD[2]^2)
    r = AgeMean[2]^2/(AgeSD[2]^2-AgeMean[2])
    nc2 = NegativeBinomial(r, p)
    
    p = 1 - (AgeSD[3]^2-AgeMean[3])/(AgeSD[3]^2)
    r = AgeMean[3]^2/(AgeSD[3]^2-AgeMean[3])
    nc3 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[4]^2-AgeMean[4])/(AgeSD[4]^2)
    r = AgeMean[4]^2/(AgeSD[4]^2-AgeMean[4])
    nc4 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[5]^2-AgeMean[5])/(AgeSD[5]^2)
    r = AgeMean[5]^2/(AgeSD[5]^2-AgeMean[5])
    nc5 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[6]^2-AgeMean[6])/(AgeSD[6]^2)
    r = AgeMean[6]^2/(AgeSD[6]^2-AgeMean[6])
    nc6 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[7]^2-AgeMean[7])/(AgeSD[7]^2)
    r = AgeMean[7]^2/(AgeSD[7]^2-AgeMean[7])
    nc7 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[8]^2-AgeMean[8])/(AgeSD[8]^2)
    r = AgeMean[8]^2/(AgeSD[8]^2-AgeMean[8])
    nc8 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[9]^2-AgeMean[9])/(AgeSD[9]^2)
    r = AgeMean[9]^2/(AgeSD[9]^2-AgeMean[9])
    nc9 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[10]^2-AgeMean[10])/(AgeSD[10]^2)
    r = AgeMean[10]^2/(AgeSD[10]^2-AgeMean[10])
    nc10 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[11]^2-AgeMean[11])/(AgeSD[11]^2)
    r = AgeMean[11]^2/(AgeSD[11]^2-AgeMean[11])
    nc11 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[12]^2-AgeMean[12])/(AgeSD[12]^2)
    r = AgeMean[12]^2/(AgeSD[12]^2-AgeMean[12])
    nc12 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[13]^2-AgeMean[13])/(AgeSD[13]^2)
    r = AgeMean[13]^2/(AgeSD[13]^2-AgeMean[13])
    nc13 = NegativeBinomial(r, p)
    

    p = 1 - (AgeSD[14]^2-AgeMean[14])/(AgeSD[14]^2)
    r = AgeMean[14]^2/(AgeSD[14]^2-AgeMean[14])
    nc14 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[15]^2-AgeMean[15])/(AgeSD[15]^2)
    r = AgeMean[15]^2/(AgeSD[15]^2-AgeMean[15])
    nc15 = NegativeBinomial(r, p)
    
    return nc1,nc2,nc3,nc4,nc5,nc6,nc7,nc8,nc9,nc10,nc11,nc12,nc13,nc14,nc15
    
end



function contact_dynamic2(h::Array{Human},P::InfluenzaParameters,Fail_Contact_Matrix,Age_group_Matrix,Number_in_age_group,Contact_Matrix_General)#,Risk_Contact,t::Int64)
    NB = N_Binomial()
    ContactMatrix = ContactMatrixFunc()
   
    for i=1:P.grid_size_human
        if h[i].health == SUSC
            h[i].daily_contacts = rand(NB[h[i].contact_group])
            for j=1:h[i].daily_contacts
                r =finding_contact2(h,i,ContactMatrix,Age_group_Matrix,Number_in_age_group)# rand(1:P.grid_size_human)#
                Contact_Matrix_General[h[i].contact_group,h[r].contact_group]+=1
                #Contact_Matrix_General[h[r].contact_group,h[i].contact_group]+=1
                if h[r].health == SYMP
                    if h[i].vaccinationStatus == 1
                        if rand()<(1-P.precaution_factorV)
                           
                            if rand() < (P.Prob_transmission*(1-h[i].vaccineEfficacy))
                                h[i].swap = LAT
                                h[i].WhoInf = r
                                break
                                
                            else
                                Fail_Contact_Matrix[h[i].contact_group,h[r].contact_group]+=1
                                h[i].NumberFails+=1
                               
                            end
                        end
                    else 
                        if rand()<(1-P.precaution_factorS)
                           
                            if rand()< P.Prob_transmission
                                h[i].swap = LAT
                                h[i].WhoInf = r
                                break
                              
                            else
                                Fail_Contact_Matrix[h[i].contact_group,h[r].contact_group]+=1
                                h[i].NumberFails+=1
                               
                            end
                        end
                    end

                elseif h[r].health == ASYMP
                    if h[i].vaccinationStatus == 1
                      
                        if rand() < (P.Prob_transmission*(1-h[i].vaccineEfficacy)*(1-P.reduction_factor))
                            h[i].swap = LAT
                            h[i].WhoInf = r
                            break
                           
                        else
                            Fail_Contact_Matrix[h[i].contact_group,h[r].contact_group]+=1
                            h[i].NumberFails+=1
                           
                        end
                    else 
                       
                        if rand()< (P.Prob_transmission*(1-P.reduction_factor))
                            h[i].swap = LAT
                            h[i].WhoInf = r
                            break
                           
                        else
                            Fail_Contact_Matrix[h[i].contact_group,h[r].contact_group]+=1
                            h[i].NumberFails+=1
                           
                        end
                    end
                end

            end

        end
    end ##close Grid human

end


function finding_contact2(h::Array{Human},index::Int64,M,Age_group_Matrix,Number_in_age_group)

    rd = rand()
    g = h[index].contact_group
    g2 = findfirst(x -> rd <= x, M[:,g])
    aux::Int64 = 0

    while aux == 0
        person1 = rand(Age_group_Matrix[g2,1:Number_in_age_group[g2]])
   
        if person1 != index
            aux = 1
            return person1
        end
    end
    

end 