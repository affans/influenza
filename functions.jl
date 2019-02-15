function contact_dynamic2(h, P::InfluenzaParameters, NB, ContactMatrix, fcm, agm, nagm, cmg,Vaccine_Strain::Array{Int8,1},rng1)
    for i=1:length(h)
        if h[i].health == SUSC
            # don't need to save this information... we might
            # h[i].daily_contacts = rand(NB[h[i].contact_group])
            _contacts = rand(NB[h[i].contact_group])
            #_contacts = rand(NB[1])
            for j=1:_contacts
                r = finding_contact2(h, i, ContactMatrix, agm, nagm)
                cmg[h[i].contact_group, h[r].contact_group] += 1
                
                ################################3
                if h[r].health == SYMP
                    available_strains = findall(x -> h[r].Vector_time[x] < (h[r].timeinstate+h[r].latenttime) && h[r].Vector_time[x] < (h[r].statetime/2.0+h[r].latenttime),1:h[r].NumberStrains)
                    VaccineEfVector = zeros(Float64,length(available_strains))

                    if h[i].vaccinationStatus == 1
                        VaccineEfVector = Calculating_Efficacy(h[r].strains_matrix[available_strains,:],length(available_strains),Vaccine_Strain,h[i].vaccineEfficacy,P)

                        if rand() < ProbOfTransmission(P.transmission_beta,VaccineEfVector)
                            TransmitingStrain = Which_One_Will_Transmit(VaccineEfVector,h[r].Vector_time[available_strains],h[r].timeinstate,h[r].latenttime,rng1)
                            #h[i].NumberStrains = h[i].NumberStrains + 1
                            #h[i].strains_matrix[h[i].NumberStrains,:] =  h[r].strains_matrix[available_strains[TransmitingStrain],:]
                            h[i].strains_matrix[1,:] =  h[r].strains_matrix[available_strains[TransmitingStrain],:]
                            
                            h[i].EfficacyVS = VaccineEfVector[TransmitingStrain]
                            h[i].swap = LAT
                            h[i].WhoInf = r
                           # break                               
                        end
                        
                    else 
                            
                        if rand()< P.transmission_beta
                            TransmitingStrain = Which_One_Will_Transmit(VaccineEfVector,h[r].Vector_time[available_strains],h[r].timeinstate,h[r].latenttime,rng1)
                            #h[i].NumberStrains = h[i].NumberStrains + 1
                            #h[i].strains_matrix[h[i].NumberStrains,:] =  h[r].strains_matrix[available_strains[TransmitingStrain],:]
                            h[i].strains_matrix[1,:] =  h[r].strains_matrix[available_strains[TransmitingStrain],:]
                            h[i].EfficacyVS = 0.0
                            h[i].swap = LAT
                            h[i].WhoInf = r
                          #  break
                        end
                    
                    end

                elseif h[r].health == ASYMP
                    available_strains = findall(x -> h[r].Vector_time[x] < (h[r].timeinstate+h[r].latenttime) && h[r].Vector_time[x] < (h[r].statetime/2.0+h[r].latenttime),1:h[r].NumberStrains)
                    VaccineEfVector = zeros(Float64,length(available_strains))

                    if h[i].vaccinationStatus == 1
                        VaccineEfVector = Calculating_Efficacy(h[r].strains_matrix[available_strains,:],length(available_strains),Vaccine_Strain,h[i].vaccineEfficacy,P)

                        if rand() < ProbOfTransmission((P.transmission_beta*(1-P.reduction_factor)),VaccineEfVector)
                            TransmitingStrain = Which_One_Will_Transmit(VaccineEfVector,h[r].Vector_time[available_strains],h[r].timeinstate,h[r].latenttime,rng1)
                            #h[i].NumberStrains = h[i].NumberStrains + 1
                            #h[i].strains_matrix[h[i].NumberStrains,:] =  h[r].strains_matrix[available_strains[TransmitingStrain],:]
                            h[i].strains_matrix[1,:] =  h[r].strains_matrix[available_strains[TransmitingStrain],:]
                            h[i].EfficacyVS = VaccineEfVector[TransmitingStrain]
                            h[i].swap = LAT
                            h[i].WhoInf = r
                           # break
                        end
                    else 
                        if rand()< (P.transmission_beta*(1-P.reduction_factor))
                            TransmitingStrain = Which_One_Will_Transmit(VaccineEfVector,h[r].Vector_time[available_strains],h[r].timeinstate,h[r].latenttime,rng1)
                            #h[i].NumberStrains = h[i].NumberStrains + 1
                            #h[i].strains_matrix[h[i].NumberStrains,:] =  h[r].strains_matrix[available_strains[TransmitingStrain],:]
                            h[i].strains_matrix[1,:] =  h[r].strains_matrix[available_strains[TransmitingStrain],:]
                            h[i].EfficacyVS = 0.0
                            h[i].swap = LAT
                            h[i].WhoInf = r
                           # break
                        end
                    end
                end

            end
        end

    end ##close Grid human
end

@inline function finding_contact2(h, idx, M, agm, nagm)
    rd = rand()
    g = h[idx].contact_group
    g2 = findfirst(x -> rd <= x, M[:,g])         
    people_to_search = agm[g2,1:nagm[g2]]
    filter!(x -> x != 18, people_to_search)
    return rand(people_to_search)
end 


