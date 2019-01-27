function contact_dynamic2(h, P::InfluenzaParameters, NB, ContactMatrix, fcm, agm, nagm, cmg)
    for i=1:length(h)
        if h[i].health == SUSC
            # don't need to save this information... we might
            # h[i].daily_contacts = rand(NB[h[i].contact_group])
            _contacts = rand(NB[h[i].contact_group])
            #_contacts = rand(NB[1])
            for j=1:_contacts
                r = finding_contact2(h, i, ContactMatrix, agm, nagm)
                cmg[h[i].contact_group, h[r].contact_group] += 1
                
                if h[r].health == SYMP
                    if rand() < (P.transmission_beta * (1 - h[i].vaccineEfficacy))
                        h[i].swap = LAT
                        h[i].WhoInf = r
                    else 
                        h[i].NumberFails += 1     
                        fcm[h[i].contact_group, h[r].contact_group] += 1
                    end
                elseif h[r].health == ASYMP                                         
                    if rand() < (P.transmission_beta * (1 - h[i].vaccineEfficacy) * (1 - P.reduction_factor))
                        h[i].swap = LAT
                        h[i].WhoInf = r                                                
                    else
                        h[i].NumberFails+=1
                        fcm[h[i].contact_group,h[r].contact_group]+=1
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


