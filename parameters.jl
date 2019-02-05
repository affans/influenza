@enum HEALTH SUSC = 1 LAT = 2 SYMP = 3 ASYMP = 4 REC = 5 UNDEF = 0


@with_kw struct InfluenzaParameters @deftype Int64
    sim_time = 200
    grid_size_human = 10000
    initial_infected = 1

    ## latent periods
    Latent_period_Max = 2
    Latent_period_Min = 1

    ## probability of going to latent -> asymptomatic/symptomatic
    ProbAsympMin::Float64 = 0.3
    ProbAsympMax::Float64 = 0.7
    reduction_factor::Float64 = 0 ## for now we set this to zero: this means prob of infection is the same whether symp or asymp

    ## if a person is vaccinated, what is their vaccine efficacy 
    ## this number eventually gets multiplied by the (calculated) frailty index.
    vaccine_efficacy::Float64 = 0.0 # change to different values 0.8, 0.6, 0.4

    transmission_beta::Float64 = 0.079 # we need to calibrate this for 1%, 4% and 8% clinical attack rate (symptomatic cases)
    log_normal_mean::Float64 = 1.0 # distribution for infectious period ## scale value
    log_normal_shape::Float64 = 0.4356
    max_infectious_period::Int64 = 15
    ####Mutation parameters

    mutation_rate::Float64 = 0.005
    matrix_strain_lines::Int64 = 100 #number of strains it is possible to save
    ##increasing the mutation rate, it is necessary to increase this number
    sequence_size::Int64 = 566 #number of strain's sites.
    number_of_states::Int8 = 20 ##number of different states of the sites
end

## age distribution discrete for humans, Canada
function distribution_age()
    ProbBirthAge = Vector{Float64}(undef, 17)
    SumProbBirthAge = Vector{Float64}(undef, 17)
    AgeMin = Vector{Int64}(undef, 17)
    AgeMax = Vector{Int64}(undef, 17)
    ProbBirthAge[1] = 0.053
    ProbBirthAge[2] = 0.055
    ProbBirthAge[3] = 0.052
    ProbBirthAge[4] = 0.056
    ProbBirthAge[5] = 0.067
    ProbBirthAge[6] = 0.070
    ProbBirthAge[7] = 0.070
    ProbBirthAge[8] = 0.068
    ProbBirthAge[9] = 0.064
    ProbBirthAge[10] = 0.066
    ProbBirthAge[11] = 0.072
    ProbBirthAge[12] = 0.073
    ProbBirthAge[13] = 0.064
    ProbBirthAge[14] = 0.054
    ProbBirthAge[15] = 0.042
    ProbBirthAge[16] = 0.029
    ProbBirthAge[17] = 0.042
    SumProbBirthAge = cumsum(ProbBirthAge)
    SumProbBirthAge[end] = 1.0
    AgeMin[1] = 0;
    AgeMax[1] = 4;

    AgeMin[2] = 5;
    AgeMax[2] = 9;

    AgeMin[3] = 10;
    AgeMax[3] = 14;

    AgeMin[4] = 15;
    AgeMax[4] = 19;

    AgeMin[5] = 20;
    AgeMax[5] = 24;

    AgeMin[6] = 25;
    AgeMax[6] = 29;

    AgeMin[7] = 30;
    AgeMax[7] = 34;

    AgeMin[8] = 35;
    AgeMax[8] = 39;

    AgeMin[9] = 40;
    AgeMax[9] = 44;

    AgeMin[10] = 45;
    AgeMax[10] = 49;

    AgeMin[11] = 50;
    AgeMax[11] = 54;

    AgeMin[12] = 55;
    AgeMax[12] = 59;

    AgeMin[13] = 60;
    AgeMax[13] = 64;

    AgeMin[14] = 65;
    AgeMax[14] = 69;

    AgeMin[15] = 70;
    AgeMax[15] = 74;

    AgeMin[16] = 75;
    AgeMax[16] = 79;

    AgeMin[17] = 80;
    AgeMax[17] = 100;


    return SumProbBirthAge,AgeMin,AgeMax
end

function ContactMatrixFunc()
    CM = zeros(Float64, 15, 15)
    CM[1, :] = [0.2287227639, 0.0578904992, 0.0147540984, 0.0100089736, 0.0190779014, 0.0375827663, 0.0554200656, 0.0437937938, 0.0462266543, 0.0236936358, 0.0184025802, 0.0231225914, 0.0311685233, 0.0212494688, 0.0201779473]
    CM[2, :] = [0.3440471891, 0.4866344605, 0.0763661202, 0.0301649755, 0.0345345345, 0.0714352327, 0.1177361029, 0.1086086086, 0.093551931, 0.0491321517, 0.0539745779, 0.0524945318, 0.0557004029, 0.0569485763, 0.0478233238]
    CM[3, :] = [0.3872637535, 0.5714170692, 0.5845628415, 0.1059570655, 0.0548489666, 0.0913923342, 0.1621394332, 0.1587420754, 0.1551593003, 0.0899990816, 0.1026370708, 0.0834288095, 0.0850912539, 0.0820229494, 0.0986653956]
    CM[4, :] = [0.4126640183, 0.5926731079, 0.6730874317, 0.6061986609, 0.1655184596, 0.1379278187, 0.1983853334, 0.1941107774, 0.2220907631, 0.1799981633, 0.1611648644, 0.1249869805, 0.1146006163, 0.1053973651, 0.1525262154]
    CM[5, :] = [0.4493800409, 0.6101449275, 0.6949453552, 0.6866155864, 0.4279279279, 0.2639186795, 0.2519552603, 0.241324658,  0.291219471,  0.2458444302, 0.2378106621, 0.1792521612, 0.1570277317, 0.1399631676, 0.1839847474]
    CM[6, :] = [0.5237751294, 0.652173913,  0.7161202186, 0.7236142749, 0.5620915033, 0.4486617551, 0.3590110167, 0.3191524858, 0.3590805375, 0.322159978, 0.3357996585, 0.2558066868, 0.2245792842, 0.2059781839, 0.2353034636]
    CM[7, :] = [0.6348862405, 0.7144122383, 0.7481557377, 0.7523296749, 0.65527292,   0.5724144363, 0.5182070473, 0.4225058392, 0.4519563931, 0.4088529709, 0.4162398027, 0.3453806895, 0.3200995497, 0.2838929027, 0.2902764538]
    CM[8, :] = [0.7340796918, 0.791626409,  0.8023907104, 0.7966452682, 0.722045575,  0.6619416208, 0.6553696073, 0.5838338338, 0.5565790586, 0.5072091101, 0.4973439575, 0.4229767733, 0.4124200047, 0.3817821221, 0.3541468065]
    CM[9, :] = [0.7954736969, 0.8639291465, 0.8740437158, 0.8515910817, 0.7824589295, 0.736174578,  0.7554452948, 0.7181348015, 0.7006676244, 0.6194324548, 0.5895465756, 0.5163003854, 0.5003555345, 0.4792463522, 0.4401016841]
    CM[10, :]= [0.8345973276, 0.8972624799, 0.9204234973, 0.9178573894, 0.8578872991, 0.8051851161, 0.8231435539, 0.8019686353, 0.805121271,  0.756359629, 0.6913299184, 0.5997291949, 0.5699217824, 0.5429947585, 0.5290753098]
    CM[11, :]= [0.8760081859, 0.9234299517, 0.9460382514, 0.9519569269, 0.9194488606, 0.8811899655, 0.8819275082, 0.8677844511, 0.8823628835, 0.8526035449, 0.8022196927, 0.7022185189, 0.6420952832, 0.6102847429, 0.5967588179]
    CM[12, :]= [0.9141687733, 0.9438808374, 0.9607240437, 0.9704562711, 0.9528351881, 0.9308962044, 0.9244806997, 0.9051551552, 0.9129552945, 0.9038479199, 0.8839878581, 0.8314758879, 0.743422612, 0.6843745573, 0.6561804894]    
    CM[13, :]= [0.9461899603, 0.9631239936, 0.9705601093, 0.979153724,  0.9685567921, 0.959153222,  0.9572786141, 0.9432766099, 0.9416884983, 0.9339700615, 0.9279074179, 0.9002187272, 0.8497274236, 0.7788638617, 0.7316491897]
    CM[14, :]= [0.9695437583, 0.9768921095, 0.9801912568, 0.9861945192, 0.9777424483, 0.9732351021, 0.9746867379, 0.9708041375, 0.9620552692, 0.9496739829, 0.9462151394, 0.9396937819, 0.9156198151, 0.8810029749, 0.8156974897]
    CM[15, :]= [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    return CM
end


function N_Binomial()

    ##Age group's mean
    AgeMean = Vector{Float64}(undef, 15)
    AgeSD = Vector{Float64}(undef, 15)

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


    nbinoms = Vector{NegativeBinomial{Float64}}(undef, 15)
    for i = 1:15
        p = 1 - (AgeSD[i]^2-AgeMean[i])/(AgeSD[i]^2)
        r = AgeMean[i]^2/(AgeSD[i]^2-AgeMean[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end

    
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
