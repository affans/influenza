@enum HEALTH SUSC = 1 LAT = 2 SYMP = 3 ASYMP = 4 REC = 5 UNDEF = 0


@with_kw struct InfluenzaParameters @deftype Int64
    sim_time = 150
    grid_size_human = 1000
    initial_infected = 1

    ProbAsympMin::Float64 = 0.3
    ProbAsympMax::Float64 = 0.7
    reduction_factor::Float64 = 0.5 # we can use 1 and calibration will adjust beta

    Latent_period_Max = 2
    Latent_period_Min = 1

    Prob_transmission::Float64 = 0.079 # we need to calibrate this for 1%, 4% and 8% clinical attack rate (symptomatic cases)
    log_normal_mean::Float64 = 1.0 # distribution for infectious period ## scale value
    log_normal_shape::Float64 = 0.4356

    precaution_factorS::Float64 = 1 # change this to 1 to ignore this parameter
    precaution_factorV::Float64 = 1 # change this to 1 to ignore this parameter

    #NumberOfContactsMin = 3
    #NumberOfContactsMax = 12

    GeneralCoverage = 1
    VaccineEfficacy::Float64 = 0.2 # change to different values 0.8, 0.6, 0.4
end

function allocation_test()
    A = Vector{Float64}(undef, 15)
    return A
end

## age distribution discrete for humans, Canada
function distribution_age()
    ProbBirthAge = Vector{Float64}(undef, 17)
    SumProbBirthAge = Vector{Float64}(undef, 17)
    AgeMin = Vector{Int64}(undef, 17)
    AgeMax = Vector{Int64}(undef, 17)


    ProbBirthAge[1] = 0.0532142275
    ProbBirthAge[2] = 0.0545793774
    ProbBirthAge[3] = 0.0523384586
    ProbBirthAge[4] = 0.0560316901
    ProbBirthAge[5] = 0.0674724602
    ProbBirthAge[6] = 0.0701439068
    ProbBirthAge[7] = 0.0695785615
    ProbBirthAge[8] = 0.0682851526
    ProbBirthAge[9] = 0.064437731
    ProbBirthAge[10] = 0.0655332187
    ProbBirthAge[11] = 0.0719434263
    ProbBirthAge[12] = 0.0731115814
    ProbBirthAge[13] = 0.064701399
    ProbBirthAge[14] = 0.0544144521
    ProbBirthAge[15] = 0.0421691092
    ProbBirthAge[16] = 0.0293566227
    ProbBirthAge[17] = 0.0426886252
    SumProbBirthAge = cumsum(ProbBirthAge)

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

# function ReductionPerGroup()
#     Reduction = Vector{Float64}(undef, 15)

#     Reduction[1] = 0.5883212623051612
#     Reduction[2] =  0.6057753194036136
#     Reduction[3] =  0.6223362260154172
#     Reduction[4] =  0.6433283687160004
#     Reduction[5] =  0.6619723535662496
#     Reduction[6] =  0.6789357497894256
#     Reduction[7] =  0.6966433712730415
#     Reduction[8] =  0.6646629167123203
#     Reduction[9] =  0.6142525576610639
#     Reduction[10] =  0.5625963990585514
#     Reduction[11] =  0.5116166856783396
#     Reduction[12] =  0.46292829585154177
#     Reduction[13] =  0.41219197569675525
#     Reduction[14] =  0.35636586929717184
#     Reduction[15] =  0.2335717105984882
#     Reduction = Reduction/0.8 ##  ??? 
#     return Reduction
# end
