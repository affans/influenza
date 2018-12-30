include("basicModelCal.jl")
################## To run this files, You must check the return of basicModelCal.jl
#######################


function run_main(P::InfluenzaParameters,numberofsims::Int64)

    resultsL = Matrix{Int64}(P.sim_time,numberofsims)
    resultsA = Matrix{Int64}(P.sim_time,numberofsims)
    resultsS = Matrix{Int64}(P.sim_time,numberofsims)
    forR0 = open(string("Data1/R0check","$(P.Prob_transmission)",".dat"),"w")
    for i=1:numberofsims
        print("$i ")
        latent,symp,asymp,b = main(i,P)
        resultsL[:,i] = latent
        resultsS[:,i] = symp
        resultsA[:,i] = asymp
        println(forR0,"$i $b")
    end
    close(forR0)
    dataprocess(resultsL,resultsS,resultsA,P)
end

function dataprocess(resultsL,resultsS,resultsA,P::InfluenzaParameters)

        writedlm("Data1/result_latent.dat",resultsL)
        writedlm("Data1/result_symp.dat",resultsS)
        writedlm("Data1/result_asymp.dat",resultsA)

end
sigma1 = 0.0
ef1 = 0.0
P=InfluenzaParameters(
    precaution_factorS = sigma1,
    precaution_factorV = 0.0,
    VaccineEfficacy = ef1,
    GeneralCoverage = 0.0,
    Prob_transmission = 0.08,
    sim_time = 365,
)

run_main(P,1000)
