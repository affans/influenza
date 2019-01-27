## this function is a simple test to learning how pmap and seeds work in julia
## STDLIB
using Random
using DelimitedFiles
using Distributed
using Base.Filesystem
## ADDED PACKAGES
using ClusterManagers
using DataFrames
using CSV
using JSON
using Statistics
using VegaLite
import REPL
using REPL.TerminalMenus

## call all the using in main package to trigger precompilation 
## the precompilation files will get shared amongst all nodes so won't clash with `@everywhere` triggering precompilation
using Parameters      ## with julia 1.1 this is now built in.
using ProgressMeter   ## can now handle parallel with progress_pmap
#using PmapProgressMeter
#using DataFrames
using Distributions
using StatsBase
using StaticArrays

addprocs(2)
@everywhere using Random
@everywhere include("Influenza.jl")
@everywhere using .InfluenzaModel

## minor tests.
@everywhere function tmp_main(simnumber)
    Random.seed!(simnumber)
    rn = rand()
    println("simnumber: $simnumber produces $rn")
    tmp_trans(simnumber)
    sleep(3)
    return nothing
end

@everywhere function tmp_trans(simnumber)
    rn = rand()
    println("... trans - simnumber: $simnumber produces $rn")
end

function tmp_sim()
    results = pmap(x -> main(x), 1:5)
end

function inf_test()
    @everywhere P = InfluenzaParameters(sim_time = 200, vaccine_efficacy = 0.0, transmission_beta=0.03)        
    results = pmap(x -> main(x, P), 1:5)
    return results
end

myfive1 = inf_test()
myfive2= inf_test()
myfive3=inf_test()

# sim 1 , latent
myfive1[1][1] == myfive2[1][1] ==  myfive3[1][1]

_t = myfive1[1][15] ## check demographics in sim 1. 
_ts1 =[sum(_t .== i) for i in unique(_t)]

_t = myfive2[1][15] ## check demographics in sim 1. 
_ts2 = [sum(_t .== i) for i in unique(_t)]
_ts1==_ts2
## module tests
function f(p)
    bar = PM.barstring(p.barlen, percentage_complete, barglyphs=p.barglyphs)
    dur = "dur"
    msg = @sprintf "%s%3u%%%s  ETA: %s" p.desc round(Int, percentage_complete) bar dur
    print(p.output, "\n" ^ (p.offset + p.numprintedvalues))
    PM.move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues)
    PM.printover(p.output, msg, p.color)
end

