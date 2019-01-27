using DataFrames
using Parameters #module
using Distributions
using StatsBase
using Profile
using BenchmarkTools


include("parameters.jl")
include("population.jl")
include("functions.jl")
P = InfluenzaParameters()
function init()
    [Human(i) for i = 1:10000]
end
 
# change default for `seconds` to 2.5
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2.50
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10

ht = @btime init(); ## good.  209.300 μs (10003 allocations: 1.45 MiB)
@btime setup_demographic($ht) ## 6 allocations? good speed. 103.300 μs (6 allocations: 1.13 KiB)

ht = init()
setup_demographic(ht)
@btime _apply_vax((0, 4), (20, 32), $ht, $P)
@btime apply_vaccination($ht, $P)
## julia> @btime _apply_vax((0, 4), (20, 32), $ht, $P)
#   44.700 μs (15 allocations: 12.78 KiB)

#  julia> @btime apply_vaccination($ht, $P)
#   563.101 μs (71 allocations: 304.34 KiB)
#this is wierd - one run 15 allocations, 12 kib  4 runs - 71 allocations 300kib

ht = init()
setup_demographic(ht)
@btime initial = setup_rand_initial_latent(ht,P) #   78.616 ns (0 allocations: 14 bytes) good


ht = init()
setup_demographic(ht)
apply_vaccination(ht, P)
function time_setup_contact(h)
    Fail_Contact_Matrix    = zeros(Int64, 15, 15)
    Contact_Matrix_General = zeros(Int64, 15, 15)
    Number_in_age_group    = zeros(Int64, 15)
    Age_group_Matrix       = zeros(Int64, 15, 10000)
    setup_contact_matrix(h, Age_group_Matrix, Number_in_age_group)
end
@btime time_setup_contact($ht) 
#614.099 μs (5 allocations: 1.15 MiB) ## the 15x10000 is 1.14 mib in size see using summarysize or sizeof

ht = init()
setup_demographic(ht)
apply_vaccination(ht, P)
function inc(_ht, _P)
    for i=1:length(_ht)
        increase_timestate(_ht[i],_P)
    end
end 
@btime inc($ht, $P)

ht = init()
setup_demographic(ht)
apply_vaccination(ht, P)
function upd(_ht, _P)
    for i=1:length(_ht)
        update_human(_ht, _P)
    end
end 
@btime upd($ht, $P)


@btime _NB = N_Binomial()
@btime _ContactMatrix = ContactMatrixFunc()


_NB = N_Binomial()


ht = init()
setup_demographic(ht)
apply_vaccination(ht, P)
Fail_Contact_Matrix    = zeros(Int64, 15, 15)
Contact_Matrix_General = zeros(Int64, 15, 15)
Number_in_age_group    = zeros(Int64, 15)
Age_group_Matrix       = zeros(Int64, 15, 10000)
setup_contact_matrix(ht, Age_group_Matrix, Number_in_age_group)
_ContactMatrix = ContactMatrixFunc()
@btime finding_contact2($ht, $1, $_ContactMatrix, Age_group_Matrix, Number_in_age_group)

function _dyn()
    _P = InfluenzaParameters()
    ht = init()
    setup_demographic(ht)
    apply_vaccination(ht, _P)
    fcm  = zeros(Int64, 15, 15)
    cmg  = zeros(Int64, 15, 15)
    nag  = zeros(Int64, 15)
    agm  = zeros(Int64, 15, 10000)
    setup_contact_matrix(ht, agm, nag)
    NB = N_Binomial()
    CM = ContactMatrixFunc()
    @btime contact_dynamic2($ht, $_P, $NB, $CM, $fcm, $agm, $nag, $cmg)
end
## trigger compilation 
_dyn()
#@btime _dyn()
Profile.clear()
@profile _dyn()
ProfileView.view()
    

## seed tests

