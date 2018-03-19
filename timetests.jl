function timetest()
  ## distinct contacts over two weeks
      PP_DD = Normal(4.6, 0.78)
      MP = ModelParameters()
      psize = MP.gridsize
      @time dpc = [Int(round(rand(PP_DD))) for i=1:psize];
      @time dpc = Int.(round.(rand(PP_DD, psize)));
      @time begin
      dpc = zeros(Int64, psize)
      for i = 1:length(dpc)
          dpc[i]=Int(round(rand(PP_DD)))
      end
  end
  end

  ## questions what happens when you start defining things outside the function.
  ## what happens when you start hardcoding parameters

function mattest()

## generate patient contact matrix, memory structure is column wise
    @time begin
    pmm = Matrix{Vector{Int64}}(10, 10)
    for i = eachindex(pmm)
        pmm[i] = zeros(Int64, 14)
    end
    end
@time reshape([zeros(Int64, 14) for i=1:100], 10, 10);
end

function mtesttwo()
    @time begin pmm = Matrix{BitArray}(10,10)
    for i = eachindex(pmm)
        pmm[i] = BitArray(14)
    end
    end
    
    @time reshape([BitArray(14) for i=1:100], 10, 10)
end

## test - memory structure of matrix of vectors
pmm = Matrix{Int64}(10, 10)
for i = eachindex(pmm)
    pmm[i] = i
end ## so column

function deltime()
    for patientid =1:100
        dpc = Int(round(rand(PP_DD)))
        elig = sample(1:100, dpc, replace=false)
        deleteat!(elig, elig .== patientid )
    end
end