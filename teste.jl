a = 1
b = 2

if a>=b
    println("igual")
else
    println("diferente")
end



a = [1 5 7 3 9]
rand(a)

d = LogNormal(1,sqrt(0.4356))
d1 = Vector{Float64}(3000)
for i = 1:length(d1)
d1[i] = rand(d)
end

writedlm("testeLogNormal.dat",d1)

arquivo = open("NB.dat","w")
nc = NegativeBinomial(r, p)

for i=1:2000
nc1 = rand(nc)
println(arquivo,"$nc1")
end
close(arquivo)

for i = 1:1
print("adsads")
end