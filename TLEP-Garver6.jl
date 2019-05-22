#Using JULIA language
#Author: Juan Luis Barberia
#Load data for 6 bus system (garver)

using CSV

branchdata  = CSV.read("branchdata.txt";  types=[Int,Int,Int,Float64,Float64,Float64])
genloaddata = CSV.read("genloaddata.txt"; types=[Int,Float64,Float64,Float64])

println(branchdata);println(genloaddata);

#Solve TLEP with DC-model, with generation rescheduling
#Using Juniper for branch and bound algorithm

using JuMP
using LinearAlgebra
using Juniper
using Ipopt
using Cbc

#Aux Variables
lenbus =length(genloaddata.bus)
lenline=length(branchdata.from)
ug     =Array(genloaddata.gen_max)
d      =Array(genloaddata.load)
un     =5
n0     =Array(branchdata.n0)
B      =Array(branchdata.b)
cost   =Array(branchdata.cost)
f_max  =Array(branchdata.f_max)
indexi=[]
indexj=[]
for i=1:length(branchdata.from); push!(indexi,branchdata.from[i]);push!(indexj,branchdata.to[i]) ;end

#Create a model
optimizer=Juniper.Optimizer
params = Dict{Symbol,Any}()
params[:nl_solver]= with_optimizer(Ipopt.Optimizer, print_level=0)
params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel=0)
params[:branch_strategy] = :Reliability
params[:traverse_strategy] = :DBFS
params[:processors] = :1
params[:strong_restart] = true
params[:log_levels] = []
params[:allow_almost_solved_integral] = false

m=Model(with_optimizer(optimizer,params))

#Def variables

@variable(m, theta[1:lenbus])
@variable(m,0<=n[1:lenline]<=un,Int)
@variable(m, 0<=g[i=1:lenbus]<=ug[i])
@variable(m, 0<=r[i=1:lenbus]<= d[i])
@variable(m,f[1:lenline])

#Def constraint

@constraint(m,f.<=(n+n0).*f_max)                                #Max flow per line
@constraint(m,f.>=(n+n0).*(-f_max))                             #Max flow per line
@constraint(m,f-B.*(n+n0).*(theta[indexi]-theta[indexj]).==0)   #KVL

#Create Y matrix for KCL

Y=Array{GenericAffExpr{Float64,VariableRef},1}()
for i=1:(lenbus^2);push!(Y,0);end
Y=reshape(Y,lenbus,lenbus)
Bn=B.*(n+n0)

for i=1:lenline
  x=branchdata.from[i]
  y=branchdata.to[i]
  Y[x,y]+=Bn[i]
  Y[y,x]+=Bn[i]
  Y[x,x]-=Bn[i]
  Y[y,y]-=Bn[i]
end

@constraint(m,Y*theta+r+g.==d)                                #KCL

#Def objective

alpha=75
@objective(m,Min,sum(n.*cost)+alpha*sum(r))

#Solve

optimize!(m)

#Print results

println("Objective:")
println(JuMP.objective_value(m))

for i=1:lenline;
  ii=indexi[i]
  ij=indexj[i]
  N=convert(Int,JuMP.value.(n[i]))
  N0=n0[i]
  println("Line $ii,$ij.........n0=$N0........nk=$N")
end
