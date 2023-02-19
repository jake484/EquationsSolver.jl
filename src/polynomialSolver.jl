struct PolyProblem <: AbstractEquationsProblem
    coff::Vector{Float64}
    PolyProblem(coff::Vector{Float64}) = new(coff)
end

function PolyProblem(coff::Vector{T}) where T<:Real
    PolyProblem(convert(Vector{Float64},coff))
end

struct Bairstov <: AbstractPolyMethod end

function solve(problem::PolyProblem,::Bairstov;p::Float64=0.0,q::Float64=0.0,maxiter::Int64=1000,abstol::Float64=1e-30,retol=1e-5)
    realRoots::Vector{Float64}=[]
    complexRoots::Vector{ComplexF64}=[]
    res::Matrix{Float64}=PolyBairstov(problem.coff,p,q,maxiter,abstol)
    k=size(res,2)
    for i=1:k-1
        temp=res[1,i]^2-4*res[2,i]
        delta=sqrt(abs(temp))
        if delta<retol*abs(res[1,i])
            temp=0
        elseif abs(res[1,i])<retol*delta
            res[1,i]=0
        end
        if temp >= 0
            push!(realRoots,(-res[1,i]+delta)/2,(-res[1,i]-delta)/2)
        else
            push!(complexRoots,(complex(-res[1,i],delta))/2,(complex(-res[1,i],-delta))/2)
        end
    end
    if isodd(length(problem.coff))
        temp=res[1,k]^2-4*res[2,k]
        delta=sqrt(abs(temp))
        if delta<retol*abs(res[1,k])
            temp=0
        elseif abs(res[1,k])<retol*delta
            res[1,k]=0
        end
        if temp >= 0
            push!(realRoots,(-res[1,k]+delta)/2,(-res[1,k]-delta)/2)
        else
            push!(complexRoots,(complex(-res[1,k],delta))/2,(complex(-res[1,k],-delta))/2)
        end
    else
        push!(realRoots,-res[2,k])
    end
    return realRoots,complexRoots
end
solve(problem::PolyProblem)=solve(problem,Bairstov())