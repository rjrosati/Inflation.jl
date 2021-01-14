import SymPy.fn_map
SymPy.fn_map["Mul"] = :__prod__

function Quote(ex)
    if ex isa Real
        return ex
    else
        return convert(Expr, ex)
    end
end

function QuoteFn(name,ex,vars=free_symbols(ex))
    quot= Expr(:function, Expr(:call, Symbol(name), map(Symbol,vars)...), Quote(ex))
    if name == "_Eh"
        println(quot)
    end
    return quot
end

function QuoteFnCSE(name,cse_ex,vars=free_symbols(cse_ex[2]))
    quot = Expr(:block)
    quot.head = :function
    push!(quot.args, :($(Symbol(name))($(map(Symbol,vars)...) )))
    fnbody = QuoteCSE(cse_ex)
    push!(quot.args,fnbody)
    return quot
end

function QuoteCSE(cse_ex)
    cse,ex = cse_ex
    ex=ex[1]
    fnbody=Expr(:block)
    for (k,v) in cse
        push!(fnbody.args, :($(Quote(k)) = $(Quote(v))) )
    end
    push!(fnbody.args,Quote(ex))
    return fnbody
end

function QuoteArrCSE(cse_ex)
    # assume we did simultaneous CSE on the array
    # cse is common list of substitutions
    # ex is array of cse'd expressions
    cse,ex = cse_ex
    fnbody=Expr(:block)
    for (k,v) in cse
        push!(fnbody.args, :($(Quote(k)) = $(Quote(v))) )
    end
    if length(ex) == 1
        ex = ex[1]
    end
    push!(fnbody.args,:([$(Quote.(ex)...)]))
    return fnbody
end
function QuoteArr(ex)
    fnbody=Expr(:block)
    push!(fnbody.args,:([$(Quote.(ex)...)]))
    return fnbody
end

QuoteFnArr(name,ex, vars=free_symbols(ex)) = Expr(:function,
                          Expr(:call, Symbol(name), map(Symbol,vars)...),
                         :([$(Quote.(ex)...)]))

function QuoteFnArrCSE(name, cse_ex, vars=free_symbols(ex))
    quot = Expr(:block)
    quot.head = :function
    push!(quot.args, :($(Symbol(name))($(map(Symbol,vars)...) )))
    fnbody = QuoteArrCSE(cse_ex)
    push!(quot.args,fnbody)
    return quot
end
