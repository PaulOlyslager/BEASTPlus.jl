# make it possible to assign a quadstrat to a specific operator or to a linear combination of operators

struct AssignQuadStrat <: BEAST.AbstractOperator
    op::BEAST.AbstractOperator
    qs
end
export AssignQuadStrat
function BEAST.assemble!(op::AssignQuadStrat,test_functions::BEAST.Space, trial_functions::BEAST.Space,
    store, threading;quadstrat=nothing)
    BEAST.assemble!(op.op,test_functions,trial_functions,store,threading; quadstrat=op.qs)
end

function BEAST.assemble(op::AssignQuadStrat, test_functions, trial_functions;
    storage_policy = Val{:bandedstorage},
    quadstrat=nothing)
    BEAST.assemble(op.op, test_functions, trial_functions; storage_policy=storage_policy, quadstrat=op.qs)
end


struct CloseQuadStrat
    sigma
    nearquadstrat
end

