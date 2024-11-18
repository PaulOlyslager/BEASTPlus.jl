struct IdentityMatrix <: BEAST.Operator end

function assemble(op::IdentityMatrix, test_functions, trial_functions;
    storage_policy = Val{:bandedstorage},
    quadstrat=nothing) 
    return sparse(I, length(test_functions), length(trial_functions))
end

# function assemble!(op::IdentityMatrix, test_functions, trial_functions;
#     storage_policy = Val{:bandedstorage},
#     quadstrat=nothing) 
#     return sparse(I, length(test_functions), length(trial_functions))
# end