using ITensors
using ITensorMPS



Id = zeros(3^3, 3^3)
for i in 1:3^3
    Id[i,i] = 1
end

mpo_A = matrix_to_mpo(Id, 3, 3)