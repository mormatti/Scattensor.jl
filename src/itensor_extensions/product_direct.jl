"""
    kron(mps1::MPS, mps2::MPS) -> MPS
    kron(mpo1::MPO, mpo2::MPO) -> MPO

Kronecker product of two MPS or MPO objects.
The resulting MPS or MPO has the length equal to the sum of the lengths of the two inputs.
It simply concatenates the tensors of the two MPS or MPO objects.

# Example
    julia> mps1 = random_mps(siteinds(3, 2))
    MPS
    [1] ((dim=3|id=322|"Site,n=1"), (dim=1|id=654|"Link,l=1"))
    [2] ((dim=1|id=654|"Link,l=1"), (dim=3|id=29|"Site,n=2"))
    julia> mps2 = random_mps(siteinds(3, 3))
    MPS
    [1] ((dim=3|id=290|"Site,n=1"), (dim=1|id=519|"Link,l=1"))
    [2] ((dim=1|id=519|"Link,l=1"), (dim=3|id=773|"Site,n=2"), (dim=1|id=498|"Link,l=2"))
    [3] ((dim=1|id=498|"Link,l=2"), (dim=3|id=337|"Site,n=3"))
    julia> mps_result = kron(mps1, mps2)
    MPS
    [1] ((dim=3|id=322|"Site,n=1"), (dim=1|id=654|"Link,l=1"))
    [2] ((dim=1|id=654|"Link,l=1"), (dim=3|id=29|"Site,n=2"))
    [3] ((dim=3|id=290|"Site,n=1"), (dim=1|id=519|"Link,l=1"))
    [4] ((dim=1|id=519|"Link,l=1"), (dim=3|id=773|"Site,n=2"), (dim=1|id=498|"Link,l=2"))
    [5] ((dim=1|id=498|"Link,l=2"), (dim=3|id=337|"Site,n=3"))
"""
function Base.kron(train1::TrainType, train2::TrainType) where {TrainType <: Union{MPS, MPO}}
    result = TrainType(length(train1) + length(train2))
    for i in eachindex(train1)
        result[i] = train1[i]
    end
    for j in eachindex(train2)
        result[length(train1) + j] = train2[j]
    end
    return result
end