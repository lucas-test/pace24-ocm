# PACE 2024 OCM

## How to use

Build with:

    make build

Run on a single graph:

    ./bin/main graphs/exact/1.gr

Run on every exact graph and creates `results.out`

    ./all.sh graphs/exact

## Lower bounds

### Between low vertices

We denote by mc(u,v) the minimal number of crossings between vertices u and v.
There are two cases: either u < v either v < u. In both cases the exact position does not influence on the number of crossings between edges adjacent to these vertices.
We denote by `cc(u,v)` the number of crossings between the edges adjacent to u and v given that u < v.
Then 

    mc(u,v) = min cc(u,v), cc(v,u)

We deduce a lower bound on the OCM

    Sum( mc(u,v) ) <= OCM(G)


