# Operational Research Course Project 2021 -- Vehicle Routing Problem

## Context
This project was in the Operational Research course scope, for the last year of the B.S. in Mathematics & Computer Science at the University of Nantes (France).

The following imaginary scenario was given: COVID-19 sick people are totally isolated at home. Drones are deployed to deliver provisions and supplies.

In each area, a unique depot (indexed `1`) stores food and medical resources. In this VRP-like _(Vehicle Routing Problem)_ problem, two constraints come into play:
1. a drone has a maximum load capacity of carrying goods;
2. every person/client/node (indexed from `2` to `n`) is visited only once.

**The aim of this problem is to minimise the distance that drones would travel.**

Two methods of solving are offered and implemented in this project, in a Julia code base.

## Exact solving
This method calls combinatorial optimisation techniques, which can however leads to very high execution time to the stakeholders, especially on very large instances.

Firstly, the set of clients (implemented as a `Vector`, and the clients are indexed from `2` to `n`) is partitionned into subsets, such that the sum of client demands does not outreach the vehicle capacity. This fulfills the first constraint (see the _Context_ section above).

Secondly, for each subset of clients, the shortest cycle is expected; this is a TSP-like _(Traveling Salesman Problem)_ problem which can be solved using the `TravelingSalesmanExact` library in Julia.

Finally, a linear program (LP) solves which cycles to eventually consider, with the objective function minimising the total traveling distance.

## Clark and Wright's algorithm as an approach

_A cycle/tour/path `1 -> x1 -> x2 -> ... -> xn -> 1` is denoted as a list/vector like: `[1,x1,x2,...,xn,1]`._

Clark and Wright's algorithm gives an approximative -- but acceptable -- solution in a much quicker time line. Its principle consists of improving iteratively an initial admissible solution by merging cycles. Here are the key steps:

1. Elementary cycles are initialised (in the form of `[1,i,1]`)
2. Given two cycles of the form `[1,...,i,1]` and `[1,j,...,1]` (where `...` can be empty), distance saved between two nodes `i` and `j` (`i` != `j`, `i,j>=2`) is computed: `s_ij := dist[i,1] + dist[1,j] - dist[i,j]`. _If this saving is positive, and if the sum of demands does not exceed the capacity limit, then these cycles can be merged._
3. All the computed savings are sorted in descending order, with the aim to target firstly the highest savings, prior to any merger.
4. The set of the ready-for-merge cycles is iterated in order to merge the cycles until there's no more merger possible. For each freshly-merged cycle of the form `[1,a,...,i,j,...,b,1]`, saving tuples `(a,b), (_,j), (i,_), (j,_), (_,i)` are ignored to minimise the number of iterations.
5. Finally, for each merger tour, the minimal distance is easily computed by summing the distances between each node/client.



## Data samples
Numerical instances are given in the `assets/{A|B}` folder as raw text in `.dat` files. Wrapper data structure and data reading function are given by the `DataVRP.jl` script, imported in both solving scripts.

For the sake of simplicity, data are integer numbers, but the project can be slightly re-adapted to handle float numbers, units, more data...

What is so special for `B` from `A`? Carriage capacity is higher in `B` data, as well as client demands. Hence, more computation to partition the clients under the capacity constraint, more client subsets, and more computational time to solve the linear program (exact solving)! Benchmarks down below have proven these statements.

## Setup
1. Make sure Julia is installed (PATH variable correctly set as well)
2. Install/Activate the environment in the local repository folder by running Julia REPL:
```
julia> #### PRESS ']' TO ENTER pkg MODE

(@v1.x) pkg> activate VRP_env

(VRP_env) pkg> instanciate    # installs the dependencies

(VRP_env) pkg> st            # should display required packages: GLPK, JuMP, TravelingSalesmanExact
```
3. From now on, activate the environment `(@v1.x) pkg> activate VRP_env` prior to running any script.

## Run
1. Enter Julia REPL
2. Activate the environment `(@v1.x) pkg> activate VRP_env` prior to running any script
3. `julia> include("{Exact|HeuristicClarkWright}Solving.jl")`

To run `ExactSolving.jl` main function: 

```
julia> timer_res_exact("assets/{A|B}/VRP{A|B}{10|..|50}.dat")
```

For `HeuristicClarkWright.jl`:
```
julia> data_then_solve_appr("assets/{A|B}/VRP{A|B}{10|..|50}.dat")
```

## Quality benchmarking

### _A_ samples

#### Exact solving

| Instance size (nbClients) | 10   | 15   | 20   | 25   | 30   | 35   | 40    | 45    | 50    |
|---------------------------|------|------|------|------|------|------|-------|-------|-------|
| Number of cycles found    | 4    | 6    | 8    | 8    | 10   | 13   | 13    | 15    | 17    |
| Value of _z_              | 3656 | 4735 | 4818 | 5932 | 7279 | 9583 | 10893 | 11889 | 11666 |

#### C&W solving

| Instance size (nbClients) | 10   | 15   | 20   | 25   | 30   | 35   | 40    | 45    | 50    |
|---------------------------|------|------|------|------|------|------|-------|-------|-------|
| Number of cycles found    | 4    | 6    | 8    | 9    | 10   | 13   | 14    | 17    | 17    |
| Value of _z_              | 3656 | 4735 | 4894 | 6289 | 7279 | 9815 | 11199 | 12879 | 11686 |
| Error margin Exact VS C&W (%)  | 0   | 0   | 1.58   | 6.02   | 0   | 2.42   | 2.81    | 8.33    | 0.17    |



### _B_ samples

#### Exact solving


| Instance size (nbClients) | 10   | 15   | 20   | 25   | 30   | 35   | 40    | 45    | 50    |
|---------------------------|------|------|------|------|------|------|-------|-------|-------|
| Number of cycles found    | 2    | 3    | 4    | 5    | 6   | _too long_   | -    | -    | -    |
| Value of _z_              | 2137 | 3080 | 3682 | 4025 | 4947 | -   | -    | -    | -    |

#### C&W solving

| Instance size (nbClients) | 10   | 15   | 20   | 25   | 30   | 35   | 40    | 45    | 50    |
|---------------------------|------|------|------|------|------|------|-------|-------|-------|
| Number of cycles found    | 3    | 3    | 4    | 5    | 6   | 6   | 8    | 10    | 11    |
| Value of _z_              | 2199 | 3080 | 3811 | 4204 | 5063 | 5794 | 7411 | 9323 | 8597 |
| Error margin Exact VS C&W (%)  | 2.90   | 0   | 3.50   | 4.45   | 2.34   | -   | -    | -    | -    |


## Time benchmarking (exact solving)

Time benchmarking is only relevant for exact solving, due to its long -- and complex -- computation, compared to C&W's algorithm. TSP functions runs far longer than other ones, especially on larger numbers of node subsets (more noticeable on *B* samples). See the measures down below.

_Configuration used: Ubuntu, Intel Core i3 @2.3GHz, 12GB RAM_

### _A_ samples
| Instance size (nbClients)         | 10 | 15 | 20 | 25 | 30  | 35 | 40 | 45  | 50 |
|-----------------------------------|--|--|--|--|--|--|--|--|--|
| Number of client subsets          | 79 | 240 | 694 | 2948 | 9126 | 4195 | 17626 | 26984 | 44145 |
| Size of the smallest subset       | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Size of the **largest** subset    | 3 | 4 | 4 | 5 | 6 | 4 | 6 | 6 | 5 |
| CPU time for partitionning        | 0.0001 | 0.0004 | 0.005 | 0.04 | 0.16 | 0.07 | 0.38 | 0.74 | 2.02 |
| CPU time for TSP solving          | 0.02 | 0.07 | 0.312 | 1.34 | 5.01 | 1.82 | 9.53 | 14.62 | 23.66 |
| CPU time for the final LP solving | 0.002 | 0.003 | 0.006 | 0.05 | 1.10 | 1.53 | 0.81 | 5.50 | 7.50 |
| **Total CPU execution time**      | **0.02** | **0.09** | **0.35** | **1.43** | **6.33** | **3.44** | **10.83** | **21.03** | **33.37** |

### _B_ samples


| Instance size (nbClients)         | 10 | 15 | 20 | 25 | 30 | 35 | 40 | 45 | 50 |
|-----------------------------------|--|--|--|--|--|--|--|--|--|
| Number of client subsets          | 288 | 4509 | 42299 | 62714 | 270307 | **_too long to run_** | - | - | - |
| Size of the smallest subset       | 1 | 1 | 1 | 1 | 1 |  |  |  |  |
| Size of the **largest** subset    | 5 | 7 | 8 | 7 | 9 |  |  |  |  |
| CPU time for partitionning        | 0.006 | 0.04 | 2.27 | 4.75 | 116.5 |  |  |  |  |
| CPU time for TSP solving          | 0.17 | 2.85 | 32.5 | 43.3 | 217.2 |  |  |  |  |
| CPU time for the final LP solving | 0.01 | 0.10 | 124.8 | 31.1 | 119.8 |  |  |  |  |
| **Total CPU execution time**      | 0.18 | 2.99 | 159.7 | 79.5 | 454.5 |  |  |  |  |
