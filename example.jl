include("main.jl")

A = [0.18 0.32 0 0 0.5 0 0 0 0 0; #A is the normalised weight matrix of the network
    0.1 0.2 0 0 0 0.2 0.3 0.1 0 0.1;
    0.2 0.2 0.8 0 0 0 0 0 0 0;
    1 0 0 0 0 0 0 0 0 0;
    0.1 0 0 0.1 0 0 0 0.1 0.6 0.1;
    0 0 0 0 0 0.6 0.3 0.1 0 0;
    0 0.3 0 0 0 0.2 0.4 0 0.1 0;
    0.9 0 0.1 0 0 0 0 0 0 0;
    0 0 0 0.7 0 0 0 0 0.3 0;
    0.1 0 0.1 0 0 0 0 0 0.8 0]
ver = ["S" "W" "A" "C" "R" "DS" "DW" "FPI" "H" "PKI"] #ver is the vector containing the names of vertices
Aanalysis(A, ver, graph=true)

C = MC(A, ver, itmc=10000, itr = 1000, var = 0.1, deb = false, mdvar=false, savef=true)
