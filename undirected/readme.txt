In undirected version, there are total three steps.
    step1: generate Influence and Susceptibility of individuals according to the unifrom distribution from 0 to 1,
              and calculate individuals f and g, finally get the estimated I, S, and estimated edge weight p_{ij}.
    step2: Run independent cascades model to disturb the edge weight p_{ij}, and get the disturbed p_{ij}, caculate disturbed f and g
    step3: take disturbed f and g as input to our IS-algorithm to get individuals' I and S.

    run the program in command line: ./main.sh

explanation of the prameters of three steps.
    ##  for program step1
    ##  parameter1: dataset name.
    ##  paraemter2: the infinitesimal to avoid divide by zero.
    ##  parameter3: the iteration times.

    ## for program step2
    ##  parameter1: dataset name.
    paraemter2: the infinitesimal
    paraemter3: realizations of Independent cascades process.

    ## for program step3
    ##  parameter1: dataset name.
    ##  realizations of Independent cascades process(here the parameter just compose the input filename)
    ##  the iteration times.

Tips:
     The effect of step2 is to simulate the noise of realistic situation.
     If the researchers do not want to implement step2, we also provide a simplified verision, which only include step3.