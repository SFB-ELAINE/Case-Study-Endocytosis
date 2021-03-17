import sessl._
import sessl.mlrules._


execute {
  new Experiment with Observation with ParallelExecution with CSVOutput with CSVInput {

    model = "../../models/M3_Wnt.mlrj"

    simulator = HybridSimulator()
    //simulator = StandardSimulator()
    parallelThreads = -1

    val stoppingTime = 61
    stopTime = stoppingTime
    replications = 1

    designFromCSV("designOfExperiment.csv")

    observe("Lrp6uPuB" ~ count("Cell/Membrane/Lrp6(uP, uB)"))
    observe("Lrp6PuB" ~ count("Cell/Membrane/Lrp6(P, uB)"))
    observe("Lrp6uPB" ~ count("Cell/Membrane/Lrp6(uP, B)"))
    observe("Lrp6PB" ~ count("Cell/Membrane/Lrp6(P, B)"))
    observe("Lrp6Axinu" ~ count("Cell/Membrane/Lrp6Axin(u)"))
    observe("Lrp6Axinp" ~ count("Cell/Membrane/Lrp6Axin(p)"))
    observe("Raft_Lrp6uPuB" ~ count("Cell/Membrane/LR/Lrp6(uP, uB)"))
    observe("Raft_Lrp6PuB" ~ count("Cell/Membrane/LR/Lrp6(P, uB)"))
    observe("Raft_Lrp6uPB" ~ count("Cell/Membrane/LR/Lrp6(uP, B)"))
    observe("Raft_Lrp6PB" ~ count("Cell/Membrane/LR/Lrp6(P, B)"))
    observe("Raft_Lrp6Axinu" ~ count("Cell/Membrane/LR/Lrp6Axin(u)"))
    observe("Raft_Lrp6Axinp" ~ count("Cell/Membrane/LR/Lrp6Axin(p)"))

    observeAt(stoppingTime)
    withRunResult(writeCSV)
  }
}
