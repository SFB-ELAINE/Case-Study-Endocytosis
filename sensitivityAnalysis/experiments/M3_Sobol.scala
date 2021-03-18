import sessl._
import sessl.mlrules._


execute {
  new Experiment with Observation with ParallelExecution with CSVOutput with CSVInput {

    model = "../../models/M3_Wnt.mlrj"

    simulator = HybridSimulator()

    parallelThreads = -1

    val stoppingTime = 61
    stopTime = stoppingTime
    replications = 1

    designFromCSV("designOfExperiment.csv")

    observe("Lrp6" ~ count("Cell/Membrane/*/Lrp6"))
    observe("Lrp6Axin" ~ count("Cell/Membrane/*/Lrp6Axin"))
    
    observeAt(stoppingTime)
    withRunResult(writeCSV)
  }
}
