import sessl._
import sessl.mlrules._
import sessl.analysis.sensitivity._


analyze((params, objective) =>
  execute(

    new Experiment with Observation with ExpressionObservation {

      model = "../models/M1_General.mlrj"
      set("nLRP6" <~ 4000)
      set("nWnt" <~ 2000)
      for ((name, value) <- params) {
        set(name <~ value)
      }

      simulator = SimpleSimulator()
      replications = 300
      stopTime = 61

      observeAt(61)
      val R = observe(count("Membrane/Lrp6 (uB) "))
      val LR = observe(count("Membrane/Lrp6 (B) "))
      val allR = observe(Expr(N => N(R) + N(LR)))

      withExperimentResult(results =>
        objective <~ results.mean(allR)
      )

    }
  )
) using new TwoLevelFullFactorialSetup {

    baseCase(
      "ke" <~ 0.2
    )
    sensitivityCase(
      "ke" <~ 0.4
    )

    withAnalysisResult(result => {
      result.printCompleteReport()
    })
}
