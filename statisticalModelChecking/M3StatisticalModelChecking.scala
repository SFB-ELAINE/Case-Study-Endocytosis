import sessl._
import sessl.mlrules._
import sessl.verification._
import sessl.verification.mitl._


object Reference {
  val ref: Trajectory[Double] = List(
    (0, 1),
    (10, 0.6408),
    (20, 0.5563),
    (30, 0.4554),
    (60, 0.3944),
    (61, 0.3944)
  )
}

execute(

  new Experiment with Observation with ExpressionObservation 
                            with StatisticalModelChecking {

    model = "../models/M3_Wnt.mlrj"
    set("nLRP6" <~ 4000)
    set("nWnt" <~ 2000)
    set("ke_raft" <~ 0.1)
    set("ke_nonraft" <~ 0.1)

    simulator = SimpleSimulator()
    stopTime = 61

    observeAt(0, 10, 20, 30, 60, 61)
    val Lrp6 = observe(count("Cell/Membrane/Lrp6(_, _)"))
    val Lrp6Axin = observe(count("Cell/Membrane/Lrp6Axin(_)"))
    val Raft_Lrp6 = observe(count("Cell/Membrane/LR/Lrp6(_, _)"))
    val Raft_Lrp6Axin = observe(count("Cell/Membrane/LR/Lrp6Axin(_)"))
    val allR = observe(Expr(N => N(Lrp6) + N(Lrp6Axin) + N(Raft_Lrp6) + N(Raft_Lrp6Axin)))

    test = SequentialProbabilityRatioTest(p = 0.8, alpha = 0.05, beta = 0.05, delta = 0.05)
    prop = MITL( G(0, 60)( 
            (OutVar(allR)/Constant(4000) - Traj(Reference.ref)).abs < Constant(0.2)) )

    withCheckResult { result =>
      println(result.satisfied)
    }
  }
)
