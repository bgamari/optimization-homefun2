import Linear
import Control.Monad
import Optimization.Constrained.ProjectedSubgradient

a = V2 1 2
b = 1
x0 = V2 0 1

constraints = [ Constr EQ 1 (V2 1 1)
              , Constr GT 0 (V2 1 0)
              , Constr GT 0 (V2 0 1)
              ]

main = do
    let steps = optimalStepSched (-1)
        proj = linearProjection constraints
    forM_ (linearProjSubgrad steps proj a b x0) print
