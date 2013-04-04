import Linear hiding (trace)
import Optimization.Constrained.ProjectedSubgradient
import Debug.Trace
import Numeric
import Data.Foldable
import Data.List

showV2 :: V2 Double -> String
showV2 = intercalate " " . toList . fmap (\x->showFFloat (Just 5) x "")

a = V2 1 2
b = 1
x0 = V2 0 1

constraints = [ Constr EQ 1 (V2 1 1)
              , Constr GT 0 (V2 1 0)
              , Constr GT 0 (V2 0 1)
              ]

main = do
    let steps = optimalStepSched 0
        proj = linearProjection constraints . \x->trace (showV2 x) x
    forM_ (linearProjSubgrad steps proj a b x0) $ \x->
        putStrLn $ showV2 x++"\t"++show (a `dot` x - b)
