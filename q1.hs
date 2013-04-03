{-# LANGUAGE RecordWildCards, StandaloneDeriving, DeriveGeneric, GeneralizedNewtypeDeriving,
             DeriveFunctor, DeriveFoldable, DeriveTraversable, RankNTypes,
             FlexibleContexts, TemplateHaskell #-}

import Control.Applicative
import Data.Function
import Data.Foldable
import Data.Foldable as F
import Data.Functor
import Data.Reflection
import Data.Traversable
import GHC.Generics
import Linear hiding (trace)
import Control.Lens
import Optimization.LineSearch

import Debug.Trace
import Numeric

data Blends a = Blends { _b1, _b2, _b3 :: a }
              deriving (Show, Read, Generic, Functor, Foldable, Traversable)
makeLenses ''Blends
instance Additive Blends where
    zero = pure 0
instance Metric Blends
instance Applicative Blends where
    pure a = Blends a a a
    Blends a1 a2 a3 <*> Blends b1 b2 b3 = Blends (a1 b1) (a2 b2) (a3 b3)

data Raws a = Raws { _r1, _r2, _r3, _r4 :: a }
            deriving (Show, Read, Generic, Functor, Foldable, Traversable)
makeLenses ''Raws
instance Additive Raws where
    zero = pure 0
instance Metric Raws
instance Applicative Raws where
    pure a = Raws a a a a
    Raws a1 a2 a3 a4 <*> Raws b1 b2 b3 b4 = Raws (a1 b1) (a2 b2) (a3 b3) (a4 b4)

-- | Available barrels per day of raw materials
bs :: Num a => Raws a
bs = Raws 4000 5050 7100 4300

-- | Maximum demand of blends
ds :: (Fractional a) => Blends a
ds = Blends 10000 (1/0) 15000

-- | Purchase price of raws
ps :: Fractional a => Raws a
ps = Raws 31.02 33.15 36.35 38.75

-- | Selling price of blends
p's :: Fractional a => Blends a
p's = Blends 45.15 42.95 40.99

-- | Selling price of excess raws
es :: Fractional a => Raws a
es = Raws lo lo hi hi
  where lo = 36.85 -- < 90 octane
        hi = 38.95 -- > 90 octane

-- | Octane content of raws
os :: Num a => Raws a
os = Raws 68 88 91 99

-- | Minimum octane content of blends
rs :: Num a => Blends a
rs = Blends 95 90 85

data Config a = Config { _xs :: Raws (Blends a)
                       , _ys :: Raws a
                       }
              deriving (Show, Read, Generic, Functor, Foldable, Traversable)
makeLenses ''Config
instance Applicative Config where
    pure a = Config (pure $ pure a) (pure a)
    Config x1 y1 <*> Config x2 y2 =
        Config ((<*>) <$> x1 <*> x2) (y1 <*> y2)
instance Additive Config where zero = pure 0
instance Metric Config

-- | To be minimized
objective :: Fractional a => Config a
objective = negated $ f2 ^+^ f3 ^-^ f1
  where f1 = set ys ps zero
        f2 = set (xs.mapped) p's zero
        f3 = set ys es zero ^-^ set xs (fmap pure es) zero

p0 :: Config Double
p0 = Config { _xs = Raws (pure 100) (pure 100) (pure 100) (pure 100)
            , _ys = Raws 1000 1000 1000 1000
            }

main = do
    let p0' = project p0
    let constant = repeat 10
        nonSummable = map (\k->0.1 / sqrt k) [1..]
        squareSummable = map (\k->1 / k) [1..]
        stepSizes = constant

    forM_ (optimize stepSizes p0') $ \p->do
      putStr $ show $ objective `dot` p
      putStr "\t"
      putStrLn $ show p

data Constraint f a = Constr Ordering a (f a)
                    deriving (Show)

constraints :: (Fractional a, Ord a) => [Constraint Config a]
constraints =
    F.toList ((\b r->Constr LT b $ set (xs.r.mapped) 1 zero)                -- (C2)
              <$> bs <*> Raws r1 r2 r3 r4)
 ++ F.toList ((\d l->Constr LT d $ set (xs.mapped.l) 1 zero)                -- (C3)
              <$> ds <*> Blends b1 b2 b3)
 ++ F.toList ((\r l l'->Constr GT 0                                         -- (C4)
                     $ (over xs (\x ->(*^)<$> os <*> x) $ set (xs.mapped.l) 1 zero)
                       ^-^ set (xs.mapped.l') r zero)
              <$> rs <*> Blends b1 b2 b3 <*> Blends b1 b2 b3)
 ++ map (\rb->Constr GT 0 $ set (xs.rb) 1 zero)                             -- x >= 0
        (do r <- [r1, r2, r3, r4]
            b <- [b1, b2, b3]
            return $ r.b)
 ++ map (\r->Constr GT 0 $ set (ys.r) 1 zero) [r1, r2, r3, r4]              -- y >= 0

project :: (Fractional a, Ord a, Show a, RealFloat a) => Config a -> Config a
project c@(Config {..}) =
    case unmet of
      []          -> c
      otherwise   -> --traceShow ( map snd $ filter (not . met c . fst) $ zip constraints [0..]
                     --          , F.concatMap (\x->showFFloat (Just 1) x "  ") $ map (`ap` c) constraints
                     --          ) $
                     project $ fixConstraint c $ maximumBy (flip compare `on` (`ap` c)) unmet
  where unmet = filter (not . met c) constraints
        ap (Constr _ b a) c = a `dot` c - b
        met c (Constr t a constr) = let y = constr `dot` c - a
                                    in case t of
                                       EQ -> abs y < 1e-1
                                       GT -> y >= 0 || abs y < 1e-1
                                       LT -> y <= 0 || abs y < 1e-1
        fixConstraint c (Constr _ b a) = c ^-^ (a `dot` c - b) *^ a ^/ quadrance a

-- | Minimize the given objective
projGradientDescent :: (Additive f, Traversable f, Metric f, Ord a, Fractional a, Show a)
                    => [a] -> (f a -> f a) -> f a -> a -> f a -> [f a]
projGradientDescent stepSizes proj a b = go stepSizes
  where go (alpha:stepSizes) x0 =
            let p = negated $ df x0
                x1 = proj $ x0 ^+^ alpha *^ p
            in x1 : go stepSizes x1
        df x = a
        f x = a `dot` x - b

optimize :: (Fractional a, Ord a, Show a, RealFloat a)
         => [a] -> Config a -> [Config a]
optimize stepSizes = projGradientDescent stepSizes project objective 0
