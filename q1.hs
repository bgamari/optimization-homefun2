{-# LANGUAGE RecordWildCards, StandaloneDeriving, DeriveGeneric, GeneralizedNewtypeDeriving,
             DeriveFunctor, DeriveFoldable, DeriveTraversable, RankNTypes,
             FlexibleContexts #-}

import Control.Applicative
import Data.Foldable
import Data.Foldable as F
import Data.Functor
import Data.Reflection
import Data.Traversable
import GHC.Generics
import Linear hiding (trace)
import Numeric.AD hiding (gradientDescent)
import Numeric.AD.Internal.Reverse
import Numeric.AD.Mode.Reverse
import Numeric.AD.Types (lowerFU, auto)
import Optimization.LineSearch

import Debug.Trace

data Blends a = Blends { b1, b2, b3 :: a }
              deriving (Show, Read, Generic, Functor, Foldable, Traversable)
instance Additive Blends where
    zero = pure 0
instance Metric Blends
instance Applicative Blends where
    pure a = Blends a a a
    Blends a1 a2 a3 <*> Blends b1 b2 b3 = Blends (a1 b1) (a2 b2) (a3 b3)

data Raws a = Raws { r1, r2, r3, r4 :: a }
            deriving (Show, Read, Generic, Functor, Foldable, Traversable)
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

data Config a = Config { xs :: Raws (Blends a)
                       , ys :: Raws a
                       , zs :: Blends a
                       }
              deriving (Show, Read, Generic, Functor, Foldable, Traversable)
instance Applicative Config where
    pure a = Config (pure $ pure a) (pure a) (pure a)
    Config x1 y1 z1 <*> Config x2 y2 z2 =
        Config ((<*>) <$> x1 <*> x2) (y1 <*> y2) (z1 <*> z2)
instance Additive Config where zero = pure 0
instance Metric Config

-- | To be minimized
objective :: Fractional a => Config a -> a
objective (Config {..}) = negate $ f2 + f3 - f1
  where f1 = ps `dot` ys
        f2 = p's `dot` zs
        f3 = es `dot` (es ^-^ fmap F.sum xs)

p0 = Config { xs = Raws (pure 100) (pure 100) (pure 100) (pure 100)
            , ys = Raws 1000 1000 1000 1000
            , zs = Blends 1000 1000 1000
            }

main = do
    print $ objective p0
    let p0' = project p0
    print $ objective p0'
    forM_ (optimize p0') $ \p->print (objective p, p)

newtype FU f a = FU (forall s. Reifies s Tape => f (Reverse a s) -> Reverse a s)

constraints :: (Fractional a, Ord a) => [FU Config a]
constraints =
    map (\(FU b)->FU $ \(Config {..})->F.sum (fmap b xs) - b zs) fubs
 ++ map (\(FU r)->FU $ \(Config {..})->min 0 $ r $ bs ^-^ fmap F.sum xs) furs
 ++ [ FU $ \(Config {..})->min 0 $ b1 $ ds ^-^ zs
    , FU $ \(Config {..})->min 0 $ b3 $ ds ^-^ zs
    ]
 ++ map (\(FU b)->FU $ \(Config {..})->min 0 $ negate $ b rs * b zs - (os `dot` fmap b xs)) fubs
 ++ map (\(FU b)->FU $ \(Config {..})->min 0 $ b $ r1 xs) fubs
 ++ map (\(FU b)->FU $ \(Config {..})->min 0 $ b $ r2 xs) fubs
 ++ map (\(FU b)->FU $ \(Config {..})->min 0 $ b $ r3 xs) fubs
 ++ map (\(FU b)->FU $ \(Config {..})->min 0 $ b $ r4 xs) fubs
 ++ map (\(FU r)->FU $ \(Config {..})->min 0 $ r ys) furs
 ++ map (\(FU b)->FU $ \(Config {..})->min 0 $ b zs) fubs
 where fubs = [FU b1, FU b2, FU b3]
       furs = [FU r1, FU r2, FU r3, FU r4]

project :: (Fractional a, Ord a, Show a) => Config a -> Config a
project c@(Config {..}) =
    case unmet of
      []          -> c
      otherwise   -> --traceShow (length unmet, c) $
                     project $ foldl' fixConstraint c unmet
  where unmet = filter (not . met c) constraints
        met c (FU constr) = let (y,_) = grad' constr c in abs y < 1
        fixConstraint c (FU constr) = head $ dropWhile (\c'->not $ met c' $ FU constr)
                                      $ gradientDescent (FU $ (^2) . constr) c

-- | Minimize the given objective
projGradientDescent :: (Additive f, Traversable f, Metric f, Ord a, Fractional a)
                    => (f a -> f a) -> FU f a -> f a -> [f a]
projGradientDescent proj (FU obj) = go
  where go x0 = let p = negated $ df x0
                    alpha0 = armijoSearch 0.1 0.2 100 f df p x0
                    f x = let (y,_) = grad' obj x in y
                    x1 = proj $ x0 ^+^ alpha0 *^ p
                in x1 : go x1
        df = grad obj

gradientDescent :: (Additive f, Traversable f, Metric f, Ord a, Fractional a)
                => FU f a -> f a -> [f a]
gradientDescent = projGradientDescent id

optimize :: (Fractional a, Ord a, Show a) => Config a -> [Config a]
optimize = projGradientDescent project (FU objective)