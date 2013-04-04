import Prelude hiding (maximum)
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import Statistics.Sample
import Data.Foldable as F
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as M
import qualified Data.Set as S
import Linear
import Data.Distributive
import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.Util as H
import qualified Data.Csv as Csv
import Data.Char (ord)
import qualified Data.ByteString.Lazy as BS
import Control.Applicative
import Data.Traversable as T
import System.Environment (getArgs)
import Debug.Trace

newtype Movie = Movie Int deriving (Show, Ord, Eq)
newtype User = User Int deriving (Show, Ord, Eq)
type Rating = Double
type Observations = M.Map (Movie,User) Rating
type Prediction = Movie -> User -> Rating

-- | Root mean squared error of prediction
rmse :: Observations -> Prediction -> Double
rmse obs pred =
    sqrt $ (F.sum $ map (\((m,u), r)->(r - pred m u)^2) $ M.assocs obs) / n
  where n = realToFrac $ M.size obs

globalMean :: Observations -> Prediction
globalMean obs = \m _ -> mu
  where mu = mean $ VU.fromList $ M.elems obs

globalMovieMean :: Observations -> Prediction
globalMovieMean obs = \m _ -> M.findWithDefault 0 m movies
  where movies = fmap (mean . VU.fromList)
                 $ M.mapKeysWith (++) (\(m,u)->m) $ fmap (:[]) obs

globalUserMean :: Observations -> Prediction
globalUserMean obs = \_ u -> M.findWithDefault 0 u users
  where users = fmap (mean . VU.fromList)
                $ M.mapKeysWith (++) (\(m,u)->u) $ fmap (:[]) obs

-- | Optimal mixture reported in Asendorf
reportedMixture = V2 0.548 0.452 :: V2 Double

mixedMean :: V2 Double -> Observations -> Prediction
mixedMean a obs =
    \m u -> let p = V2 (M.findWithDefault 0 m movies) (M.findWithDefault 0 u users)
            in a `dot` p
  where movies = fmap (mean . VU.fromList)
                 $ M.mapKeysWith (++) (\(m,u)->m) $ fmap (:[]) obs
        users = fmap (mean . VU.fromList)
                $ M.mapKeysWith (++) (\(m,u)->u) $ fmap (:[]) obs

obsToHMatrix :: Observations -> H.Matrix Double
obsToHMatrix obs = H.buildMatrix (nMovies+1) (nUsers+1) f
  where Movie nMovies = S.findMax $ S.map fst $ M.keysSet obs
        User nUsers   = S.findMax $ S.map snd $ M.keysSet obs
        f (m,u) = M.findWithDefault 0 (Movie m, User u) obs

froebeniusNorm :: H.Matrix Double -> Double
froebeniusNorm = H.sumElements . H.mapMatrix (^2)

svdThresh :: Double -> Double -> [Double] -> Observations -> IO Prediction
svdThresh tol tau deltas t = do
    rStar <- go $ svdThresh' tau deltas t
    return $ predict rStar
  where go (r:rs) = let err = relError r
                    in if err > tol
                           then putStrLn ("SVT error = "++show err) >> go rs
                           else return r
        t' = obsToHMatrix t
        normTProj = froebeniusNorm (proj t t')
        relError r = froebeniusNorm (proj t (t' `H.sub` r)) / normTProj
        predict r (Movie m) (User u) = r H.@@> (m,u)

-- | Generate a list of recommendation iterates @R_q@
svdThresh' :: Double -> [Double] -> Observations -> [H.Matrix Double]
svdThresh' tau deltas t = go deltas (H.zeros (nMovies+1) (nUsers+1))
  where t' = obsToHMatrix t
        go :: [Double] -> H.Matrix Double -> [H.Matrix Double]
        go (d:deltas) y0 =
            let r1 = shrink tau y0
                y1 = y0 `H.add` (d `H.scale` proj t (t' `H.sub` r1))
            in r1 : go deltas y1
        Movie nMovies = S.findMax $ S.map fst $ M.keysSet t
        User nUsers   = S.findMax $ S.map snd $ M.keysSet t

-- | Shrink operator
shrink :: Double -> H.Matrix Double -> H.Matrix Double
shrink tau x = let (u,s,v) = H.thinSVD x
                   s' = H.diag $ H.mapVector (\x->max (x-tau) 0) s
               in u `H.mXm` s' `H.mXm` H.trans v

-- | Project recommendations onto observations
proj :: Observations -> H.Matrix Double -> H.Matrix Double
proj t r =
    obsToHMatrix $ M.mapWithKey (\(Movie m, User u) _ -> r H.@@> (m,u)) t

nuclearNorm :: H.Matrix Double -> Double
nuclearNorm = H.sumElements . H.singularValues

robustCompletion :: Double -> Double -> Double -> Observations -> IO Prediction
robustCompletion mu0 rho lambda t =
    let go (a:as) = let error = rmse t (predict a)
                    in do print error
                          go as
        predict a = \(Movie m) (User u) -> a H.@@> (m,u)
    in go $ robustCompletion' mu0 rho lambda t

robustCompletion' :: Double -> Double -> Double -> Observations -> [H.Matrix Double]
robustCompletion' mu0 rho lambda t =
    let go (mu:mus) b0 e0 y0 =
            let l1 = shrink (1/mu) (r `H.sub` b0 `H.sub` e0 `H.sub` H.scale (1/mu) y0)
                d0 = r `H.sub` l1 `H.sub` e0 `H.sub` H.scale (1/mu) y0
                b1 = H.mapMatrix (\d->if abs d <= lambda/mu
                                        then 0 else d - signum d * lambda/mu) d0
                e1 = projComp t (r `H.sub` b1 `H.sub` l1 `H.sub` (1/mu `H.scale` y0))
                y1 = y0 `H.add` H.scale mu (l1 `H.add` b1 `H.add` e1 `H.sub` r)
            in l1 : go mus b1 e1 y1
        r = obsToHMatrix t
        Movie nMovies = S.findMax $ S.map fst $ M.keysSet t
        User nUsers   = S.findMax $ S.map snd $ M.keysSet t
        zeros = H.zeros (nMovies+1) (nUsers+1)
    in go (iterate (*rho) mu0) zeros zeros zeros

-- | Project into the complement of t
projComp :: Observations -> H.Matrix Double -> H.Matrix Double
projComp t =
    H.mapMatrixWithIndex (\(m,u) r->if (Movie m, User u) `M.member` t
                                      then 0 else r)

compToLinear :: Observations -> Map Movie (Map User Rating)
compToLinear =
    M.unionsWith M.union . map (\((m,u),r)->M.singleton m (M.singleton u r)) . M.assocs

-- | Non-negative matrix factorization with non-increasing Euclidean distance
nmf :: (H.Product a, Fractional a, Ord a)
    => Int                        -- ^ The rank of the resulting decomposition
    -> H.Matrix a                 -- ^ Initial W
    -> H.Matrix a                 -- ^ Initial H
    -> H.Matrix a                 -- ^ The matrix to be factored
    -> [(H.Matrix a, H.Matrix a)]
nmf r w0 h0 v = go w0 h0
  where go w0 h0 = let wv  = H.trans w0 `H.mXm` v
                       vh  = v  `H.mXm` H.trans h0
                       eps = 1e-5
                       wwh = H.mapMatrix (max eps) $ H.trans w0 `H.mXm` w0 `H.mXm` h0
                       whh = H.mapMatrix (max eps) $ w0 `H.mXm` h0 `H.mXm` H.trans h0
                       w1  = H.mapMatrixWithIndex (\i w->w * vh H.@@> i / whh H.@@> i) w0
                       h1  = H.mapMatrixWithIndex (\i h->h * wv H.@@> i / wwh H.@@> i) h0
                   in (w1,h1) : go w1 h1

nmfPredict :: Int -> H.Matrix Double -> [Prediction]
nmfPredict r v = map predict $ nmf r w0 h0 v
  where predict (w,h) = let v = w `H.mXm` h
                        in \(Movie m) (User u) -> v H.@@> (m,u)
        w0 = H.diagRect 0 (H.buildVector r (const 1)) n r
        h0 = H.diagRect 0 (H.buildVector r (const 1)) r m
        (n,m) = (H.rows v, H.cols v)

tau = 6000
deltas = repeat 1.9

main = do
    train:test:_ <- getArgs
    trainD' <- readMovieLens train
    testD' <- readMovieLens test
    let maps = buildMaps $ trainD' `M.union` testD'
        trainD = remapIds maps trainD'
        testD = remapIds maps testD'
    putStrLn "Training:" >> describeObservations trainD
    putStrLn "Test:" >> describeObservations testD

    rSvd <- svdThresh 0.2 tau deltas trainD
    --rRobust <- robustCompletion 1 0.2 0.1 trainD
    let nmf50 = head $ drop 10 $ nmfPredict 100 (obsToHMatrix trainD)
    let preds = [ ("global mean",       globalMean trainD)
                , ("movie mean",        globalMovieMean trainD)
                , ("user mean",         globalUserMean trainD)
                , ("reported mixture",  mixedMean reportedMixture trainD)
                -- , ("robust completion", rRobust)
                , ("NMF 50",            nmf50)
                 , ("SVD threshold",     rSvd)
                ]
    forM_ preds $ \(name,pred) -> do
        putStrLn $ (take 30 $ name++repeat ' ')++show (rmse testD pred)

decodeOpts = Csv.defaultDecodeOptions { Csv.decDelimiter = fromIntegral $ ord '\t' }

describeObservations :: Observations -> IO ()
describeObservations obs = do
    putStrLn $ "Movies : "++show (S.size $ S.map fst $ M.keysSet obs)
    putStrLn $ "Users  : "++show (S.size $ S.map snd $ M.keysSet obs)
    putStrLn $ "Observations: "++show (M.size obs)
    putStrLn ""

-- | Parse MovieLens data
-- Format described at <http://www.grouplens.org/system/files/ml-100k-README.txt>
readMovieLens :: FilePath -> IO Observations
readMovieLens fname = do
    recs <- either error id . Csv.decodeWith decodeOpts False <$> BS.readFile fname
    return $ M.unions $ map f $ V.toList recs
  where f :: (Int, Int, Int, Int) -> Observations
        f (user, movie, rating, _) =
            M.singleton (Movie movie, User user) (realToFrac rating)

-- | Remap movie and user identifiers to domain starting with zero
buildMaps :: Observations -> (M.Map Movie Movie, M.Map User User)
buildMaps obs = (movieMap, userMap)
  where mkMapping :: Ord a => [b] -> S.Set a -> M.Map a b
        mkMapping ids keys = M.fromList $ zip (S.toList keys) ids
        movieMap = mkMapping [Movie i | i <- [0..]] $ S.map fst $ M.keysSet obs
        userMap = mkMapping [User i | i <- [0..]] $ S.map snd $ M.keysSet obs

remapIds :: (M.Map Movie Movie, M.Map User User) -> Observations -> Observations
remapIds (movieMap, userMap) = M.mapKeys (\(m,u)->(movieMap M.! m, userMap M.! u))
