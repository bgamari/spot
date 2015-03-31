{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE PatternSynonyms #-}

import Data.Array.Repa.FFTW
import qualified Data.Array.Repa as A
import Data.Array.Repa (Z(..), DIM1, DIM3, Array, (:.)(..))
import Data.Array.Repa.Eval (Elt(..))
import Data.Array.Repa.Repr.ForeignPtr (F)
import Data.Complex
import Control.Parallel.Strategies
--import Numeric.SpecFunctions.Bessel
import Numeric.GSL.Special.Bessel

import Data.Monoid
import qualified Data.ByteString.Builder as B
import System.IO (stdout)

data Pos a = Pos { posR, posPhi, posZ :: !a }
           deriving (Show, Eq, Ord, Functor)

showPos :: Pos Double -> B.Builder
showPos (Pos r phi z) = B.doubleDec r <> "\t" <> B.doubleDec phi <> "\t" <> B.doubleDec z

j0 = bessel_J0
j1 = bessel_J1
j2 = bessel_Jn 2

energyDensity :: Int    -- ^ number of thetas to sample in integration
              -> Double -- ^ alpha: the half-angle of objective opening in radians
              -> Double -- ^ beta: the ratio of the back-aperture radius to the beam radius
              -> Double -- ^ wavenumber
              -> Pos Double -- ^ position
              -> Double
energyDensity nThetas alpha beta k (Pos r phi z) = 
    1 / 16 / pi * ((magnitude psi0)^2
                   + 4 * (magnitude psi1)^2 * (cos phi)^2
                   + (magnitude psi2)^2
                   + 2 * cos (2*phi) * realPart (psi0 * conjugate psi2)
                  )
  where
    dtheta = alpha / realToFrac nThetas

    -- optical coordinates
    v = k * r * sin alpha
    u = k * z * (sin alpha)^2

    psis (Z :. i) =
        let theta, q, vq :: Double
            theta = realToFrac i * dtheta
            q = sin theta / sin alpha
            vq = v * q

            qq :: Complex Double
            qq = cis (u * cos theta / (sin alpha)^2)

            a :: Complex Double
            a = realToFrac $ sqrt (cos theta * exp (-2 * beta^2 * q^2))

            sinTh, cosTh :: Complex Double
            sinTh = realToFrac $ sin theta
            cosTh = realToFrac $ cos theta

            !psi0 = a * sinTh * (1 + cosTh) * realToFrac (j0 vq) * qq
            !psi1 = a * sinTh^2 * realToFrac (j1 vq) * qq
            !psi2 = a * sinTh * (1 - cosTh) * realToFrac (j2 vq) * qq
        in (psi0, psi1, psi2)
    psi0, psi1, psi2 :: Complex Double
    (psi0, psi1, psi2) = A.foldAllS addTerms (0,0,0) (A.fromFunction (A.ix1 nThetas) psis) `divBy` realToFrac dtheta
    addTerms (a,b,c) (x,y,z) = (a+x, b+y, c+z)
    divBy (a,b,c) s = (a/s, b/s, c/s)

instance (Num a, Elt a) => Data.Array.Repa.Eval.Elt (Complex a) where
    touch (a :+ b) = touch a >> touch b
    zero = 0 :+ 0
    one = 1 :+ 0

pattern Ix3 x y z = Z :. x :. y :. z

excitationDensity :: A.Array A.D DIM3 Double
excitationDensity = A.map (energyDensity 100 alpha beta k) coords
  where
    alpha = 1
    beta = 1
    k = 2*pi/514

-- lengths in nanometers
-- times in microseconds
coords :: A.Array A.D DIM3 (Pos Double)
coords = A.fromFunction sh (\(Ix3 r phi z) -> Pos (dr*realToFrac r) (dphi*realToFrac phi) (dz*realToFrac z))
  where
    n = 128
    sh@(Ix3 nr nphi nz) = A.ix3 n n n
    dr = 1000 / realToFrac nr
    dphi = 2*pi / realToFrac nphi
    dz = 4000 / realToFrac nz

test = do
     let slice = A.Z :. A.All :. A.All :. (0::Int)
     let exc = A.traverse2 (A.map showPos coords) (A.map B.doubleDec excitationDensity) (const id) (\f g i->f i<>"\t"<>g i<>"\n")
     B.hPutBuilder stdout $ foldl (<>) mempty $ A.toList $ A.slice exc slice

main = doCorr

doCorr = do
     a <- A.computeUnboxedP excitationDensity
     let d = 1
         taus = map (10**) $ linspace (-4) 7 128
         concs = map (\t->A.map (\(Pos r _ _)->diffGreen d t r) coords) taus
         corr = correlate a concs
     putStrLn $ unlines $ zipWith (\t c->show t <> "\t" <> show c) taus corr

linspace :: RealFrac a => a -> a -> Int -> [a]
linspace a b n = [a + d * realToFrac i | i <- [1..n-1]]
  where
    d = (b - a) / realToFrac n

diffGreen :: Double -- ^ diffusivity
          -> Double -- ^ time
          -> Double -- ^ radial position
          -> Double
diffGreen d t r = 1 / (4*pi*d*t)**1.5 * exp (-r^2 / (4*d*t))

correlate :: A.Array A.U DIM3 Double    -- ^ observation volume power density field
          -> [A.Array A.D DIM3 Double]  -- ^ concentration field at discrete time steps
          -> [Double]
correlate obs concs = withStrategy (parList rseq) $ map go concs
  where
    o = fft3d (A.computeS $ A.map realToFrac obs)
    go :: A.Array A.D DIM3 Double -> Double
    go conc = A.sumAllS $ a A.*^ obs
      where
        conc' = A.map realToFrac conc
        a = A.map realPart $ ifft3d $ A.computeS $ conc' A.*^ o
