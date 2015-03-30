{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE PatternSynonyms #-}

import Data.Array.Repa.FFTW
import qualified Data.Array.Repa as A
import Data.Array.Repa (Z(..), DIM1, DIM3, Array, (:.)(..))
import Data.Array.Repa.Eval (Elt(..))
import Data.Complex
import Numeric.SpecFunctions.Bessel
import Numeric.GSL.Special.Bessel

data Pos a = Pos { posR, posPhi, posZ :: !a }
           deriving (Show, Eq, Ord, Functor)

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

j2 = bessel_Jn 2

pattern Ix3 x y z = Z :. x :. y :. z

-- lengths in nanometers

main =
    let sh@(Ix3 nr nphi nz) = A.ix3 128 128 128
        dr = 10 / realToFrac nr
        dphi = 2*pi / realToFrac nphi
        dz = 20 / realToFrac nz
        coords = A.fromFunction sh (\(Ix3 r phi z) -> fmap realToFrac $ Pos (dr*realToFrac r) (dphi*realToFrac phi) (dz*realToFrac z))
        alpha = 1
        beta = 1
        k = 2*pi/514
     in do a <- A.computeUnboxedP $ A.map (energyDensity 100 alpha beta k) coords
           print a

