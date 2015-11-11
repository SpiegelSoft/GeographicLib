﻿namespace GeographicLib

open System.Runtime.InteropServices
open System

module MathLib =
    let degreesPerRadian = 1.0<deg/rad> * 180.0/Math.PI
    let radiansPerDegree = 1.0<rad/deg> * Math.PI/180.0
    let radians (x : float<deg>) = x * radiansPerDegree
    let degrees (x : float<deg>) = x * degreesPerRadian
    let sin (x : float<deg>) = x |> radians |> float |> Math.Sin
    let cos (x : float<deg>) = x |> radians |> float |> Math.Cos
    let hypot (x, y) =
        let mutable x = abs x
        let y = abs y
        let mutable t = min x y
        x <- max x y
        t <- t / x
        x * sqrt (1.0 + t * t)
    let norm (x, y) =
        let h = hypot (x, y)
        (x / h, y / h)

    [<Struct>]
    type ValueAndError(value : float<deg>, error : float<deg>) =
        member this.Value = value
        member this.Error = error
    let liftToKeepUnits f (arg : float<'u>) : float<'u> =
        float arg |> f |> LanguagePrimitives.FloatWithMeasure
    let liftToKeepUnits2 f (arg1 : float<'u>, arg2 : float<'u>) : float<'u> =
        (float arg1, float arg2) |> f |> LanguagePrimitives.FloatWithMeasure
    let remainder = Math.IEEERemainder |> liftToKeepUnits2
    let sum(u : float<deg>, v : float<deg>) =
        let s = u + v;
        let mutable up = s - v;
        let mutable vpp = s - up;
        up <- up - u;
        vpp <- vpp - v;
        new ValueAndError(s, -(up + vpp))

    let atanh(x) = (log(1.0 + x) - log(1.0 - x))/2.0
    let rec getMachineEpsilon eps = if 1.0 + eps = 1.0 then eps else getMachineEpsilon ((abs eps) / 2.0)
    let eatanhe(x, es) = if es >= 0.0 then es * atanh(es * x) else -es * atan(es * x)
    let sq(x : float<'a>) = x * x
    let minFloat = Double.Epsilon
    let epsilon = getMachineEpsilon 0.1
    let epsilonSquared = sq(epsilon)
    let epsilonToTheFourth = sq(epsilonSquared)
    let tolerance = sqrt epsilon
    let tolerance0 = tolerance * epsilon |> sqrt |> sqrt
    let piOver2 = Math.PI / 2.0
    let polyval(N, p : float[], x) =
        let mutable y = if N < 0 then 0.0 else p.[0]
        let mutable N = N
        let mutable i = 1
        while N > 0 do
            N <- N - 1
            y <- y * x + p.[i]
            i <- i + 1
        y
    let z = 1.0<deg>/16.0
    let angRound(x : float<deg>) =
        let mutable y = abs(x)
        y <- if y < z then z - (z - y) else y
        match Constants.Precision with
        | 4 -> if x <= 0.0<deg> then 0.0<deg> - y else y
        | 5 -> if x < 0.0<deg> then 0.0<deg> - y else y
        | _ -> if x < 0.0<deg> then 0.0<deg> - y else y

    let angNormalise(x : float<deg>) =
        let reverse180 angle = if angle = 180.0<deg> then -180.0<deg> else angle 
        reverse180(remainder(x, 360.0<deg>))
    let angDiff(x : float<deg>, y : float<deg>) =
        let differenceOfRemainders = sum(remainder(x, 360.0<deg>), remainder(-y, 360.0<deg>))
        let normalisedDiff = angNormalise differenceOfRemainders.Value
        (if normalisedDiff = 180.0<deg> & differenceOfRemainders.Error < 0.0<deg> then -180.0<deg> else normalisedDiff) - differenceOfRemainders.Error

    // Evaluate the sine and cosine function with the argument in degrees
    //
    // @tparam T the type of the arguments.
    // @param[in] x in degrees.
    // @param[out] sinx sin(<i>x</i>).
    // @param[out] cosx cos(<i>x</i>).
    //
    // The results obey exactly the elementary properties of the trigonometric
    // functions, e.g., sin 9&deg; = cos 81&deg; = &minus; sin 123456789&deg;.
    //*********************************************************************/
    let sincos(x : float<deg>) =
        // In order to minimize round-off errors, this function exactly reduces
        // the argument to the range [-45, 45] before converting it to radians.
        let mutable r = remainder(x, 360.0<deg>)
        let q = Math.Floor(r / 90.0<deg> + 0.5) |> int
        r <- r - (90.0<deg> * float q)
        match (q &&& 3) with
        | 0 -> (sin r, cos r)
        | 1 -> (cos r, -sin r)
        | 2 -> (-sin r, -cos r)
        | _ -> (-cos r, sin r)
