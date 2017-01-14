namespace GeographicLib

open System.Runtime.InteropServices
open System

module internal Utilities =
    let inline swap<'a> (x : 'a byref) (y : 'a byref) =
        let temp = x
        x <- y
        y <- temp

module MathLib =
    let degrees (x : float) = x*1.0<rad> |> UnitConversion.degrees
    let sin (x : float<deg>) = x |> UnitConversion.radians |> float |> Math.Sin
    let cos (x : float<deg>) = x |> UnitConversion.radians |> float |> Math.Cos
    let hypot (x, y) =
        let mutable x = abs x
        let y = abs y
        let mutable t = min x y
        x <- max x y
        t <- t / x
        x * sqrt (1.0 + t * t)
    let norm (x : float byref) (y : float byref) =
        let h = hypot (x, y)
        x <- x / h
        y <- y / h

    [<Struct>]
    type ValueAndError(value : float<deg>, error : float<deg>) =
        member __.Value = value
        member __.Error = error
    let private liftToKeepUnits f (arg : float<'u>) : float<'u> =
        float arg |> f |> LanguagePrimitives.FloatWithMeasure
    let private liftToKeepUnits2 f (arg1 : float<'u>, arg2 : float<'u>) : float<'u> =
        (float arg1, float arg2) |> f |> LanguagePrimitives.FloatWithMeasure
    let private remainder = Math.IEEERemainder |> liftToKeepUnits2
    let sum(u : float<deg>, v : float<deg>) =
        let s = u + v;
        let mutable up = s - v;
        let mutable vpp = s - up;
        up <- up - u;
        vpp <- vpp - v;
        new ValueAndError(s, -(up + vpp))

    let atanh(x) = (log(1.0 + x) - log(1.0 - x))/2.0
    let rec getMachineEpsilon eps = if 1.0 + eps = 1.0 then eps else getMachineEpsilon ((abs eps) / 2.0)
    let eatanhe(x, es) = if es > 0.0 then es * atanh(es * x) else -es * atan(es * x)
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
        (if normalisedDiff = 180.0<deg> && differenceOfRemainders.Error < 0.0<deg> then -180.0<deg> else normalisedDiff) - differenceOfRemainders.Error

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

    let sinCosSeries(sinp, sinx, cosx, (c : float[]), n) =
        // Evaluate
        // y = sinp ? sum(c[i] * sin(2*i * x), i, 1, n) : sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
        // using Clenshaw summation.  N.B. c[0] is unused for sin series
        // Approx operation count = (n + 5) mult and (2 * n + 2) add
        let currentIndex = ref (n + (if sinp then 1 else 0))
        let ar = 2.0 * (cosx - sinx) * (cosx + sinx) // 2 * cos(2 * x)
        let mutable y0 = 0.0
        let mutable y1 = 0.0
        if n &&& 1 > 0 then
            decr currentIndex
            y0 <- c.[currentIndex.Value]
        let step = ref (currentIndex.Value/2)
        while step.Value > 0 do
            decr step
            decr currentIndex
            y1 <- ar * y0 - y1 + c.[currentIndex.Value]
            decr currentIndex
            y0 <- ar * y1 - y0 + c.[currentIndex.Value]
        if sinp then 
            2.0 * sinx * cosx * y0 // sin(2 * x) * y0
        else
            cosx * (y0 - y1)

    let atan2d y x =
        let mutable y, x = y, x
        // In order to minimize round-off errors, this function rearranges the
        // arguments so that result of atan2 is in the range [-pi/4, pi/4] before
        // converting it to degrees and mapping the result to the correct
        // quadrant.
        let mutable q = 0
        if abs(y) > abs(x) then 
            Utilities.swap &x &y
            q <- 2
        if x < 0.0 then 
            x <- -x
            q <- q + 1
        // here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4]
        let mutable ang = atan2 y x |> degrees
        match q with
        // Note that atan2d(-0.0, 1.0) will return -0.  However, we expect that
        // atan2d will not be called with y = -0.  If need be, include
        //
        //   | 0: ang = 0 + ang; break;
        //
        // and handle mpfr as in AngRound.
        | 1 -> ang <- (if y > 0.0 then 180.0<deg> else -180.0<deg>) - ang
        | 2 -> ang <- 90.0<deg> - ang
        | 3 -> ang <- -90.0<deg> + ang
        | _ -> raise <| new ArgumentException()
        ang
