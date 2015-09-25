namespace GeographicLib

open System

type LowToHighRatio =
    struct
        val Ratio : float
        new(ratio : float) = { Ratio = if ratio <= 1.0 then ratio else 1.0 / ratio }
    end

[<Measure>] type rad
[<Measure>] type deg
[<Measure>] type km
[<Measure>] type m

module Constants = 
    [<Literal>]
    let WGS84_a = 6378137.0<m>
    let WGS84_f = LowToHighRatio(1.0 / ( 298257223563.0 / 1000000000.0 ))

module MathLib =
    let atanh(x : float) = MathNet.Numerics.Trig.Atanh x //(log(1.0 + x) - log(1.0 - x))/2.0
    let rec getMachineEpsilon eps = if 1.0 + eps = 1.0 then eps else getMachineEpsilon ((abs eps) / 2.0)
    let eatanhe(x, es) = if es >= 0.0 then es * atanh(es * x) else -es * atan(es * x)
    let sq(x) = x * x
    let epsilon = getMachineEpsilon 0.1
    let epsilonSquared = sq(epsilon)
    let epsilonToTheFourth = sq(epsilonSquared)
    let tolerance = sqrt epsilon
    let tolerance0 = tolerance * epsilon |> sqrt |> sqrt
    let piOver2 = Math.PI / 2.0

type EllipticFunction(k2 : float, ?alpha2 : float, ?kp2 : float, ?alphap2 : float) =
    let alpha2 = defaultArg alpha2 0.0
    let kp2 = defaultArg kp2 1.0 - k2
    let alphap2 = defaultArg alphap2 1.0 - alpha2

    let RF(x, y) =
        let tolRG0 = 2.7 * sqrt(MathLib.epsilon * 0.01)
        let (xn, yn) = (sqrt x, sqrt y)
        let mutable (xn, yn) = (max xn yn, min xn yn)
        while abs(xn-yn) > tolRG0 * xn do
            let t = (xn + yn) / 2.0
            yn <- sqrt (xn * yn)
            xn <- t
        Math.PI / (xn + yn)

    let RG(x, y) =
        let tolRG0 = 2.7 * sqrt(MathLib.epsilon * 0.01)
        let x0 = sqrt(max x y)
        let y0 = sqrt(min x y)
        let mutable xn = x0
        let mutable yn = y0
        let mutable s = 0.0
        let mutable mul = 0.25
        while abs(xn - yn) > tolRG0 * xn do
            let mutable t = (xn + yn) / 2.0
            yn <- sqrt(xn * yn)
            xn <- t
            mul <- 2.0 * mul
            t <- xn - yn
            s <- s + mul * t * t
        (MathLib.sq( (x0 + y0) / 2.0 ) - s) * Math.PI / (2.0 * (xn + yn))

    let RD(x, y, z) =
        let tolRD = Math.Pow(0.2 * MathLib.epsilon * 0.01, 1.0 / 8.0)
        let A0 = (x + y + 3.0*z) / 5.0
        let mutable An = A0
        let Q = max(max (abs(A0-x)) (abs(A0-y))) (abs(A0-z)) / tolRD
        let mutable x0 = x
        let mutable y0 = y
        let mutable z0 = z
        let mutable mul = 1.0
        let mutable s = 0.0
        while Q >= mul * abs(An) do
            let lam = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0)
            s <- s + 1.0/(mul * sqrt(z0) * (z0 + lam))
            An <- (An + lam)/4.0
            x0 <- (x0 + lam)/4.0
            y0 <- (y0 + lam)/4.0
            z0 <- (z0 + lam)/4.0
            mul <- mul * 4.0

        let X = (A0 - x) / (mul * An)
        let Y = (A0 - y) / (mul * An)
        let Z = -(X + Y) / 3.0
        let E2 = X*Y - 6.0*Z*Z
        let E3 = (3.0*X*Y - 8.0*Z*Z)*Z
        let E4 = 3.0 * (X*Y - Z*Z) * Z*Z
        let E5 = X*Y*Z*Z*Z
        // http://dlmf.nist.gov/19.36.E2
        // Polynomial is
        // (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
        //    - E2^3/16 + 3*E3^2/40 + 3*E2*E4/20 + 45*E2^2*E3/272
        //    - 9*(E3*E4+E2*E5)/68)
        ((471240.0 - 540540.0 * E2) * E5 + (612612.0 * E2 - 540540.0 * E3 - 556920.0) * E4 +
            E3 * (306306.0 * E3 + E2 * (675675.0 * E2 - 706860.0) + 680680.0) +
            E2 * ((417690.0 - 255255.0 * E2) * E2 - 875160.0) + 4084080.0) / (4084080.0 * mul * An * sqrt(An)) + 3.0 * s

    let eps = k2 / MathLib.sq(sqrt(kp2) + 1.0)
    let (Kc, Ec, Dc) = if k2 = 0.0 then (MathLib.piOver2, MathLib.piOver2, MathLib.piOver2 / 2.0) else ((if kp2 = 0.0 then Double.PositiveInfinity else RF(kp2, 1.0)), (if kp2 = 0.0 then 1.0 else 2.0 * RG(kp2, 1.0)), (if kp2 = 0.0 then Double.PositiveInfinity else RD(0.0, kp2, 1.0) / 3.0)) 


type AlbersEqualArea(semiMajorAxis : float<m>, flattening : LowToHighRatio, sinLat1 : float, cosLat1 : float, sinLat2 : float, cosLat2 : float, k1 : float) =
    let a = semiMajorAxis
    let f = flattening.Ratio
    let fm = 1.0 - f
    let e2 = f * (2.0 - f)
    let e = e2 |> abs |> sqrt
    let atanhee(x) = if f > 0.0 then MathLib.atanh(e * x) / e else if f < 0.0 then (atan2 (e * abs x) (if x < 0.0 then -1.0 else 1.0))/e else x
    let e2m = 1.0 - e2
    let qZ = 1.0 + e2m * atanhee(1.0)
    let qx = qZ / (2.0 * e2m)

type TransverseMercator(semiMajorAxis : float<m>, flattening : LowToHighRatio, scale : float) =
    let a = semiMajorAxis
    let f = flattening.Ratio
    let e2 = f * (2.0 - f)
    let es = (if f < 0.0 then -1.0 else 1.0) * sqrt(abs(e2))
    let e2m = 1.0 - e2
    let c =  sqrt(e2m) * exp(MathLib.eatanhe(1.0, es)) 
    let n = f / (2.0  - f)

type Ellipsoid(semiMajorAxis : float<m>, flattening : LowToHighRatio) =
    let a = semiMajorAxis
    let f = flattening.Ratio
    let f1 = 1.0 - f
    let f12 = sqrt f1
    let e2 = f * (2.0 - f)
    let es = (if f < 0.0 then -1.0 else 1.0) * sqrt(abs(e2))
    let e12 = e2 / (1.0 - e2)
    let n = f / (2.0  - f)
    let b = a * f1
    let tm = TransverseMercator(a, flattening, 1.0)
    let e11 = EllipticFunction(-e12)
    let au = AlbersEqualArea(a, flattening, 0.0, 1.0, 0.0, 1.0, 1.0)
    static member WGS84 = Ellipsoid(Constants.WGS84_a, Constants.WGS84_f)