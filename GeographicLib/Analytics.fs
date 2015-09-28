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
    // The precision of floating point numbers used in %GeographicLib.  1 means
    // float (single precision); 2 (the default) means double; 3 means long double;
    // 4 is reserved for quadruple precision.  Nearly all the testing has been
    // carried out with doubles and that's the recommended configuration.  In order
    // for long double to be used, GEOGRAPHICLIB_HAVE_LONG_DOUBLE needs to be
    // defined.  Note that with Microsoft Visual Studio, long double is the same as
    // double.
    [<Literal>]
    let Precision = 2

module MathLib =
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

module GeodesicCoefficients =
    let GeodesicOrder =
        match Constants.Precision with
        | 2 -> 6
        | 1 -> 3
        | 3 -> 7
        | _ -> 8
    let nA1 = GeodesicOrder
    let nC1 = GeodesicOrder
    let nC1p = GeodesicOrder
    let nA2 = GeodesicOrder
    let nC2 = GeodesicOrder
    let nA3 = GeodesicOrder
    let nA3x = nA3
    let nC3 = GeodesicOrder
    let nC3x = (nC3 * (nC3 - 1)) / 2;
    let nC4 = GeodesicOrder
    let nC4x = (nC4 * (nC4 + 1)) / 2;
    // Size for temporary array
    // nC = max(max(nC1_, nC1p_, nC2_) + 1, max(nC3_, nC4_))
    let nC = GeodesicOrder + 1
    [<Literal>]
    let maxit1 = 20
    let A3Coeff =
        match GeodesicOrder with
        | 3 -> 
            [| 
                -1.0; 4.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 4 -> 
            [|
                -1.0; 16.0; // A3 coeff of eps^3
                -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 5 -> 
            [|
                -3.0; 64.0; // A3 coeff of eps^4
                -3.0; -1.0; 16.0; // A3 coeff of eps^3
                3.0; -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 6 -> 
            [|
                -3.0; 128.0; // A3 coeff of eps^5
                -2.0; -3.0; 64.0; // A3 coeff of eps^4
                -1.0; -3.0; -1.0; 16.0; // A3 coeff of eps^3
                3.0; -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 7 -> 
            [|
                -5.0; 256.0; // A3 coeff of eps^6
                -5.0; -3.0; 128.0; // A3 coeff of eps^5
                -10.0; -2.0; -3.0; 64.0; // A3 coeff of eps^4
                5.0; -1.0; -3.0; -1.0; 16.0; // A3 coeff of eps^3
                3.0; -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 8 -> 
            [|
                -25.0; 2048.0; // A3 coeff of eps^7
                -15.0; -20.0; 1024.0; // A3 coeff of eps^6
                -5.0; -10.0; -6.0; 256.0; // A3 coeff of eps^5
                -5.0; -20.0; -4.0; -6.0; 128.0; // A3 coeff of eps^4
                5.0; -1.0; -3.0; -1.0; 16.0; // A3 coeff of eps^3
                3.0; -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | _ -> raise <| new ArgumentException()


//   \brief %Geodesic calculations
//   
//   The shortest path between two points on a ellipsoid at (\e lat1, \e lon1)
//   and (\e lat2, \e lon2) is called the geodesic.  Its length is \e s12 and
//   the geodesic from point 1 to point 2 has azimuths \e azi1 and \e azi2 at
//   the two end points.  (The azimuth is the heading measured clockwise from
//   north.  \e azi2 is the "forward" azimuth, i.e., the heading that takes you
//   beyond point 2 not back to point 1.)  In the figure below, latitude if
//   labeled &phi;, longitude &lambda; (with &lambda;<sub>12</sub> =
//   &lambda;<sub>2</sub> &minus; &lambda;<sub>1</sub>), and azimuth &alpha;.
//   
//   <img src="http://upload.wikimedia.org/wikipedia/commons/c/cb/Geodesic_problem_on_an_ellipsoid.svg" width=250 alt="spheroidal triangle">
//   
//   Given \e lat1, \e lon1, \e azi1, and \e s12, we can determine \e lat2, \e
//   lon2, and \e azi2.  This is the \e direct geodesic problem and its
//   solution is given by the function Geodesic::Direct.  (If \e s12 is
//   sufficiently large that the geodesic wraps more than halfway around the
//   earth, there will be another geodesic between the points with a smaller \e
//   s12.)
//   
//   Given \e lat1, \e lon1, \e lat2, and \e lon2, we can determine \e azi1, \e
//   azi2, and \e s12.  This is the \e inverse geodesic problem, whose solution
//   is given by Geodesic::Inverse.  Usually, the solution to the inverse
//   problem is unique.  In cases where there are multiple solutions (all with
//   the same \e s12, of course), all the solutions can be easily generated
//   once a particular solution is provided.
//   
//   The standard way of specifying the direct problem is the specify the
//   distance \e s12 to the second point.  However it is sometimes useful
//   instead to specify the arc length \e a12 (in degrees) on the auxiliary
//   sphere.  This is a mathematical construct used in solving the geodesic
//   problems.  The solution of the direct problem in this form is provided by
//   Geodesic::ArcDirect.  An arc length in excess of 180&deg; indicates that
//   the geodesic is not a shortest path.  In addition, the arc length between
//   an equatorial crossing and the next extremum of latitude for a geodesic is
//   90&deg;.
//   
//   This class can also calculate several other quantities related to
//   geodesics.  These are:
//   - <i>reduced length</i>.  If we fix the first point and increase \e azi1
//     by \e dazi1 (radians), the second point is displaced \e m12 \e dazi1 in
//     the direction \e azi2 + 90&deg;.  The quantity \e m12 is called
//     the "reduced length" and is symmetric under interchange of the two
//     points.  On a curved surface the reduced length obeys a symmetry
//     relation, \e m12 + \e m21 = 0.  On a flat surface, we have \e m12 = \e
//     s12.  The ratio <i>s12</i>/\e m12 gives the azimuthal scale for an
//     azimuthal equidistant projection.
//   - <i>geodesic scale</i>.  Consider a reference geodesic and a second
//     geodesic parallel to this one at point 1 and separated by a small
//     distance \e dt.  The separation of the two geodesics at point 2 is \e
//     M12 \e dt where \e M12 is called the "geodesic scale".  \e M21 is
//     defined similarly (with the geodesics being parallel at point 2).  On a
//     flat surface, we have \e M12 = \e M21 = 1.  The quantity 1/\e M12 gives
//     the scale of the Cassini-Soldner projection.
//   - <i>area</i>.  The area between the geodesic from point 1 to point 2 and
//     the equation is represented by \e S12; it is the area, measured
//     counter-clockwise, of the geodesic quadrilateral with corners
//     (<i>lat1</i>,<i>lon1</i>), (0,<i>lon1</i>), (0,<i>lon2</i>), and
//     (<i>lat2</i>,<i>lon2</i>).  It can be used to compute the area of any
//     simple geodesic polygon.
//   
//   Overloaded versions of Geodesic::Direct, Geodesic::ArcDirect, and
//   Geodesic::Inverse allow these quantities to be returned.  In addition
//   there are general functions Geodesic::GenDirect, and Geodesic::GenInverse
//   which allow an arbitrary set of results to be computed.  The quantities \e
//   m12, \e M12, \e M21 which all specify the behavior of nearby geodesics
//   obey addition rules.  If points 1, 2, and 3 all lie on a single geodesic,
//   then the following rules hold:
//   - \e s13 = \e s12 + \e s23
//   - \e a13 = \e a12 + \e a23
//   - \e S13 = \e S12 + \e S23
//   - \e m13 = \e m12 \e M23 + \e m23 \e M21
//   - \e M13 = \e M12 \e M23 &minus; (1 &minus; \e M12 \e M21) \e m23 / \e m12
//   - \e M31 = \e M32 \e M21 &minus; (1 &minus; \e M23 \e M32) \e m12 / \e m23
//   
//   Additional functionality is provided by the GeodesicLine class, which
//   allows a sequence of points along a geodesic to be computed.
//   
//   The shortest distance returned by the solution of the inverse problem is
//   (obviously) uniquely defined.  However, in a few special cases there are
//   multiple azimuths which yield the same shortest distance.  Here is a
//   catalog of those cases:
//   - \e lat1 = &minus;\e lat2 (with neither point at a pole).  If \e azi1 =
//     \e azi2, the geodesic is unique.  Otherwise there are two geodesics and
//     the second one is obtained by setting [\e azi1, \e azi2] = [\e azi2, \e
//     azi1], [\e M12, \e M21] = [\e M21, \e M12], \e S12 = &minus;\e S12.
//     (This occurs when the longitude difference is near &plusmn;180&deg; for
//     oblate ellipsoids.)
//   - \e lon2 = \e lon1 &plusmn; 180&deg; (with neither point at a pole).  If
//     \e azi1 = 0&deg; or &plusmn;180&deg;, the geodesic is unique.  Otherwise
//     there are two geodesics and the second one is obtained by setting [\e
//     azi1, \e azi2] = [&minus;\e azi1, &minus;\e azi2], \e S12 = &minus;\e
//     S12.  (This occurs when \e lat2 is near &minus;\e lat1 for prolate
//     ellipsoids.)
//   - Points 1 and 2 at opposite poles.  There are infinitely many geodesics
//     which can be generated by setting [\e azi1, \e azi2] = [\e azi1, \e
//     azi2] + [\e d, &minus;\e d], for arbitrary \e d.  (For spheres, this
//     prescription applies when points 1 and 2 are antipodal.)
//   - s12 = 0 (coincident points).  There are infinitely many geodesics which
//     can be generated by setting [\e azi1, \e azi2] = [\e azi1, \e azi2] +
//     [\e d, \e d], for arbitrary \e d.
//   
//   The calculations are accurate to better than 15 nm (15 nanometers) for the
//   WGS84 ellipsoid.  See Sec. 9 of
//   <a href="http://arxiv.org/abs/1102.1215v1">arXiv:1102.1215v1</a> for
//   details.  The algorithms used by this class are based on series expansions
//   using the flattening \e f as a small parameter.  These are only accurate
//   for |<i>f</i>| &lt; 0.02; however reasonably accurate results will be
//   obtained for |<i>f</i>| &lt; 0.2.  Here is a table of the approximate
//   maximum error (expressed as a distance) for an ellipsoid with the same
//   major radius as the WGS84 ellipsoid and different values of the
//   flattening.<pre>
//       |f|      error
//       0.01     25 nm
//       0.02     30 nm
//       0.05     10 um
//       0.1     1.5 mm
//       0.2     300 mm
//   </pre>
//   For very eccentric ellipsoids, use GeodesicExact instead.
//   
//   The algorithms are described in
//   - C. F. F. Karney,
//     <a href="https://dx.doi.org/10.1007/s00190-012-0578-z">
//     Algorithms for geodesics</a>,
//     J. Geodesy <b>87</b>, 43--55 (2013);
//     DOI: <a href="https://dx.doi.org/10.1007/s00190-012-0578-z">
//     10.1007/s00190-012-0578-z</a>;
//     addenda: <a href="http://geographiclib.sf.net/geod-addenda.html">
//     geod-addenda.html</a>.
//
//   For more information on geodesics see \ref geodesic.
//   
//   Example of use:
//   \include example-Geodesic.cpp
//   
//   <a href="GeodSolve.1.html">GeodSolve</a> is a command-line utility
//   providing access to the functionality of Geodesic and GeodesicLine.
type Geodesic(semiMajorAxis : float<m>, flattening : LowToHighRatio) =
    let tiny = sqrt MathLib.minFloat
    let tol0 = MathLib.epsilon
    let tol1 = 200.0 * tol0
    let tol2 = sqrt tol0
    let tolb = tol0 * tol2
    let xthresh = 1000.0 * tol2
    let a = semiMajorAxis
    let f = flattening.Ratio
    let f1 = 1.0 - f
    let e2 = f * (2.0 - f)
    let ep2 = e2 / MathLib.sq(f1)
    let n = f / (2.0 - f)
    let b = a * f1
    let c2 = (MathLib.sq(a) + MathLib.sq(b) * MathLib.eatanhe(1.0, (if f < 0.0 then -1.0 else 1.0) * sqrt(abs(e2))) / e2) / 2.0
    let etol2 = 0.1 * tol2 / (sqrt(max 0.001 (abs f) * (min 1.0 1.0 - f/2.0) / 2.0))
    let mutable o = 0
    let A3x = [|0..GeodesicCoefficients.nA3 - 1|] |> Array.rev |> Array.mapi (fun k j ->
        let m = min (GeodesicCoefficients.nA3 - j - 1) j
        let value = MathLib.polyval(m, GeodesicCoefficients.A3Coeff.[o..], n) / GeodesicCoefficients.A3Coeff.[o + m + 1]
        o <- o + m + 2
        value)

    member this.A3 with get() = A3x