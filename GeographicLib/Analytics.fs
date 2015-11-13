namespace GeographicLib

open System

[<Struct>]
type GeodesicLocation(latitude : float<deg>, longitude : float<deg>) =
    member this.Latitude = if abs latitude > 90.0<deg> then Double.NaN |> LanguagePrimitives.FloatWithMeasure else latitude
    member this.Longitude = longitude

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

[<System.FlagsAttribute>]
type PermissionFlags = 
    CapNone         = 0b0000000000000000
    | CapC1         = 0b0000000000000001
    | CapC1p        = 0b0000000000000010     
    | CapC2         = 0b0000000000000100
    | CapC3         = 0b0000000000001000
    | CapC4         = 0b0000000000010000
    | CapAll        = 0b0000000000011111
    | OutAll        = 0b0111111110000000
    | OutMask       = 0b1111111110000000
                    
type Mask =         
    None            = 0b0000000000000000
    | Latitude      = 0b0000000010000000
    | Longitude     = 0b0000000100001000
    | Azimuth       = 0b0000001000000000
    | Distance      = 0b0000010000000001
    | DistanceIn    = 0b0000100000000011
    | ReducedLength = 0b0001000000000101
    | GeodesicScale = 0b0010000000000101
    | Area          = 0b0100000000010000
    | LongUnroll    = 0b1000000000010000
    | All           = 0b0111111110011111

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
    let generatePolynomial nCoeff (coeff : float[]) =
        o <- 0
        [|0..nCoeff - 1|] |> Array.rev |> Array.mapi (fun k j ->
            let m = min (nCoeff - j - 1) j
            let value = MathLib.polyval(m, coeff.[o..], n) / coeff.[o + m + 1]
            o <- o + m + 2
            value)
        
    let A3x = generatePolynomial GeodesicCoefficients.nA3 GeodesicCoefficients.A3Coeff
    let C3x = generatePolynomial GeodesicCoefficients.nC3 GeodesicCoefficients.C3Coeff
    let C4x = generatePolynomial GeodesicCoefficients.nC4 GeodesicCoefficients.C4Coeff

    let Lengths (eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2) outMask (s12b : float byref, m12b : float byref, m0 : float byref, gs12 : float byref, gs21 : float byref) =
        let outMask = outMask &&& int PermissionFlags.OutMask
        let mutable m0x = 0.0
        let mutable J12 = 0.0 
        let mutable A1 = 0.0
        let mutable A2 = 0.0
        let mutable Ca = [||]
        let mutable Cb = [||]
        if outMask &&& int (Mask.Distance ||| Mask.ReducedLength ||| Mask.GeodesicScale) > 0 then
            Ca <- GeodesicCoefficients.C1Fourier eps
            let mutable Cb = [||]
            A1 <- GeodesicCoefficients.A1m1f(eps)
            if outMask &&& int (Mask.ReducedLength ||| Mask.GeodesicScale) > 0 then
                A2 <- GeodesicCoefficients.A2m1f(eps);
                Cb <- GeodesicCoefficients.C2Fourier eps;
                m0x <- A1 - A2;
                A2 <- 1.0 + A2;
            A1 <- 1.0 + A1;

        if outMask &&& int (Mask.Distance) > 0 then
            let B1 = MathLib.sinCosSeries(true, ssig2, csig2, Ca, GeodesicCoefficients.nC1) - MathLib.sinCosSeries(true, ssig1, csig1, Ca, GeodesicCoefficients.nC1)
            s12b <- A1 * (sig12 + B1)
            if outMask &&& int (Mask.ReducedLength ||| Mask.GeodesicScale) > 0 then
                let B2 = MathLib.sinCosSeries(true, ssig2, csig2, Cb, GeodesicCoefficients.nC2) - MathLib.sinCosSeries(true, ssig1, csig1, Cb, GeodesicCoefficients.nC2)
                J12 <- m0x * sig12 + (A1 * B1) - (A2 * B2)
        else if outMask &&& int (Mask.ReducedLength ||| Mask.GeodesicScale) > 0 then
            for l in 1..GeodesicCoefficients.nC2 do
                Cb.[l] <- A1 * Ca.[l] - A2 * Cb.[l]
            J12 <- m0x * sig12 + MathLib.sinCosSeries(true, ssig2, csig2, Cb, GeodesicCoefficients.nC2) - MathLib.sinCosSeries(true, ssig1, csig1, Cb, GeodesicCoefficients.nC2)

        if outMask &&& int (Mask.ReducedLength) > 0 then
            m0 <- m0x
            // Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
            // accurate cancellation in the case of coincident points.
            m12b <- dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12

        if outMask &&& int (Mask.GeodesicScale) > 0 then
            let csig12 = csig1 * csig2 + ssig1 * ssig2
            let t = ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
            gs12 <- csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1
            gs21 <- csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2

    let GenInverse (location1 : GeodesicLocation, location2 : GeodesicLocation) outMask (s12 : float<m> byref, azi1 : float<deg> byref, azi2 : float<deg> byref, m12 : float<m> byref, gs12 : float byref, gs21 : float byref, ga12 : float<m^2> byref) =
        let outMask = outMask &&& PermissionFlags.OutMask
        // Compute longitude difference (AngDiff does this carefully).  Result is
        // in [-180, 180] but -180 is only for west-going geodesics.  180 is for
        // east-going and meridional geodesics.
        // If very close to being on the same half-meridian, then make it so.
        let mutable lon12 = MathLib.angDiff(location1.Longitude, location2.Longitude) |> MathLib.angRound
        let mutable lonsign = if lon12 >= 0.0<deg> then 1.0 else -1.0
        lon12 <- lon12 * lonsign;
        // If really close to the equator, treat as on equator.
        let mutable lat1 = MathLib.angRound(location1.Latitude)
        let mutable lat2 = MathLib.angRound(location2.Latitude)
        // Swap points so that point with higher (abs) latitude is point 1
        let swapp = if abs(lat1) >= abs(lat2) then 1.0 else -1.0
        if swapp < 0.0 then
            let (l2, l1) = (lat1, lat2)
            lat1 <- l1
            lat2 <- l2
        // Make lat1 <= 0
        let latsign = if lat1 < 0.0<deg> then 1.0 else -1.0;
        lat1 <- lat1 * latsign;
        lat2 <- lat2 * latsign;

        let normalisedSinCos lat =
            let mutable sbet, cbet = MathLib.sincos lat
            sbet <- sbet * f1
            let norm = MathLib.norm (sbet, cbet)
            sbet <- fst norm
            cbet <- snd norm
            cbet <- max tiny cbet
            (sbet, cbet)
        // Now we have
        //
        //     0 <= lon12 <= 180
        //     -90 <= lat1 <= 0
        //     lat1 <= lat2 <= -lat1
        //
        // longsign, swapp, latsign register the transformation to bring the
        // coordinates to this canonical form.  In all cases, 1 means no change was
        // made.  We make these transformations so that there are few cases to
        // check, e.g., on verifying quadrants in atan2.  In addition, this
        // enforces some symmetries in the results returned.
        let mutable sbet1, cbet1 = normalisedSinCos lat1
        let mutable sbet2, cbet2 = normalisedSinCos lat2

        // If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
        // |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
        // a better measure.  This logic is used in assigning calp2 in Lambda12.
        // Sometimes these quantities vanish and in that case we force bet2 = +/-
        // bet1 exactly.  An example where is is necessary is the inverse problem
        // 48.522876735459 0 -48.52287673545898293 179.599720456223079643
        // which failed with Visual Studio 10 (Release and Debug)
        if cbet1 < -sbet1 then
            if cbet2 = cbet1 then sbet2 <- if sbet2 < 0.0 then sbet1 else -sbet1
        else
            if -sbet1 = abs sbet2 then cbet2 <- cbet1

        let dn1 = sqrt (1.0 + ep2 * MathLib.sq(sbet1))
        let dn2 = sqrt (1.0 + ep2 * MathLib.sq(sbet2))

        let lam12 = MathLib.radians lon12
        let slam12, clam12 = MathLib.sincos lon12
        let mutable a12, sig12, calp1, salp1, calp2, salp2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

        let mutable meridian = (lat1 = -90.0<deg> || slam12 = 0.0)
        if meridian then
            // Endpoints are on a single full meridian, so the geodesic might lie on
            // a meridian.
            calp1 <- clam12 
            salp1 <- slam12; // Head to the target longitude
            calp2 <- 1.0
            salp2 <- 0.0  // At the target we're heading north
            let ssig1, csig1 = sbet1, calp1 * cbet1
            let ssig2, csig2 = sbet2, calp2 * cbet2
            sig12 <- Math.Atan2(max (csig1 * ssig2 - ssig1 * csig2) 0.0, csig1 * csig2 + ssig1 * ssig2)
        0.0

    static member WGS84 = Geodesic(Constants.WGS84_a, Constants.WGS84_f)

    member this.A3 with get() = A3x    
    member this.C3 with get() = C3x    
    member this.C4 with get() = C4x