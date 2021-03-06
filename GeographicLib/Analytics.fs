﻿namespace GeographicLib

open System
open Newtonsoft.Json
open Newtonsoft.Json.Linq

type GeodesicLocationSerialiser() =
    inherit JsonConverter()
    override __.CanConvert(objectType) = objectType.Equals(typeof<GeodesicLocation>)
    override __.WriteJson(_, _, _) = raise <| new NotImplementedException()
    override __.CanWrite = false
    override __.ReadJson(reader, _, _, _) =
        let jsonObject = JObject.Load(reader)
        let latitude = jsonObject.GetValue("Latitude").Value<float>()
        let longitude = jsonObject.GetValue("Longitude").Value<float>()
        new GeodesicLocation(latitude * 1.0<deg>, longitude * 1.0<deg>) :> obj

and [<Struct; StructuralEquality; NoComparison; JsonConverter(typeof<GeodesicLocationSerialiser>)>] GeodesicLocation(latitude : float<deg>, longitude : float<deg>) =
    member __.Latitude = if abs latitude > 90.0<deg> then Double.NaN |> LanguagePrimitives.FloatWithMeasure else latitude
    member __.Longitude = longitude

module internal Analytics =

    let generateA3Polynomial n nCoeff (coeff : float[]) nArray =
        let mutable k, o = 0, 0
        let array = Array.zeroCreate nArray
        for j in (nCoeff - 1)..(-1)..0 do
            let m = min (nCoeff - j - 1) j
            array.[k] <- MathLib.polyval(m, coeff.[o..], n) / coeff.[o + m + 1]
            k <- k + 1
            o <- o + m + 2
        array

    let generateC3Polynomial n nCoeff (coeff : float[]) nArray =
        let mutable o, k = 0, 0
        let array = Array.zeroCreate nArray
        for l in [|1..nCoeff - 1|] do
            for j in Array.rev [|l..nCoeff - 1|] do
                let m = min (nCoeff - j - 1) j
                array.[k] <- MathLib.polyval(m, coeff.[o..], n) / coeff.[o + m + 1]
                k <- k + 1
                o <- o + m + 2
        array

    let generateC4Polynomial n nCoeff (coeff : float[]) nArray =
        let mutable o, k = 0, 0
        let array = Array.zeroCreate nArray
        for l in [|0..nCoeff - 1|] do
            for j in Array.rev [|l..nCoeff - 1|] do
                let m = nCoeff - j - 1
                array.[k] <- MathLib.polyval(m, coeff.[o..], n) / coeff.[o + m + 1]
                k <- k + 1
                o <- o + m + 2
        array

    let A3x n = generateA3Polynomial n GeodesicCoefficients.nA3 GeodesicCoefficients.A3Coeff GeodesicCoefficients.nA3x
    let C3x n = generateC3Polynomial n GeodesicCoefficients.nC3 GeodesicCoefficients.C3Coeff GeodesicCoefficients.nC3x
    let C4x n = generateC4Polynomial n GeodesicCoefficients.nC4 GeodesicCoefficients.C4Coeff GeodesicCoefficients.nC4x
    let A3f n eps = MathLib.polyval(GeodesicCoefficients.nA3 - 1, A3x n, eps)
    let C3f n eps (c : float[]) =
        // Evaluate C3 coeffs
        // Elements c[1] thru c[nC3_ - 1] are set
        let mutable mult = 1.0
        let mutable o = 0
        for l in [|1..(GeodesicCoefficients.nC3 - 1)|] do
            let m = GeodesicCoefficients.nC3 - l - 1 // order of polynomial in eps
            mult <- mult * eps
            c.[l] <- mult * MathLib.polyval(m, (C3x n).[o..], eps)
            o <- o + m + 1

    [<System.FlagsAttribute>]
    type internal PermissionFlags = 
        CapNone         = 0b0000000000000000
        | CapC1         = 0b0000000000000001
        | CapC1p        = 0b0000000000000010     
        | CapC2         = 0b0000000000000100
        | CapC3         = 0b0000000000001000
        | CapC4         = 0b0000000000010000
        | CapAll        = 0b0000000000011111
        | OutAll        = 0b0111111110000000
        | OutMask       = 0b1111111110000000
                    
    type internal Mask =         
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

open Analytics
type internal GeodesicLine(geodesic: Geodesic, location: GeodesicLocation, azimuth: float<deg>) = 
    let tiny = geodesic.Tiny
    let lat1, lon1 = location.Latitude, location.Longitude
    let azi1 = MathLib.angNormalise azimuth
    let salp1, calp1 = azi1 |> MathLib.angRound |> MathLib.sincos
    let f = geodesic.FlatteningRatio
    let b, f1, ep2 = geodesic.Parameters
    let mutable sbet1, cbet1 = lat1 |> MathLib.angRound |> MathLib.sincos
    do 
        sbet1 <- sbet1 * f1
        MathLib.norm &sbet1 &cbet1
        cbet1 <- max tiny cbet1
    let salp0, calp0 = salp1 * cbet1, MathLib.hypot(calp1, salp1 * sbet1)
    // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    // sig = 0 is nearest northward crossing of equator.
    // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
    // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    let somg1 = salp0 * sbet1
    let comg1 = if sbet1 <> 0.0 || calp1 <> 0.0 then cbet1 * calp1 else 1.0
    let mutable ssig1, csig1 = sbet1, comg1
    do MathLib.norm &ssig1 &csig1
    let k2 = calp0 * calp0 * ep2
    let eps = k2 / (2.0 * (1.0 + sqrt(1.0 + k2)) + k2)
    let A1m1 = GeodesicCoefficients.A1m1f(eps)
    let C1a = Array.create(GeodesicCoefficients.nC1 + 1) 0.0
    do GeodesicCoefficients.C1Fourier eps C1a
    let B11 = MathLib.sinCosSeries(true, ssig1, csig1, C1a, GeodesicCoefficients.nC1)
    let s, c = sin B11, cos B11
    let stau1 = ssig1 * c + csig1 * s
    let ctau1 = csig1 * c - ssig1 * s
    let C3a = Array.create GeodesicCoefficients.nC3 0.0
    let n = f / (2.0 - f)
    do C3f n eps C3a
    let A3c = -f * salp0 * A3f n eps
    let B31 = MathLib.sinCosSeries(true, ssig1, csig1, C3a, GeodesicCoefficients.nC3 - 1)
    member __.Position (distance: float<m>, lat2: byref<float<deg>>, lon2: byref<float<deg>>, azi2: byref<float<deg>>) =
        let tau12 = distance / (b * (1.0 + A1m1))
        let s, c = sin(tau12), cos(tau12)
        let C1pa = Array.create(GeodesicCoefficients.nC1p + 1) 0.0
        GeodesicCoefficients.C1pFourier eps C1pa
        let mutable B12 = - MathLib.sinCosSeries(true, stau1 * c + ctau1 * s, ctau1 * c - stau1 * s, C1pa, GeodesicCoefficients.nC1p)
        let mutable sig12 = tau12 - (B12 - B11)
        let mutable ssig12, csig12 = sin sig12, cos sig12
        if abs f > 0.01 then
            // Reverted distance series is inaccurate for |f| > 1/100, so correct
            // sig12 with 1 Newton iteration.  The following table shows the
            // approximate maximum error for a = WGS_a() and various f relative to
            // GeodesicExact.
            //     erri = the error in the inverse solution (nm)
            //     errd = the error in the direct solution (series only) (nm)
            //     errda = the error in the direct solution (series + 1 Newton) (nm)
            //
            //       f     erri  errd errda
            //     -1/5    12e6 1.2e9  69e6
            //     -1/10  123e3  12e6 765e3
            //     -1/20   1110 108e3  7155
            //     -1/50  18.63 200.9 27.12
            //     -1/100 18.63 23.78 23.37
            //     -1/150 18.63 21.05 20.26
            //      1/150 22.35 24.73 25.83
            //      1/100 22.35 25.03 25.31
            //      1/50  29.80 231.9 30.44
            //      1/20   5376 146e3  10e3
            //      1/10  829e3  22e6 1.5e6
            //      1/5   157e6 3.8e9 280e6
            let ssig2 = ssig1 * csig12 + csig1 * ssig12
            let csig2 = csig1 * csig12 - ssig1 * ssig12
            B12 <- MathLib.sinCosSeries(true, ssig2, csig2, C1a, GeodesicCoefficients.nC1)
            let serr = (1.0 + A1m1) * (sig12 + (B12 - B11)) - (distance / b)
            sig12 <- sig12 - serr / sqrt(1.0 + k2 * ssig2 * ssig2)
            ssig12 <- sin(sig12); csig12 <- cos(sig12)
        let ssig2 = ssig1 * csig12 + csig1 * ssig12
        let mutable csig2 = csig1 * csig12 - ssig1 * ssig12
        let sbet2 = calp0 * ssig2
        let mutable cbet2 = MathLib.hypot(salp0, calp0 * csig2)
        if cbet2 = 0.0 then 
            cbet2 <- tiny; csig2 <- tiny
        let salp2, calp2 = salp0, calp0 * csig2
        let somg2, comg2 = salp0 * ssig2, csig2
        let E = MathLib.copySign(1.0, salp0)
        let omg12 = E * (sig12 - ((atan2 ssig2 csig2) - (atan2 ssig1 csig1)) + (atan2(E * somg2) (comg2) - atan2(E * somg1) (comg1)))
        let lam12 = omg12 + A3c * (sig12 + (MathLib.sinCosSeries(true, ssig2, csig2, C3a, GeodesicCoefficients.nC3 - 1) - B31))
        let lon12 = lam12 * 1.0<rad> |> UnitConversion.degrees
        lon2 <- lon1 + lon12
        lat2 <- 1.0<rad> * atan2 sbet2 (f1 * cbet2) |> UnitConversion.degrees
        azi2 <- 1.0<rad> * atan2 salp2 calp2 |> UnitConversion.degrees

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
and Geodesic(semiMajorAxis : float<m>, flattening : LowToHighRatio) =
 
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
       
    let C4f eps (c : float[]) =
        // Evaluate C4 coeffs
        // Elements c[0] thru c[nC4_ - 1] are set
        let mutable mult = 1.0
        let mutable o = 0
        for l in [|0..(GeodesicCoefficients.nC4 - 1)|] do// l is index of C4[l]
            let m = GeodesicCoefficients.nC4 - l - 1 // order of polynomial in eps
            c.[l] <- mult * MathLib.polyval(m, (C4x n).[o..], eps)
            o <- o + m + 1
            mult <- mult * eps

    let Lengths (eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2) outMask (s12b : float byref) (m12b : float byref) (m0 : float byref) (gs12 : float byref) (gs21 : float byref) (Ca : float[]) =
        let outMask = outMask &&& int PermissionFlags.OutMask
        let mutable m0x = 0.0
        let mutable J12 = 0.0
        let mutable A1 = 0.0
        let mutable A2 = 0.0
        let Cb = Array.create (GeodesicCoefficients.nC2 + 1) 0.0
        if outMask &&& int (Mask.Distance ||| Mask.ReducedLength ||| Mask.GeodesicScale) > 0 then
            A1 <- GeodesicCoefficients.A1m1f(eps)
            GeodesicCoefficients.C1Fourier eps Ca
            if outMask &&& int (Mask.ReducedLength ||| Mask.GeodesicScale) > 0 then
                A2 <- GeodesicCoefficients.A2m1f(eps);
                GeodesicCoefficients.C2Fourier eps Cb;
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

    let Astroid x y =
        // Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
        // This solution is adapted from Geocentric::Reverse.
        let mutable k = 0.0
        let p, q = MathLib.sq(x), MathLib.sq(y)
        let r = (p + q - 1.0) / 6.0
        if not (q = 0.0 && r <= 0.0) then
            // Avoid possible division by zero when r = 0 by multiplying equations
            // for s and t by r^3 and r, resp.
            let S = p * q / 4.0 // S = r^3 * s
            let r2 = MathLib.sq(r)
            let r3 = r * r2
            // The discriminant of the quadratic equation for T3.  This is zero on
            // the evolute curve p^(1/3)+q^(1/3) = 1
            let disc = S * (S + 2.0 * r3)
            let mutable u = r
            if disc >= 0.0 then
                let mutable T3 = S + r3
                // Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
                // of precision due to cancellation.  The result is unchanged because
                // of the way the T is used in definition of u.
                T3 <- T3 + if T3 < 0.0 then -sqrt(disc) else sqrt(disc) // T3 = (r * t)^3
                // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
                let T = T3 ** (1.0 / 3.0) // T = r * t
                // T can be zero; but then r2 / T -> 0.
                u <- u + T + (if T <> 0.0 then r2 / T else 0.0)
            else
                // T is complex, but the way u is defined the result is real.
                let ang = atan2 (sqrt(-disc)) -(S + r3)
                // There are three possible cube roots.  We choose the root which
                // avoids cancellation.  Note that disc < 0 implies that r < 0.
                u <- u + 2.0 * r * cos(ang / 3.0);

            let v = sqrt(MathLib.sq(u) + q) // guaranteed positive
            // Avoid loss of accuracy when u < 0.
            let uv = if u < 0.0 then q / (v - u) else u + v // u+v, guaranteed positive
            let w = (uv - q) / (2.0 * v) // positive?
            // Rearrange expression for k to avoid loss of accuracy due to
            // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
            k <- uv / (sqrt(uv + MathLib.sq(w)) + w) // guaranteed positive
        else // q == 0 && r <= 0
            // y = 0 with |x| <= 1.  Handle this case directly.
            // for y small, positive root is k = abs(y)/sqrt(1-x^2)
            k <- 0.0
        k

    let InverseStart (sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12) (salp1 : float byref) (calp1 : float byref) (salp2 : float byref) (calp2 : float byref) (dnm : float byref) (Ca : float[]) =
        let mutable sig12 = -1.0
        let sbet12 = sbet2 * cbet1 - cbet2 * sbet1
        let cbet12 = cbet2 * cbet1 + sbet2 * sbet1
        let sbet12a = sbet2 * cbet1 + cbet2 * sbet1
        let shortline = cbet12 >= 0.0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5
        let mutable omg12 = lam12
        if shortline then
            let mutable sbetm2 = MathLib.sq(sbet1 + sbet2)
            // sin((bet1+bet2)/2)^2
            // =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
            sbetm2 <- sbetm2 / (sbetm2 + MathLib.sq(cbet1 + cbet2))
            dnm <- sqrt (1.0 + ep2 * sbetm2)
            omg12 <- omg12 / (f1 * dnm)

        let mutable somg12 = sin omg12
        let mutable comg12 = cos omg12
        salp1 <- cbet2 * somg12
        calp1 <- if comg12 >= 0.0 then sbet12 + (cbet2 * sbet1 * MathLib.sq somg12 / (1.0 + comg12)) else sbet12a - cbet2 * sbet1 * MathLib.sq somg12 / (1.0 - comg12)
        let ssig12 = MathLib.hypot(salp1, calp1)
        let csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12

        if shortline && ssig12 < etol2 then
            // really short lines
            salp2 <- cbet1 * somg12
            calp2 <- sbet12 - cbet1 * sbet2 * (if comg12 >= 0.0 then MathLib.sq(somg12) / (1.0 + comg12) else 1.0 - comg12)
            MathLib.norm &salp2 &calp2
            // Set return value
            sig12 <- atan2 ssig12  csig12
        else if abs(n) > 0.1 || csig12 >= 0.0 || ssig12 >= 6.0 * abs(n) * Math.PI * MathLib.sq(cbet1) then
            0 |> ignore // Nothing to do, zeroth order spherical approximation is OK
        else
            // Scale lam12 and bet2 to x, y coordinate system where antipodal point
            // is at origin and singular point is at y = 0, x = -1.
            let mutable y = 0.0
            let mutable lamscale = 0.0
            let mutable betscale = 0.0
            let mutable x = 0.0
            if f >= 0.0 then // In fact f == 0 does not get here
                // x = dlong, y = dlat
                let k2 = MathLib.sq(sbet1) * ep2
                let eps = k2 / (2.0 * (1.0 + sqrt(1.0 + k2)) + k2)
                lamscale <- f * cbet1 * (A3f n eps) * Math.PI
                betscale <- lamscale * cbet1
                x <- (lam12 - Math.PI) / lamscale
                y <- sbet12a / betscale
            else                  // _f < 0
                // x = dlat, y = dlong
                let cbet12a = cbet2 * cbet1 - sbet2 * sbet1
                let bet12a = atan2 sbet12a cbet12a
                let mutable m12b = 0.0
                let mutable m0 = 0.0
                let mutable dummy = 0.0
                // In the case of lon12 = 180, this repeats a calculation made in Inverse.
                Lengths (n, Math.PI + bet12a, sbet1, -cbet1, dn1, sbet2, cbet2, dn2, cbet1, cbet2) (int Mask.ReducedLength) &dummy &m12b &m0 &dummy &dummy Ca
                x <- -1.0 + m12b / (cbet1 * cbet2 * m0 * Math.PI)
                betscale <- if x < -0.01 then sbet12a / x else -f * MathLib.sq(cbet1) * Math.PI
                lamscale <- betscale / cbet1
                y <- (lam12 - Math.PI) / lamscale

            if y > -tol1 && x > -1.0 - xthresh then
                // strip near cut
                // Need real(x) here to cast away the volatility of x for min/max
                if f >= 0.0 then
                    salp1 <- min 1.0 -x
                    calp1 <- - sqrt(1.0 - MathLib.sq(salp1))
                else
                    calp1 <- max (if x > -tol1 then 0.0 else -1.0) x
                    salp1 <- sqrt(1.0 - MathLib.sq(calp1))
            else
                // Estimate alp1, by solving the astroid problem.
                //
                // Could estimate alpha1 = theta + pi/2, directly, i.e.,
                //   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
                //   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
                //
                // However, it's better to estimate omg12 from astroid and use
                // spherical formula to compute alp1.  This reduces the mean number of
                // Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
                // (min 0 max 5).  The changes in the number of iterations are as
                // follows:
                //
                // change percent
                //    1       5
                //    0      78
                //   -1      16
                //   -2       0.6
                //   -3       0.04
                //   -4       0.002
                //
                // The histogram of iterations is (m = number of iterations estimating
                // alp1 directly, n = number of iterations estimating via omg12, total
                // number of trials = 148605):
                //
                //  iter    m      n
                //    0   148    186
                //    1 13046  13845
                //    2 93315 102225
                //    3 36189  32341
                //    4  5396      7
                //    5   455      1
                //    6    56      0
                //
                // Because omg12 is near pi, estimate work with omg12a = pi - omg12
                let k = Astroid x y
                let omg12a = lamscale * (if f >= 0.0 then -x * k/(1.0 + k) else -y * (1.0 + k)/k)
                somg12 <- sin(omg12a)
                comg12 <- -cos(omg12a)
                // Update spherical estimate of alp1 using omg12 instead of lam12
                salp1 <- cbet2 * somg12
                calp1 <- sbet12a - cbet2 * sbet1 * MathLib.sq(somg12) / (1.0 - comg12)
        if not (salp1 <= 0.0) then 
            MathLib.norm &salp1 &calp1
        else
            salp1 <- 1.0
            calp1 <- 0.0

        sig12

    let Lambda12 (sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1) (salp2 : float byref) (calp2 : float byref) (sig12 : float byref) (ssig1 : float byref) (csig1 : float byref) (ssig2 : float byref) (csig2 : float byref) (eps : float byref) (domg12 : float byref) diffp (dlam12 : float byref) (Ca : float[]) =
        let calp1 = if sbet1 = 0.0 && calp1 = 0.0 then -tiny else calp1
        let salp0 = salp1 * cbet1
        let calp0 = MathLib.hypot(calp1, salp1 * sbet1) // calp0 > 0

        let mutable somg1, comg1, somg2, comg2, omg12, lam12 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        // tan(bet1) = tan(sig1) * cos(alp1)
        // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
        ssig1 <- sbet1
        somg1 <- salp0 * sbet1
        comg1 <- calp1 * cbet1
        csig1 <- comg1
        MathLib.norm &ssig1 &csig1
        // Math::norm(somg1, comg1); -- don't need to normalize!

        // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
        // about this case, since this can yield singularities in the Newton
        // iteration.
        // sin(alp2) * cos(bet2) = sin(alp0)
        salp2 <- if cbet2 <> cbet1 then salp0 / cbet2 else salp1
        // calp2 = sqrt(1 - sq(salp2))
        //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
        // and subst for calp0 and rearrange to give (choose positive sqrt
        // to give alp2 in [0, pi/2]).
        calp2 <- 
            if cbet2 <> cbet1 || abs sbet2 <> -sbet1 then
                sqrt(MathLib.sq(calp1 * cbet1) + (if cbet1 < -sbet1 then (cbet2 - cbet1) * (cbet1 + cbet2) else (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2
            else
                abs(calp1);
        // tan(bet2) = tan(sig2) * cos(alp2)
        // tan(omg2) = sin(alp0) * tan(sig2).
        ssig2 <- sbet2
        somg2 <- salp0 * sbet2
        comg2 <- calp2 * cbet2
        csig2 <- comg2
        MathLib.norm &ssig2 &csig2
        // Math::norm(somg2, comg2); -- don't need to normalize!

        // sig12 = sig2 - sig1, limit to [0, pi]
        sig12 <- atan2 (max(csig1 * ssig2 - ssig1 * csig2) 0.0) (csig1 * csig2 + ssig1 * ssig2)

        // omg12 = omg2 - omg1, limit to [0, pi]
        omg12 <- atan2(max(comg1 * somg2 - somg1 * comg2) 0.0) (comg1 * comg2 + somg1 * somg2)
        let mutable B312, h0 = 0.0, 0.0
        let k2 = MathLib.sq(calp0) * ep2
        eps <- k2 / (2.0 * (1.0 + sqrt(1.0 + k2)) + k2)
        C3f n eps Ca
        B312 <- (MathLib.sinCosSeries(true, ssig2, csig2, Ca, GeodesicCoefficients.nC3 - 1) - MathLib.sinCosSeries(true, ssig1, csig1, Ca, GeodesicCoefficients.nC3 - 1))
        h0 <- -f * A3f n eps
        domg12 <- salp0 * h0 * (sig12 + B312)
        lam12 <- omg12 + domg12

        if diffp then
            if calp2 = 0.0 then
                dlam12 <- -2.0 * f1 * dn1 / sbet1
            else
                let mutable dummy = 0.0
                Lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2) (int Mask.ReducedLength) &dummy &dlam12 &dummy &dummy &dummy Ca
                dlam12 <- dlam12 * f1 / (calp2 * cbet2)

        lam12

    let GenInverse (location1 : GeodesicLocation, location2 : GeodesicLocation) outMask (s12 : float<m> byref) (azi1 : float<deg> byref) (azi2 : float<deg> byref) (m12 : float<m> byref) (gs12 : float byref) (gs21 : float byref) (ga12 : float<m^2> byref) =
        let maxit1 = GeodesicCoefficients.maxit1
        let maxit2 = maxit1 + 64 + 10
        let outMask = outMask &&& int PermissionFlags.OutMask
        // Compute longitude difference (AngDiff does this carefully).  Result is
        // in [-180, 180] but -180 is only for west-going geodesics.  180 is for
        // east-going and meridional geodesics.
        // If very close to being on the same half-meridian, then make it so.
        let mutable lon12 = MathLib.angDiff(location1.Longitude, location2.Longitude) |> MathLib.angRound
        let mutable lon12s = lon12
        let mutable lonsign = if lon12 >= 0.0<deg> then 1.0 else -1.0
        lon12 <- lon12 * lonsign;
        // If very close to being on the same half-meridian, then make it so.
        lon12s <- MathLib.angRound((180.0<deg> - lon12) - lonsign * lon12s)
        let lam12 = UnitConversion.radians lon12
        let slam12, clam12 = 
            match lon12 > 90.0<deg> with
            | true -> 
                let sincos = MathLib.sincos lon12s; 
                (fst sincos, -snd sincos)
            | false -> MathLib.sincos lon12
        // If really close to the equator, treat as on equator.
        let mutable lat1 = MathLib.angRound(location1.Latitude)
        let mutable lat2 = MathLib.angRound(location2.Latitude)
        // Swap points so that point with higher (abs) latitude is point 1
        let swapp = if abs(lat1) >= abs(lat2) then 1.0 else -1.0
        if swapp < 0.0 then 
            Utilities.swap &lat1 &lat2
            lonsign <- -1.0 * lonsign

        // Make lat1 <= 0
        let latsign = if lat1 < 0.0<deg> then 1.0 else -1.0;
        lat1 <- lat1 * latsign;
        lat2 <- lat2 * latsign;

        let normalisedSinCos lat =
            let mutable sbet, cbet = MathLib.sincos lat
            sbet <- sbet * f1
            MathLib.norm &sbet &cbet
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
        let mutable s12b = 0.0
        let mutable m12b = 0.0
        let mutable s12x = 0.0<m>
        let mutable m12x = 0.0<m>

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

        let mutable a12, sig12, calp1, salp1, calp2, salp2 = 0.0<deg>, 0.0, 0.0, 0.0, 0.0, 0.0

        let mutable meridian = (lat1 = -90.0<deg> || slam12 = 0.0)
        let Ca = Array.create GeodesicCoefficients.nC 0.0
        if meridian then
            // Endpoints are on a single full meridian, so the geodesic might lie on
            // a meridian.
            calp1 <- clam12 
            salp1 <- slam12; // Head to the target longitude
            calp2 <- 1.0
            salp2 <- 0.0  // At the target we're heading north
            let ssig1, csig1 = sbet1, calp1 * cbet1
            let ssig2, csig2 = sbet2, calp2 * cbet2
            sig12 <- atan2 (max (csig1 * ssig2 - ssig1 * csig2) 0.0) csig1 * csig2 + ssig1 * ssig2
            let mutable dummy = 0.0
            Lengths (n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2) (outMask ||| int (Mask.Distance ||| Mask.ReducedLength)) &s12b &m12b &dummy &gs12 &gs21 Ca
            if sig12 < 1.0 || m12x >= 0.0<m> then
                // Need at least 2, to handle 90 0 90 180
                if sig12 < 3.0 * tiny then
                    s12x <- 0.0<m>
                    m12x <- 0.0<m>
                    sig12 <- 0.0
                m12x <- m12x * b / 1.0<m>
                s12x <- s12x * b / 1.0<m>
                a12 <- MathLib.degrees sig12
            else
                // m12 < 0, i.e., prolate and too close to anti-podal
                meridian <- false
        
        let mutable omg12 = 0.0
        if not meridian && sbet1 = 0.0 && (f <= 0.0 || lam12 <= Math.PI - f * Math.PI) then
            calp2 <- 0.0
            calp1 <- 0.0
            salp1 <- 1.0
            salp2 <- 1.0
            s12x <- a * lam12
            omg12 <- lam12 / f1
            sig12 <- omg12
            m12x <- b * sin sig12
            if outMask &&& int Mask.GeodesicScale > 0 then
                gs21 <- cos sig12
                gs21 <- gs12
            a12 <- lon12 / f1
        else if not meridian then
            let mutable dnm = 0.0
            sig12 <- InverseStart (sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12) &salp1 &calp1 &salp2 &calp2 &dnm Ca
            if sig12 >= 0.0 then
                // Short lines (InverseStart sets salp2, calp2, dnm)
                s12x <- sig12 * b * dnm
                m12x <- MathLib.sq(dnm) * b * sin(sig12 / dnm)
                if outMask &&& int Mask.GeodesicScale > 0 then
                    gs21 <- cos(sig12 / dnm)
                    gs12 <- gs21
                a12 <- MathLib.degrees sig12
                omg12 <- lam12 / (f1 * dnm)
            else
                // Newton's method.  This is a straightforward solution of f(alp1) =
                // lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
                // root in the interval (0, pi) and its derivative is positive at the
                // root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
                // alp1.  During the course of the iteration, a range (alp1a, alp1b) is
                // maintained which brackets the root and with each evaluation of
                // f(alp) the range is shrunk, if possible.  Newton's method is
                // restarted whenever the derivative of f is negative (because the new
                // value of alp1 is then further from the solution) or if the new
                // estimate of alp1 lies outside (0,pi); in this case, the new starting
                // guess is taken to be (alp1a + alp1b) / 2.
                //
                // initial values to suppress warnings (if loop is executed 0 times)
                let mutable ssig1, csig1, ssig2, csig2, eps = 0.0, 0.0, 0.0, 0.0, 0.0
                // Bracketing range
                let mutable salp1a, calp1a, salp1b, calp1b = tiny, 1.0, tiny, -1.0
                let mutable tripn, tripb = false, false
                let rec newton numit converged =
                    if not converged then
                        if numit > maxit2 then raise <| new ArithmeticException("Newton's method failed to converge.")
                        // the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
                        // WGS84 and random input: mean = 2.85, sd = 0.60
                        let mutable dv = 0.0
                        let v = (Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1) &salp2 &calp2 &sig12 &ssig1 &csig1 &ssig2 &csig2 &eps &omg12 (numit < maxit1) &dv Ca) - lam12
                        // 2 * tol0 is approximately 1 ulp for a number in [0, pi].
                        // Reversed test to allow escape with NaNs
                        let mutable useNextMidpoint = true
                        if tripb || not (abs(v) >= (if tripn then 8.0 else 2.0) * tol0) then 
                            newton (numit + 1) true
                        else
                            // Update bracketing values
                            if v > 0.0 && (numit > maxit1 || calp1/salp1 > calp1b/salp1b) then
                                salp1b <- salp1
                                calp1b <- calp1
                            else if v < 0.0 && (numit > maxit1 || calp1/salp1 < calp1a/salp1a) then
                                salp1a <- salp1
                                calp1a <- calp1
                            if numit < maxit1 && dv > 0.0 then
                                let dalp1 = -v/dv;
                                let sdalp1, cdalp1 = sin(dalp1), cos(dalp1)
                                let nsalp1 = salp1 * cdalp1 + calp1 * sdalp1
                                if nsalp1 > 0.0 && abs(dalp1) < Math.PI then
                                    calp1 <- calp1 * cdalp1 - salp1 * sdalp1
                                    salp1 <- nsalp1;
                                    MathLib.norm &salp1 &calp1
                                    // In some regimes we don't get quadratic convergence because
                                    // slope -> 0.  So use convergence conditions based on epsilon
                                    // instead of sqrt(epsilon).
                                    tripn <- abs(v) <= 16.0 * tol0
                                    useNextMidpoint <- false
                                    newton (numit + 1) false

                            if useNextMidpoint then
                                // Either dv was not postive or updated value was outside legal
                                // range.  Use the midpoint of the bracket as the next estimate.
                                // This mechanism is not needed for the WGS84 ellipsoid, but it does
                                // catch problems with more eccentric ellipsoids.  Its efficacy is
                                // such for the WGS84 test set with the starting guess set to alp1 =
                                // 90deg:
                                // the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
                                // WGS84 and random input: mean = 4.74, sd = 0.99
                                salp1 <- (salp1a + salp1b)/2.0
                                calp1 <- (calp1a + calp1b)/2.0
                                MathLib.norm &salp1 &calp1
                                tripn <- false
                                tripb <- abs(salp1a - salp1) + (calp1a - calp1) < tolb || abs(salp1 - salp1b) + (calp1 - calp1b) < tolb
                                newton (numit + 1) false
                    else
                        0 |> ignore
                
                newton 1 false
                let mutable dummy = 0.0
                // Ensure that the reduced length and geodesic scale are computed in
                // a "canonical" way, with the I2 integral.
                let lengthMask = outMask ||| (if outMask &&& int (Mask.ReducedLength ||| Mask.GeodesicScale) > 0 then int Mask.Distance else int Mask.None)
                Lengths (eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2) lengthMask &s12b &m12b &dummy &gs12 &gs21 Ca
                m12x <- m12b * b
                s12x <- s12b * b;
                a12 <- MathLib.degrees sig12
                omg12 <- lam12 - omg12
        
        if outMask &&& int Mask.Distance > 0 then
            s12 <- 0.0<m> + s12x // Convert -0 to 0

        if outMask &&& int Mask.ReducedLength > 0 then
            m12 <- 0.0<m> + m12x // Convert -0 to 0

        if outMask &&& int Mask.Area > 0 then
            let salp0 = salp1 * cbet1
            let calp0 = MathLib.hypot(calp1, salp1 * sbet1) // calp0 > 0
            let mutable alp12 = 0.0
            if (calp0 <> 0.0 && salp0 <> 0.0) then
                // From Lambda12: tan(bet) = tan(sig) * cos(alp)
                let mutable ssig1 = sbet1
                let mutable csig1 = calp1 * cbet1
                let mutable ssig2 = sbet2
                let mutable csig2 = calp2 * cbet2
                let k2 = MathLib.sq(calp0) * ep2
                let eps = k2 / (2.0 * (1.0 + sqrt(1.0 + k2)) + k2)
                // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
                let A4 = MathLib.sq(a) * calp0 * salp0 * e2
                MathLib.norm &ssig1 &csig1
                MathLib.norm &ssig2 &csig2
                C4f eps Ca
                let B41 = MathLib.sinCosSeries(false, ssig1, csig1, Ca, GeodesicCoefficients.nC4)
                let B42 = MathLib.sinCosSeries(false, ssig2, csig2, Ca, GeodesicCoefficients.nC4)
                ga12 <- A4 * (B42 - B41)
            else
                // Avoid problems with indeterminate sig1, sig2 on equator
                ga12 <- 0.0<m^2>

            if not meridian && omg12 < 0.75 * Math.PI && sbet2 - sbet1 < 1.75 then // Lat difference too big
                // Use tan(Gamma/2) = tan(omg12/2)
                // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
                // with tan(x/2) = sin(x)/(1+cos(x))
                let somg12, domg12 = sin(omg12), 1.0 + cos(omg12)
                let dbet1, dbet2 = 1.0 + cbet1, 1.0 + cbet2
                alp12 <- 2.0 * atan2 (somg12 * (sbet1 * dbet2 + sbet2 * dbet1 )) (domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ))
            else
                // alp12 = alp2 - alp1, used in atan2 so no need to normalize
                let mutable salp12 = salp2 * calp1 - calp2 * salp1
                let mutable calp12 = calp2 * calp1 + salp2 * salp1
                // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
                // salp12 = -0 and alp12 = -180.  However this depends on the sign
                // being attached to 0 correctly.  The following ensures the correct
                // behavior.
                if salp12 = 0.0 && calp12 < 0.0 then
                    salp12 <- tiny * calp1
                    calp12 <- -1.0

                alp12 <- atan2 salp12 calp12
            ga12 <- ga12 + c2 * alp12;
            ga12 <- ga12 * swapp * lonsign * latsign;
            // Convert -0 to 0
            ga12 <- ga12 + 0.0<m^2>

        // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
        if swapp < 0.0 then
            Utilities.swap &salp1 &salp2
            Utilities.swap &calp1 &calp2
            if outMask &&& int Mask.GeodesicScale > 0 then
                Utilities.swap &gs12 &gs21

        salp1 <- salp1 * swapp * lonsign
        calp1 <- calp1 * swapp * latsign
        salp2 <- salp2 * swapp * lonsign
        calp2 <- calp2 * swapp * latsign

        if outMask &&& int Mask.Azimuth > 0 then
            azi1 <- MathLib.atan2d salp1 calp1
            azi2 <- MathLib.atan2d salp2 calp2

        // Returned value in [0, 180]
        a12

    static member WGS84 = Geodesic(Constants.WGS84_a, Constants.WGS84_f)

    member val internal FlatteningRatio = f
    member val internal Parameters: float<m> * float * float = b, f1, ep2
    member val internal Tiny = tiny
    member __.Distance location1 location2 =
        let mutable t, azi1, azi2, tm, s12, a = 0.0, 0.0<deg>, 0.0<deg>, 0.0<m>, 0.0<m>, 0.0<m^2>
        let _ = GenInverse (location1, location2) (int Mask.Distance) &s12 &azi1 &azi2 &tm &t &t &a
        s12
    member __.Azimuths location1 location2 =
        let mutable t, azi1, azi2, tm, s12, a = 0.0, 0.0<deg>, 0.0<deg>, 0.0<m>, 0.0<m>, 0.0<m^2>
        let _ = GenInverse (location1, location2) (int Mask.Distance) &s12 &azi1 &azi2 &tm &t &t &a
        azi1, azi2
    member this.Location (location: GeodesicLocation) (azimuth: float<deg>) (distance: float<m>) =
        let geodesicLine = new GeodesicLine(this, location, azimuth)
        let mutable lat2, lon2, azi2 = 0.0<deg>, 0.0<deg>, 0.0<deg>
        geodesicLine.Position(distance, &lat2, &lon2, &azi2)
        new GeodesicLocation(lat2, lon2)
