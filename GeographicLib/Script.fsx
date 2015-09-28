// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#load "UnitsOfMeasure.fs"
#load "Constants.fs"
#load "MathLib.fs"
#load "GeodesicCoefficients.fs"
#load "Analytics.fs"

open GeographicLib

let wgs84Ellipsoid = Ellipsoid.WGS84
// Test a polynomial
// (x + 1)(x + 2) = x^2 + 3x + 2
// so p = [|2; 3; 1|], m = 2 and x = 1 gives a result of 6.
let shouldEqual6 = MathLib.polyval(2, [|2.0; 3.0; 1.0|], 1.0)
// (x + 1)(x + 2)(x + 3) = x^3 + 6x^2 + 11x + 6
// so p = [|6; 11; 6; 1|], m = 3 and x = 1 gives a result of 24.
let shouldEqual24 = MathLib.polyval(3, [|6.0; 11.0; 6.0; 1.0|], 1.0)
    
let geo = new Geodesic(Constants.WGS84_a, Constants.WGS84_f)



    
