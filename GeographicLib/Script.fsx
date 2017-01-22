// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "..\packages\Newtonsoft.Json.9.0.1\lib\portable-net45+wp80+win8+wpa81\Newtonsoft.Json.dll"
#load "UnitsOfMeasure.fs"
#load "Constants.fs"
#load "MathLib.fs"
#load "GeodesicCoefficients.fs"
#load "Analytics.fs"

open GeographicLib

// Test a polynomial
// (x + 1)(x + 2) = x^2 + 3x + 2
// so p = [|2; 3; 1|], m = 2 and x = 1 gives a result of 6.
let shouldEqual6 = MathLib.polyval(2, [|2.0; 3.0; 1.0|], 1.0)
// (x + 1)(x + 2)(x + 3) = x^3 + 6x^2 + 11x + 6
// so p = [|6; 11; 6; 1|], m = 3 and x = 1 gives a result of 24.
let shouldEqual24 = MathLib.polyval(3, [|6.0; 11.0; 6.0; 1.0|], 1.0)
    
let geo = new Geodesic(Constants.WGS84_a, Constants.WGS84_f)
let s12 = geo.Distance (new GeodesicLocation(23.0<deg>, 117.0<deg>)) (new GeodesicLocation(34.0<deg>, 99.1<deg>))

let res = geo.Location (new GeodesicLocation(54.1<deg>, -0.02<deg>)) -90.0<deg> 1000.0<m>

printf "%s, %s" (res.Latitude.ToString()) (res.Longitude.ToString())

    
