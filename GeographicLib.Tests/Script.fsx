// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.

#r "bin/Release/GeographicLib.dll"

open GeographicLib

let geo = Geodesic.WGS84
let landsEndLocation = new GeodesicLocation(50.0657<deg>, -5.7113<deg>)
let johnOGroatsLocation = new GeodesicLocation(58.6373<deg>, 3.0689<deg>)
let distance = geo.Distance landsEndLocation johnOGroatsLocation




