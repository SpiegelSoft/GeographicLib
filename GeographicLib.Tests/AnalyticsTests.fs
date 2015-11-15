namespace GeographicLib.Tests

open GeographicLib
open Xunit
open Xunit.Extensions
open System

type ``Geodesic Distances``()=

    [<Fact>]
    member test.``Example distance``() =
        let geo = new Geodesic(Constants.WGS84_a, Constants.WGS84_f)
        let s12 = geo.Inverse (new GeodesicLocation(23.0<deg>, 117.0<deg>), new GeodesicLocation(34.0<deg>, 99.1<deg>))
        s12 |> ignore
