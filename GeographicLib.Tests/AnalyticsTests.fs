namespace GeographicLib.Tests

open GeographicLib
open Xunit
open Xunit.Extensions
open System

type ``Geodesic Distances``()=

    [<Theory>]
    [<InlineData(36.530042355041<deg>, -48.164270779097768864<deg>, 5.762344694676510456<deg>, 9398502.0434687<m>)>]
    [<InlineData(62.912967770911<deg>, -23.186512533703375483<deg>, 68.567247430960081525<deg>, 11254667.2007643<m>)>]
    [<InlineData(2.881248229541<deg>, 53.997072295385487038<deg>, 44.520619105667619923<deg>, 6958264.1576889<m>)>]
    member test.``Example distance``(lat1, lat2, long2, expected) =
        let geo = new Geodesic(Constants.WGS84_a, Constants.WGS84_f)
        let actual = geo.Inverse (new GeodesicLocation(lat1, 0.0<deg>), new GeodesicLocation(lat2, long2))
        Assert.Equal(expected/1.0<m>, actual/1.0<m>, 6)
