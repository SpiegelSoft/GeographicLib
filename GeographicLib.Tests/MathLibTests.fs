namespace GeographicLib.Tests

open GeographicLib
open Xunit
open Xunit.Extensions
open System

type ``Trig functions with degree arguments``()=

    static member Angles0To360 = [| 0.0<deg>; 30.0<deg>; 45.0<deg>; 60.0<deg>; 90.0<deg>; 120.0<deg>; 135.0<deg>; 180.0<deg>; 210.0<deg>; 225.0<deg>; 240.0<deg>; 270.0<deg>; 300.0<deg>; 315.0<deg>; 330.0<deg>; 360.0<deg> |] |> Seq.map (fun a -> [| a :> Object |])

    [<Fact>]
    member test.``(x + 1)(x + 2)``() =
        Assert.Equal(MathLib.polyval(2, [|2.0; 3.0; 1.0|], 1.0), 6.0)

    [<Fact>]
    member test.``(x + 1)(x + 2)(x + 3)``() =
        Assert.Equal(MathLib.polyval(3, [|6.0; 11.0; 6.0; 1.0|], 1.0), 24.0)

    [<Theory>]
    [<MemberData("Angles0To360")>]
    member test.``cos(x) = sin(90 - x)``(x : float<deg>) =
        let _, cosx = MathLib.sincos(x)
        let sin_90_x, _ = MathLib.sincos(90.0<deg> - x)
        Assert.Equal(cosx, sin_90_x, 14)

    [<Theory>]
    [<MemberData("Angles0To360")>]
    member test.``sin(x) = cos(90 - x)``(x : float<deg>) =
        let sinx, _ = MathLib.sincos(x)
        let _, cos_90_x = MathLib.sincos(90.0<deg> - x)
        Assert.Equal(sinx, cos_90_x, 14)

    [<Theory>]
    [<MemberData("Angles0To360")>]
    member test.``sincos(x) = sincos(360 + x)``(x : float<deg>) =
        let sinx, cosx = MathLib.sincos(x)
        let sin_360_x, cos_360_x = MathLib.sincos(360.0<deg> + x)
        Assert.Equal(cosx, cos_360_x, 15)
        Assert.Equal(sinx, sin_360_x, 15)

    [<Theory>]
    [<InlineData(0.0<deg>, 0.0, 1.0)>]
    [<InlineData(90.0<deg>, 1.0, 0.0)>]
    [<InlineData(180.0<deg>, 0.0, -1.0)>]
    [<InlineData(270.0<deg>, -1.0, 0.0)>]
    [<InlineData(360.0<deg>, 0.0, 1.0)>]
    member test.``Canonical sin values``(x : float<deg>, expectedsin, expectedcos : float) =
        let sinx, cosx = MathLib.sincos(x)
        Assert.Equal(sinx, expectedsin, 15)
        Assert.Equal(cosx, expectedcos, 15)
