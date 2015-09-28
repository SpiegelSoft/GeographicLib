namespace GeographicLib.Tests

open GeographicLib
open Xunit
open Xunit.Extensions

type ``Polynomial Evaluation``()=

    [<Fact>]
    member test.``(x + 1)(x + 2)``()=
        Assert.Equal(MathLib.polyval(2, [|2.0; 3.0; 1.0|], 1.0), 6.0)
