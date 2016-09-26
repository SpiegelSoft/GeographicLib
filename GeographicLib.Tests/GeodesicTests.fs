namespace GeographicLib.Tests

open GeographicLib
open Xunit
open Xunit.Extensions
open System
open System.IO

type ``Geodesic Distances``()=

    [<Theory>]
    [<MemberData("TestData")>]
    member test.``Loaded test data``(lat1, lat2, long2, expected) =
        let geo = Geodesic.WGS84
        let actual = geo.Distance (new GeodesicLocation(lat1, 0.0<deg>)) (new GeodesicLocation(lat2, long2))
        Assert.Equal(expected/1.0<m>, actual/1.0<m>, 4)

    static member TestData
        with get() =
            let readLines = seq {
                use ms = new MemoryStream(System.Text.Encoding.UTF8.GetBytes(GeographicLib.Tests.TestData.data))
                use sr = new StreamReader (ms)
                while not sr.EndOfStream do
                yield sr.ReadLine ()
            }
            seq { for line in (readLines) do 
                    let array = line.Split(' ')
                    yield [|array.[0] :> Object; array.[3] :> Object; array.[4] :> Object; array.[6] :> Object |]
            }
            |> Seq.toArray
        
