namespace GeographicLib.Tests

open GeographicLib
open Xunit
open System
open System.IO
open Newtonsoft.Json

type ``GeodesicLocation serialisation``()=
    [<Fact>]
    member __.``Forwards and backwards serialisation of GeodesicLocation struct``() =
        let location = new GeodesicLocation(51.5194894<deg>, -0.0524782<deg>)
        let serialisedLocation = JsonConvert.SerializeObject(location)
        let deserialisedLocation = JsonConvert.DeserializeObject<GeodesicLocation>(serialisedLocation)
        Assert.Equal(location, deserialisedLocation)

type ``Geodesic Distances``() =
    [<Theory>]
    [<MemberData("TestData")>]
    member __.``Zero distances``(lat, _: float, long, _: float) =
        let location = new GeodesicLocation(lat * 1.0<deg> - 0.001<deg>, long * 1.0<deg> + 0.0001<deg>)
        let result = Geodesic.WGS84.Distance location location
        Assert.Equal(result/1.0<m>, 0.0, 8)

    [<Theory>]
    [<MemberData("TestData")>]
    member __.``Loaded test data`` (lat1, lat2, long2, expected) =
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
        
