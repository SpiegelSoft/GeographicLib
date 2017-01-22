namespace GeographicLib.Tests

open GeographicLib

open Xunit

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
    static let testData = seq {
        use ms = new MemoryStream(System.Text.Encoding.UTF8.GetBytes(GeographicLib.Tests.TestData.data))
        use sr = new StreamReader (ms)
        while not sr.EndOfStream do
        yield sr.ReadLine ()
    }

    [<Theory>]
    [<MemberData("InverseTestData")>]
    member __.``Zero distances``(lat, _: float, long, _: float) =
        let location = new GeodesicLocation(lat * 1.0<deg>, long * 1.0<deg>)
        let result = Geodesic.WGS84.Distance location location
        Assert.Equal(result/1.0<m>, 0.0, 8)

    [<Theory>]
    [<MemberData("InverseTestData")>]
    member __.``Inverse Problem`` (lat1, lat2, long2, expected) =
        let actual = Geodesic.WGS84.Distance (new GeodesicLocation(lat1, 0.0<deg>)) (new GeodesicLocation(lat2, long2))
        Assert.Equal(expected/1.0<m>, actual/1.0<m>, 4)

    [<Theory>]
    [<MemberData("ForwardLocationTestData")>]
    member __.``Forward Problem (Latitude, Longitude)`` (lat1, long1, azi1, distance, expectedLat, expectedLong) =
        let actual = Geodesic.WGS84.Location (new GeodesicLocation(lat1, long1)) (azi1 * 1.0<deg>) (distance * 1.0<m>)
        Assert.Equal(actual.Latitude/1.0<deg>, expectedLat, 10)
        Assert.Equal(actual.Longitude/1.0<deg>, expectedLong, 10)

    static member InverseTestData
        with get() =
            seq { for line in testData do 
                    let array = line.Split(' ')
                    yield [|array.[0] :> obj; array.[3] :> obj; array.[4] :> obj; array.[6] :> obj |]
            } |> Seq.toArray

    static member ForwardLocationTestData
        with get() =
            seq { for line in testData do 
                    let array = line.Split(' ')
                    yield [|array.[0] :> obj; array.[1] :> obj; array.[2] :> obj; array.[6] :> obj; array.[3] :> obj; array.[4] :> obj |]
            } |> Seq.toArray
        
