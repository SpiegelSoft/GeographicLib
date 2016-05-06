# GeographicLib
A port of GeographicLib into F#

My goodness, this package has had more downloads than I first expected.  I hope those of you who have downloaded it have had the wherewithal to look at the tests in the source code to work out how to actually use the functionality.

The geodesic function we expose at present is the Distance function.  If you want another function then please feel free to log an issue.

You can use this package in a C# project, but then you won't reap the fabulous benefits of units of measure.

How to calculate distance between two locations:

```FSharp
open GeographicLib

let geo = Geodesic.WGS84
let landsEndLocation = new GeodesicLocation(50.0775475<deg>, -5.6352355<deg>)
let johnOGroatsLocation = new GeodesicLocation(58.6366688<deg>, -3.0827024<deg>)
let distance = geo.Distance landsEndLocation johnOGroatsLocation
```
