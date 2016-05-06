# GeographicLib
A port of GeographicLib into F#

My goodness, this package has had more downloads than I first expected.  I hope those of you who have downloaded it have had the wherewithal to look at the source code for the tests to work out how to actually use the functionality.

The only function this exposes at present is the Distance function.  If you want another function then please feel free to log an issue.

How to calculate distance between two locations:

```FSharp
open GeographicLib

let geo = Geodesic.WGS84
let landsEndLocation = new GeodesicLocation(50.0775475<deg>, -5.6352355<deg>)
let johnOGroatsLocation = new GeodesicLocation(58.6366688<deg>, -3.0827024<deg>)
let distance = geo.Distance landsEndLocation johnOGroatsLocation
```
