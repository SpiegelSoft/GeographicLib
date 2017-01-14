namespace GeographicLib

open System

[<Measure>] type rad
[<Measure>] type deg
[<Measure>] type km
[<Measure>] type m

module UnitConversion =
    let metres (distance: float<km>) = distance / 1.0<km> * 1000.0<m>
    let kilometres (distance: float<m>) = distance / 1000.0<m> * 1.0<km>
    let private degreesPerRadian = 1.0<deg/rad> * 180.0/Math.PI
    let private radiansPerDegree = 1.0</deg> * Math.PI/180.0
    let radians (x : float<deg>) = x * radiansPerDegree
    let degrees (x : float<rad>) = x * degreesPerRadian
