namespace GeographicLib

open System

module MathLib =
    [<Struct>]
    type ValueAndError(value : float<deg>, error : float<deg>) =
        member this.Value = sum
        member this.Error = error

    let liftToKeepUnits f (arg : float<'u>) : float<'u> =
        float arg |> f |> LanguagePrimitives.FloatWithMeasure
    let liftToKeepUnits2 f (arg1 : float<'u>, arg2 : float<'u>) : float<'u> =
        (float arg1, float arg2) |> f |> LanguagePrimitives.FloatWithMeasure
    let remainder = Math.IEEERemainder |> liftToKeepUnits2    
    let sum(u : float<deg>, v : float<deg>) =
      let s = u + v;
      let mutable up = s - v;
      let mutable vpp = s - up;
      up <- up - u;
      vpp <- vpp - v;
      new ValueAndError(s, -(up + vpp))

    let atanh(x) = (log(1.0 + x) - log(1.0 - x))/2.0
    let rec getMachineEpsilon eps = if 1.0 + eps = 1.0 then eps else getMachineEpsilon ((abs eps) / 2.0)
    let eatanhe(x, es) = if es >= 0.0 then es * atanh(es * x) else -es * atan(es * x)
    let sq(x : float<'a>) = x * x
    let minFloat = Double.Epsilon
    let epsilon = getMachineEpsilon 0.1
    let epsilonSquared = sq(epsilon)
    let epsilonToTheFourth = sq(epsilonSquared)
    let tolerance = sqrt epsilon
    let tolerance0 = tolerance * epsilon |> sqrt |> sqrt
    let piOver2 = Math.PI / 2.0
    let polyval(N, p : float[], x) =
        let mutable y = if N < 0 then 0.0 else p.[0]
        let mutable N = N
        let mutable i = 1
        while N > 0 do
            N <- N - 1
            y <- y * x + p.[i]
            i <- i + 1
        y
    let z = 1.0/16.0
    let angRound(x) =
        let mutable y = abs x
        y <- if y < z then z - (z - y) else y
        match Constants.Precision with
        | 4 -> if x <= 0.0 then 0.0 - y else y
        | 5 -> if x < 0.0 then 0.0 - y else y
        | _ -> if x < 0.0 then 0.0 - y else y

    let angNormalise(x : float<deg>) =
        let reverse180 angle = if angle = 180.0<deg> then -180.0<deg> else angle 
        reverse180(remainder(x, 360.0<deg>))
    let angDiff(x : float<deg>, y : float<deg>) =
        let differenceOfRemainders = sum(remainder(x, 360.0<deg>), remainder(-y, 360.0<deg>))
        let normalisedDiff = angNormalise differenceOfRemainders.Value
        (if normalisedDiff = 180.0<deg> && differenceOfRemainders.Error < 0.0<deg> then -180.0<deg> else normalisedDiff) - differenceOfRemainders.Error



